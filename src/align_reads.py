"""
align_reads.py
==============
Módulo para el alineamiento de lecturas RNA-seq curadas contra el genoma 
de referencia utilizando HISAT2, y procesamiento de alineamientos con samtools.

Funciones principales:
    build_index             Construye el índice de HISAT2 (omite si ya existe).
    convert_sam_to_bam      Convierte formato SAM a BAM ordenado e indexado (.bai), 
                            realizando la limpieza del archivo SAM temporal.
    setup_logging           Configura el logging del módulo.
    my_parser               Define los argumentos de línea de comandos.
    main                    Ejecuta el flujo completo.

Argumentos:
    -r, --reference         Ruta al genoma de referencia (FASTA).
    -i, --index             Ruta base y prefijo para el índice HISAT2.
    -gtf, --annotation      Ruta al archivo de anotación GTF (opcional).
    -d, --data_dir          Directorio de entrada con las lecturas curadas.
    -o, --outdir            Directorio de salida para los BAM.
    -w, --workers           Hilos de procesamiento (default: 6).
    -u, --umbral            Umbral mínimo de tasa de alineamiento (default 70.0)
    -v, --verbose           Muestra mensajes DEBUG.
    --log_file              Guarda el archivo de log.
"""

import argparse as ap
import glob
import logging
import os
import re
import subprocess
import sys

logger = logging.getLogger(__name__)

def build_index(fasta_path: str = "data/genome/genome_Anopheles_gambiae.fasta", 
                index_prefix: str = "data/index/hisat2/genome", 
                threads: int = 6, 
                gtf_path: str | None = "data/genome/genome_Anopheles_gambiae.gtf") -> None:
    """
    Construye el índice de HISAT2 a partir del genoma de referencia de Anopheles gambiae.
    Si se proporciona un archivo GTF, extrae los sitios de empalme (splice sites) y 
    exones para optimizar el alineamiento de lecturas de ARN.

    Args:
        fasta_path (str): Ruta al archivo FASTA del genoma (default: data/genome/genome_Anopheles_gambiae.fasta).
        index_prefix (str): Directorio y prefijo del índice (default: data/index/hisat2/genome).
        threads (int): Número de hilos para construir el índice (default: 6).
        gtf_path (str | None): Ruta opcional al archivo GTF de anotación (default: data/genome/genome_Anopheles_gambiae.gtf).

    Returns:
        None
    """
    if os.path.exists(f"{index_prefix}.1.ht2"):
        logger.info("El índice de HISAT2 ya existe en %s. Omitiendo construcción.", index_prefix)
        return

    if not os.path.exists(fasta_path):
        logger.warning("No se encontró el genoma FASTA en: %s", fasta_path)
        raise FileNotFoundError(f"Archivo FASTA faltante: {fasta_path}")

    os.makedirs(os.path.dirname(index_prefix), exist_ok=True)

    # se contruye en una lista el comando minimo que ocupa hisat2-build
    cmd = ["hisat2-build", "-p", str(threads)]

    if gtf_path and os.path.exists(gtf_path):
        logger.info("Archivo GTF detectado. Extrayendo sitios de empalme y exones...")
        ss_file = f"{index_prefix}.ss"
        exon_file = f"{index_prefix}.exon"

        try:
            with open(ss_file, "w") as ss_out:
                subprocess.run(["hisat2_extract_splice_sites.py", gtf_path], 
                               check=True, stdout=ss_out, stderr=subprocess.PIPE, text=True)
            with open(exon_file, "w") as exon_out:
                subprocess.run(["hisat2_extract_exons.py", gtf_path], 
                               check=True, stdout=exon_out, stderr=subprocess.PIPE, text=True)
                
            logger.info("Sitios de empalme y exones extraídos correctamente.")

            # Añadimos los flags al comando principal
            cmd.extend(["--ss", ss_file, "--exon", exon_file])

        except subprocess.CalledProcessError as e:
            logger.error("Error al procesar el GTF con utilidades de HISAT2: %s", e.stderr)
            raise RuntimeError(f"Fallo en extracción del GTF: {e.stderr}") from e
    elif gtf_path:
        logger.warning("Se proporcionó la ruta del GTF pero el archivo no existe: %s. Se omitirá.", gtf_path)

    # Se completa el comando del hisat2-build con el FASTA y el prefijo de salida
    cmd.extend([fasta_path, index_prefix])

    try:
        subprocess.run(cmd, check=True, stderr=subprocess.PIPE, text=True)
        logger.info("Índice de HISAT2 construido exitosamente en %s", index_prefix)
    except subprocess.CalledProcessError as e:
        logger.error("Error al construir el índice con hisat2-build: %s", e.stderr)
        raise RuntimeError(f"Fallo en hisat2-build: {e.stderr}") from e

def convert_sam_to_bam(sam_path: str, bam_path: str, threads: int = 6) -> None:
    """
    Convierte un archivo SAM a un archivo BAM ordenado por coordenadas e 
    inmediatamente genera su índice (.bai) usando samtools.

    Args:
        sam_path (str): Ruta del archivo SAM de entrada.
        bam_path (str): Ruta del archivo BAM de salida.
        threads (int): Número de hilos para el ordenamiento.

    Returns:
        None
    """
    if not os.path.exists(sam_path):
        logger.warning("No se encontró el archivo SAM para convertir: %s", sam_path)
        raise FileNotFoundError(f"Archivo SAM faltante: {sam_path}")

    os.makedirs(os.path.dirname(bam_path), exist_ok=True)
    logger.info("Convirtiendo SAM a BAM ordenado: %s", bam_path)

    # Comando para ordenar directamente desde SAM a BAM
    cmd_sort = ["samtools", "sort", "-@", str(threads), "-o", bam_path, sam_path]
    # Comando para indexar el BAM
    cmd_index = ["samtools", "index", "-@", str(threads), bam_path]
    
    try:
        subprocess.run(cmd_sort, check=True, stderr=subprocess.PIPE, text=True)
        logger.info("BAM ordenado generado con éxito.")
        
        logger.info("Generando índice para %s (.bai)", bam_path)
        subprocess.run(cmd_index, check=True, stderr=subprocess.PIPE, text=True)
        logger.info("Índice .bai generado con éxito.")

        # Limpieza: Borra el SAM original para evitar una saturación de espacio en donde se ejecuta el pipeline
        logger.info("Eliminando archivo SAM temporal: %s", sam_path)
        os.remove(sam_path)

    except subprocess.CalledProcessError as e:
        logger.error("Error en la ejecución de samtools: %s", e.stderr)
        raise RuntimeError(f"samtools falló durante la conversión: {e.stderr}") from e

def get_alignment_rate(summary_dir: str, suffix: str = ".summary.txt", umbral: float = 70.0) -> list:
    """
    Lee los resúmenes de HISAT2 de un directorio y extrae la tasa de alineamiento global
    de cada muestra. Posteriormente, mostrará una advertencia sí la tasa de alineamiento
    de la muestra es menor al umbral seleccionado (por defecto 70.0)

    Args:
        summary_dir (str) -> Directorio con los archivos de resumen (uno por muestra,
            generados con --summary-file de hisat2).
        suffix (str) -> Sufijo de los archivos de resumen; se quita para obtener el
            nombre de la muestra.
        umbral (int) -> Umbral mínimo para el filtro por tasa de alineamiento de cada muestra

    Returns:
        tasas (dict[str, float]) -> Diccionario {muestra: tasa de alineamiento en %}.
    """
    pattern = re.compile(r"([\d.]+)% overall alignment rate")
    tasas = {}

    for ruta in sorted(glob.glob(os.path.join(summary_dir, f"*{suffix}"))):
        muestra = os.path.basename(ruta).removesuffix(suffix)

        try:
            with open(ruta, encoding="utf-8") as f:
                contenido = f.read()
        except OSError as e:
            logger.error("No se pudo leer el resumen %s: %s", ruta, e)
            raise RuntimeError(f"Fallo al leer {ruta}: {e}") from e

        m = pattern.search(contenido)
        if m is None:
            logger.warning("No se encontró la tasa de alineamiento en %s", ruta)
            raise ValueError(f"No se encontró la tasa de alineamiento en {ruta}")
        try:
            tasas[muestra] = float(m.group(1))
        except ValueError as e:
            logger.error("Tasa inválida en %s: %s", ruta, m.group(1))
            raise ValueError(f"Tasa de alineamiento no numérica en {ruta}: {m.group(1)}") from e

    if not tasas:
        logger.warning("No hay archivos *%s en %s", suffix, summary_dir)
        raise FileNotFoundError(f"No hay archivos *{suffix} en {summary_dir}")

    for muestra, t in tasas.items():
        if t < umbral:
            logger.warning("Muestra %s con tasa de alineamiento baja: %.2f%% (umbral %.1f%%)",
                           muestra, t, umbral)

    return tasas

def setup_logging(verbose: bool = False, log_file: str | None = None) -> None:
    """
    Configura el logging del módulo: consola (stdout) y, opcionalmente, archivo.

    Args:
        verbose (bool) -> Si es True, muestra mensajes DEBUG
        log_file (str | None) -> Ruta del archivo de log (None para no guardarlo)

    Returns:
        None
    """
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)
        handlers.append(logging.FileHandler(log_file, encoding="utf-8"))

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=handlers,
    )

def my_parser() -> ap.Namespace:
    """
    Define y procesa los argumentos de línea de comandos requeridos para el alineamiento.

    Returns:
        ap.Namespace: Argumentos indicados por el usuario.
    """
    parser = ap.ArgumentParser(description="Pipeline de alineamiento RNA-seq con HISAT2")
    parser.add_argument("-r", "--reference", required=True,
                        help="Ruta al genoma de referencia (FASTA)")
    parser.add_argument("-i", "--index", required=True, 
                        help="Prefijo para el índice de HISAT2")
    parser.add_argument("-gtf", "--annotation", required=False, 
                        default= None, 
                        help="Ruta al archivo de anotación GTF (opcional para el índice)")
    parser.add_argument("-d", "--data_dir", required=False, 
                        default="data/curated/", 
                        help="Directorio raíz de las lecturas curadas")
    parser.add_argument("-o", "--outdir", required=False, 
                        default="data/aligned/", 
                        help="Directorio de salida para los archivos BAM")
    parser.add_argument("-w", "--workers", type=int, required= False, default= None, 
                        help="Número de hilos para hisat2 y samtools (default: 6)")
    parser.add_argument("-u", "--umbral", type=float, default=70.0,
                        help="Tasa de alineamiento mínima de chequeo de las muestras (default: 70.0)")
    parser.add_argument("-v", "--verbose", action="store_true", 
                        help="Muestra mensajes DEBUG en consola")
    parser.add_argument("--log_file", default="logs/align_reads.log", 
                        help="Ruta para guardar el archivo de log")
    
    return parser.parse_args()
    

def main():
    """
    Ejecuta el flujo de alineamiento de lecturas RNA-seq contra el genoma de referencia
    con HISAT2 y procesa los alineamientos con samtools.

    Pasos, en orden:
        1. Construye el índice de HISAT2 a partir del genoma de referencia (omite la
           construcción si el índice ya existe). Si se proporciona un GTF, extrae los
           sitios de empalme y exones para optimizar el alineamiento.
        2. Alinea con HISAT2 cada par de lecturas curadas (_val_1.fq, _val_2.fq) y
           guarda el resumen de alineamiento de cada muestra.
        3. Convierte cada archivo SAM a BAM ordenado e indexado (.bai) y elimina el SAM
           temporal.
        4. Lee los resúmenes de HISAT2 y registra una advertencia por cada muestra cuya
           tasa de alineamiento esté por debajo del umbral.
    """

    args = my_parser()
    setup_logging(args.verbose, args.log_file)
    logger.info("Inicializando módulo de alineamiento de lecturas")

    try:
        # Construir el indice
        build_index(args.reference, args.index, args.workers, args.annotation)

        # Aquí se colocará la función "run_hisat2"

        # Convierte el archivo SAM a formato BAM
        # convert_sam_to_bam(sam_path, bam_path)

        # Revisa las tasas de alineamiento de las muestras
        get_alignment_rate(summary_dir=args.outdir, umbral=args.umbral)
    except (FileNotFoundError, RuntimeError) as e:
        logger.error("El pipeline se detuvo por un error crítico en la indexación: %s.", e)
        sys.exit(1)

if __name__ == "__main__":
    main()