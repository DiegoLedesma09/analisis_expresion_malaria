"""
align_reads.py
==============
Módulo para el alineamiento de lecturas RNA-seq curadas contra el genoma 
de referencia utilizando HISAT2, y procesamiento de alineamientos con samtools.

Requisitos:
    - HISAT2 (hisat2, hisat2-build disponibles en el PATH)
    - samtools (disponible en el PATH)
    - Lecturas pareadas limpias (_val_1.fq, _val_2.fq)

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
    -w, --workers           Hilos de procesamiento (default: 8).
    -v, --verbose           Muestra mensajes DEBUG.
    --log_file              Guarda el archivo de log.
"""

import argparse as ap
import logging
import os
import subprocess
import sys

logger = logging.getLogger(__name__)

def build_index(fasta_path: str, index_prefix: str, threads: int = 6, gtf_path: str | None = None) -> None:
    """
    Construye el índice de HISAT2 a partir del genoma de referencia de Anopheles gambiae.
    Si se proporciona un archivo GTF, extrae los sitios de empalme (splice sites) y 
    exones para optimizar el alineamiento de lecturas de ARN.

    Args:
        fasta_path (str): Ruta al archivo FASTA del genoma.
        index_prefix (str): Directorio y prefijo del índice (ej. data/index/hisat2/genome).
        threads (int): Número de hilos para construir el índice.
        gtf_path (str | None): Ruta opcional al archivo GTF de anotación.

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

def my_parse() -> ap.Namespace:
    """
    Define y procesa los argumentos de línea de comandos requeridos para el alineamiento.

    Returns:
        ap.Namespace: Argumentos indicados por el usuario.
    """
    parser = ap.ArgumentParser(description="Pipeline de alineamiento RNA-seq con HISAT2")
    parser.add_argument("-r", "--reference", required=False, 
                        default="data/genome/genome_Anopheles_gambiae.fasta", 
                        help="Ruta al genoma de referencia (FASTA)")
    parser.add_argument("-i", "--index", required=False, 
                        default="data/index/hisat2/genome", 
                        help="Prefijo para el índice de HISAT2")
    parser.add_argument("-gtf", "--annotation", required=False, 
                        default="data/genome/genome_Anopheles_gambiae.gtf", 
                        help="Ruta al archivo de anotación GTF (opcional para el índice)")
    parser.add_argument("-d", "--data_dir", required=False, 
                        default="data/curated/", 
                        help="Directorio raíz de las lecturas curadas")
    parser.add_argument("-o", "--outdir", required=False, 
                        default="data/aligned/", 
                        help="Directorio de salida para los archivos BAM")
    parser.add_argument("-w", "--workers", type=int, default=6, 
                        help="Número de hilos para hisat2 y samtools (default: 6)")
    parser.add_argument("-v", "--verbose", action="store_true", 
                        help="Muestra mensajes DEBUG en consola")
    parser.add_argument("--log_file", default="logs/align_reads.log", 
                        help="Ruta para guardar el archivo de log")
    
    return parser.parse_args()
    

def main():
    print("Inicializando módulo / Esto se va a borrar")
    args = my_parse()
    setup_logging(args.verbose, args.log_file)
    logger.info("Inicializando módulo de alineamiento de lecturas")


    # 1. Construir el indice
    try:
        build_index(args.reference, args.index, args.workers, args.annotation)
    except (FileNotFoundError, RuntimeError) as e:
        logger.error("El pipeline se detuvo por un error crítico en la indexación.")
        sys.exit(1)

if __name__ == "__main__":
    main()