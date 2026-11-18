"""
counts_features.py
==============

Módulo para la creación de la tabla de conteos de los genes por muestra a partir de los
archivos .bam generados por el módulo de alineamiento (align_reads.py) creados por HISAT2.

Usa featureCounts para hacer la matriz de conteos y procede a limpiar el resultado para
regresar una matriz de genes X muestras.

Finalmente, regresa un plot que representa cuántas muestras asociadas tiene cada gen.

Requisitos:
    - Archivos .bam generados por HISAT2 (align_reads.py module).
    - Archivo .gtf extraído desde NCBI (download_data.py module).
    - Ambiente de Conda (malaria_env) activado.

Funciones principales:
    find_bams           Localiza los BAM por grupo.
    run_featurecounts   Ejecuta featureCounts para la creación de la matriz de conteos
    parse_summary       Lee el .summary y devuelve el porcentaje de lecturas asignadas por muestra.
    clean_counts        Deja una matriz limpia (genes x muestras)
    plot_assigment      Gráfica opcional de lecturas asignadas por muestra.
    my_parser           Define los argumentos de línea de comandos.
    setup_logging       Configura el logging del módulo.
    main                Ejecuta el flujo completo.

# Agregar archivos generados aquí

Argumentos:
    -w,  --workers      Número de hilos para featureCounts (default: 6).
    -o,  --outdir       Directorio de salida (default: results).
    -v,  --verbose      Muestra mensajes DEBUG.
    --log_file          Ruta del archivo de log.
    -g,  --graphic      Incluye el gráfico comparativo del trimming (default).
         --no-graphic   Omite el gráfico comparativo del trimming.

# Agregar uso aquí
"""

import argparse as ap
import logging
import os
import sys
import glob
import subprocess
import pandas as pd

logger = logging.getLogger(__name__)

def find_bams(data_dir:str, grupo: str) -> list:
    """
    Localiza los archivos BAM ordenados e indexados de un grupo específico.

    Args: 
        data_dir (str): Directorio donde están los archivos BAM (ej. "data/aligned/"). 
        grupo (str): Nombre del grupo a buscar.
    
    Returns:
        list: Lista con las rutas de los archivos .bam encontrados. 
    """
    logger.info("Buscando archivos BAM en el directorio %s para el grupo %s", data_dir, grupo)

    dir_grupo = os.path.join(data_dir, grupo)

    if not os.path.isdir(dir_grupo):
        raise FileNotFoundError(
        f"No se encontró la carpeta del grupo '{grupo}' en {data_dir}"
    )
    bam_files = sorted(glob.glob(os.path.join(dir_grupo, "*.bam")))
    if not bams:
        logger.warning("No se encontraron archivos .bam en %s", dir_grupo)
    else:
        logger.info("Se encontraron %d archivos BAM en el grupo '%s'", len(bams), grupo)
        
    return bam_files

def run_featureCounts(bam_files: list, gtf_path: str, outdir: str, workers: int) -> None:
    """
    Ejecuta featureCounts para generar la matriz de conteos.

    Args:
        bam_files (list): Lista de rutas a los archivos .bam.
        gtf_path (str): Ruta al archivo de anotación .gtf.
        outdir (str): Directorio de salida.
        workers (int): Número de hilos.

    Returns:
        None
    """
    if not bam_files:
            logger.error("La lista de BAMs está vacía. No se puede ejecutar featureCounts.")
            return
    
    if not os.path.exists(gtf_path):
            logger.error("No se encontró el archivo GTF en: %s", gtf_path)
            raise FileNotFoundError(f"Falta el GTF: {gtf_path}")
    
    command = [
        "featureCounts",
        "-T", str(workers),
        "-p",
        "-B",
        "-C",
        "-t", "exon",
        "-g", "gene_id",
        "-a", gtf_path,
        "-o", outdir,
        *bam_files
    ]
    
    logger.info("Ejecutando featureCounts con %d hilos", workers)
    
    try:
        subprocess.run(command, check=True, stderr=subprocess.PIPE, text=True)
        logger.info("featureCounts completado. Resultados guardados en %s", outdir)
    except subprocess.CalledProcessError as e:
        logger.error("Error al ejecutar featureCounts: %s", e)
        raise

def parse_summary():
    print("Luisosoma")

def clean_counts(path_matrix:str, outdir:str) -> pd.DataFrame:
    """
    Limpia la matriz de conteos realizada por featureCounts, dejando una matriz (genes x muestras)

    Args:
        path_matrix (str) -> Ruta a la matriz de conteos realizada por featureCounts
        outdir (str) -> Ruta hacia el directorio de salida

    Returns:
        clean_matrix (pd.DataFrame) -> Matriz de conteos limpia (genes x muestras)
    """

    logger.info("Iniciando limpieza de la matriz de conteos generada por featureCounts")

    try:
        df = pd.read_csv(path_matrix, sep="\t", comment="#")

        sample_columns = [col for col in df.columns if col not in ["Geneid", "Chr", "Start", "End", "Strand", "Length"]]
        clean_matrix = df.set_index("Geneid")[sample_columns]

        clean_matrix.columns = [col.split('/')[-1].replace('.bam', '') for col in clean_matrix.columns]

        os.makedirs(outdir, exist_ok=True)
        output_file = os.path.join(outdir, "clean_matrix.tsv")
        clean_matrix.to_csv(output_file, sep="\t")

        logger.info("Se ha realizado la limpieza de la matriz realizada por featureCounts y se ha guardado en %s", output_file)

        return clean_matrix

    except FileNotFoundError as e:
        logger.warning("No se ha encontrado el archivo en la ruta %s", e)
    except ValueError as e:
        logger.warning("La limpieza de la matriz a fallado %s", e)
        raise ValueError(f"Se detecto el error {e}")

def plot_assigment():
    print("Luisosoma")

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
    parser = ap.ArgumentParser(description="Pipeline para la creación de la matriz de conteos con featureCounts")
    parser.add_argument("-o", "--outdir", required=False, default="data/counts/", 
                        help="Directorio de salida")
    parser.add_argument("-w", "--workers", type=int, required= False, default= 6, 
                        help="Número de hilos")
    parser.add_argument("-v", "--verbose", action="store_true", 
                        help="Muestra mensajes DEBUG en consola")
    parser.add_argument("-g", "--graphic", action=ap.BooleanOptionalAction, default=True,
                    help="Incluye el grafico de la distribución de muestras asignadas")
    parser.add_argument("--log_file", default="logs/feature_counts.log", 
                        help="Ruta para guardar el archivo de log")
    
    return parser.parse_args()

def main():
    """
    Ejecuta el flujo completo de creación de la matriz de conteos.

    Pasos, en orden:
        1. Encuentra los archivos .bam asociados a las lecturas
        2. Corre featureCounts y crea la matriz de conteos.
        3. Lee el archivo .summary y devuelve el porcentaje de lecturas asociadas por muestra.
        4. Limpia la matriz de conteos regresando una matriz de genes x muestras.
        5. Si --graphic está activo, grafica la distribución de lecturas asociadas por muestra

    Los argumentos se leen de la línea de comandos mediante my_parser().
    """
    args = my_parser()
    setup_logging(args.verbose, args.log_file)
    logger.info("Inicializando módulo de creación de tabla de conteos")

    find_bams()

    run_featureCounts()

    parse_summary()

    clean_counts(outdir=args.outdir)

    if args.graphic:
        plot_assigment()