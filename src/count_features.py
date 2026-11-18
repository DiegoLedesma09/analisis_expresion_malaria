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
    plot_assignment      Gráfica opcional de lecturas asignadas por muestra.
    my_parser           Define los argumentos de línea de comandos.
    setup_logging       Configura el logging del módulo.
    main                Ejecuta el flujo completo.

Archivos generados:
    data/counts/raw_counts.txt                  Salida cruda de featureCounts.
    data/counts/raw_counts.txt.summary          Resumen estadístico de asignación.
    data/counts/counts_matrix.tsv               Matriz limpia de conteos (genes x muestras).
    data/counts/samples_metadata.tsv            Tabla de metadatos (muestra y condición).
    results/images/assignment_proportion.png    Gráfica de lecturas asignadas por muestra.
    logs/count_features.log                     Archivo de registro de ejecución.

Argumentos:
    -d,   --data_dir    Directorio raíz con los BAM por condición (default: data/aligned/).
    -gtf, --annotation  Ruta al archivo de anotación GTF.
    -c1,  --condition1  Condición 1 (default: resistant).
    -c2,  --condition2  Condición 2 (default: susceptible).
    -s,   --strand      Especificidad de cadena: 0 (unstranded), 1 (forward), 2 (reverse) (default: 0).
    -w,  --workers      Número de hilos para featureCounts (default: 6).
    -o,  --outdir       Directorio de salida (default: results).
    -u,   --umbral      Umbral mínimo de lecturas asignadas en % (default: 70.0).
    -v,  --verbose      Muestra mensajes DEBUG.
    --log_file          Ruta del archivo de log.
    -g,  --graphic      Incluye el gráfico comparativo del conteo (default).
         --no-graphic   Omite el gráfico comparativo del conteo.

Uso:
    python src/count_features.py -gtf data/genome/genome_Anopheles_gambiae.gtf -v
"""

import argparse as ap
import glob
import logging
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
        logger.error("No se encontró la carpeta del grupo %s en %s", grupo, data_dir)
        raise FileNotFoundError(f"No se encontró la carpeta del grupo '{grupo}' en {data_dir}")

    bam_files = sorted(glob.glob(os.path.join(dir_grupo, "*.bam")))
    if not bam_files:
        logger.error("No se encontraron archivos .bam en %s", dir_grupo)
        raise FileNotFoundError(f"No se encontraron archivos .bam en {dir_grupo}")

    else:
        logger.info("Se encontraron %d archivos BAM en el grupo '%s'", len(bam_files), grupo)
        
    return bam_files

def run_featureCounts(bam_files: list, gtf_path: str, output_file: str, workers: int, strand: int = 0) -> None:
    """
    Ejecuta featureCounts para generar la matriz de conteos.

    Args:
        bam_files (list): Lista de rutas a los archivos .bam.
        gtf_path (str): Ruta al archivo de anotación .gtf.
        output_file (str): Ruta del archivo de salida de featureCounts.
        workers (int): Número de hilos.
        strand (int): Especificidad de cadena: 0 (unstranded), 1 (forward), 2 (reverse).

    Returns:
        None
    """
    if not bam_files:
            logger.error("La lista de BAMs está vacía. No se puede ejecutar featureCounts.")
            raise ValueError("No se proporcionaron archivos BAM para contar.")
    
    if not os.path.exists(gtf_path):
            logger.error("No se encontró el archivo GTF en: %s", gtf_path)
            raise FileNotFoundError(f"Falta el GTF: {gtf_path}")
    
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    
    command = [
        "featureCounts",
        "-T", str(workers),
        "-p", "--countReadPairs",
        "-B",
        "-C",
        "-s", str(strand),
        "-t", "exon",
        "-g", "gene_id",
        "-a", gtf_path,
        "-o", output_file,
        *bam_files
    ]
    
    logger.info("Ejecutando featureCounts con %d hilos", workers)
    
    try:
        subprocess.run(command, check=True, stderr=subprocess.PIPE, text=True)
        logger.info("featureCounts completado. Resultados guardados en %s", output_file)
    except subprocess.CalledProcessError as e:
        logger.error("Error al ejecutar featureCounts: %s", e.stderr)
        raise RuntimeError(f"Fallo en featureCounts: {e.stderr}") from e

def parse_summary(summary_file: str, umbral: float = 70.0) -> dict[str, float]:
    """
    Lee el archivo .summary generado por featureCounts, calcula el porcentaje de
    lecturas asignadas por muestra y emite una advertencia si alguna queda por
    debajo del umbral establecido.

    Args:
        summary_file (str): Ruta al archivo .summary generado por featureCounts.
        umbral (float): Umbral mínimo esperado de lecturas asignadas en % (default: 70.0).

    Returns:
        dict[str, float]: Diccionario {muestra: porcentaje_asignado}.
    """
    if not os.path.exists(summary_file):
        logger.error("No se encontró el resumen: %s", summary_file)
        raise FileNotFoundError(f"Resumen faltante: {summary_file}")

    try:
        df = pd.read_csv(summary_file, sep="\t", index_col=0)
    except (OSError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        logger.error("No se pudo leer el resumen %s: %s", summary_file, e)
        raise RuntimeError(f"Fallo al leer el archivo de resumen {summary_file}: {e}") from e

    if df.empty or "Assigned" not in df.index:
        logger.error("El archivo %s está vacío o no contiene la fila 'Assigned'", summary_file)
        raise ValueError(f"Formato inválido en el resumen de featureCounts: {summary_file}")

    assignment_rates: dict[str, float] = {}

    for col in df.columns:
        sample_name = os.path.basename(col).removesuffix(".bam")
        assigned = df.loc["Assigned", col]
        total = df[col].sum()

        if total <= 0:
            logger.error("La muestra %s tiene 0 lecturas totales en %s", sample_name, summary_file)
            raise ValueError(f"Conteo total igual a cero para la muestra {sample_name}")

        pct = float((assigned / total) * 100)
        assignment_rates[sample_name] = pct

        logger.info(
            "Muestra %s: %.2f%% de lecturas asignadas (%d / %d)",
            sample_name,
            pct,
            assigned,
            total,
        )

        if pct < umbral:
            logger.warning(
                "Muestra %s con porcentaje de asignación bajo: %.2f%% (umbral %.1f%%)",
                sample_name,
                pct,
                umbral,
            )

    return assignment_rates 

def clean_counts(raw_counts_file: str, outdir: str = "data/counts/clean") -> pd.DataFrame:
    """
    Limpia la matriz de conteos realizada por featureCounts, dejando una matriz (genes x muestras)

    Args:
        raw_counts_file (str) -> Ruta a la matriz de conteos cruda realizada por featureCounts
        outdir (str) -> Ruta hacia el directorio de salida

    Returns:
        clean_matrix (pd.DataFrame) -> Matriz de conteos limpia (genes x muestras)
    """

    logger.info("Iniciando limpieza de la matriz de conteos generada por featureCounts")

    if not os.path.exists(raw_counts_file):
        logger.error("No se encontró la matriz cruda: %s", raw_counts_file)
        raise FileNotFoundError(f"Matriz faltante: {raw_counts_file}")

    df = pd.read_csv(raw_counts_file, sep="\t", comment="#").set_index("Geneid")
    clean_matrix = df.drop(columns=["Chr", "Start", "End", "Strand", "Length"])
    clean_matrix.columns = [os.path.basename(c).removesuffix(".bam") for c in clean_matrix.columns]

    os.makedirs(outdir, exist_ok=True)
    output_file = os.path.join(outdir, "clean_matrix.tsv")
    clean_matrix.to_csv(output_file, sep="\t")

    logger.info("Se ha realizado la limpieza de la matriz realizada por featureCounts y se ha guardado en %s", output_file)

    return clean_matrix

def plot_assignment(assignment_rates: dict[str, float], outdir: str = "results/images/", umbral: float = 70.0) -> None:
    """
    Genera y guarda una gráfica de barras con el porcentaje de lecturas asignadas
    por muestra, resaltando aquellas que estén por debajo del umbral.

    Args:
        assignment_rates (dict[str, float]): Diccionario {muestra: porcentaje_asignado}.
        outdir (str): Directorio de salida de la imagen (default: results/images/).
        umbral (float): Umbral mínimo de referencia en % (default: 70.0).

    Returns:
        None
    """
    if not assignment_rates:
        logger.warning("No hay tasas de asignación disponibles para graficar.")
        raise ValueError("El diccionario de tasas de asignación está vacío.")

    os.makedirs(outdir, exist_ok=True)
    output_path = os.path.join(outdir, "assignment_proportion.png")

    samples = list(assignment_rates.keys())
    rates = list(assignment_rates.values())
    colors = ["#1b9e77" if r >= umbral else "#d95f02" for r in rates]

    fig, ax = plt.subplots(figsize=(max(6, len(samples) * 1.4), 4.5))
    try:
        bars = ax.bar(samples, rates, color=colors, width=0.55)

        for bar, r in zip(bars, rates):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 1.5,
                f"{r:.1f}%",
                ha="center",
                va="bottom",
                fontsize=9,
            )

        ax.axhline(
            umbral,
            color="#7570b3",
            linestyle="--",
            linewidth=1.2,
            label=f"Umbral mínimo ({umbral:.1f}%)",
        )

        ax.set_ylabel("Lecturas asignadas (%)")
        ax.set_xlabel("Muestra")
        ax.set_title("Porcentaje de lecturas asignadas por muestra (featureCounts)")
        ax.set_ylim(0, 112)
        ax.tick_params(axis="x", rotation=20)
        ax.legend(loc="upper right")

        fig.tight_layout()
        fig.savefig(output_path, dpi=300)
        logger.info("Gráfica de asignación guardada en: %s", output_path)

    except OSError as e:
        logger.error("No se pudo guardar la gráfica en %s: %s", output_path, e)
        raise RuntimeError(f"Fallo al guardar la gráfica en {output_path}: {e}") from e
    finally:
        plt.close(fig)

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
    parser.add_argument("-d", "--data_dir", required=False, default="data/aligned/",
                        help="Directorio raíz con los BAM organizados por condición (default: data/aligned/)")
    parser.add_argument("-gtf", "--annotation", required=False,
                        default="data/genome/genome_Anopheles_gambiae.gtf",
                        help="Ruta al archivo de anotación GTF")
    parser.add_argument("-c1", "--condition1", required=False, default="resistant",
                        help="Condición 1 de los datos de expresión (default: resistant)")
    parser.add_argument("-c2", "--condition2", required=False, default="susceptible",
                        help="Condición 2 de los datos de expresión (default: susceptible)")
    parser.add_argument("-s", "--strand", type=int, choices=[0, 1, 2], default=0,
                        help="Especificidad de cadena: 0 (unstranded), 1 (forward), 2 (reverse) (default: 0)")
    parser.add_argument("-o", "--outdir", required=False, default="data/counts/", 
                        help="Directorio de salida")
    parser.add_argument("-w", "--workers", type=int, required= False, default= 6, 
                        help="Número de hilos")
    parser.add_argument("-u", "--umbral", type=float, default=70.0,
                        help="Umbral mínimo de lecturas asignadas en %% (default: 70.0)")
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

    raw_counts_path = os.path.join(args.outdir, "raw_counts.txt")
    summary_path = f"{raw_counts_path}.summary"

    try:
        grupos = [args.condition1, args.condition2]

        bam_files = []
        for condition in grupos:
            bam_files += find_bams(data_dir=args.data_dir, grupo=condition)

        run_featureCounts(bam_files=bam_files, gtf_path=args.annotation, output_file=raw_counts_path,
                          workers=args.workers, strand=args.strand)

        assignment_rates = parse_summary(summary_path, umbral=args.umbral)

        clean_counts(raw_counts_path, outdir=os.path.join(args.outdir, "clean"))

        if args.graphic:
            plot_assignment(assignment_rates, outdir="results/images/", umbral=args.umbral)

        logger.info("Módulo de conteo finalizado correctamente")

    except (FileNotFoundError, RuntimeError, ValueError) as e:
        logger.error("El pipeline de conteo se detuvo por un error crítico: %s", e)
        sys.exit(1)

if __name__ == "__main__":
    main()