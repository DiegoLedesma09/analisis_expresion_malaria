"""
qc_trimming.py
==============

Módulo para el control de calidad y la limpieza de lecturas paired-end de
Anopheles gambiae (grupos resistente y susceptible).

Corre FastQC y MultiQC sobre las lecturas crudas (solo diagnóstico), recorta
adaptadores y bases de baja calidad con Trim Galore, vuelve a correr FastQC y
MultiQC sobre las lecturas recortadas y, con ese segundo reporte, decide qué
muestras entran al análisis. Opcionalmente grafica cuántas muestras pasan el
filtro antes y después del trimming.

Requisitos:
    - FastQC y Trim Galore (con Cutadapt) disponibles en el PATH.
    - multiqc, pandas y matplotlib.
    - Lecturas crudas en data/raw/<grupo>/sra_runs/ (por ejemplo, las que
      genera download_data.py), con nombres <SRR>_1.fastq y <SRR>_2.fastq.

Funciones principales:
    run_fastqc          Corre FastQC sobre las lecturas crudas o recortadas.
    run_multiqc         Genera un reporte de MultiQC por grupo.
    read_tsv            Lee multiqc_fastqc.txt de cada grupo en un DataFrame.
    classify_sample     Asigna pass/warn/fail a cada muestra según OVERALL_FLAGS.
    delete_fails        Elimina las muestras clasificadas como fail.
    pair_samples        Conserva solo los SRR con ambas lecturas (_1 y _2).
    run_trim_galore     Recorta los pares con Trim Galore (modo --paired).
    plot_qc_comparison  Grafica las muestras que pasan el filtro antes y después
                        del trimming.
    my_parser           Define los argumentos de línea de comandos.
    setup_logging       Configura el logging del módulo.
    main                Ejecuta el flujo completo.

Archivos generados (rutas relativas al directorio de salida, default results/):
    results/fastqc/<raw|trimmed>/<grupo>/    Reportes de FastQC
    results/qc/<raw|trimmed>/<grupo>/        Reportes de MultiQC
    results/trimmed/<grupo>/                 Lecturas recortadas (<SRR>_val_1.fq y
                                             <SRR>_val_2.fq) y reportes de trimming
    results/images/qc_comparison.png         Gráfica antes/después del trimming

Argumentos:
    -w,  --workers      Número de hilos para FastQC y Trim Galore (default: 8).
    -o,  --outdir       Directorio de salida (default: results).
    -v,  --verbose      Muestra mensajes DEBUG.
    --log_file          Ruta del archivo de log.
    -g,  --graphic      Incluye el gráfico comparativo del trimming (default).
         --no-graphic   Omite el gráfico comparativo del trimming.

Uso:
    python src/qc_trimming.py -w 8 -o results -v
    python src/qc_trimming.py --no-graphic

Nota:
    - Ejecutar desde la raíz del proyecto, ya que las rutas de entrada
      (data/raw/) son relativas.
    - Ejecútese después del módulo download_data.py.
"""

import argparse as ap
import logging
import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import multiqc
import pandas as pd
from multiqc.core.update_config import ClConfig

logger = logging.getLogger(__name__)

# Voy a definir una lista de todos los campos que evalúa FASTQC
FASTQC_FLAG_COLUMNS = [
    "per_base_sequence_quality",
    "per_tile_sequence_quality",
    "per_sequence_quality_scores",
    "per_base_sequence_content",
    "per_sequence_gc_content",
    "per_base_n_content",
    "sequence_length_distribution",
    "sequence_duplication_levels",
    "overrepresented_sequences",
    "adapter_content",
]

# Y los campos que nos interesan para nuestro filtrado de lecturas
OVERALL_FLAGS = [
    "adapter_content",
    "per_base_n_content",
    "per_base_sequence_content",
]

FASTQ_SUFFIX_FWD = "_1.fastq"
FASTQ_SUFFIX_REV = "_2.fastq"

FASTQ_SUFFIX_FWD = "_1.fastq"
FASTQ_SUFFIX_REV = "_2.fastq"

# Parámetros de Trim Galore
TRIM_QUALITY = 20      # Calidad Phred mínima para recortar el extremo 3
TRIM_MIN_LENGTH = 36   # Longitud mínima de lectura tras el recorte

def run_fastqc(w:int= 4, outdir:str = "results/", stage:str = "raw") -> list: 
    """
    run_fastqc corre FastQC sobre los .fastq de cada grupo de expresión del organismo Anopheles gambiae
    y deja los reportes en <outdir>/fastqc/<stage>/<grupo>.

    Args:
        w (int) -> Número de hilos que usará FastQC (por defecto 8)
        outdir (str) -> Ruta del directorio de salida (por defecto results/)
        stage (str) -> "raw" para las lecturas crudas (data/raw/<grupo>/sra_runs) o "trimmed" para las
                       lecturas recortadas por Trim Galore (<outdir>/trimmed/<grupo>, archivos *_val_*.fq)

    Returns:
        fastqc_dirs (list) -> Lista de rutas de salida de FastQC, una por grupo
    """
    fastqc_dirs = []

    for grupo in ["resistant", "susceptible"]:
        if stage == "raw":
            entrada = Path("data/raw") / grupo / "sra_runs"
            fastqs = sorted(entrada.glob("*.fastq*"))
        else:
            entrada = Path(outdir) / "trimmed" / grupo
            fastqs = sorted(entrada.glob("*_val_*.fq"))

        salida = Path(f"{outdir}/fastqc") / stage / grupo
        salida.mkdir(parents=True, exist_ok=True)

        if not fastqs:
            logger.error("[%s] No hay .fastq en %s", stage, entrada)
            raise FileNotFoundError(f"No hay .fastq en {entrada}")

        logger.info("[%s] FastQC sobre %d archivos del grupo '%s' (hilos: %d)", stage, len(fastqs), grupo, w)

        try:
            subprocess.run(
                ["fastqc", "-t", str(w), "--outdir", str(salida), *map(str, fastqs)],
                check=True,
            )
        except FileNotFoundError:
            logger.error("No se encontró el ejecutable 'fastqc'. ¿Está instalado y en el PATH?")
            raise
        except subprocess.CalledProcessError as e:
            logger.error("[%s] FastQC falló en el grupo '%s' (código %d)", stage, grupo, e.returncode)
            raise

        logger.info("[%s] FastQC terminado para '%s', reportes en %s", stage, grupo, salida)
        fastqc_dirs.append(salida)

    return fastqc_dirs

def run_multiqc(fastqc_dirs:list, outdir:str = "results/", stage:str = "raw") -> list:
    """
    run_multiqc se encargará de recorrer los directorios con los reportes de FastQC de cada grupo
    (resistente y susceptible) y generará un reporte de MultiQC por grupo, que se guardará en la
    carpeta results/qc.

    Args:
        fastqc_dirs (list) -> Lista de rutas de salida generadas por run_fastqc, una por grupo
        outdir (str) -> Ruta al directorio de salida (por defecto results/)

    Returns:
         output_dirs (list) -> Lista de rutas de salida generadas por MultiQC, una por grupo

    """
    output_dirs = []

    for i in fastqc_dirs:
        output_dir = f"{outdir}/qc/{stage}/{i.name}"
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        logger.info("[%s] MultiQC sobre %s", stage, i)
        
        try:
            multiqc.run(str(i), cfg=ClConfig(output_dir=output_dir), return_html=True)
        except Exception:
            logger.exception("[%s] MultiQC falló en %s", stage, i)
            raise

        logger.info("[%s] Reporte MultiQC en %s", stage, output_dir)
        output_dirs.append(output_dir)

    return output_dirs

def read_tsv(output_dirs:list) -> pd.DataFrame:
    """ 
    read_tsv se encarga de leer el tsv de FastQC, su única responsabilidad es parsear el archivo y devolver un DataFrame por
    muestra, tal y como vienen de MultiQC
    
    Args: 
         output_dirs (list) -> Lista de rutas de salida generadas por run_multiqc, una por grupo
    
    Returns:
         df_qc (DataFrame) -> DataFrame por muestra de MultiQC, parseado y por el orden del MultiQC
    """

    dfs = []

    for output_dir in output_dirs:
        tsv_path = Path(output_dir) /"multiqc_data"/"multiqc_fastqc.txt"

        try:
            df_grupo = pd.read_csv(tsv_path, sep="\t")
        except FileNotFoundError:
            logger.warning("No existe %s, se omite este grupo", tsv_path)
            continue
        except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
            logger.warning("No se pudo parsear %s (%s), se omite este grupo", tsv_path, e)
            continue

        df_grupo["grupo"] = Path(output_dir).name
        logger.info("Leídas %d filas de %s", len(df_grupo), tsv_path)
        dfs.append(df_grupo)

    if not dfs:
        logger.error("No se pudo leer ningún multiqc_fastqc.txt")
        raise ValueError("No se pudo leer ningún multiqc_fastqc.txt")

    df_qc = pd.concat(dfs, ignore_index=True)
    logger.info("DataFrame de QC con %d muestras en total", len(df_qc))

    return df_qc

def classify_sample(df_qc:pd.DataFrame) -> dict:
    """
    classify_sample se encarga de construir una estructura de diccionarios, en por cada muestra se asocia la calidad, dependendiendo
    principalmente de la presencia de adaptadores en la muestra, presencia de n y contenido de bases. De ésta manera se clasifican 
    las muestras sobre las que realmente se van a usar en el análisis.
    
    Args: 
        df_qc (DataFrame) -> Estructura de datos creada por la función read_tsv(). DataFrame de la información de las lecturas por el 
                        MultiQC
    
    Returns:
        dict_qc (dict) -> Diccionario en el que se le asocia la calidad de la muestra con el identificador de la lectura
    """
    dict_qc = {}

    for i, row in df_qc.iterrows():
        try:
            sample = row["Sample"]

            banderas = {col: row[col] for col in FASTQC_FLAG_COLUMNS}
        except KeyError as e:
            logger.warning("Fila %s sin la columna %s, se omite la muestra", i, e)
            continue

        estados = [banderas[col] for col in OVERALL_FLAGS]
        if "fail" in estados:
            banderas["overall"] = "fail"
        elif "warn" in estados:
            banderas["overall"] = "warn"
        else:
            banderas["overall"] = "pass"

        banderas["grupo"] = row["grupo"]

        logger.debug("Muestra %s (%s): overall=%s", sample, banderas["grupo"], banderas["overall"])

        dict_qc[sample] = banderas

    logger.info("Clasificadas %d muestras", len(dict_qc))

    return dict_qc

def delete_fails(dict_qc: dict) -> dict:
    """ 
    delete_fails detecta los SRRs que son clasificados como fails, los elimina de los archivos que se usarán en el análisis.
    
    Args: 
        dict_qc (dict) ->  Diccionario creado por la función classify_sample

    Returns:
        dict_qc_filtered (dict) -> Diccionario que asocia la calidad de la muestra con el identificador de la lectura clasificada
                                          por la calidad de la lectura eliminando los "fails"
    """
    dict_qc_filtered = {
        sample: banderas
        for sample, banderas in dict_qc.items()
        if banderas["overall"] != "fail"
    }

    descartadas = sorted(set(dict_qc) - set(dict_qc_filtered))
    if descartadas:
        logger.warning("Se descartan %d muestras con 'fail': %s", len(descartadas), ", ".join(descartadas))
    else:
        logger.info("Ninguna muestra clasificada como 'fail'")

    if not dict_qc_filtered:
        logger.warning("Ninguna muestra pasó el filtro. Revisa OVERALL_FLAGS")

    return dict_qc_filtered

def pair_samples(dict_qc_filtered:dict) -> dict:
    """
    pair_samples agrupa las lecturas _1 y _2 por SRR y se queda solo con los pares completos,
    es decir, aquellos en los que ambas lecturas pasaron el control de calidad. Funciona con los
    nombres de las lecturas crudas (SRR_1) y de las recortadas por Trim Galore (SRR_val_1).

    Args:
        dict_qc_filtered (dict) -> Diccionario creado por delete_fails(), con "grupo" incluido por muestra

    Returns:
        pares (dict) -> {SRR: {"grupo": grupo}}, el formato que espera run_trim_galore
    """
    mates = {}

    for sample, info in dict_qc_filtered.items():
        sample = sample.replace("_val_", "_")
        srr, _, mate = sample.rpartition("_")
        entrada = mates.setdefault(srr, {"grupo": info["grupo"], "mates": set()})
        entrada["mates"].add(mate)

    pares = {
        srr: {"grupo": d["grupo"]}
        for srr, d in mates.items()
        if d["mates"] == {"1", "2"}
    }

    incompletos = sorted(set(mates) - set(pares))
    if incompletos:
        logger.warning("Se descartan %d SRR con par incompleto: %s", len(incompletos), ", ".join(incompletos))

    logger.info("Pares completos: %d de %d SRR", len(pares), len(mates))

    return pares


def run_trim_galore(pares:dict, raw_dir:str = "data/raw", output_dir:str = "results/trimmed", w:int = 4) -> dict:
    """
    run_trim_galore se encarga de correr Trim Galore (modo --paired) sobre cada par de lecturas, usando el
    campo "grupo" para localizar los fastq crudos. Trim Galore detecta el adaptador automáticamente y recorta
    las bases de baja calidad. Por cada SRR se generan <SRR>_val_1.fq y <SRR>_val_2.fq, que contienen solo las
    lecturas que siguen teniendo pareja, además de los reportes de trimming.

    Args:
        pares (dict) -> Diccionario creado por pair_samples(), {SRR: {"grupo": grupo}}
        raw_dir (str) -> Ruta al directorio donde están los fastq crudos, organizados por grupo/sra_runs (por defecto data/raw)
        output_dir (str) -> Ruta al directorio donde Trim Galore escribirá las lecturas ya limpias (por defecto results/trimmed)
        w (int) -> Número de cores que usará Trim Galore (por defecto 8)

    Returns:
        output_paths (dict) -> Diccionario que asocia cada SRR con la carpeta de salida de su trimming
    """
    output_paths = {}
    fallidos = []

    for n, (srr, info) in enumerate(pares.items(), start=1):
        grupo = info["grupo"]

        fwd = f"{raw_dir}/{grupo}/sra_runs/{srr}{FASTQ_SUFFIX_FWD}"
        rev = f"{raw_dir}/{grupo}/sra_runs/{srr}{FASTQ_SUFFIX_REV}"

        if not (Path(fwd).exists() and Path(rev).exists()):
            logger.warning("[%d/%d] %s: no se encontraron %s y/o %s, se omite", n, len(pares), srr, fwd, rev)
            fallidos.append(srr)
            continue

        sample_output_dir = Path(output_dir) / grupo
        sample_output_dir.mkdir(parents=True, exist_ok=True)

        logger.info("[%d/%d] Trim Galore sobre %s (grupo '%s')", n, len(pares), srr, grupo)

        try:
            subprocess.run([
                "trim_galore", "--paired",
                "--cores", str(w),
                "--phred33",
                "--quality", str(TRIM_QUALITY),
                "--length", str(TRIM_MIN_LENGTH),
                "--basename", srr,
                "--output_dir", str(sample_output_dir),
                fwd, rev,
            ], check=True)

        except FileNotFoundError:
            logger.error("No se encontró el ejecutable 'trim_galore'. ¿Está instalado y en el PATH?")
            raise
        except subprocess.CalledProcessError as e:
            logger.warning("[%d/%d] Trim Galore falló en %s (código %d), se omite", n, len(pares), srr, e.returncode)
            fallidos.append(srr)
            continue

        output_paths[srr] = str(sample_output_dir)

    if fallidos:
        logger.warning("Trimming fallido u omitido en %d SRR: %s", len(fallidos), ", ".join(fallidos))

    if not output_paths:
        logger.error("Trim Galore no produjo ninguna salida")
        raise RuntimeError("Trim Galore no produjo ninguna salida")

    logger.info("Trimming completado en %d de %d SRR", len(output_paths), len(pares))

    return output_paths

def plot_qc_comparison(dict_qc_raw: dict, dict_qc_trim: dict, outdir: str = "results/") -> Path:
    """
    plot_qc_comparison grafica, por grupo (resistente y susceptible), cuántas muestras pasan el filtro
    (overall distinto de "fail") antes y después del trimming, para evaluar el efecto de Trimmomatic.

    Args:
        dict_qc_raw (dict) -> Diccionario de classify_sample() sobre las lecturas crudas
        dict_qc_trim (dict) -> Diccionario de classify_sample() sobre las lecturas recortadas
        outdir (str) -> Ruta al directorio de salida (por defecto results/)

    Returns:
        output_path (Path) -> Ruta de la gráfica generada
    """
    grupos = sorted({b["grupo"] for b in (*dict_qc_raw.values(), *dict_qc_trim.values())})

    def contar(dict_qc: dict, grupo: str) -> tuple:
        total = sum(1 for b in dict_qc.values() if b["grupo"] == grupo)
        pasan = sum(1 for b in dict_qc.values() if b["grupo"] == grupo and b["overall"] != "fail")
        return pasan, total

    conteos_raw = [contar(dict_qc_raw, g) for g in grupos]
    conteos_trim = [contar(dict_qc_trim, g) for g in grupos]

    ancho = 0.35
    posiciones = range(len(grupos))

    fig, ax = plt.subplots(figsize=(7, 5))
    for offset, conteos, etiqueta in [
        (-ancho / 2, conteos_raw, "Antes del trimming"),
        (ancho / 2, conteos_trim, "Después del trimming"),
    ]:
        barras = ax.bar([p + offset for p in posiciones], [c[0] for c in conteos], ancho, label=etiqueta)
        for barra, (pasan, total) in zip(barras, conteos):
            ax.text(barra.get_x() + barra.get_width() / 2, barra.get_height(),
                    f"{pasan}/{total}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(list(posiciones))
    ax.set_xticklabels(grupos)
    ax.set_ylabel("Muestras que pasan el filtro")
    ax.set_title("Muestras que pasan el filtro de calidad antes y después del trimming")
    ax.legend()

    output_path = Path(outdir) / "images" / "qc_comparison.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)

    logger.info("Gráfica de comparación guardada en %s", output_path)

    return output_path

def my_parser() -> ap.Namespace:
    """
    Retiene los argumentos deseados para el módulo desde línea de comandos.

    Args:
        None
    
    Returns:
        ap.Namespace: Argumentos indicados por el usuario
    """    
    parser = ap.ArgumentParser(
    description="Pipeline para la limpieza de las lecturas mediante el uso de FastQC, MultiQC y Trim Galore")
    parser.add_argument("-w", "--workers", default=4, type=int, help="No. de workers para trabajar")
    parser.add_argument("-o", "--outdir", default="results/", help="Ruta al directorio de salida")
    parser.add_argument("-v", "--verbose", action="store_true", help="Muestra mensajes DEBUG")
    parser.add_argument("--log_file", default="logs/qc_trimming.log", help="Archivo de log")
    parser.add_argument("-g", "--graphic", action=ap.BooleanOptionalAction, default=True,
                        help="Incluye el grafico comparativo del uso de trimming (usa --no-graphic para omitirlo)")
    return parser.parse_args()

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


def main():
    """
    Ejecuta el flujo completo de control de calidad y trimming.

    Pasos, en orden:
        1. Corre FastQC y MultiQC sobre las lecturas crudas y clasifica las
           muestras (solo diagnóstico, no se descarta nada).
        2. Empareja las lecturas _1 y _2 por SRR y las recorta con Trim Galore.
        3. Corre FastQC y MultiQC sobre las lecturas recortadas y las clasifica.
        4. Si --graphic está activo, grafica cuántas muestras pasan el filtro
           antes y después del trimming.
        5. Elimina las muestras con fail y conserva solo los pares completos,
           que son los que entran al análisis.

    Los argumentos se leen de la línea de comandos mediante my_parser().

    Returns:
        pares_finales (dict) -> {SRR: {"grupo": grupo}} de los pares que
                                pasan el filtro final
    """
    args = my_parser()
    setup_logging(args.verbose, args.log_file)

    logger.info("Inicio del pipeline (workers=%d, outdir=%s)", args.workers, args.outdir)

    try:
        fastqc_dirs = run_fastqc(w=args.workers, outdir=args.outdir, stage="raw")
        output_dirs = run_multiqc(fastqc_dirs, outdir=args.outdir, stage="raw")
        df_qc = read_tsv(output_dirs)
        dict_qc = classify_sample(df_qc)

        logger.info("Muestras crudas totales: %d", len(dict_qc))
        logger.info("Muestras crudas que pasarían el filtro (solo informativo): %d", len(delete_fails(dict_qc)))

        pares = pair_samples(dict_qc)
        logger.info("Pares completos a recortar: %d", len(pares))

        run_trim_galore(pares, output_dir=f"{args.outdir}/trimmed", w=args.workers)

        fastqc_dirs_trim = run_fastqc(w=args.workers, outdir=args.outdir, stage="trimmed")
        output_dirs_trim = run_multiqc(fastqc_dirs_trim, outdir=args.outdir, stage="trimmed")
        df_qc_trim = read_tsv(output_dirs_trim)
        dict_qc_trim = classify_sample(df_qc_trim)

        if args.graphic:
            plot_qc_comparison(dict_qc, dict_qc_trim, outdir=args.outdir)

        dict_qc_filtered = delete_fails(dict_qc_trim)
        pares_finales = pair_samples(dict_qc_filtered)

        logger.info("Muestras recortadas totales: %d", len(dict_qc_trim))
        logger.info("Muestras recortadas que pasan el filtro: %d", len(dict_qc_filtered))
        logger.info("Pares completos que pasan el filtro final: %d", len(pares_finales))

        if not pares_finales:
            logger.warning("Ningún par pasó el filtro final. Revisa OVERALL_FLAGS o los nombres en 'Sample'")

    except Exception as e:
        logger.error("%s", e)
        logger.debug("Detalle del error", exc_info=True)
        sys.exit(1)

    logger.info("Pipeline terminado")
    return pares_finales

if __name__ == "__main__":
    main()