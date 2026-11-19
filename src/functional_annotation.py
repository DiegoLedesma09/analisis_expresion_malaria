"""
functional_annotation.py
================

Módulo para el análisis de enriquecimiento funcional de Gene Ontology (GO) a
partir de los resultados de expresión diferencial entre muestras resistentes y
susceptibles de Anopheles gambiae.

Utiliza la prueba exacta de Fisher (vía goatools) para identificar términos GO
sobrerrepresentados en las listas de genes sobreexpresados y subexpresados,
aplicando corrección por pruebas múltiples (Benjamini-Hochberg).
Genera tablas con los resultados y gráficos de puntos (dotplots) y barras 
(barplots) de los términos más significativos.

Requisitos:
    - Python 3.12
    - Librerías: pandas, numpy, matplotlib, goatools
    - Archivos de datos: 
        * Ontología en data/go/go-basic.obo
        * Anotaciones gene2go en data/go/gene2go (descargado de NCBI)
    - Tabla de resultados de differential_expression.py.

Funciones principales:
    load_results        Busca y lee la tabla de resultados de expresión diferencial.
    select_genes        Extrae las listas de genes sobreexpresados y subexpresados.
    load_go_resources   Carga la ontología (.obo) y la anotación gene -> GO.
    build_background    Define el conjunto de fondo (universo) con los genes expresados.
    run_go_enrichment   Ejecuta la prueba exacta de Fisher con corrección múltiple. 
    plot_enrichment     Genera la gráfica opcional de los términos más enriquecidos.
    setup_logging       Configura el logging del módulo.
    my_parser           Define y procesa los argumentos de línea de comandos.
    main                Orquesta el flujo completo de enriquecimiento GO.

    Archivos generados:
    results/go/go_enrichment_results.tsv          Tabla con todos los términos GO evaluados.
    results/images/go_barplot_<grupo>.png             Gráfico de barras por grupo (up/down).
    results/images/go_dotplot_<grupo>.png             Gráfico de puntos por grupo (up/down).
    logs/go_enrichment.log                        Archivo de registro.

Argumentos:
    -w,   --workdir       Directorio raíz del proyecto (default: .).
    -i,   --input         Ruta a la tabla de expresión diferencial.
    --obo                 Ruta al archivo go-basic.obo.
    --annotation          Ruta al archivo gene2go.
    -o,   --outdir        Directorio de salida (default: results/go).
    -t,   --taxid         TaxID de NCBI para A. gambiae (default: 7165).
    --pvalue              Umbral de significancia ajustado (default: 0.05).
    --method              Corrección múltiple: 'fdr_bh' o 'bonferroni' (default: fdr_bh).
    --aspectos            Aspectos de GO a evaluar (BP, MF, CC).
"""
import argparse as ap
import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.obo_parser import GODag

logger = logging.getLogger(__name__)

def load_results(ruta_resultados: str) -> pd.DataFrame:
    """
    Busca y lee la tabla de resultados de expresión diferencial (separada por comas)

    Args:
        ruta_resultados (str) -> Ruta hacia el archivo de la tabla de resultados realizada
            por el módulo de diferential_expression.py

    Returns:
        df (pd.DataFrame) -> Dataframe representativo de los resultados
    """
    if not os.path.exists(ruta_resultados):
        logger.error("No se encontró el archivo de resultados %s", ruta_resultados)
        raise FileNotFoundError(f"No se encontró el archivo de resultados: {ruta_resultados}")
    df_res = pd.read_csv(ruta_resultados, sep=",", index_col=0)

    # gene2go usa el GeneID de NCBI como entero: LOC1276231 -> 1276231
    es_loc = df_res.index.astype(str).str.fullmatch(r"LOC\d+")
    n_sin_loc = int((~es_loc).sum())
    if n_sin_loc:
        ejemplos = df_res.index[~es_loc][:5].tolist()
        logger.warning("%d genes no tienen formato LOC<número> y no se pueden mapear a GeneID "
                       "(ejemplos: %s); quedarán sin anotación GO.", n_sin_loc, ejemplos)
    df_res.index = [int(str(g)[3:]) if loc else g for g, loc in zip(df_res.index, es_loc)]
    return df_res


def select_genes(df_res: pd.DataFrame, condition1: str = "resistant", condition2: str = "susceptible"):
    """
    Extrae las listas de genes sobreexpresados en cada condición a partir de la
    columna 'categoria' generada por diferential_expression.py.

    Args:
        df_res (pd.DataFrame) -> Dataframe representativo de los resultados
        condition1 (str) -> Condición cuya categoría es 'up_{condition1}'
        condition2 (str) -> Condición cuya categoría es 'up_{condition2}'

    Returns:
        genes_up1 (list) -> Genes sobreexpresados en condition1
        genes_up2 (list) -> Genes sobreexpresados en condition2
    """
    if "categoria" not in df_res.columns:
        raise ValueError("La tabla no tiene la columna 'categoria'")

    if "padj" not in df_res.columns or df_res["padj"].isna().all():
        logger.warning("Tabla sin padj: los genes se seleccionaron solo por corte de LFC; "
                       "el enriquecimiento GO es exploratorio.")

    cat1, cat2 = f"up_{condition1}", f"up_{condition2}"
    categorias = set(df_res["categoria"].dropna().unique())
    for cat in (cat1, cat2):
        if cat not in categorias:
            logger.warning("La categoría '%s' no aparece en la tabla (categorías: %s)",
                           cat, sorted(categorias))

    genes_up1 = df_res.index[df_res["categoria"] == cat1].tolist()
    genes_up2 = df_res.index[df_res["categoria"] == cat2].tolist()
    return genes_up1, genes_up2

def load_go_resources(obo_path: str, gene2go_path: str, taxid:int):
    """
    Carga la ontología (.obo) y la anotación gene -> GO.

    Args:
        obo_path (str) -> Ruta al archivo descargado en el repositorio .obo
        gene2go_path (str) -> Ruta al archivo gene2go en el repositorio
        taxid (int) -> Identificador de organismo según su taxid de NCBI

    Returns:
        godag (str) -> Archivo .obo de Gene Ontology
        gene2gos (str) -> Archivo gene2gos de Gene Ontology

    """
    logger.info("Cargando ontología desde %s", obo_path)
    godag = GODag(obo_path, prt=sys.stdout)
    
    logger.info("Cargando anotaciones desde %s", gene2go_path)
    # Gene2GoReader para leer anotaciones formato gene2go (taxid 7165 para A. gambiae)
    anno_reader = Gene2GoReader(gene2go_path, taxid=taxid)
    gene2gos = anno_reader.get_ns2assc()
    return godag, gene2gos

def build_background(df_res: pd.DataFrame) -> set:
    """
    Define el conjunto de fondo (universo) con los genes expresados evaluados.

    Args:
        df_res (pd.DataFrame) -> Dataframe representativo de los resultados

    Returns:
        background (set) -> Set del universo sobre el que va a funcionar GO
    """
    return set(df_res.index.tolist())

def run_go_enrichment(genes_estudio: list[str] | set[str],
    background: set,
    godag: GODag,
    gene2gos: dict[str, dict[str, set]],
    nombre_lista: str = "estudio",
    alpha: float = 0.05,
    method: str = "fdr_bh",
    aspectos: list[str] | None = None,) -> pd.DataFrame:
    """
    Ejecuta el análisis de sobrerrepresentación (prueba exacta de Fisher con
    corrección por pruebas múltiples) en los aspectos de Gene Ontology acordados
    y devuelve los resultados enriquecidos en un DataFrame de pandas.

    Args:
        genes_estudio (list[str] | set[str]) -> Lista o conjunto de genes
            diferencialmente expresados (ej. genes_up1 o genes_up2).
        background (set) -> Conjunto de fondo (universo) devuelto por
            build_background.
        godag (GODag) -> Ontología GO cargada desde go-basic.obo.
        gene2gos (dict[str, dict[str, set]]) -> Diccionario por aspecto
            ('BP', 'MF', 'CC') devuelto por anno_reader.get_ns2gene2gos().
        nombre_lista (str) -> Etiqueta del grupo analizado (ej. 'up_resistant' o
            'up_susceptible') para registro en log y columna de salida.
        alpha (float) -> Umbral de significancia para p ajustado (default: 0.05).
        method (str) -> Método de corrección por pruebas múltiples: 'fdr_bh'
            (Benjamini-Hochberg) o 'bonferroni' (default: 'fdr_bh').
        aspectos (list[str] | None) -> Aspectos de GO a conservar ('BP', 'MF', 'CC').
            Si es None, evalúa los tres aspectos.

    Returns:
        df_go (pd.DataFrame) -> Tabla de términos sobrerrepresentados ordenada por
            p ajustado, con columnas: grupo, go_id, aspecto, termino, k_estudio,
            n_estudio, K_universo, N_universo, gene_ratio, p_value, padj,
            significativo, profundidad y genes.
    """

    aspectos_validos = set(aspectos) if aspectos else {"BP", "MF", "CC"}
    ns2assoc_filtrado = {
        ns: assoc for ns, assoc in gene2gos.items() if ns in aspectos_validos
    }

    genes_anotados_totales = set().union(*ns2assoc_filtrado.values()) if ns2assoc_filtrado else set()
    universo_anotado = set(background) & genes_anotados_totales
    estudio_set = set(genes_estudio)
    estudio_valido = estudio_set & universo_anotado
    sin_anotacion = estudio_set - universo_anotado

    logger.info(
        "[%s] Genes de entrada: %d | Con anotación en universo: %d | Sin anotación/equivalencia: %d",
        nombre_lista,
        len(estudio_set),
        len(estudio_valido),
        len(sin_anotacion),
    )

    if len(estudio_valido) < 3:
        logger.warning(
            "[%s] Lista vacía o demasiado corta (%d gen(es) anotados); se omite el análisis GO.",
            nombre_lista,
            len(estudio_valido),
        )
        return pd.DataFrame()

    goea = GOEnrichmentStudyNS(
        pop=universo_anotado,
        ns2assoc=ns2assoc_filtrado,
        godag=godag,
        propagate_counts=True,
        alpha=alpha,
        methods=[method],
        prt=None,
    )

    resultados = goea.run_study(estudio_valido, prt=None)
    col_padj = f"p_{method}"

    registros = []
    for r in resultados:
        if r.enrichment != "e":
            continue

        k, n = r.ratio_in_study
        K, N = r.ratio_in_pop
        p_ajustado = getattr(r, col_padj, r.p_uncorrected)

        registros.append(
            {
                "grupo": nombre_lista,
                "go_id": r.GO,
                "aspecto": r.NS,
                "termino": r.name,
                "k_estudio": k,
                "n_estudio": n,
                "K_universo": K,
                "N_universo": N,
                "gene_ratio": round(k / n, 6) if n > 0 else 0.0,
                "p_value": r.p_uncorrected,
                "padj": p_ajustado,
                "significativo": bool(p_ajustado < alpha),
                "profundidad": r.depth,
                "genes": ";".join(sorted(str(g) for g in r.study_items)),
            }
        )

    if not registros:
        logger.warning("[%s] No se encontraron términos GO sobrerrepresentados.", nombre_lista)
        return pd.DataFrame()

    df_go = pd.DataFrame(registros).sort_values(["padj", "p_value"], ascending=[True, True])
    logger.info(
        "[%s] Términos sobrerrepresentados: %d | Significativos (padj < %.3f): %d",
        nombre_lista,
        len(df_go),
        alpha,
        int(df_go["significativo"].sum()),
    )
    return df_go

def filter_and_sort(df_go: pd.DataFrame, grupo: str, top_n: int = 10) -> pd.DataFrame:
    """
    Filtra los términos significativos (si existen) y devuelve los principales
    ordenados de menor a mayor padj, listos para graficar de arriba hacia abajo.

    Args:
        df_go (pd.DataFrame) -> DataFrame con resultados de enriquecimiento GO.
        grupo (str) -> Nombre del grupo a filtrar.
        top_n (int) -> Número máximo de términos a conservar (default: 10).

    Returns:
        df_top (pd.DataFrame) -> DataFrame filtrado y ordenado para el eje Y.
    """
    df_g = df_go[df_go["grupo"] == grupo].copy()
    if df_g.empty:
        logger.debug("Se recibió un DataFrame vacío para filtrar y ordenar.")
        return df_g

    if "significativo" in df_g.columns and df_g["significativo"].any():
        df_g = df_g[df_g["significativo"]]
        logger.debug("Se filtraron %d términos significativos para graficar.", len(df_g))
    else:
        logger.debug("No hay términos significativos o falta la columna.")

    return df_g.sort_values("padj", ascending=True).head(top_n).iloc[::-1]

def barplot(df_top: pd.DataFrame, grupo: str, outdir: str) -> None:
    """
    Genera y guarda un barplot horizontal de enriquecimiento GO (X: k_estudio,
    Y: termino, color: padj).

    Args:
        df_top (pd.DataFrame) -> Tabla con los términos GO seleccionados.
        grupo (str) -> Nombre del grupo a filtrar.
        outdir (str) -> Directorio donde se guardará la figura.

    Returns:
        None
    """
    logger.debug("Generando barplot '%s' con %d términos en %s", grupo, len(df_top), outdir)
    fig, ax = plt.subplots(figsize=(8, 6))
    y = np.arange(len(df_top))
    padj_seguro = np.clip(df_top["padj"].astype(float), 1e-300, 1.0)

    norm = plt.Normalize(vmin=padj_seguro.min(), vmax=padj_seguro.max())
    cmap = plt.get_cmap("RdYlBu")

    ax.barh(
        y,
        df_top["k_estudio"].astype(float),
        color=cmap(norm(padj_seguro)),
        edgecolor="black",
        linewidth=0.5,
    )

    etiquetas = [f"{row['termino']} ({row['aspecto']})" for _, row in df_top.iterrows()]
    ax.set_yticks(y)
    ax.set_yticklabels(etiquetas, fontsize=6)
    ax.set_xlabel("Count (k_estudio)")
    ax.set_title(f"GO Barplot del {grupo}")
    ax.grid(True, axis="x", linestyle="--", alpha=0.4)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("padj")

    fig.tight_layout()
    ruta = os.path.join(outdir, f"go_barplot_{grupo}.png")
    fig.savefig(ruta, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("Barplot guardado en %s", ruta)

def dotplot(df_top: pd.DataFrame, grupo: str, outdir: str) -> None:
    """
    Genera y guarda un dotplot de enriquecimiento GO (X: gene_ratio, Y: termino,
    tamaño: k_estudio, color: padj).

    Args:
        df_top (pd.DataFrame) -> Tabla con los términos GO seleccionados.
        grupo (str) -> Nombre del grupo a filtrar.
        outdir (str) -> Directorio donde se guardará la figura.

    Returns:
        None
    """
    logger.debug("Generando dotplot '%s' con %d términos en %s", grupo, len(df_top), outdir)
    fig, ax = plt.subplots(figsize=(8, 6))
    y = np.arange(len(df_top))
    padj_seguro = np.clip(df_top["padj"].astype(float), 1e-300, 1.0)

    sc = ax.scatter(
        df_top["gene_ratio"].astype(float),
        y,
        s=df_top["k_estudio"].astype(float) * 25,
        c=padj_seguro,
        cmap="GnBu",
        edgecolors="black",
        linewidths=0.6,
    )

    etiquetas = [f"{row['termino']} ({row['aspecto']})" for _, row in df_top.iterrows()]
    ax.set_yticks(y)
    ax.set_yticklabels(etiquetas, fontsize=6)
    ax.set_xlabel("GeneRatio")
    ax.set_title(f"GO Dotplot del {grupo}")
    ax.grid(True, linestyle="--", alpha=0.4)

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("padj")

    fig.tight_layout()
    ruta = os.path.join(outdir, f"go_dotplot_{grupo}.png")
    fig.savefig(ruta, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("Dotplot guardado en %s", ruta)

def plot_enrichment(df_go: pd.DataFrame, outdir: str, top_n: int = 10) -> None:
    """
    Genera dotplots y barplots por grupo experimental.

    Args:
        df_go (pd.DataFrame): Tabla combinada devuelta por el análisis de enriquecimiento GO.
        outdir (str): Directorio de salida donde se guardarán las figuras.
        top_n (int): Número máximo de términos a mostrar por gráfica (default: 10).
    """
    if df_go.empty:
        logger.warning("No hay resultados de GO para graficar.")
        return

    os.makedirs(outdir, exist_ok=True)
    logger.info("Iniciando generación de gráficas GO (top_n=%d)", top_n)
    
    for grupo in df_go["grupo"].dropna().unique():
        df_top = filter_and_sort(df_go, grupo=grupo, top_n=top_n)
        if df_top.empty:
            logger.warning("[%s] Sin términos para graficar.", grupo)
            continue
        dotplot(df_top, grupo=grupo, outdir=outdir)
        barplot(df_top, grupo=grupo, outdir=outdir)

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
    parser = ap.ArgumentParser(description="Pipeline para la creación de la anotación funcional con GO")
    parser.add_argument("-w", "--workdir", required=False, default=".",
                        help="Directorio raíz del proyecto (default: .)")
    parser.add_argument("-i", "--input", required=False, default="results/exp_dif/resultados_completos.csv",
                        help="Ruta a la tabla de expresión diferencial")
    parser.add_argument("--obo", required=False, default="data/go/go-basic.obo",
                        help="Ruta al archivo go-basic.obo")
    parser.add_argument("--annotation", required=False, default="data/go/gene2go",
                        help="Ruta al archivo gene2go")
    parser.add_argument("-o", "--outdir", required=False, default="results/go",
                        help="Directorio de salida (default: results/go)")
    parser.add_argument("-t", "--taxid", type=int, required=False, default=7165)
    parser.add_argument("--condition1", default="resistant",
                        help="Primera condición; se busca la categoría 'up_<condition1>' (default: resistant)")
    parser.add_argument("--condition2", default="susceptible",
                        help="Segunda condición; se busca la categoría 'up_<condition2>' (default: susceptible)")
    parser.add_argument("--pvalue", type=float, required=False, default=0.05,
                        help="Umbral de significancia para p ajustado (default: 0.05)")
    parser.add_argument("--method", choices=["fdr_bh", "bonferroni"], default="fdr_bh",
                        help="Corrección por pruebas múltiples (default: fdr_bh)")
    parser.add_argument("--aspectos", nargs="+", choices=["BP", "MF", "CC"], default=["BP", "MF", "CC"],
                        help="Aspectos de GO a evaluar (default: BP MF CC)")
    parser.add_argument("-v", "--verbose", action="store_true", 
                        help="Muestra mensajes DEBUG en consola")
    parser.add_argument("-g", "--graphic", action=ap.BooleanOptionalAction, default=True,
                    help="Incluye las gráficas de enriquecimiento GO")
    parser.add_argument("--log_file", default="logs/go_enrichment.log", 
                        help="Ruta para guardar el archivo de log")
    parser.add_argument( "--top-n", type=int, required=False, default=10,
        help="Número máximo de términos GO a mostrar en cada gráfica (default: 10)")
    
    return parser.parse_args()

def main() -> None:
    """
    Orquesta el flujo completo de enriquecimiento funcional con Gene Ontology.

    Pasos, en orden:
        1. Parsea los argumentos de línea de comandos y configura el registro de logs.
        2. Resuelve las rutas de entrada y salida (relativas a workdir o absolutas)
           y carga los resultados de expresión diferencial.
        3. Selecciona los genes regulados al alza por condición, carga los recursos
           GO (.obo y anotación) y construye el universo de genes (background).
        4. Ejecuta el enriquecimiento GO para cada condición y consolida los
           resultados no vacíos en un único DataFrame.
        5. Guarda la tabla consolidada en go_enrichment_results.tsv dentro de outdir
           (con el esquema de columnas definido si no hubo términos enriquecidos).
        6. Si --graphic está activo, genera las gráficas de enriquecimiento; en caso
           contrario, registra que las gráficas están desactivadas.

    Los argumentos se leen de la línea de comandos mediante my_parser().

    Returns:
        None
    """
    args = my_parser()
    log_path = args.log_file if os.path.isabs(args.log_file) else os.path.join(args.workdir, args.log_file)
    setup_logging(args.verbose, log_path)
    logger.info("Inicio del módulo de enriquecimiento de Gene Ontology")

    try:
        input_path = args.input if os.path.isabs(args.input) else os.path.join(args.workdir, args.input)
        obo_path = args.obo if os.path.isabs(args.obo) else os.path.join(args.workdir, args.obo)
        anno_path = args.annotation if os.path.isabs(args.annotation) else os.path.join(args.workdir, args.annotation)
        outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(args.workdir, args.outdir)

        df_res = load_results(input_path)
        genes_up1, genes_up2 = select_genes(df_res, args.condition1, args.condition2)
        godag, gene2gos = load_go_resources(obo_path, anno_path, taxid=args.taxid)
        background = build_background(df_res)

        tablas = []
        for nombre, lista in [(f"up_{args.condition1}", genes_up1),
                              (f"up_{args.condition2}", genes_up2)]:
            df_go = run_go_enrichment(
                genes_estudio=lista,
                background=background,
                godag=godag,
                gene2gos=gene2gos,
                nombre_lista=nombre,
                alpha=args.pvalue,
                method=args.method,
                aspectos=args.aspectos,
            )
            if not df_go.empty:
                tablas.append(df_go)

        os.makedirs(outdir, exist_ok=True)
        ruta_salida = os.path.join(outdir, "go_enrichment_results.tsv")

        if tablas:
            df_final = pd.concat(tablas, ignore_index=True)
        else:
            df_final = pd.DataFrame(
                columns=[
                    "grupo", "go_id", "aspecto", "termino", "k_estudio", "n_estudio",
                    "K_universo", "N_universo", "gene_ratio", "p_value", "padj",
                    "significativo", "profundidad", "genes",
                ]
            )

        df_final.to_csv(ruta_salida, sep="\t", index=False)
        logger.info("Resultados guardados en %s (%d filas)", ruta_salida, len(df_final))

        if args.graphic:
            logger.info("Gráficas activadas (--graphic).")
            plot_enrichment(
                df_go=df_final,
                outdir=outdir,
                top_n=args.top_n
            )
        else:
            logger.info("Gráficas desactivadas (--no-graphic).")

    except (FileNotFoundError, ValueError, RuntimeError) as e:
        logger.error("Error crítico en el análisis GO: %s", e)
        sys.exit(1)

if __name__ == "__main__":
    main()
