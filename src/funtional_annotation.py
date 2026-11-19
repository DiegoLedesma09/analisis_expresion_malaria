"""
funtional_annotation.py
================

Módulo para el análisis de enriquecimiento funcional de Gene Ontology (GO) a
partir de los resultados de expresión diferencial entre muestras resistentes y
susceptibles de Anopheles gambiae.

Funciones principales:
    load_results        Busca y lee la tabla de resultados de expresión diferencial.
    select_genes        Extrae las listas de genes sobreexpresados y subexpresados.
    load_go_resources   Carga la ontología (.obo) y la anotación gene -> GO, y
                        registra sus versiones en el log.
    build_background    Define el conjunto de fondo (universo) con los genes expresados.
    run_go_enrichment   Ejecuta la prueba exacta de Fisher con corrección múltiple
                        en los aspectos de GO y devuelve un DataFrame enriquecido.
    plot_enrichment     Genera la gráfica opcional de los términos más enriquecidos.
    setup_logging       Configura el logging del módulo.
    my_parser           Define y procesa los argumentos de línea de comandos.
    main                Orquesta el flujo completo de enriquecimiento GO.
"""
import argparse as ap
import logging
import os
import sys

import pandas as pd
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.go_enrichment import GOEnrichmentStudyNS
from goatools.obo_parser import GODag

logger = logging.getLogger(__name__)

def load_results(ruta_resultados: str) -> pd.DataFrame:
    """
    Busca y lee la tabla de resultados de expresión diferencial

    Args:
        ruta_resultados (str) -> Ruta hacia el archivo de la tabla de resultados realizada
            por el módulo de diferential_expression.py

    Returns:
        df (pd.DataFrame) -> Dataframe representativo de los resultados
    """
    if not os.path.exists(ruta_resultados):
        logger.error("No se encontró el archivo de resultados %s", ruta_resultados)
        raise FileNotFoundError(f"No se encontró el archivo de resultados: {ruta_resultados}")
    df_res = pd.read_csv(ruta_resultados, sep="\t", index_col=0)
    return df_res


def select_genes(df_res: pd.DataFrame, lfc_umbral: float = 1.0, padj_umbral: float = 0.05):
    """
    Extrae las listas de genes sobreexpresados (up) y subexpresados (down) del DataFrame
    representativo de los resultados generados por el módulo diferential_expression.py

    Args:
        df_res (pd.DataFrame) -> Dataframe representativo de los resultados
        lfc_unmbral (float) -> Valor de corte mínimo de Log2FoldChange para determinar genes
            sobreexpresados y subexpresados
        padj_umbral (float) -> Valor de corte mínimo del padj para determinar genes sobrexpresados
            y subexpresados

    Returns:
        genes_up (pd.DataFrame) -> DataFrame con los genes sobreexpresados
        genes_down (pd.DataFrame) -> DataFrame con los genes subexpresados
    """
    if 'padj' in df_res.columns and df_res['padj'].notna().any():
        sig = df_res[(df_res['padj'].fillna(1) < padj_umbral) & (df_res['log2FoldChange'].abs() >= lfc_umbral)]
    else:
        sig = df_res[df_res['log2FoldChange'].abs() >= lfc_umbral]
    
    genes_up = sig[sig['log2FoldChange'] > 0].index.tolist()
    genes_down = sig[sig['log2FoldChange'] < 0].index.tolist()

    return genes_up, genes_down

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
    gene2gos = anno_reader.get_ns2gene2gos()
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
            diferencialmente expresados (ej. genes_up o genes_down).
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
    
    return parser.parse_args()


def main() -> None:
    """
    Orquesta el flujo completo de enriquecimiento funcional con Gene Ontology.
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
        genes_up, genes_down = select_genes(df_res)
        godag, gene2gos = load_go_resources(obo_path, anno_path, taxid=args.taxid)
        background = build_background(df_res)

        tablas = []
        for nombre, lista in [("up_resistant", genes_up), ("up_susceptible", genes_down)]:
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

    except (FileNotFoundError, ValueError, RuntimeError) as e:
        logger.error("Error crítico en el análisis GO: %s", e)
        sys.exit(1)

if __name__ == "__main__":
    main()