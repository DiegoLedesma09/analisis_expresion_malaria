import argparse as ap
import logging
import os
import sys

import pandas as pd
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.go_enrichment import GOEnrichmentStudy
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
    parser.add_argument("-o", "--outdir", required=False, default="data/counts/", 
                        help="Directorio de salida")
    parser.add_argument("-w", "--workers", type=int, required= False, default= 6, 
                        help="Número de hilos")
    parser.add_argument("-t", "--taxid", type=int, required=False, default=7165)
    parser.add_argument("-v", "--verbose", action="store_true", 
                        help="Muestra mensajes DEBUG en consola")
    parser.add_argument("-g", "--graphic", action=ap.BooleanOptionalAction, default=True,
                    help="Incluye el grafico de la distribución de muestras asignadas")
    parser.add_argument("--log_file", default="logs/feature_counts.log", 
                        help="Ruta para guardar el archivo de log")
    
    return parser.parse_args()


def main():
    load_results()
if __name__ == "__main__":
    main()