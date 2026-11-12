"""
differential_expression.py
==============
"""

import argparse as ap
import logging
import os
import sys

import pandas as pd

logger = logging.getLogger(__name__)

def load_counts(ruta_matriz:str) -> pd.DataFrame:
    """
    lee la matriz de conteos de featureCounts y la deja limpia (genes × muestras, 
    sin las columnas de anotación). Indexa en la columna "Geneid" y quita el prefijo .bam
    de las muestras

    Args:
        ruta_matriz (str) -> Cadena de texto que sea la ruta global hacia la matriz de conteos 
            generada por FeatureCounts

    Returns:
        matriz_conteos (pd.DataFrame) -> Matriz de conteos parseada sin columnas innecesarias, indexada 
            en "Geneid" y con los nombres de las muestras correctos.
    """

    df = pd.read_csv(ruta_matriz, sep="\t", comment="#")
    df = df.set_index("Geneid")

    col_del = ["Chr", "Start", "End", "Strand", "Length"]
    matriz_conteos = df.drop(columns=col_del)
    matriz_conteos.columns = [os.path.basename(c).removesuffix(".bam") for c in matriz_conteos.columns]

    assert matriz_conteos.index.is_unique
    matriz_conteos = matriz_conteos.astype(int)

    return matriz_conteos


def load_metadata(matriz_conteos:pd.DataFrame) -> pd.DataFrame:
    """
    Crea la tabla de metadatos que usará InMoose a partir de las dos diferentes poblaciones para
    poder hacer la comparación.

    Args:
        matriz_conteos (pd.DataFrame) -> Matriz de conteos parseada sin columnas innecesarias, indexada 
            en "Geneid" y con los nombres de las muestras correctos.
    
    Returns:
        meta (pd.DataFrame) -> Matriz de metadatos asociada a cada columna según la condición (resistant, 
            susceptible) de cada lectura.
    """

    meta = pd.DataFrame({"condition": ["resistant", "resistant", "susceptible", "susceptible"]},
        index=matriz_conteos.columns,
    )
    return meta

def filter_genes():
    """
    Elimina genes con baja expresión
    """
    

    return

def run_differential_expression():
    return

def classify_genes():
    return

def plot_pca():
    return

def plot_volcano():
    return

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
    Retiene los argumentos deseados para el módulo desde línea de comandos.

    Args:
        None
    
    Returns:
        ap.Namespace: Argumentos indicados por el usuario
    """    
    parser = ap.ArgumentParser(
    description="Pipeline para del análisis de expresión diferencial con el uso de InMoose")
    parser.add_argument("-o", "--outdir", required=False, default="results/exp_dif", help="Directorio de salida")
    parser.add_argument("-w", "--workers", required=False, type=int, default=8, help="Workers para InMoose")
    parser.add_argument("-pca", "--pcaplot", action=ap.BooleanOptionalAction, default=True,
                        help="Incluye el grafico de PCA para la comparación de genes (usa --no-pcaplot para omitirlo)")
    parser.add_argument("-vol", "--volcanoplot", action=ap.BooleanOptionalAction, default=True,
                        help="Incluye el grafico de volcano para la comparación de genes (usa --no-volcanoplot para omitirlo)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Muestra mensajes DEBUG")
    parser.add_argument("--log_file", default="logs/download_data.log", help="Archivo de log")
    return parser.parse_args()


def main():
    args = my_parser()
    setup_logging(args.verbose, args.log_file)
    logger.info("Inicio del pipeline para el Análisis de Expresión diferencial")

if __name__ == "__main__":
    main()