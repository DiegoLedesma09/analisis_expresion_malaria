"""
differential_expression.py
==========================

Módulo para el análisis de expresión diferencial entre muestras resistentes y
susceptibles, a partir de la matriz de conteos generada por featureCounts.

Lee y limpia la matriz de conteos, arma la tabla de metadatos, filtra los genes
de baja expresión y corre el análisis con InMoose (DESeq2), tomando
susceptible como nivel de referencia (un log2FoldChange positivo indica mayor
expresión en resistentes). Después clasifica los genes (up_resistant,
up_susceptible o ns), guarda los resultados y, opcionalmente, genera el PCA de
las muestras y el volcano plot.

Requisitos:
    - InMoose (inmoose.deseq2)
    - pandas, numpy, matplotlib y scikit-learn
    - Matriz de conteos de featureCounts (genes × muestras, un BAM por muestra).

Funciones principales:
    load_counts                 Lee la matriz de featureCounts y la deja limpia
                                (genes × muestras, sin columnas de anotación).
    load_metadata               Crea la tabla de metadatos con la condición de
                                cada muestra.
    filter_genes                Elimina genes de baja expresión
                                (>= 10 conteos en al menos 2 muestras).
    run_differential_expression Ejecuta el análisis de expresión diferencial
                                con InMoose (DESeq2).
    classify_genes              Clasifica los genes en up_resistant,
                                up_susceptible o ns
                                (padj < 0.05 y |log2FoldChange| >= 1).
    save_results                Guarda la tabla completa y las listas de genes
                                por categoría.
    plot_pca                    Grafica el PCA de las muestras (log2(CPM+1)).
    plot_volcano                Grafica el volcano plot de los genes.
    setup_logging               Configura el logging del módulo.
    my_parser                   Define los argumentos de línea de comandos.
    main                        Ejecuta el flujo completo.

Archivos generados (rutas relativas al directorio donde se ejecuta):
    results/exp_dif/resultados_completos.csv   Tabla completa con la columna
                                               "categoria".
    results/exp_dif/up_resistant.txt           Genes más expresados en resistentes.
    results/exp_dif/up_susceptible.txt         Genes más expresados en susceptibles.
    results/images/pca.png
    results/images/volcano.png
    logs/differential_expression.log

Argumentos:
    -m,   --matrix        Ruta a la matriz de conteos de featureCounts (obligatorio).
    -o,   --outdir        Directorio de salida de los resultados (default: results/exp_dif).
    -pca, --pcaplot       Incluye el gráfico de PCA (default).
          --no-pcaplot    Omite el gráfico de PCA.
    -vol, --volcanoplot   Incluye el volcano plot (default).
          --no-volcanoplot Omite el volcano plot.
    -v,   --verbose       Muestra mensajes DEBUG.
    --log_file            Archivo de log (default: logs/differential_expression.log).

Uso:
    python src/differential_expression.py -m results/matriz_conteos.txt -v
    python src/differential_expression.py -m results/matriz_conteos.txt --no-pcaplot --no-volcanoplot

Nota:
    Ejecutar desde la raíz del proyecto, ya que las rutas de salida por defecto
    son relativas. La tabla de metadatos asume 4 muestras en este orden de
    columnas: 2 resistentes y después 2 susceptibles.
"""

import argparse as ap
import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from inmoose.deseq2 import DESeq, DESeqDataSet
from sklearn.decomposition import PCA

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
    meta["condition"] = pd.Categorical(meta["condition"],
                                   categories=["susceptible", "resistant"])
    return meta

def filter_genes(matriz_conteos:pd.DataFrame) -> pd.DataFrame:
    """
    Elimina genes cuyo umbral de conteos es >= 10 y con un mínimo de muestras >= 2 con la finalidad
    de aumentar la potencia estadística al eliminar genes que no tienen información suficiente, al igual
    que reducir el ruido computacional al eliminar genes no informativos.

    Args:
        matriz_conteos (pd.DataFrame) -> Matriz de conteos parseada sin columnas innecesarias, indexada 
            en "Geneid" y con los nombres de las muestras correctos.
        
    Returns:
        matriz_filtrada (pd.DataFrame) -> Matriz de conteos filtrada sin genes no informativos.
    """
    UMBRAL_CONTEOS = 10
    MINIMO_MUESTRAS = 2
    
    mask = (matriz_conteos >= UMBRAL_CONTEOS).sum(axis=1) >= MINIMO_MUESTRAS
    matriz_filtrada = matriz_conteos[mask]

    return matriz_filtrada

def run_differential_expression(matriz_filtrada:pd.DataFrame, meta:pd.DataFrame) -> pd.DataFrame:
    """
    Ejecuta el análisis de expresión diferencial con InMoose (DESeq2) comparando las condiciones
    resistant vs susceptible, con susceptible como nivel de referencia. Un log2FoldChange positivo
    indica mayor expresión en resistentes.

    Args:
        matriz_filtrada (pd.DataFrame) -> Matriz de conteos crudos (enteros) genes × muestras, sin
            genes de baja expresión.
        meta (pd.DataFrame) -> Tabla de metadatos indexada por muestra, con la columna "condition"
            como categórica (referencia: susceptible).

    Returns:
        df_res (pd.DataFrame) -> Tabla de resultados por gen (baseMean, log2FoldChange, pvalue, padj,
            etc.) ordenada por padj de menor a mayor.
    """

    dds = DESeqDataSet(countData = matriz_filtrada,
                        clinicalData = meta,
                        design = "~ condition"
                        )
    
    dds = DESeq(dds)
    res = dds.results()
    df_res = pd.DataFrame(res).sort_values("padj")
    return df_res

def classify_genes(df_res: pd.DataFrame) -> pd.DataFrame:
    """
    Clasifica cada gen según su padj y su log2FoldChange, agregando la columna "categoria".
    No elimina genes: devuelve la tabla completa para poder usarla en el volcano.

    Args:
        df_res (pd.DataFrame) -> Tabla de resultados de expresión diferencial, con las columnas
            "padj" y "log2FoldChange".

    Returns:
        df_clasificado (pd.DataFrame) -> Copia de df_res con la columna "categoria":
            "up_resistant", "up_susceptible" o "ns" (no significativo).
    """
    df_clasificado = df_res.copy()

    sig = df_clasificado["padj"] < 0.05

    condiciones = [
        sig & (df_clasificado["log2FoldChange"] >= 1.0),
        sig & (df_clasificado["log2FoldChange"] <= -1.0),
    ]
    df_clasificado["categoria"] = np.select(
        condiciones, ["up_resistant", "up_susceptible"], default="ns"
    )

    return df_clasificado

def save_results(df_clasificado: pd.DataFrame, outdir: str = "results/exp_dif") -> None:
    """
    Guarda los DataFrames generados antes y después del análisis de expresión
    diferencial, además de separarlos según la condición de muestreo (susceptible o resistente)
    en sus respectivos archivos (resultados_completos.csv, up_resistant.txt y up_susceptible.txt)
    en la carpeta de "results/exp_dif"

    Args:
        df_clasificado (pd.DataFrame) -> Copia de df_res con la columna "categoria":
            "up_resistant", "up_susceptible" o "ns" (no significativo).
        
        outdir (str) -> Directorio donde se guardan los resultados (por defecto results/exp_dif/).
    """
    os.makedirs(outdir, exist_ok=True)

    df_clasificado.to_csv(os.path.join(outdir, "resultados_completos.csv"))

    for cat in ["up_resistant", "up_susceptible"]:
        genes = df_clasificado.index[df_clasificado["categoria"] == cat]
        pd.Series(genes).to_csv(os.path.join(outdir, f"{cat}.txt"), index=False, header=False)

def plot_pca(matriz_filtrada: pd.DataFrame, meta: pd.DataFrame, outdir: str = "results/images/") -> None:
    """
    Genera el gráfico de PCA de las muestras y lo guarda como imagen. Los conteos se normalizan
    a log2(CPM+1) para corregir por profundidad de secuenciación y estabilizar la varianza, y el
    PCA se calcula sobre las muestras (genes como variables). Cada punto es una muestra,
    etiquetada con su nombre y coloreada según su condición (rojo: resistant, azul: susceptible).
    Los ejes indican el porcentaje de varianza explicada por PC1 y PC2.

    Args:
        matriz_filtrada (pd.DataFrame) -> Matriz de conteos crudos (enteros) genes × muestras, sin
            genes de baja expresión.
        meta (pd.DataFrame) -> Tabla de metadatos indexada por muestra, con la columna "condition"
            ("resistant" o "susceptible").
        outdir (str) -> Directorio donde se guarda la figura (por defecto results/images/).

    Returns:
        None. Guarda la figura en <outdir>/pca.png (300 dpi).
    """

    cpm = matriz_filtrada / matriz_filtrada.sum(axis=0) * 1e6
    log_cpm = np.log2(cpm + 1)

    pca = PCA(n_components=2)
    coords = pca.fit_transform(log_cpm.T)
    var = pca.explained_variance_ratio_ * 100

    COLORES = {"resistant": "#d62728", "susceptible": "#1f77b4"}

    fig, ax = plt.subplots(figsize=(6, 5))
    for i, muestra in enumerate(log_cpm.columns):
        cond = meta.loc[muestra, "condition"]
        ax.scatter(coords[i, 0], coords[i, 1], c=COLORES[cond], s=80,
                   label=cond if cond not in ax.get_legend_handles_labels()[1] else None)
        ax.annotate(muestra, (coords[i, 0], coords[i, 1]),
                    textcoords="offset points", xytext=(6, 6), fontsize=8)

    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.legend(frameon=False)

    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, "pca.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

def plot_volcano(df_clasificado: pd.DataFrame, outdir: str = "results/images/") -> None:
    """
    Genera el volcano plot de la expresión diferencial (resistant vs susceptible) y lo guarda como
    imagen. Cada punto es un gen: el eje X es el log2FoldChange y el eje Y es -log10(padj). Los
    genes se colorean según la columna "categoria" (rojo: up_resistant, azul: up_susceptible,
    gris: ns). Los genes con padj = NaN se excluyen del gráfico porque no se pueden ubicar en el eje Y.

    Args:
        df_clasificado (pd.DataFrame) -> Tabla de resultados devuelta por classify_genes, con las
            columnas "log2FoldChange", "padj" y "categoria".
        outdir (str) -> Directorio donde se guarda la figura (se crea si no existe).

    Returns:
        None.
    """

    COLORES = {
        "up_resistant": "#d62728",
        "up_susceptible": "#1f77b4",
        "ns": "#b0b0b0"
    }

    df = df_clasificado.dropna(subset=["padj"]).copy()
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(7, 6))
    for cat, color in COLORES.items():
        sub = df[df["categoria"] == cat]
        ax.scatter(sub["log2FoldChange"], sub["neg_log10_padj"],
                   c=color, s=10, alpha=0.7, label=f"{cat} (n={len(sub)})")

    ax.set_xlabel("log2FoldChange (resistant vs susceptible)")
    ax.set_ylabel("-log10(padj)")
    ax.legend(frameon=False)

    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, "volcano.png"), dpi=300, bbox_inches="tight")
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
    Retiene los argumentos deseados para el módulo desde línea de comandos.

    Args:
        None
    
    Returns:
        ap.Namespace: Argumentos indicados por el usuario
    """    
    parser = ap.ArgumentParser(
    description="Pipeline para del análisis de expresión diferencial con el uso de InMoose")
    parser.add_argument("-m", "--matrix", required=True, help="Ruta a la matriz de conteo generada por featureCounts")
    parser.add_argument("-o", "--outdir", required=False, default="results/exp_dif", help="Directorio de salida")
    parser.add_argument("-pca", "--pcaplot", action=ap.BooleanOptionalAction, default=True,
                        help="Incluye el grafico de PCA para la comparación de genes (usa --no-pcaplot para omitirlo)")
    parser.add_argument("-vol", "--volcanoplot", action=ap.BooleanOptionalAction, default=True,
                        help="Incluye el grafico de volcano para la comparación de genes (usa --no-volcanoplot para omitirlo)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Muestra mensajes DEBUG")
    parser.add_argument("--log_file", default="logs/differential_expression.log", help="Archivo de log")
    return parser.parse_args()


def main():
    """
    Ejecuta el flujo completo de expresión diferencial con InMoose (DESeq2) para la comparación
    resistant vs susceptible.

    Pasos, en orden:
        1. Configura el logging (consola y archivo).
        2. Lee la matriz de conteos de featureCounts y la deja limpia (genes × muestras).
        3. Crea la tabla de metadatos con la condición de cada muestra.
        4. Filtra los genes de baja expresión.
        5. Corre el análisis de expresión diferencial.
        6. Clasifica los genes en up_resistant, up_susceptible o ns.
        7. Guarda en --outdir la tabla completa, el universo de genes probados y las listas de
           genes por categoría (insumos para el análisis de GO).
        8. Si --pcaplot está activo, genera el PCA de las muestras.
        9. Si --volcanoplot está activo, genera el volcano plot.

    Los argumentos se leen de la línea de comandos mediante my_parser().

    Returns:
        None
    """

    args = my_parser()
    setup_logging(args.verbose, args.log_file)
    logger.info("Inicio del pipeline para el Análisis de Expresión diferencial")

    # Lectura y parseo de la matriz de conteos generada por featureCounts
    matriz_conteos = load_counts(args.matrix)

    # Creación de la matriz de metadatos asociada a las columnas de la matriz de conteos
    meta = load_metadata(matriz_conteos)

    # Filtramos la matriz por incidencia de genes
    matriz_filtrada = filter_genes(matriz_conteos)

    # Corremos InMoos para el análisis de expresión diferencial
    df_res = run_differential_expression(matriz_filtrada, meta)

    # Clasificamos como sobreexpresados y subexpresados a los genes
    df_clasificado = classify_genes(df_res)

    # Se guardan los resultados en la carpeta "results/exp_dif"
    save_results(df_clasificado)

    # Se plotea el PCA plot
    if args.pcaplot:
        plot_pca(matriz_filtrada, meta)

    # Se plotea el volcano plot
    if args.volcanoplot:
        plot_volcano(df_clasificado)

if __name__ == "__main__":
    main()