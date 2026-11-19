"""
differential_expression.py
==========================

Módulo para el análisis de expresión diferencial entre muestras resistentes y
susceptibles, a partir de la matriz de conteos generada por featureCounts.

Lee y limpia la matriz de conteos, arma la tabla de metadatos, filtra los genes
de baja expresión y corre el análisis con InMoose (DESeq2) si existen réplicas
suficientes (>= 3 muestras por condición), tomando susceptible como nivel de 
referencia (un log2FoldChange positivo indica mayor expresión en resistentes). 
Cuando no hay réplicas suficientes (ej. experimentos con 2 muestras), calcula
automáticamente el log2FoldChange directo sobre log2(CPM + 1). Después clasifica
los genes (up_resistant, up_susceptible o ns), guarda los resultados y, 
opcionalmente, genera los gráficos correspondientes (PCA y volcano plot o volcano 
sin p-valor).

Requisitos:
    - InMoose (inmoose.deseq2)
    - pandas, numpy, matplotlib y scikit-learn
    - Matriz de conteos de featureCounts (genes × muestras, un BAM por muestra).

Funciones principales:
    load_counts                 Lee la matriz de featureCounts y la deja limpia
                                (genes × muestras, sin columnas de anotación).
    infer_conditions            Obtiene la condición de cada muestra a partir de la
                                carpeta del BAM en el encabezado de featureCounts.
    load_metadata               Crea la tabla de metadatos con la condición de
                                cada muestra.
    filter_genes                Elimina genes de baja expresión
                                (>= 10 conteos en al menos 2 muestras, o 1 si hay <= 2).
    run_differential_expression Ejecuta el análisis de expresión diferencial
                                con InMoose (DESeq2) si hay réplicas suficientes.
    classify_genes              Clasifica los genes en up_resistant,
                                up_susceptible o ns
                                (padj < 0.05 y |log2FoldChange| >= 1).
    save_results                Guarda la tabla completa y las listas de genes
                                por categoría.
    plot_pca                    Grafica el PCA de las muestras (log2(CPM+1)) si
                                hay >= 3 muestras en total; con <= 2 se omite.
    plot_volcano                Grafica el volcano plot de los genes.
    plot_volcano_sin_p          Grafica el volcano plot basado en log2FoldChange
                                cuando no existen p-valores ajustados.
    setup_logging               Configura el logging del módulo.
    my_parser                   Define los argumentos de línea de comandos.
    main                        Ejecuta el flujo completo.

Archivos generados (rutas relativas al directorio donde se ejecuta):
    results/exp_dif/resultados_completos.csv   Tabla completa con la columna
                                               "categoria".
    results/exp_dif/up_resistant.txt           Genes más expresados en resistentes.
    results/exp_dif/up_susceptible.txt         Genes más expresados en susceptibles.
    results/images/pca.png                     (si hay > 2 muestras en total)
    results/images/volcano.png                 (con réplicas / padj)
    results/images/volcano_sin_p.png           (sin réplicas / solo log2FoldChange)
    logs/differential_expression.log

Argumentos:
    -m,   --matrix        Ruta a la matriz cruda de featureCounts (default: data/counts/raw_counts.txt).
    -md,  --metadata      Tabla de metadatos con columnas "sample" y "condition" (opcional).
    -o,   --outdir        Directorio de salida de los resultados (default: results/exp_dif).
    -c1,  --condition1    Condición de prueba (default: resistant).
    -c2,  --condition2    Condición de referencia (default: susceptible).
    -pca, --pcaplot       Incluye el gráfico de PCA (default).
          --no-pcaplot    Omite el gráfico de PCA.
    -vol, --volcanoplot   Incluye el volcano plot (default).
          --no-volcanoplot Omite el volcano plot.
    -v,   --verbose       Muestra mensajes DEBUG.
    --log_file            Archivo de log (default: logs/differential_expression.log).

Uso:
    python src/differential_expression.py -v
    python src/differential_expression.py -m data/counts/raw_counts.txt --no-pcaplot --no-volcanoplot

Nota:
    Ejecutar desde la raíz del proyecto, ya que las rutas de salida por defecto
    son relativas. Si no se indica --metadata, la condición de cada muestra se
    toma de la carpeta de su BAM en el encabezado de la matriz cruda
    (data/aligned/<condición>/<muestra>.bam), por lo que no sirve la matriz ya limpia.
"""

import argparse as ap
import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from inmoose.deseq2 import DESeq, DESeqDataSet
from sklearn.decomposition import PCA

logger = logging.getLogger(__name__)

COLUMNAS_ANOTACION = ["Chr", "Start", "End", "Strand", "Length"]

def load_counts(ruta_matriz: str) -> pd.DataFrame:
    """
    Lee la matriz de conteos de featureCounts y la deja limpia (genes × muestras, 
    sin las columnas de anotación). Indexa en la columna "Geneid" y quita el prefijo .bam
    de las muestras.

    Args:
        ruta_matriz (str) -> Cadena de texto que sea la ruta hacia la matriz de conteos.

    Returns:
        matriz_conteos (pd.DataFrame) -> Matriz de conteos parseada sin columnas innecesarias, indexada 
            en "Geneid" y con los nombres de las muestras correctos.
    """
    if not os.path.exists(ruta_matriz):
        logger.error("No se encontró la matriz de conteos en: %s", ruta_matriz)
        raise FileNotFoundError(f"Matriz de conteos faltante: {ruta_matriz}")

    df = pd.read_csv(ruta_matriz, sep="\t", comment="#")
    df = df.set_index("Geneid")

    col_presentes = [c for c in COLUMNAS_ANOTACION if c in df.columns]
    matriz_conteos = df.drop(columns=col_presentes)
    
    matriz_conteos.columns = [os.path.basename(c).removesuffix(".bam") for c in matriz_conteos.columns]

    assert matriz_conteos.index.is_unique
    matriz_conteos = matriz_conteos.astype(int)

    return matriz_conteos

def infer_conditions(ruta_matriz: str) -> dict[str, str]:
    """
    Obtiene la condición de cada muestra a partir de la carpeta del BAM que featureCounts
    deja en el encabezado de la matriz cruda (por ejemplo, data/aligned/resistant/SRR1.bam
    corresponde a la condición "resistant").

    Args:
        ruta_matriz (str) -> Ruta a la matriz cruda de featureCounts (con las rutas de los BAM
            en el encabezado).

    Returns:
        condiciones (dict[str, str]) -> Diccionario {muestra: condición}.
    """
    encabezado = pd.read_csv(ruta_matriz, sep="\t", comment="#", nrows=0)
    columnas = [c for c in encabezado.columns if c != "Geneid" and c not in COLUMNAS_ANOTACION]

    condiciones = {}
    for col in columnas:
        carpeta = os.path.basename(os.path.dirname(col))
        if not carpeta:
            raise ValueError(
                f"La columna '{col}' no trae la ruta del BAM, no se puede inferir la condición. "
                "Usa la matriz cruda de featureCounts o indica un archivo con --metadata."
            )
        condiciones[os.path.basename(col).removesuffix(".bam")] = carpeta

    return condiciones

def load_metadata(matriz_conteos: pd.DataFrame, ruta_matriz: str, metadata_path: str | None = None,
                  condition1: str = "resistant", condition2: str = "susceptible") -> pd.DataFrame:
    """
    Crea la tabla de metadatos que usará InMoose. Si se indica un archivo de metadatos, lo lee;
    de lo contrario, toma la condición de cada muestra de la carpeta de su BAM en la matriz cruda.
    La condición de referencia es condition2, por lo que un log2FoldChange positivo indica mayor
    expresión en condition1.

    Args:
        matriz_conteos (pd.DataFrame) -> Matriz de conteos limpia (genes × muestras).
        ruta_matriz (str) -> Ruta a la matriz cruda de featureCounts.
        metadata_path (str | None) -> Ruta a un archivo de metadatos con las columnas "sample"
            y "condition" (opcional).
        condition1 (str) -> Condición de prueba (default: resistant).
        condition2 (str) -> Condición de referencia (default: susceptible).

    Returns:
        meta (pd.DataFrame) -> Tabla de metadatos indexada por muestra, con la columna "condition"
            como categórica (referencia: condition2) y en el mismo orden que las columnas de la matriz.
    """
    if metadata_path:
        if not os.path.exists(metadata_path):
            logger.error("No se encontró el archivo de metadatos: %s", metadata_path)
            raise FileNotFoundError(f"Metadatos faltantes: {metadata_path}")

        meta = pd.read_csv(metadata_path, sep=None, engine="python")
        if "sample" in meta.columns:
            meta = meta.set_index("sample")
        else:
            meta = meta.set_index(meta.columns[0])

        faltantes = [m for m in matriz_conteos.columns if m not in meta.index]
        if faltantes:
            raise ValueError(f"Muestras sin condición en {metadata_path}: {faltantes}")

        meta = meta.loc[matriz_conteos.columns, ["condition"]]
    else:
        condiciones = infer_conditions(ruta_matriz)

        faltantes = [m for m in matriz_conteos.columns if m not in condiciones]
        if faltantes:
            raise ValueError(f"Muestras sin condición inferida: {faltantes}")

        meta = pd.DataFrame({"condition": [condiciones[m] for m in matriz_conteos.columns]},
                            index=matriz_conteos.columns)

    desconocidas = set(meta["condition"]) - {condition1, condition2}
    if desconocidas:
        raise ValueError(f"Condiciones no reconocidas {sorted(desconocidas)}; se esperaban "
                         f"'{condition1}' y '{condition2}'")

    meta["condition"] = pd.Categorical(meta["condition"], categories=[condition2, condition1])

    conteo_por_condicion = meta["condition"].value_counts()
    for condicion, n in conteo_por_condicion.items():
        if n < 2:
            logger.warning("La condición '%s' tiene %d muestra(s); DESeq2 necesita réplicas para "
                           "estimar la dispersión", condicion, n)

    return meta

def filter_genes(matriz_conteos:pd.DataFrame) -> pd.DataFrame:
    """
    Elimina genes de baja expresión (>= 10 conteos en al menos 2 muestras, o en al menos
    1 muestra si el experimento solo cuenta con 2 muestras en total). Con la finalidad
    de aumentar la potencia estadística al eliminar genes que no tienen información suficiente, al igual
    que reducir el ruido computacional al eliminar genes no informativos.

    Args:
        matriz_conteos (pd.DataFrame) -> Matriz de conteos parseada sin columnas innecesarias, indexada 
            en "Geneid" y con los nombres de las muestras correctos.
        
    Returns:
        matriz_filtrada (pd.DataFrame) -> Matriz de conteos filtrada sin genes no informativos.
    """
    UMBRAL_CONTEOS = 10
    MINIMO_MUESTRAS = 1 if matriz_conteos.shape[1] <= 2 else 2
    
    mask = (matriz_conteos >= UMBRAL_CONTEOS).sum(axis=1) >= MINIMO_MUESTRAS
    matriz_filtrada = matriz_conteos[mask]

    logger.info("Filtrado de genes: %d de %d genes retenidos (>= %d conteos en >= %d muestra(s))",
                len(matriz_filtrada), len(matriz_conteos), UMBRAL_CONTEOS, MINIMO_MUESTRAS)

    return matriz_filtrada

def run_differential_expression(matriz_filtrada:pd.DataFrame, meta:pd.DataFrame) -> pd.DataFrame:
    """
    Ejecuta el análisis de expresión diferencial comparando condition1 vs condition2
    (con condition2 como referencia; log2FoldChange positivo indica mayor expresión en condition1).

    - Si hay >= 3 muestras por condición: ejecuta InMoose (DESeq2).
    - Si hay < 3 muestras en alguna condición (ej. solo 2 muestras): calcula directamente
      el log2FoldChange a partir de log2(CPM + 1), dejando pvalue y padj como NaN.
    

    Args:
        matriz_filtrada (pd.DataFrame) -> Matriz de conteos crudos (enteros) genes × muestras, sin
            genes de baja expresión.
        meta (pd.DataFrame) -> Tabla de metadatos indexada por muestra, con la columna "condition"
            como categórica (referencia: susceptible).

    Returns:
        df_res (pd.DataFrame) -> Tabla de resultados por gen (baseMean, log2FoldChange, pvalue, padj,
            etc.) ordenada por padj de menor a mayor.
    """
    conteo_por_condicion = meta["condition"].value_counts()
    tiene_replicas = (conteo_por_condicion >= 3).all() and len(conteo_por_condicion) >= 2

    if tiene_replicas:  
        logger.info("Réplicas suficientes detectadas (%s). Ejecutando DESeq2 con InMoose.",
                    conteo_por_condicion.to_dict())

        dds = DESeqDataSet(countData = matriz_filtrada.T,
                            clinicalData = meta,
                            design = "~ condition"
                            )
    
        dds = DESeq(dds)
        res = dds.results()
        df_res = pd.DataFrame(res).sort_values("padj")
        return df_res

    logger.warning(
        "Réplicas insuficientes (%s). Calculando log2FoldChange directo sobre log2(CPM + 1).",
        conteo_por_condicion.to_dict(),
    )
    
    cond_ref, cond_test = list(meta["condition"].cat.categories)
    muestras_test = meta.index[meta["condition"] == cond_test]
    muestras_ref = meta.index[meta["condition"] == cond_ref]

    cpm = matriz_filtrada.div(matriz_filtrada.sum(axis=0), axis=1) * 1e6
    log_cpm = np.log2(cpm + 1.0)

    mean_test = log_cpm[muestras_test].mean(axis=1)
    mean_ref = log_cpm[muestras_ref].mean(axis=1)

    df_res = pd.DataFrame(
        {
            "baseMean": cpm.mean(axis=1),
            f"log2CPM_{cond_test}": mean_test,
            f"log2CPM_{cond_ref}": mean_ref,
            "log2FoldChange": mean_test - mean_ref,
            "pvalue": np.nan,
            "padj": np.nan,
        },
        index=matriz_filtrada.index,
    )

    df_res = df_res.reindex(df_res["log2FoldChange"].abs().sort_values(ascending=False).index)
    return df_res

def classify_genes(df_res: pd.DataFrame) -> pd.DataFrame:
    """
    Clasifica cada gen en 'up_resistant', 'up_susceptible' o 'ns'. Cuando existen
    p-valores ajustados utiliza padj < 0.05 y |log2FoldChange| >= 1.0; cuando no
    hay p-valores (análisis sin réplicas suficientes), clasifica únicamente por
    |log2FoldChange| >= 1.0.

    Args:
        df_res (pd.DataFrame) -> Tabla de resultados devuelta por run_differential_expression,
            con las columnas "log2FoldChange" y "padj".

    Returns:
        df_clasificado (pd.DataFrame) -> Copia de df_res con la columna adicional "categoria".
    """
    df_clasificado = df_res.copy()

    if df_clasificado["padj"].isna().all():
        logger.warning("No hay p-valores ajustados. Clasificando solo por |log2FoldChange| >= 1.0")
        condiciones = [
            (df_clasificado["log2FoldChange"] >= 1.0),
            (df_clasificado["log2FoldChange"] <= -1.0),
        ]
    else:
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
    etiquetada con su nombre y coloreada según su condición (azul: condición de referencia,
    rojo: condición de prueba). Los ejes indican el porcentaje de varianza explicada por PC1 y PC2.

    Args:
        matriz_filtrada (pd.DataFrame) -> Matriz de conteos crudos (enteros) genes × muestras, sin
            genes de baja expresión.
        meta (pd.DataFrame) -> Tabla de metadatos indexada por muestra, con la columna "condition"
            como categórica (la primera categoría es la referencia).
        outdir (str) -> Directorio donde se guarda la figura (por defecto results/images/).

    Returns:
        None. Guarda la figura en <outdir>/pca.png (300 dpi).
    """

    cpm = matriz_filtrada / matriz_filtrada.sum(axis=0) * 1e6
    log_cpm = np.log2(cpm + 1)

    pca = PCA(n_components=2)
    coords = pca.fit_transform(log_cpm.T)
    var = pca.explained_variance_ratio_ * 100

    categorias = list(meta["condition"].cat.categories)
    COLORES = {categorias[0]: "#1f77b4", categorias[1]: "#d62728"}

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

def plot_volcano(df_clasificado: pd.DataFrame, outdir: str = "results/images/",
                 condition1: str = "resistant", condition2: str = "susceptible") -> None:
    """
    Genera el volcano plot de la expresión diferencial ({condition1} vs {condition2}) y lo guarda como
    imagen. Cada punto es un gen: el eje X es el log2FoldChange y el eje Y es -log10(padj). Los
    genes se colorean según la columna "categoria" (rojo: up_resistant, azul: up_susceptible,
    gris: ns). Los genes con padj = NaN se excluyen del gráfico porque no se pueden ubicar en el eje Y.

    Args:
        df_clasificado (pd.DataFrame) -> Tabla de resultados devuelta por classify_genes, con las
            columnas "log2FoldChange", "padj" y "categoria".
        outdir (str) -> Directorio donde se guarda la figura (se crea si no existe).
        condition1 (str) -> Condición de prueba (default: resistant).
        condition2 (str) -> Condición de referencia (default: susceptible).

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

    ax.set_xlabel(f"log2FoldChange ({condition1} vs {condition2})")
    ax.set_ylabel("-log10(padj)")
    ax.legend(frameon=False)

    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, "volcano.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    
def plot_volcano_sin_p(df_raise: pd.DataFrame, outdir: str = "results/images/",
                       condition1: str = "resistant", condition2: str = "susceptible") -> None:
    """
    Genera un volcano plot basado únicamente en log2FoldChange cuando no se cuenta
    con p-valores ajustados (por ejemplo, análisis sin réplicas biológicas).
    El eje X representa log2FoldChange y el eje Y muestra la expresión promedio
    log2(baseMean + 1) para visualizar la magnitud del cambio por nivel de expresión.

    Args:
        df_raise (pd.DataFrame) -> Tabla de resultados clasificada con las columnas
            "log2FoldChange", "baseMean" y "categoria".
        outdir (str) -> Directorio donde se guarda la figura (default: results/images/).
        condition1 (str) -> Condición de prueba (default: resistant).
        condition2 (str) -> Condición de referencia (default: susceptible).

    Returns:
        None. Guarda la figura en <outdir>/volcano_sin_p.png (300 dpi).
    """
    COLORES = {
        "up_resistant": "#d62728",
        "up_susceptible": "#1f77b4",
        "ns": "#b0b0b0",
    }

    df = df_raise.dropna(subset=["log2FoldChange"]).copy()
    y_vals = (
        np.log2(df["baseMean"] + 1.0)
        if "baseMean" in df.columns
        else pd.Series(np.zeros(len(df)), index=df.index)
    )
    df["y_expr"] = y_vals

    fig, ax = plt.subplots(figsize=(7, 6))
    for cat in ["ns", "up_susceptible", "up_resistant"]:
        sub = df[df["categoria"] == cat]
        ax.scatter(
            sub["log2FoldChange"],
            sub["y_expr"],
            c=COLORES[cat],
            s=10,
            alpha=0.7,
            label=f"{cat} (n={len(sub)})",
        )

    ax.axvline(1.0, color="#333333", linestyle="--", linewidth=1)
    ax.axvline(-1.0, color="#333333", linestyle="--", linewidth=1)
    ax.axvline(0.0, color="#666666", linestyle="-", linewidth=0.8)

    ax.set_xlabel(f"log2FoldChange ({condition1} vs {condition2})")
    ax.set_ylabel("log2(baseMean + 1)" if "baseMean" in df.columns else "Genes")
    ax.set_title("Volcano Plot sin p-valor (basado en log2FoldChange)")
    ax.legend(frameon=False)

    os.makedirs(outdir, exist_ok=True)
    ruta_salida = os.path.join(outdir, "volcano_sin_p.png")
    fig.savefig(ruta_salida, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("Volcano plot sin p-valor guardado en: %s", ruta_salida)

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

def plot_heatmap(matriz_filtrada: pd.DataFrame, df_clasificado: pd.DataFrame, outdir: str = "results/images/", top_n: int = 25) -> None: 
    """
    Genera un Heatmap de los genes más significativos (Top N up y Top N down) 
    utilizando Seaborn. 

    Args:
        matriz_filtrada (pd.DataFrame): Matriz de conteos filtrada.
        df_clasificado (pd.DataFrame): Resultados con la columna "categoria".
        outdir (str): Directorio de salida para la imagen.
        top_n (int): Número de genes top a graficar.
    
    Returns:
        None
    """
    
    up = df_clasificado[df_clasificado["categoria"] == "up_resistant"].head(top_n).index
    down = df_clasificado[df_clasificado["categoria"] == "up_susceptible"].head(top_n).index
    genes = list(up) + list(down)
    
    if not genes:
        logger.warning("No hay genes significativos para graficar en el Heatmap.")
        return 
    
    sub_matriz = matriz_filtrada.loc[genes]
    cpm = sub_matriz / matriz_filtrada.sum(axis=0) * 1e6
    log_cpm = np.log2(cpm + 1)

    std_dev = log_cpm.std(axis=1).replace(0, 1)
    z_scores = log_cpm.sub(log_cpm.mean(axis=1), axis=0).div(std_dev, axis=0)

    fig, ax = plt.subplots(figsize=(8, max(6, len(genes) * 0.35)))
    
    sns.heatmap(
        z_scores, 
        cmap="RdBu_r",
        center=0,
        vmin=-2, vmax=2,
        annot=False,
        linewidths=0.5,
        cbar_kws={'label': 'Z-score'},
        ax=ax
    )
    
    # 5. Personalización
    ax.set_title(f"Expresión de los Top {top_n} genes significativos\n(Sobreexpresados y Subexpresados)", fontsize=12, pad=15)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)

    os.makedirs(outdir, exist_ok=True)
    output_path = os.path.join(outdir, "heatmap.png")
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    
    logger.info("Heatmap guardado exitosamente en %s", output_path)


def my_parser() -> ap.Namespace:
    """
    Retiene los argumentos deseados para el módulo desde línea de comandos.

    Args:
        None
    
    Returns:
        ap.Namespace: Argumentos indicados por el usuario
    """    
    parser = ap.ArgumentParser(
    description="Pipeline para el análisis de expresión diferencial con el uso de InMoose")
    parser.add_argument("-m", "--matrix", required=False, default="data/counts/raw_counts.txt",
                        help="Ruta a la matriz cruda generada por featureCounts")
    parser.add_argument("-md", "--metadata", required=False, default=None,
                        help="Tabla de metadatos con las columnas 'sample' y 'condition' (opcional)")
    parser.add_argument("-o", "--outdir", required=False, default="results/exp_dif", help="Directorio de salida")
    parser.add_argument("-pca", "--pcaplot", action=ap.BooleanOptionalAction, default=True,
                        help="Incluye el grafico de PCA para la comparación de genes (usa --no-pcaplot para omitirlo)")
    parser.add_argument("-c1", "--condition1", required=False, default="resistant",
                        help="Condición de prueba (default: resistant)")
    parser.add_argument("-c2", "--condition2", required=False, default="susceptible",
                        help="Condición de referencia (default: susceptible)")
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
        5. Corre el análisis de expresión diferencial: corre InMoose (DESeq2) si
           ambas condiciones cuentan con >= 3 réplicas, o calcula log2FoldChange
           directo sobre log2(CPM + 1) cuando no hay réplicas suficientes.
        6. Clasifica los genes en up_resistant, up_susceptible o ns usando
           padj < 0.05 y |log2FoldChange| >= 1.0 (o solo |log2FoldChange| >= 1.0 si
           no existen p-valores ajustados).
        7. Guarda en --outdir la tabla completa y las listas de genes por categoría.
        8. Si --pcaplot está activo existen más de 2 muestras en total, genera el
           gráfico de PCA en results/images/pca.png; con <= 2 muestras lo omite.
        9. Si --volcanoplot está activo, genera el volcano plot si hay padj, o 
           el volcano plot sin p-valor (volcano_sin_p.png) si el análisis fue 
           exploratorio.

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
    meta = load_metadata(matriz_conteos, args.matrix, args.metadata,
                         condition1=args.condition1, condition2=args.condition2)

    # Filtramos la matriz por incidencia de genes
    matriz_filtrada = filter_genes(matriz_conteos)

    # Corremos InMoos para el análisis de expresión diferencial
    df_res = run_differential_expression(matriz_filtrada, meta)

    # Clasificamos como sobreexpresados y subexpresados a los genes
    df_clasificado = classify_genes(df_res)

    # Se guardan los resultados en la carpeta "results/exp_dif"
    save_results(df_clasificado, args.outdir)

    n_muestras = len(matriz_filtrada.columns)
    # Se plotea el PCA plot si hay más de 2 muestras en total
    if args.pcaplot:
        if n_muestras <= 2:
            logger.warning("Se detectaron %d muestra(s); se requieren al menos 3 muestras en total para generar el PCA plot.",
                           n_muestras,)
        else:
            plot_pca(matriz_filtrada, meta)

    # Se plotea el volcano plot
    if args.volcanoplot:
        if df_clasificado["padj"].isna().all():
            plot_volcano_sin_p(df_clasificado, condition1=args.condition1, condition2=args.condition2)
        else:
            plot_volcano(df_clasificado, condition1=args.condition1, condition2=args.condition2)
    logger.info("Análisis de expresión diferencial finalizado correctamente")

if __name__ == "__main__":
    main()