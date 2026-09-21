import multiqc
from multiqc.core.update_config import ClConfig
from pathlib import Path
import pandas as pd

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

def run_multiqc():
    # ===============================================================================================================================
    # run_multiqc se encargará de recorrer una lista de los directorios que contengan los archivos .fasta de datos de 
    # expresión del organismo Anopheles gambiae. Así como identificará la calidad de las lecturas y escribirá un dataframe
    # que se guardará en la carpeta .results, en dicho Dataframe se identificará cada lectura y se dará un resumen de la calidad. 
    #
    # Args:
    #
    # Returns:
    #      output_dirs (List) -> Lista de rutas de salida generadas por MultiQC, una por grupo
    #
    # ===============================================================================================================================
    dirs_fasta = []
    dirs_fasta.append(Path("../data/raw/resistent"))
    dirs_fasta.append(Path("../data/raw/susceptible"))

    output_dirs = []

    for i in dirs_fasta:
        output_dir = f"../results/qc/{i.name}"
        multiqc.run(i, cfg=ClConfig(output_dir=output_dir), return_html=True)
        output_dirs.append(output_dir)

    return output_dirs

def read_tsv(output_dirs):
    # ===============================================================================================================================
    # read_tsv se encarga de leer el tsv de FastQC, su única responsabilidad es parsear el archivo y devolver un DataFrame por
    # muestra, tal y como vienen de MultiQC
    #
    # Args: 
    #      output_dirs (List) -> Lista de rutas de salida generadas por run_multiqc, una por grupo
    #
    # Returns:
    #      df_qc (DataFrame) -> DataFrame por muestra de MultiQC, parseado y por el orden del MultiQC
    #
    # ===============================================================================================================================

    dfs = []

    for output_dir in output_dirs:
        tsv_path = Path(output_dir) /"multiqc_data"/"multiqc_fastqc.txt"
        df_grupo = pd.read_csv(tsv_path, sep="\t")
        df_grupo["grupo"] = Path(output_dir).name
        dfs.append(df_grupo)

    df_qc = pd.concat(dfs, ignore_index=True)

    return df_qc

def classify_sample(df_qc):
    # ===============================================================================================================================
    # classify_sample se encarga de construir una estructura de diccionarios, en por cada muestra se asocia la calidad, dependendiendo
    # principalmente de la presencia de adaptadores en la muestra. De ésta manera se clasifican las muestras sobre las que realmente se
    # van a usar en el análisis
    # 
    # Args: 
    #     df_qc (DataFrame) -> Estructura de datos creada por la función read_tsv(). DataFrame de la información de las lecturas por el 
    #                     MultiQC
    #
    # Returns:
    #     dict_qc (Diccionario) -> Diccionario en el que se le asocia la calidad de la muestra con el identificador de la lectura
    # ===============================================================================================================================
    dict_qc = {}

    for i, row in df_qc.iterrows():
        sample = row["Sample"]

        banderas = {col: row[col] for col in FASTQC_FLAG_COLUMNS}
        banderas["overall"] = banderas["adapter_content"]
        banderas["grupo"] = row["grupo"]

        dict_qc[sample] = banderas

    return dict_qc

def delete_fails(dict_qc):
    # ===============================================================================================================================
    # delete_fails detecta los SRRs que son clasificados como fails, los elimina de los archivos que se usarán en el análisis.
    # 
    # Args: 
    #     dict_qc (Diccionario) ->  Diccionario creado por la función classify_sample
    #                     #
    # Returns:
    #     dict_qc_filtered (Diccionario) -> Diccionario que asocia la calidad de la muestra con el identificador de la lectura clasificada
    #                                       por la calidad de la lectura eliminando los "fails"
    # ===============================================================================================================================
    dict_qc_filtered = {
        sample: banderas
        for sample, banderas in dict_qc.items()
        if banderas["overall"] != "fail"
    }

    return dict_qc_filtered