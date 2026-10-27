import subprocess
from pathlib import Path

import multiqc
import pandas as pd
from multiqc.core.update_config import ClConfig

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

def run_fastqc(threads=4):
    """
    run_fastqc corre FastQC sobre los .fastq de cada grupo de expresión del organismo Anopheles gambiae
    y deja los reportes en results/fastqc/<grupo>.

    Args:
        threads (int) -> Número de hilos que usará FastQC

    Returns:
        fastqc_dirs (List) -> Lista de rutas de salida de FastQC, una por grupo
    """
    fastqc_dirs = []

    for grupo in ["resistent", "susceptible"]:
        entrada = Path("data/raw") / grupo / "sra_runs"
        salida = Path("results/fastqc") / grupo
        salida.mkdir(parents=True, exist_ok=True)

        fastqs = sorted(entrada.glob("*.fastq*"))
        if not fastqs:
            raise FileNotFoundError(f"No hay .fastq en {entrada}")

        subprocess.run(
            ["fastqc", "-t", str(threads), "-o", str(salida), *map(str, fastqs)],
            check=True,
        )
        fastqc_dirs.append(salida)

    return fastqc_dirs

def run_multiqc(fastqc_dirs):
    """
    run_multiqc se encargará de recorrer los directorios con los reportes de FastQC de cada grupo
    (resistente y susceptible) y generará un reporte de MultiQC por grupo, que se guardará en la
    carpeta results/qc.

    Args:
        fastqc_dirs (List) -> Lista de rutas de salida generadas por run_fastqc, una por grupo

    Returns:
         output_dirs (List) -> Lista de rutas de salida generadas por MultiQC, una por grupo

    """
    output_dirs = []

    for i in fastqc_dirs:
        output_dir = f"results/qc/{i.name}"
        multiqc.run(str(i), cfg=ClConfig(output_dir=output_dir), return_html=True)
        output_dirs.append(output_dir)

    return output_dirs

def read_tsv(output_dirs):
    """ 
    read_tsv se encarga de leer el tsv de FastQC, su única responsabilidad es parsear el archivo y devolver un DataFrame por
    muestra, tal y como vienen de MultiQC
    
    Args: 
         output_dirs (List) -> Lista de rutas de salida generadas por run_multiqc, una por grupo
    
    Returns:
         df_qc (DataFrame) -> DataFrame por muestra de MultiQC, parseado y por el orden del MultiQC
    """

    dfs = []

    for output_dir in output_dirs:
        tsv_path = Path(output_dir) /"multiqc_data"/"multiqc_fastqc.txt"
        df_grupo = pd.read_csv(tsv_path, sep="\t")
        df_grupo["grupo"] = Path(output_dir).name
        dfs.append(df_grupo)

    df_qc = pd.concat(dfs, ignore_index=True)

    return df_qc

def classify_sample(df_qc):
    """
    classify_sample se encarga de construir una estructura de diccionarios, en por cada muestra se asocia la calidad, dependendiendo
    principalmente de la presencia de adaptadores en la muestra. De ésta manera se clasifican las muestras sobre las que realmente se
    van a usar en el análisis
    
    Args: 
        df_qc (DataFrame) -> Estructura de datos creada por la función read_tsv(). DataFrame de la información de las lecturas por el 
                        MultiQC
    
    Returns:
        dict_qc (Diccionario) -> Diccionario en el que se le asocia la calidad de la muestra con el identificador de la lectura
    """
    dict_qc = {}

    for i, row in df_qc.iterrows():
        sample = row["Sample"]

        banderas = {col: row[col] for col in FASTQC_FLAG_COLUMNS}
        banderas["overall"] = banderas["adapter_content"]
        banderas["grupo"] = row["grupo"]

        dict_qc[sample] = banderas

    return dict_qc

def delete_fails(dict_qc):
    """ 
    delete_fails detecta los SRRs que son clasificados como fails, los elimina de los archivos que se usarán en el análisis.
    
    Args: 
        dict_qc (Diccionario) ->  Diccionario creado por la función classify_sample
                        #
    Returns:
        dict_qc_filtered (Diccionario) -> Diccionario que asocia la calidad de la muestra con el identificador de la lectura clasificada
                                          por la calidad de la lectura eliminando los "fails"
    """
    dict_qc_filtered = {
        sample: banderas
        for sample, banderas in dict_qc.items()
        if banderas["overall"] != "fail"
    }

    return dict_qc_filtered

def main():
    fastqc_dirs = run_fastqc()
    output_dirs = run_multiqc(fastqc_dirs)
    df_qc = read_tsv(output_dirs)
    dict_qc = classify_sample(df_qc)
    dict_qc_filtered = delete_fails(dict_qc)

    print(f"Muestras totales: {len(dict_qc)}")
    print(f"Muestras que pasan el filtro: {len(dict_qc_filtered)}")

    return dict_qc_filtered

if __name__ == "__main__":
    main()