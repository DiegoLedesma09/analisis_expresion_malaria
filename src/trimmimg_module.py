import subprocess

FASTQ_SUFFIX_FWD = "_1.fastq.gz"
FASTQ_SUFFIX_REV = "_2.fastq.gz"


def run_trim_galore(dict_qc_filtered, raw_dir="../data/raw", output_dir="../results/trimmed"):
    # ===============================================================================================================================
    # run_trim_galore se encarga de correr Trim Galore sobre cada muestra que pasó el control de calidad (dict_qc_filtered),
    # usando el campo "grupo" de cada muestra para localizar el fastq crudo correspondiente. Las lecturas "unpaired" se
    # conservan (--retain_unpaired) por si se necesitan después.
    #
    # Args:
    #     dict_qc_filtered (Diccionario) -> Diccionario creado por delete_fails(), con "grupo" incluido por muestra
    #     raw_dir (str) -> Ruta base donde están los fastq crudos, organizados en subcarpetas por grupo
    #     output_dir (str) -> Ruta base donde Trim Galore escribirá las lecturas ya limpias
    #
    # Returns:
    #     output_paths (Diccionario) -> Diccionario que asocia cada muestra con la ruta de salida usada para su trimming
    # ===============================================================================================================================
    output_paths = {}

    for sample, info in dict_qc_filtered.items():
        grupo = info["grupo"]

        fwd = f"{raw_dir}/{grupo}/{sample}{FASTQ_SUFFIX_FWD}"
        rev = f"{raw_dir}/{grupo}/{sample}{FASTQ_SUFFIX_REV}"
        sample_output_dir = f"{output_dir}/{grupo}"

        subprocess.run([
            "trim_galore",
            "--paired",
            "--retain_unpaired",
            "--output_dir", sample_output_dir,
            fwd, rev,
        ], check=True)

        output_paths[sample] = sample_output_dir

    return output_paths