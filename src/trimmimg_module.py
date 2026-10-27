import subprocess
from pathlib import Path

from multiqc_module import main as run_qc

FASTQ_SUFFIX_FWD = "_1.fastq"
FASTQ_SUFFIX_REV = "_2.fastq"

# Archivo de adaptadores que usa Trimmomatic (verifica que corresponda a tu plataforma)
ADAPTERS = "TruSeq3-PE-2.fa"

def pair_samples(dict_qc_filtered):
    """
    pair_samples agrupa las lecturas _1 y _2 por SRR y se queda solo con los pares completos,
    es decir, aquellos en los que ambas lecturas pasaron el control de calidad.

    Args:
        dict_qc_filtered (Diccionario) -> Diccionario creado por delete_fails(), con "grupo" incluido por muestra

    Returns:
        pares (Diccionario) -> {SRR: {"grupo": grupo}}, el formato que espera run_trimmomatic
    """
    mates = {}

    for sample, info in dict_qc_filtered.items():
        srr, _, mate = sample.rpartition("_")
        entrada = mates.setdefault(srr, {"grupo": info["grupo"], "mates": set()})
        entrada["mates"].add(mate)

    return {
        srr: {"grupo": d["grupo"]}
        for srr, d in mates.items()
        if d["mates"] == {"1", "2"}
    }


def run_trimmomatic(pares, raw_dir="data/raw", output_dir="results/trimmed", threads=4):
    """
    run_trimmomatic se encarga de correr Trimmomatic (modo PE) sobre cada par de lecturas que pasó el control de
    calidad, usando el campo "grupo" para localizar los fastq crudos. Por cada SRR se generan cuatro archivos:
    paired y unpaired para cada lectura (las unpaired se conservan por si se necesitan después).

    Args:
        pares (Diccionario) -> Diccionario creado por pair_samples(), {SRR: {"grupo": grupo}}
        raw_dir (str) -> Ruta base donde están los fastq crudos, organizados por grupo/sra_runs
        output_dir (str) -> Ruta base donde Trimmomatic escribirá las lecturas ya limpias
        threads (int) -> Número de hilos que usará Trimmomatic

    Returns:
        output_paths (Diccionario) -> Diccionario que asocia cada SRR con la carpeta de salida de su trimming
    """
    output_paths = {}

    for srr, info in pares.items():
        grupo = info["grupo"]

        fwd = f"{raw_dir}/{grupo}/sra_runs/{srr}{FASTQ_SUFFIX_FWD}"
        rev = f"{raw_dir}/{grupo}/sra_runs/{srr}{FASTQ_SUFFIX_REV}"

        sample_output_dir = Path(output_dir) / grupo
        sample_output_dir.mkdir(parents=True, exist_ok=True)

        fwd_paired = sample_output_dir / f"{srr}_1_paired.fastq"
        fwd_unpaired = sample_output_dir / f"{srr}_1_unpaired.fastq"
        rev_paired = sample_output_dir / f"{srr}_2_paired.fastq"
        rev_unpaired = sample_output_dir / f"{srr}_2_unpaired.fastq"

        subprocess.run([
            "trimmomatic", "PE",
            "-threads", str(threads),
            "-phred33",
            fwd, rev,
            str(fwd_paired), str(fwd_unpaired),
            str(rev_paired), str(rev_unpaired),
            f"ILLUMINACLIP:{ADAPTERS}:2:30:10",
            "LEADING:3",
            "TRAILING:3",
            "SLIDINGWINDOW:4:15",
            "MINLEN:36",
        ], check=True)

        output_paths[srr] = str(sample_output_dir)

    return output_paths


def main():
    dict_qc_filtered = run_qc()
    pares = pair_samples(dict_qc_filtered)

    print(f"Pares completos que pasan el filtro: {len(pares)}")

    output_paths = run_trimmomatic(pares)
    return output_paths


if __name__ == "__main__":
    main()