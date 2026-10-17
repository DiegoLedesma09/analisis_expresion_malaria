import glob
import gzip
import os
import re

import matplotlib.pyplot as plt


def count_reads(fastq_path: str) -> int:
    """
    Cuenta las lecturas de un archivo FASTQ (cada lectura ocupa 4 líneas).
    Acepta archivos comprimidos (.gz).

    Args:
        fastq_path (str) -> Ruta del archivo FASTQ

    Returns:
        int -> Número de lecturas
    """
    opener = gzip.open if fastq_path.endswith(".gz") else open
    n_lines = 0
    with opener(fastq_path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            n_lines += chunk.count(b"\n")
    return n_lines // 4

def plot_reads_proportion(data_dir: str = "data/", outdir: str = "results/") -> None:
    """
    Cuenta las lecturas de las muestras resistentes y susceptibles y grafica
    la proporción de cada grupo. Guarda la figura en results/reads_proportion.png.

    En datos pareados cuenta solo el archivo _1, para no duplicar fragmentos.

    Args:
        data_dir (str) -> Directorio de datos (por defecto 'data/')
        outdir (str) -> Directorio de salida de la figura (por defecto 'results/')

    Returns:
        dict[str, int] -> Lecturas por fenotipo
    """
    counts = {}
    for phenotype in ("resistant", "susceptible"):
        runs_dir = os.path.join(data_dir, "raw", phenotype, "sra_runs")
        if not os.path.isdir(runs_dir):
            raise FileNotFoundError(f"No existe el directorio: {runs_dir}")

        fastqs = [f for f in glob.glob(os.path.join(runs_dir, "*.fastq*"))
                  if not re.search(r"_2\.fastq", f)]
        counts[phenotype] = sum(count_reads(f) for f in fastqs)

    total = sum(counts.values())
    if total == 0:
        raise ValueError("No se encontraron lecturas en los directorios de datos")

    os.makedirs(outdir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(5, 4))
    bars = ax.bar(counts.keys(), counts.values(), color=["#d95f02", "#1b9e77"])
    for bar, n in zip(bars, counts.values()):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{n / total:.1%}\n({n:,})", ha="center", va="bottom")

    ax.set_ylabel("Número de lecturas")
    ax.set_title("Proporción de lecturas por fenotipo")
    ax.set_ylim(0, max(counts.values()) * 1.2)

    output_path = os.path.join(outdir, "reads_proportion.png")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main():
    plot_reads_proportion(data_dir="data/", outdir="results/")

