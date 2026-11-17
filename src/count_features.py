import argparse as ap
import logging
import os
import sys

logger = logging.getLogger(__name__)

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
    parser = ap.ArgumentParser(description="Pipeline para la creación de la matriz de conteos con featureCounts")
    parser.add_argument("-o", "--outdir", required=False, default="data/aligned/", 
                        help="Directorio de salida")
    parser.add_argument("-w", "--workers", type=int, required= False, default= 6, 
                        help="Número de hilos")
    parser.add_argument("-v", "--verbose", action="store_true", 
                        help="Muestra mensajes DEBUG en consola")
    parser.add_argument("--log_file", default="logs/align_reads.log", 
                        help="Ruta para guardar el archivo de log")
    
    return parser.parse_args()

def main():
    """
    Ejecuta el flujo de alineamiento de lecturas RNA-seq contra el genoma de referencia
    con HISAT2 y procesa los alineamientos con samtools.
    Pasos, en orden:
        1. Construye el índice de HISAT2 a partir del genoma de referencia (omite la
           construcción si el índice ya existe). Si se proporciona un GTF, extrae los
           sitios de empalme y exones para optimizar el alineamiento.
        2. Alinea con HISAT2 cada par de lecturas curadas (_val_1.fq, _val_2.fq) y
           guarda el resumen de alineamiento de cada muestra.
        3. Convierte cada archivo SAM a BAM ordenado e indexado (.bai) y elimina el SAM
           temporal.
        4. Lee los resúmenes de HISAT2 y registra una advertencia por cada muestra cuya
           tasa de alineamiento esté por debajo del umbral.
    """
    args = my_parser()
    setup_logging(args.verbose, args.log_file)
    logger.info("Inicializando módulo de alineamiento de lecturas")
