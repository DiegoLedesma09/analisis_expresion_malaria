"""
En este módulo, nos encargaremos de la descarga de los datos. De igual manera, haremos una función para cada una de las siguientes cosas
con el fin de hacerlo más reproducible. Las instrucciones de cada uno de los módulos son los siguientes.

Funciones:
    - Una función que permita introducir el email y la api.key
    - Una función que obtenga toda la información del BioProject del proyecto (Metadatos) más importante.
    Guardar el resultado en un archivo denominado bioproject_metadata y guardarlo en docs/. 
    - Una función que obtenga todos los biosamples y SRA asociados al correspondiente Bioproject. Regresar dos listas.
    Una con los IDs de Biosamples asociados y otra con los SRA asociados al BioProject (Imprimir uno debajo de
    otro).
    - Dos funciones. De la misma manera guardar los mismos metadatos de los BioSamples asociados al BioProject y de los SRA.
     Guardarlos en los archivos correspondientes biosample_metadata y sra_metadata.
    - Una función que se encargue de buscar y descargar los archivos FASTA asociados a los SRA. Guardar las lecturas obtenidas en
    la carpeta de data/raw/resistent o data/raw/susceptible
    - Una función que guarde los metadatos del Genoma del organismo del BioProject, desde la plataforma de Genome.
    Guardar el resultado en un archivo denominado genome_metadata en docs/
    - Una función que se encargue de buscar y descargar el genoma del organismo del BioProject desde la plataforma
    de Genome. Guardar el resultado en un archivo denominado genome_{organismo}.fasta en la data/genome/
    - Una función argparse con los argumentos de -p --project (Una string que de el BioProject), -e --email (String),
    -o --outdir (path), -c --concat (Un booleano que indique sí quiere concatenar), -r --reference (Un booleano que indique sí se quiere el genoma o no),
""" 

import logging

from Bio import Entrez


def configurate_entrez(email:str, api_key:str):
    """
    Ésta función se encarga de configurar Entrez.

    Args:
        email (str): Email del usuario
        api_key (str): Llave de NCBI para agilizar descargas
    """

    if not email:
        raise ValueError("Se requiere de una dirección Email para NCBI Entrez")
    
    Entrez.email = email
    Entrez.api_key = api_key if api_key else None
    logging.info(f"Se ha configurado Entrez con el email {email}, API key: {'Entregada' if api_key else 'None'}")  # noqa: LOG015
    print(f"Se ha configurado Entrez con el email: {email}")