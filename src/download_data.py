<<<<<<< HEAD
<<<<<<< HEAD
"""
En este módulo, nos encargaremos de la descarga de los datos. De igual manera, haremos una función para cada una de las siguientes cosas
con el fin de hacerlo más reproducible. Las instrucciones de cada uno de los módulos son los siguientes.

Funciones:
Y    - Una función que permita introducir el email y la api.key
Y    - Una función que obtenga toda la información del BioProject del proyecto (Metadatos) más importante.
    Guardar el resultado en un archivo denominado bioproject_metadata y guardarlo en docs/. 
Y    - Una función que obtenga todos los biosamples y SRA asociados al correspondiente Bioproject. Regresar dos listas.
    Una con los IDs de Biosamples asociados y otra con los SRA asociados al BioProject (Imprimir uno debajo de
    otro).
Y    - Dos funciones. De la misma manera guardar los mismos metadatos de los BioSamples asociados al BioProject y de los SRA.
     Guardarlos en los archivos correspondientes biosample_metadata y sra_metadata.
    - Una función que se encargue de buscar y descargar los archivos FASTA asociados a los SRA. Guardar las lecturas obtenidas en
    la carpeta de data/raw/resistent o data/raw/susceptible
    - Una función que guarde los metadatos del Genoma del organismo del BioProject, usando el database assembly.
    Guardar el resultado en un archivo denominado genome_metadata en docs/
    - Una función que se encargue de buscar y descargar el genoma del organismo del BioProject desde la plataforma
    de Genome. Guardar el resultado en un archivo denominado genome_{organismo}.fasta en la data/genome/
    - Una función argparse con los argumentos de -p --project (Una string que de el BioProject), -e --email (String),
    -o --outdir (path), -c --concat (Un booleano que indique sí quiere concatenar), -r --reference (Un booleano que indique sí se quiere el genoma o no),
""" 

import logging
import os

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

def get_bioproject_uid(bioproject_id: str) -> str:
    """
    Función auxiliar para obtener el UID interno numérico de un BioProject.

    Args:
        bioproject_id (str) -> Identificador del Bioproject.

    Returns:
        id_list[0] (str) -> Identificador UID numérico del Bioproject.
    """
    handle = Entrez.esearch(db="bioproject", term=bioproject_id)
    record = Entrez.read(handle)
    handle.close()
    
    id_list = record.get("IdList", [])
    if not id_list:
        raise ValueError(f"No se encontró el BioProject: {bioproject_id}")
        
    return id_list[0]


def get_bioproject_metadata(bioproject_id: str, outdir = "docs"):
    """
    Obtiene la información y los metadatos más importantes de un BioProject
        y los guarda en docs/bioproject_metadata
    Args:
        bioproject_id (str): El identificador del BioProject (ej. 'PRJNAxxxxxx')
        outdir (str): Carpeta de destino (por defecto 'docs/')
    
    Returns:
        None
    """

    # Asegurar que la carpeta del destino exista
    os.makedirs(outdir, exist_ok=True)


    project_uid = get_bioproject_uid(bioproject_id)

    resumen_handle = Entrez.esummary(db="bioproject", id=project_uid)
    resumen_record = Entrez.read(resumen_handle)
    resumen_handle.close()

    # Guardar el resultado en un archivo denominado bioproject_metadata en docs/
    output_path = os.path.join(outdir, "bioproject_metadata")
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f"======Metadatos del BioProject: {bioproject_id}======\n")

        doc = resumen_record["DocumentSummary"][0]
        f.writelines(f"{key} \t {value}\n" for key, value in doc.items())

    print(f"Metadatos guardados correctamente en: {output_path}")

    return resumen_record

def get_biosamples_and_sra(bioproject_id: str) -> tuple[list[str], list[str]]:
    """
    Obtiene todos los Biosamples y SRA asociados al BioProject.
    Imprime los IDs uno debajo del otro y regresa dos listas.

    Args:
        bioproject_id (str) -> Identicador del BioProject

    Returns:
        biosample_ids, sra_ids (tuple[list[str], list[str]]) -> Tupla de listas:
            biosample_ids: Lista con los ids de los biosamples asociados al bioproject.
            sra_ids: Lista con los ids de los sras asociados al bioproject. 
    """

    # Se llama la funcion para obtener el UID interno del BioProject
    
    bioproject_uid = get_bioproject_uid(bioproject_id)
    
    # Buscar los BioSamples asociados
    
    bs_link_handle = Entrez.elink(dbfrom="bioproject", db="biosample", id=bioproject_uid)
    bs_link_record = Entrez.read(bs_link_handle)
    bs_link_handle.close()

    biosample_ids = []
    if bs_link_record[0].get("LinkSetDb"):
        for enlace in bs_link_record[0]["LinkSetDb"][0]["Link"]:
            biosample_ids.append(enlace["Id"])


    # Buscar los SRA vinculados usando elink
    sra_link_handle = Entrez.elink(dbfrom="bioproject", db="sra", id=bioproject_uid)
    sra_link_record = Entrez.read(sra_link_handle)
    sra_link_handle.close()

    sra_ids = []
    if sra_link_record[0].get("LinkSetDb"):
        for enlace in sra_link_record[0]["LinkSetDb"][0]["Link"]:
            sra_ids.append(enlace["Id"])

    # Imprimir uno debajo del otro
    print(f"\n======BioSamples asociados a {bioproject_id}======")
    for bs in biosample_ids:
        print(bs)

    print(f"\n======SRA asociados a {bioproject_id}======")
    for sra in sra_ids:
        print(sra)

    return biosample_ids, sra_ids

def save_biosample_metadata(biosample_ids: list[str], bioproject_id: str, outdir: str = "docs"):
    """
    Obtiene los metadatos de una lista de BioSamples y los guarda en 
    docs/biosample_metadata extrayendo únicamente los campos de interés.

    Args:
        biosample_ids (list[str]) -> Lista con los ids de los biosamples asociados al bioproject
        bioproject_id (str) -> Identificador del Bioproject
        outdir (str) -> Directorio de salida (por defecto '/docs')

    Returns:
        None
    """
    if not biosample_ids:
        print("No hay BioSamples para guardar")
        return

    os.makedirs(outdir, exist_ok=True)

    # Unir la lista de Ids con comas para descargarlos todos juntos
    ids_str = ",".join(biosample_ids)

    metadata_bs_handle = Entrez.esummary(db="biosample", id=ids_str)
    metadata_bs_record = Entrez.read(metadata_bs_handle)
    metadata_bs_handle.close()

    output_path= os.path.join(outdir, "biosample_metadata")
    with open(output_path, "w", encoding="utf-8") as f:

        f.write(f"====== Metadatos de BioSamples asociados a {bioproject_id} ======\n")

        documentos = metadata_bs_record.get("DocumentSummarySet",{}).get("DocumentSummary",[])

        for doc in documentos:
            f.write(f"\n--- BioSample ID: {doc.get('Id', 'Desconocido')} ---\n")
            f.write(f"Accession \t {doc.get('Accession', 'N/A')}\n")
            f.write(f"Id \t {doc.get('Id', 'N/A')}\n")
            f.write(f"Title \t {doc.get('Title', 'N/A')}\n")
            f.write(f"Organism \t {doc.get('Organism', 'N/A')}\n")
            f.write(f"ProjectAccession \t {doc.get('ProjectAccession', 'N/A')}\n")

    print(f"Metadatos de BioSamples guardados correctamente en: {output_path}")

def save_sra_metadata(sra_ids: list[str], bioproject_id: str, outdir: str = "docs"):
    """
    Obtiene los metadatos de una lista de SRA y los guarda en 
    docs/sra_metadata extrayendo únicamente los campos de interés.

    Args:
        sra_ids (list[str]) -> Lista con los ids de los sra asociados al bioproject
        bioproject_id (str) -> Identificador del Bioproject
        outdir (str) -> Directorio de salida (por defecto '/docs')

    Returns:
        None
    """
    if not sra_ids:
        print("No hay SRA para guardar.")
        return
    
    os.makedirs(outdir, exist_ok=True)

    ids_str_sra = ",".join(sra_ids)

    handle_sra = Entrez.esummary(db="sra", id=ids_str_sra)
    record_sra = Entrez.read(handle_sra)
    handle_sra.close()

    output_path = os.path.join(outdir, "sra_metadata")
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f"======Metadatos de SRA asociados al BioProject: {bioproject_id}======\n")

        documentos_sra = record_sra.get("DocumentSummary",[])

        for doc in documentos_sra:  
            f.write(f"\n--- SRA ID: {doc.get('Id', 'Desconocido')} ---\n")  
            # Campos directos
            f.write(f"Experiment Accession (SRX) \t {doc.get('ExtId', 'N/A')}\n")
            f.write(f"Id \t {doc.get('Id', 'N/A')}\n")
            f.write(f"Title \t {doc.get('Title', 'N/A')}\n")

            ## Extraccion del numero SRR(Run)
            runs_xml = doc.get('Runs', '')
            srr_accession = "N/A"
            if 'acc="' in runs_xml:
                partes = runs_xml.split('acc="')
                srr_accession = partes[1].split('"')[0]
            f.write(f"Run Accession (SRR) \t {srr_accession}\n")

    print(f"Metadatos de SRA guardados correctamente en: {output_path}")

def download_sra_lectures(bioproject_id:str, outdir:str = "data/"):
    """
    Descarga las lecturas correspondientes a los SRA asociados al Bioproject,
    las guarda dentro de nuevos directorios correspondientes dependiendo de 
    la condición del muestreo.

    Args:
        bioproject_id (str) -> Identificador del Bioproject
        outdir (str) -> Directorio de salida (por defecto '/data')
    """
    _ , sra_ids = get_biosamples_and_sra(bioproject_id)

    os.makedirs(f"{outdir}/raw/susceptible/sra_runs", exist_ok=True)
    os.makedirs(f"{outdir}/raw/resistent/sra_runs", exist_ok=True)

    for sra_id in sra_ids:
        handle_sra = Entrez.efetch(db="sra", id=sra_id)
        record = Entrez.read(handle_sra)
        print(record)
        handle_sra.close()

    def get_genome_metadata(bioproject_id: str, outdir: str = "docs/"):
        """
        Decarga los metadatos del genoma de referencia asociado al organismo del Bioproject desde la base de datos 'assembly' de 
        NCBI y los guarda en el archivo docs/genome_metadata.
        Args:
            bioproject_id (str) -> Identificador del BioProject. 
            outdir (str) -> Directorio de salida (por defecto '/data').
        """
        # Asegurar que la carpeta del destino exista
        os.makedirs(outdir, exist_ok=True)

        project_uid = get_bioproject_uid(bioproject_id)
        handle_bp = Entrez.esummary(db="bioproject", id=project_uid)
        record_bp = Entrez.read(handle_bp)
        handle_bp.close()

        doc_bp = record_bp["DocumentSummarySet"]["DocumentSummary"][0]
        organism = doc_bp.get("Organism_Name")
        if not organism:
        print(f"No se pudo determinar el organismo para {bioproject_id}")
        return

        print(f"Buscando genoma de referencia para: {organism}")

        handle_search = Entrez.esearch(db="assembly", term=f'{organism}[Organism] AND "reference genome"[filter]')
        record_search = Entrez.read(handle_search)
        handle_search.close()

        assembly_uids = record_search["IdList"]
        if not assembly_uids:
            print(f"No se encontró ensamblaje de referencia para {organism}")
            return
        
        assembly_uid = assembly_uids[0]

        handle_genome = Entrez.esummary(db="assembly", id=assembly_uid, report="full")
        record_genome = Entrez.read(handle_genome)
        handle_genome.close()

        output_path = os.path.join(outdir, "genome_metadata")
        
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(f"====== Metadatos del Genoma asociado al BioProject: {bioproject_id} ======\n\n")

            doc = record_genome["DocumentSummarySet"]["DocumentSummary"][0]
        campos_interes = [
            "AssemblyAccession",
            "AssemblyName",
            "AssemblyStatus",
            "Organism",
            "RefSeq_category",
            "FtpPath_RefSeq",
            "FtpPath_GenBank",
        ]
        for campo in campos_interes:
            f.write(f"{campo} \t {doc.get(campo, 'N/A')}\n")

    print(f"Metadatos del genoma guardados en: {output_path}")







        
        


=======
=======
>>>>>>> b04a8ea (Resolving the merge conflicts in download_data.py)
"""
download_data.py
================

Módulo para la descarga y consulta de datos de un BioProject de NCBI.

A partir de un identificador de BioProject (por ejemplo, PRJNA170440), obtiene
sus metadatos, los BioSamples y SRA asociados, descarga las lecturas
clasificadas por fenotipo (susceptible o resistente) y, opcionalmente,
consulta el genoma de referencia del organismo.

Requisitos:
    - Biopython (Bio.Entrez)
    - sra-tools (fasterq-dump disponible en el PATH). En Windows debe
      ejecutarse desde WSL (Ubuntu), ya que sra-tools no existe para win-64.
    - Conexión a internet y un correo electrónico para NCBI Entrez.

Funciones principales:
    configurate_entrez          Configura el email y la API key de Entrez.
    get_bioproject_uid          Obtiene el UID numérico interno de un BioProject.
    get_bioproject_metadata     Guarda los metadatos del BioProject.
    get_biosamples_and_sra      Regresa los IDs de BioSamples y SRA asociados.
    save_biosample_metadata     Guarda los metadatos de los BioSamples.
    save_sra_metadata           Guarda los metadatos de los SRA.
    download_sra_lectures       Descarga las lecturas con fasterq-dump y las
                                separa por fenotipo.
    get_bioproject_organism     Regresa el organismo asociado al Bioproject.
    get_genome_metadata         Guarda los metadatos del genoma de referencia
                                (base de datos assembly).
    get_reference_assembly      Busca el ensamble de referencia del organismo
                                del Bioproject.
    count_reads                 Cuenta las lecturas de los archivos.
    plot_reads_proportion       Cuenta las lecturas de las muestras resistentes y
                                susceptibles y grafica la proporción por grupo.
    download_reference_genome   Descarga del genoma de referencia.
    setup_logging               Configura el logging del módulo.            
    my_parser                   Define los argumentos de línea de comandos.
    main                        Ejecuta el flujo completo.

Archivos generados (rutas relativas al directorio donde se ejecuta):
    docs/bioproject_metadata
    docs/biosample_metadata
    docs/sra_metadata
    docs/genome_metadata
    data/raw/susceptible/sra_runs/     Lecturas de la cepa susceptible
    data/raw/resistent/sra_runs/       Lecturas de la cepa resistente
    results/reads_proportion

Argumentos:
    -p,  --project      BioProject ID (obligatorio).
    -e,  --email        Email para NCBI Entrez (obligatorio).
    -a,  --api_key      API key de NCBI para agilizar las consultas.
    -om, --outdirm      Directorio de salida de los metadatos (default: docs/).
    -od, --outdird      Directorio de salida de los datos (default: data/).
    --thr_fq            Workers para fasterq-dump (default: 8).
    -r,  --reference    Incluye el genoma de referencia (default).
         --no-reference Omite el genoma de referencia.
    -g,  --graphic      Incluye el gráfico comparativo de proporción de lecturas (default).
         --no-graphic   Omite el gráfico comparativo de proporción de lecturas.
    -v,  --verbose      Muestra mensajes DEBUG.
    --log_file          Guarda el archivo de log.

Uso:
    python src/download_data.py -p PRJNA170440 -e correo@ejemplo.com -v
    python src/download_data.py -p PRJNA170440 -e correo@ejemplo.com --no-reference --no-graphic -log_file

Nota:
    Ejecutar desde la raíz del proyecto, ya que las rutas de salida por defecto
    son relativas.
"""

import argparse as ap
import glob
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
import urllib.request
from urllib.error import HTTPError, URLError

import matplotlib.pyplot as plt
from Bio import Entrez

logger = logging.getLogger(__name__)

def configurate_entrez(email:str, api_key:str) -> None:
    """
    Ésta función se encarga de configurar Entrez.

    Args:
        email (str): Email del usuario
        api_key (str): Llave de NCBI para agilizar descargas

    Returns:
        None
    """

    if not email:
        raise ValueError("Se requiere de una dirección Email para NCBI Entrez")
    
    Entrez.email = email
    Entrez.api_key = api_key if api_key else None
    logger.info("Entrez configurado con el email: %s (API key: %s)", email, "entregada" if api_key else "no")

def get_bioproject_uid(bioproject_id: str) -> str:
    """
    Función auxiliar para obtener el UID interno numérico de un BioProject.

    Args:
        bioproject_id (str) -> Identificador del Bioproject.

    Returns:
        id_list[0] (str) -> Identificador UID numérico del Bioproject.
    """

    try:
        logger.debug("Consultando UID de %s", bioproject_id)
        with Entrez.esearch(db="bioproject", term=bioproject_id) as handle:
            record = Entrez.read(handle)
            id_list = record.get("IdList", [])
            if not id_list:
                raise ValueError(f"No se encontró el BioProject: {bioproject_id}")

        logger.debug("UID de %s: %s", bioproject_id, id_list[0])   
        return id_list[0]

    except HTTPError as e:
        raise ConnectionError(f"NCBI respondió con error {e.code} para {bioproject_id}") from e
    except (URLError, RuntimeError) as e:
        raise ConnectionError(f"No se pudo consultar NCBI para {bioproject_id}: {e}") from e

def get_bioproject_metadata(bioproject_id: str, outdir = "docs") -> None:
    """
    Obtiene la información y los metadatos más importantes de un BioProject
        y los guarda en docs/bioproject_metadata
    Args:
        bioproject_id (str): El identificador del BioProject (ej. 'PRJNAxxxxxx')
        outdir (str): Carpeta de destino (por defecto 'docs/')
    
    Returns:
        None
    """

    # Asegurar que la carpeta del destino exista
    os.makedirs(outdir, exist_ok=True)

    project_uid = get_bioproject_uid(bioproject_id)

    resumen_handle = Entrez.esummary(db="bioproject", id=project_uid)
    resumen_record = Entrez.read(resumen_handle)
    resumen_handle.close()

    # Guardar el resultado en un archivo denominado bioproject_metadata en docs/
    output_path = os.path.join(outdir, "bioproject_metadata")
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f"======Metadatos del BioProject: {bioproject_id}======\n")

        doc = resumen_record["DocumentSummarySet"]["DocumentSummary"][0]
        f.writelines(f"{key} \t {value}\n" for key, value in doc.items())

    logger.info("... guardados en: %s", output_path)

def get_biosamples_and_sra(bioproject_id: str) -> tuple[list[str], list[str]]:
    """
    Obtiene todos los Biosamples y SRA asociados al BioProject.
    Imprime los IDs uno debajo del otro y regresa dos listas.

    Args:
        bioproject_id (str) -> Identicador del BioProject

    Returns:
        biosample_ids, sra_ids (tuple[list[str], list[str]]) -> Tupla de listas:
            biosample_ids: Lista con los ids de los biosamples asociados al bioproject.
            sra_ids: Lista con los ids de los sras asociados al bioproject. 
    """

    # Se llama la funcion para obtener el UID interno del BioProject
    
    bioproject_uid = get_bioproject_uid(bioproject_id)
    
    # Buscar los BioSamples asociados
    
    bs_link_handle = Entrez.elink(dbfrom="bioproject", db="biosample", id=bioproject_uid)
    bs_link_record = Entrez.read(bs_link_handle)
    bs_link_handle.close()

    biosample_ids = []
    if bs_link_record[0].get("LinkSetDb"):
        for enlace in bs_link_record[0]["LinkSetDb"][0]["Link"]:
            biosample_ids.append(enlace["Id"])


    # Buscar los SRA vinculados usando elink
    sra_link_handle = Entrez.elink(dbfrom="bioproject", db="sra", id=bioproject_uid)
    sra_link_record = Entrez.read(sra_link_handle)
    sra_link_handle.close()

    sra_ids = []
    if sra_link_record[0].get("LinkSetDb"):
        for enlace in sra_link_record[0]["LinkSetDb"][0]["Link"]:
            sra_ids.append(enlace["Id"])

    # Imprimir uno debajo del otro
    print(f"\n======BioSamples asociados a {bioproject_id}======")
    for bs in biosample_ids:
        logger.info("%s", bs)

    print(f"\n======SRA asociados a {bioproject_id}======")
    for sra in sra_ids:
        logger.info("%s", sra)

    logger.info("%d BioSamples y %d SRA asociados a %s", len(biosample_ids), len(sra_ids), bioproject_id)

    return biosample_ids, sra_ids

def save_biosample_metadata(biosample_ids: list[str], bioproject_id: str, outdir: str = "docs") -> None:
    """
    Obtiene los metadatos de una lista de BioSamples y los guarda en 
    docs/biosample_metadata extrayendo únicamente los campos de interés.

    Args:
        biosample_ids (list[str]) -> Lista con los ids de los biosamples asociados al bioproject
        bioproject_id (str) -> Identificador del Bioproject
        outdir (str) -> Directorio de salida (por defecto '/docs')

    Returns:
        None
    """
    if not biosample_ids:
        logger.warning("No hay BioSamples para guardar")
        return

    os.makedirs(outdir, exist_ok=True)

    # Unir la lista de Ids con comas para descargarlos todos juntos
    ids_str = ",".join(biosample_ids)

    metadata_bs_handle = Entrez.esummary(db="biosample", id=ids_str)
    metadata_bs_record = Entrez.read(metadata_bs_handle)
    metadata_bs_handle.close()

    output_path= os.path.join(outdir, "biosample_metadata")
    with open(output_path, "w", encoding="utf-8") as f:

        f.write(f"====== Metadatos de BioSamples asociados a {bioproject_id} ======\n")

        documentos = metadata_bs_record.get("DocumentSummarySet",{}).get("DocumentSummary",[])

        for doc in documentos:
            f.write(f"\n--- BioSample ID: {doc.get('Id', 'Desconocido')} ---\n")
            f.write(f"Accession \t {doc.get('Accession', 'N/A')}\n")
            f.write(f"Id \t {doc.get('Id', 'N/A')}\n")
            f.write(f"Title \t {doc.get('Title', 'N/A')}\n")
            f.write(f"Organism \t {doc.get('Organism', 'N/A')}\n")
            f.write(f"ProjectAccession \t {doc.get('ProjectAccession', 'N/A')}\n")

    logger.info("... guardados en: %s", output_path)

def save_sra_metadata(sra_ids: list[str], bioproject_id: str, outdir: str = "docs") -> None:
    """
    Obtiene los metadatos de una lista de SRA y los guarda en 
    docs/sra_metadata extrayendo únicamente los campos de interés.

    Args:
        sra_ids (list[str]) -> Lista con los ids de los sra asociados al bioproject
        bioproject_id (str) -> Identificador del Bioproject
        outdir (str) -> Directorio de salida (por defecto '/docs')

    Returns:
        None
    """
    if not sra_ids:
        logger.warning("No hay SRA para guardar.")
        return
    
    os.makedirs(outdir, exist_ok=True)

    ids_str_sra = ",".join(sra_ids)

    handle_sra = Entrez.esummary(db="sra", id=ids_str_sra)
    record_sra = Entrez.read(handle_sra)
    handle_sra.close()

    output_path = os.path.join(outdir, "sra_metadata")
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f"======Metadatos de SRA asociados al BioProject: {bioproject_id}======\n")

        for doc in record_sra:  
            f.write(f"\n--- SRA ID: {doc.get('Id', 'Desconocido')} ---\n")  
            # Campos directos
            f.write(f"Experiment Accession (SRX) \t {doc.get('ExtId', 'N/A')}\n")
            f.write(f"Id \t {doc.get('Id', 'N/A')}\n")
            f.write(f"Title \t {doc.get('Title', 'N/A')}\n")

            ## Extraccion del numero SRR(Run)
            runs_xml = doc.get('Runs', '')
            srr_accession = "N/A"
            if 'acc="' in runs_xml:
                partes = runs_xml.split('acc="')
                srr_accession = partes[1].split('"')[0]
            f.write(f"Run Accession (SRR) \t {srr_accession}\n")

    logger.info("... guardados en: %s", output_path)

def download_sra_lectures(bioproject_id:str, outdir:str = "data/", w:int = 8) -> None:
    """
    Descarga las lecturas correspondientes a los SRA asociados al Bioproject,
    las guarda dentro de nuevos directorios correspondientes dependiendo de 
    la condición del muestreo.

    Args:
        bioproject_id (str) -> Identificador del Bioproject
        outdir (str) -> Directorio de salida (por defecto '/data')
        w (int) -> Número de workers para la descarga de datos (por defecto: 8)
    
    Returns:
        None
    """
    _ , sra_ids = get_biosamples_and_sra(bioproject_id)

    # Creo los directorios de salida
    susceptible_dir = os.path.join(outdir, "raw", "susceptible", "sra_runs")
    resistent_dir = os.path.join(outdir, "raw", "resistent", "sra_runs")

    # Creo los directorios de salida en caso de que no existan
    os.makedirs(susceptible_dir, exist_ok=True)
    os.makedirs(resistent_dir, exist_ok=True)

    # Busco en Sra el sra_id buscando sí es la cepa Susceptible o Resistente al insecticida
    for i, sra_id in enumerate(sra_ids, start=1):
        handle_sra = Entrez.esummary(db="sra", id=sra_id)
        record = Entrez.read(handle_sra)
        exp_xml = record[0]["ExpXml"]
        handle_sra.close()

        # Mediante una función regex, busco el campo que me diga el fenotipo
        match = re.search(r':\s*([a-zA-Z]+);', exp_xml)

        # Si encuentra un fenotipo asociado, lo guarda
        if match:
            fenotipo = match.group(1)
        else:
            raise ValueError(f"No se encontró fenotipo asociado al {sra_id}")
        
        # Le asocio el nombre del fenotipo al directorio de salida
        file_domain = fenotipo.lower()
        if "susceptible" in file_domain:
            domain_dir = susceptible_dir
        elif "resistent" in file_domain:
            domain_dir = resistent_dir
        else:
            domain_dir = os.path.join(outdir, "raw", file_domain, "sra_runs")
            os.makedirs(domain_dir, exist_ok=True)
        
        logger.info("[%d/%d] Descargando %s (%s) hacia %s", i, len(sra_ids), sra_id, fenotipo, domain_dir)
        
        # Mediante subprocess, uso la libreria de fasterq-dump para descargar las lecturas
        # desde SRA, separando las dos lecturas pareadas y uso 8 workers
        cmd = ["fasterq-dump", 
            sra_id, 
            "-O", domain_dir, 
            "--split-files", 
            "-e", str(w)
            ]
        
        try:
            subprocess.run(cmd, check=True)
            logger.info("Descarga exitosa: %s", sra_id)
        except subprocess.CalledProcessError as e:
            raise ValueError(f"Error al descargar {sra_id} con fasterq-dump {e}")

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
    
>>>>>>> 7ba079d (Updating download_sra_runs functions in download_data.py module)
=======
=======
def get_genome_metadata(bioproject_id: str, outdir: str = "docs/"):
=======
def get_bioproject_organism(bioproject_id:str) -> str:
    """
    Encuentra al organismo muestra del Bioproject asociado

    Args:
        bioproject_id (str) -> Identificador del Bioproject
    
    Returns:
        organism (str) -> Organismo asociado al Bioproject
    """
    try:
        project_uid = get_bioproject_uid(bioproject_id)
        with Entrez.esummary(db="bioproject", id=project_uid) as handle_bp:
            record_bp = Entrez.read(handle_bp)

            doc_bp = record_bp["DocumentSummarySet"]["DocumentSummary"][0]
            organism = doc_bp.get("Organism_Name")

            if not organism:
                raise ValueError(f"No se encontró organismo asociado al Bioproject {bioproject_id}")

            return organism

    except HTTPError as e:
        raise ConnectionError(f"NCBI respondió con error {e.code} para {bioproject_id}") from e
    except (URLError, RuntimeError) as e:
        raise ConnectionError(f"No se pudo consultar NCBI para {bioproject_id}: {e}") from e

def get_genome_metadata(bioproject_id:str, outdir:str="docs/"):
>>>>>>> e6c0d04 (Adding functions and documentation at download_data.py module)
        """
        Decarga los metadatos del genoma de referencia asociado al 
        organismo del Bioproject desde la base de datos 'assembly' de 
        NCBI y los guarda en el archivo docs/genome_metadata.

        Args:
            bioproject_id (str) -> Identificador del BioProject. 
            outdir (str) -> Directorio de salida (por defecto '/docs').
        """
        # Asegurar que la carpeta del destino exista
        os.makedirs(outdir, exist_ok=True)

        organism = get_bioproject_organism(bioproject_id)

        logger.info("Buscando genoma de referencia para: %s", organism)

        handle_search = Entrez.esearch(db="assembly", term=f'{organism}[Organism] AND "reference genome"[filter]')
        record_search = Entrez.read(handle_search)
        handle_search.close()

        assembly_uids = record_search["IdList"]
        if not assembly_uids:
            logger.warning("No se encontró ensamblaje de referencia para %s", organism)
            return
        
        assembly_uid = assembly_uids[0]

        handle_genome = Entrez.esummary(db="assembly", id=assembly_uid, report="full")
        record_genome = Entrez.read(handle_genome)
        handle_genome.close()

        output_path = os.path.join(outdir, "genome_metadata")
        
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(f"====== Metadatos del Genoma asociado al BioProject: {bioproject_id} ======\n\n")

            doc = record_genome["DocumentSummarySet"]["DocumentSummary"][0]

            campos_interes = [
                "AssemblyAccession",
                "AssemblyName",
                "AssemblyStatus",
                "Organism",
                "RefSeq_category",
                "FtpPath_RefSeq",
                "FtpPath_GenBank",
            ]
            f.writelines(f"{campo} \t {doc.get(campo, 'N/A')}\n" for campo in campos_interes)

<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> b04a8ea (Resolving the merge conflicts in download_data.py)
def download_reference_genome(bioproject_id:str, outdir:str = "data/", w:int = 8) -> None:
=======
            print(f"Metadatos del genoma guardados en: {output_path}")
=======
            logger.info("Metadatos del genoma guardados en: %s", output_path)
>>>>>>> 7866ad0 (Adding logging funtions at download_data.py module to see mistakes)

def get_reference_assembly(organism: str) -> dict:
>>>>>>> e6c0d04 (Adding functions and documentation at download_data.py module)
    """
    Busca en la base de datos 'assembly' el ensamblaje de referencia de un
    organismo y regresa su documento de metadatos.

    Args:
        organism (str) -> Nombre del organismo

    Returns:
        dict -> Documento del ensamblaje (incluye campos como FtpPath_RefSeq)
    """
    try:
        term = f'{organism}[Organism] AND "reference genome"[filter]'
        with Entrez.esearch(db="assembly", term=term) as handle:
            record_search = Entrez.read(handle)

        assembly_uids = record_search["IdList"]
        if not assembly_uids:
            raise ValueError(f"No se encontró ensamblaje de referencia para {organism}")

        with Entrez.esummary(db="assembly", id=assembly_uids[0], report="full") as handle:
            record_genome = Entrez.read(handle)

        return record_genome["DocumentSummarySet"]["DocumentSummary"][0]

    except HTTPError as e:
        raise ConnectionError(f"NCBI respondió con error {e.code} para {organism}") from e
    except (URLError, RuntimeError) as e:
        raise ConnectionError(f"No se pudo consultar NCBI para {organism}: {e}") from e


def download_reference_genome(bioproject_id: str, outdir: str = "data/") -> None:
    """
    Descarga el genoma de referencia del organismo del BioProject y lo guarda
    descomprimido en data/genome/genome_{organismo}.fasta.

    Args:
        bioproject_id (str) -> Identificador del Bioproject
        outdir (str) -> Directorio de salida (por defecto 'data/')

    Returns:
        None
    """
    organism = get_bioproject_organism(bioproject_id)
    assembly = get_reference_assembly(organism)

    # Prefiere RefSeq y usa GenBank como respaldo
    ftp_path = assembly.get("FtpPath_RefSeq") or assembly.get("FtpPath_GenBank")
    if not ftp_path:
        raise ValueError(f"El ensamblaje de {organism} no tiene ruta de descarga")

    # La ruta apunta a una carpeta; el archivo repite su último segmento
    base_url = ftp_path.replace("ftp://", "https://", 1).rstrip("/")
    file_name = f"{base_url.split('/')[-1]}_genomic.fna.gz"
    url = f"{base_url}/{file_name}"

    genome_dir = os.path.join(outdir, "genome")
    os.makedirs(genome_dir, exist_ok=True)

    safe_name = re.sub(r"\W+", "_", organism).strip("_")
    fasta_path = os.path.join(genome_dir, f"genome_{safe_name}.fasta")

    if os.path.exists(fasta_path):
        logger.info("El genoma ya existe, se omite la descarga: %s", fasta_path)
        return

    # Archivos temporales: solo se renombran al terminar bien
    gz_tmp = fasta_path + ".gz.part"
    fasta_tmp = fasta_path + ".part"

    (f"Descargando genoma de {organism} desde {url}")
    try:
        urllib.request.urlretrieve(url, gz_tmp)

        with gzip.open(gz_tmp, "rb") as f_in, open(fasta_tmp, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        os.replace(fasta_tmp, fasta_path)
        print(f"Genoma guardado en: {fasta_path}")

    except HTTPError as e:
        raise ConnectionError(f"No se pudo descargar {url}: error {e.code}") from e
    except URLError as e:
        raise ConnectionError(f"No se pudo descargar {url}: {e}") from e
    except (gzip.BadGzipFile, EOFError) as e:
        raise ValueError(f"El archivo descargado está corrupto o incompleto: {e}") from e
    finally:
        for tmp in (gz_tmp, fasta_tmp):
            if os.path.exists(tmp):
                os.remove(tmp)

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
    for phenotype in ("resistent", "susceptible"):
        runs_dir = os.path.join(data_dir, "raw", phenotype, "sra_runs")
        if not os.path.isdir(runs_dir):
            raise FileNotFoundError(f"No existe el directorio: {runs_dir}")

        fastqs = [f for f in glob.glob(os.path.join(runs_dir, "*.fastq*"))
                  if not re.search(r"_2\.fastq", f)]
        counts[phenotype] = sum(count_reads(f) for f in fastqs)
        logger.info("%s: %d lecturas", phenotype, counts[phenotype])

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

    logger.info("Figura guardada en: %s", output_path)

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
    description="Pipeline para la descarga y consulta de lecturas de SRA y genoma dereferencia asocidos a un BioProject")
    parser.add_argument("-p", "--project", required=True, help="BioProject ID (e.g., PRJNA877658)")
    parser.add_argument("-e", "--email", required=True, help="Email para NCBI Entrez")
    parser.add_argument("-a", "--api_key", required=False, help="API key para NCBI Entrez")
    parser.add_argument("-om", "--outdirm", required=False, default="docs/", help="Directorio de salida para metadatos")
    parser.add_argument("-od", "--outdird", required=False, default="data/", help="Directorio de salida para datos")
    parser.add_argument("--thr_fq", required=False, type=int, default=8, help="Workers para faster-dump")
<<<<<<< HEAD
    parser.add_argument("-r", "--reference", required=False, default=True, help="Descarga del genoma de referencia")
<<<<<<< HEAD
    return parser.parse_args()
>>>>>>> 1d019bd (Adding funtions structure download_genome and improving documentation in download_data.py module)
=======
    return parser.parse_args()
>>>>>>> b04a8ea (Resolving the merge conflicts in download_data.py)
=======
    parser.add_argument("-r", "--reference", action=ap.BooleanOptionalAction, default=True, 
                        help="Descarga del genoma de referencia (usa --no-reference para omitirlo)")
    parser.add_argument("-g", "--graphic", action=ap.BooleanOptionalAction, default=True,
                        help="Incluye el grafico comparativo de proporcion de lecturas (usa --no-graphic para omitirlo)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Muestra mensajes DEBUG")
    parser.add_argument("--log_file", default="logs/download_data.log", help="Archivo de log")
    return parser.parse_args()

def main() -> None:
    """
    Ejecuta el flujo completo de descarga para un BioProject.

    Pasos, en orden:
        1. Configura Entrez con el email y la API key.
        2. Guarda los metadatos del BioProject.
        3. Obtiene los BioSamples y SRA asociados, y guarda sus metadatos.
        4. Descarga las lecturas de los SRA, separadas por fenotipo.
        5. Si --reference está activo, guarda los metadatos del genoma de
           referencia y lo descarga.
        6. Si --graphic está activo, crea el grafico comparativo de proporcion
            de lecturas

    Los argumentos se leen de la línea de comandos mediante my_parser().

    Returns:
        None
    """
    args = my_parser()
    setup_logging(args.verbose, args.log_file)
    logger.info("Inicio del pipeline para el BioProject %s", args.project)

    try:
        configurate_entrez(args.email, args.api_key)

        # Metadatos del BioProject
        get_bioproject_metadata(args.project, outdir=args.outdirm)

        # BioSamples y SRA (se consultan una vez y se reutilizan)
        biosample_ids, sra_ids = get_biosamples_and_sra(args.project)
        save_biosample_metadata(biosample_ids, args.project, outdir=args.outdirm)
        save_sra_metadata(sra_ids, args.project, outdir=args.outdirm)

        # Lecturas
        download_sra_lectures(args.project, outdir=args.outdird, w=args.thr_fq)

        # Genoma de referencia (opcional)
        if args.reference:
            get_genome_metadata(args.project, outdir=args.outdirm)
            download_reference_genome(args.project, outdir=args.outdird)
        
        # Grafico comparativo (opcional)
        if args.graphic:
            plot_reads_proportion(data_dir="data/", outdir="results/")

        logger.info("Pipeline terminado correctamente")

    except (ConnectionError, ValueError) as e:
        logger.error("%s", e)
        logger.debug("Detalle del error", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
>>>>>>> e6c0d04 (Adding functions and documentation at download_data.py module)
