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
    - Una función que guarde los metadatos del Genoma del organismo del BioProject, desde la plataforma de Genome.
    Guardar el resultado en un archivo denominado genome_metadata en docs/
    - Una función que se encargue de buscar y descargar el genoma del organismo del BioProject desde la plataforma
    de Genome. Guardar el resultado en un archivo denominado genome_{organismo}.fasta en la data/genome/
    - Una función argparse con los argumentos de -p --project (Una string que de el BioProject), -e --email (String),
    -o --outdir (path), -c --concat (Un booleano que indique sí quiere concatenar), -r --reference (Un booleano que indique sí se quiere el genoma o no),
""" 

import os
import re
import subprocess

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

def download_sra_lectures(bioproject_id:str, outdir:str = "data/", w:int = "8"):
    """
    Descarga las lecturas correspondientes a los SRA asociados al Bioproject,
    las guarda dentro de nuevos directorios correspondientes dependiendo de 
    la condición del muestreo.

    Args:
        bioproject_id (str) -> Identificador del Bioproject
        outdir (str) -> Directorio de salida (por defecto '/data')
        w (int) -> Número de workers para la descarga de datos (por defecto: 8)
    """
    _ , sra_ids = get_biosamples_and_sra(bioproject_id)


    # Creo los directorios de salida
    susceptible_dir = os.path.join(outdir, "raw", "susceptible", "sra_runs")
    resistent_dir = os.path.join(outdir, "raw", "resistent", "sra_runs")

    # Creo los directorios de salida en caso de que no existan
    os.makedirs(susceptible_dir, exist_ok=True)
    os.makedirs(resistent_dir, exist_ok=True)

    # Busco en Sra el sra_id buscando sí es la cepa Susceptible o Resistente al insecticida
    for sra_id in sra_ids:
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
            domain_dir = os.path.join("raw", file_domain, "sra_runs")
            os.makedirs(domain_dir, exist_ok=True)
        
        print(f"Descargando la corrida {sra_id} ({fenotipo}, hacia {domain_dir})")
        
        # Mediante subprocess, uso la libreria de fasterq-dump para descargar las lecturas
        # desde SRA, separando las dos lecturas pareadas y uso 8 workers
        cmd = ["fasterq-dump", 
            sra_id, 
            "-O", domain_dir, 
            "--split-files", 
            "-e", w
            ]
        
        try:
            subprocess.run(cmd, check=True)
            print(f"Descarga exitosa {sra_id}")
        except subprocess.CalledProcessError as e:
            raise ValueError(f"Error al descargar {sra_id} con fasterq-dump {e}")

    