
<div align="left">
<img width="439" height="246" alt="Mosquito_cp" src="https://github.com/user-attachments/assets/80d26d3f-81de-4460-ae09-b20f52c08f29" />

<br>
<br>
<img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white" height="48" alt="Python" />
<img src="https://img.shields.io/badge/Conda-44A833?style=for-the-badge&logo=anaconda&logoColor=white" height="48" alt="Conda" />
<img src="https://img.shields.io/badge/WSL-0078D4?style=for-the-badge&logo=linux&logoColor=white" height="48" alt="WSL" />
<br>
<br>
</div>

# MOSQUITO DGE
## Expresión diferencial entre organismos de *Anopheles gambiae*, mosquito vector de la malaria, resistentes y susceptibles a insecticidas

Se analizarán los datos disponibles en la base de datos del NCBI para comparar los perfiles de expresión génica (transcriptomas) de mosquitos *Anopheles gambiae* resistentes a deltametrina frente a los susceptibles. El objetivo es identificar qué genes se expresan de manera diferencial y que, por lo tanto, podrían estar implicados en la resistencia.

Para conocer el planteamiento, las preguntas de investigación, la metodología y
los resultados, consulta el [reporte de proyecto](docs/reporte-proyecto.md).

## Funcionalidades

- Descarga lecturas (FASTQ), genomas y metadatos desde NCBI utilizando un BioProject.
- Clasifica automáticamente las muestras por fenotipo (resistente/susceptible).
- Evalúa la calidad de las lecturas crudas y recortadas (FastQC y MultiQC).
- Recorta adaptadores y filtra lecturas de baja calidad (Trim Galore).
- Construye índices de referencia y alinea lecturas paired-end al genoma (HISAT2).
- Ordena e indexa archivos de alineamiento y evalúa la tasa de éxito (samtools).
- Cuantifica la expresión génica generando una matriz de conteos (featureCounts).
- Ejecuta análisis de expresión diferencial con InMoose (DESeq2) o cálculo directo de Log2FC.
- Clasifica genes sobreexpresados y subexpresados entre condiciones.
- Realiza enriquecimiento funcional de Gene Ontology (GO) con pruebas estadísticas (goatools).
- Genera visualizaciones automáticas: Volcano plots, PCA, Heatmaps y gráficos de enriquecimiento GO.

## Estructura del Proyecto

```text
analisis_expresion_malaria/
├── data/                   # Almacenamiento de datos
│   ├── raw/                # Lecturas crudas descargadas del SRA
│   ├── curated/            # Lecturas limpias (post-trimming) listas para alinear
│   ├── genome/             # Genoma de referencia (FASTA) y anotación (GTF)
│   ├── index/              # Índices construidos por HISAT2
│   ├── aligned/            # Archivos BAM ordenados e indexados por muestra
│   ├── counts/             # Matrices de conteos generadas por featureCounts
│   └── go/                 # Archivos de referencia para GO (go-basic.obo, gene2go)
├── docs/                   # Documentación del proyecto
│   ├── docs_modules/       # Documentación detallada de cada script (.md)
│   ├── metadata/           # Metadatos descargados de NCBI (BioProject, SRA, etc.)
│   └── reporte-proyecto.md # Reporte final del proyecto
├── envs/                   # Archivos YAML para crear los entornos de Conda
├── logs/                   # Archivos .log generados por cada módulo
├── results/                # Resultados finales para el usuario
│   ├── exp_dif/            # Tablas de expresión diferencial (CSV, TXT)
│   ├── go/                 # Tablas de enriquecimiento funcional (TSV)
│   ├── fastqc/             # Reportes HTML de calidad (FastQC)
│   └── images/             # Gráficas finales 
├── src/                    # Código fuente (los 7 módulos del pipeline en Python)
├── tests/                  # Pruebas
├── CITATION.cff            
├── LICENSE
├── codemeta.json
└── README.md
```

## Requisitos

- Python 3.12 o posterior.
- Git.
- Entorno Conda (ver envs/malaria_env_linux.yml) o las siguientes herramientas instaladas:
  - FastQC y Trim Galore.
  - HISAT2 y samtools.
  - Subread (featureCounts).
- Conexión a internet para la descarga de datos

## Datos

<!--
AYUDA:
Indiquen la fuente, el comando o enlace de descarga y la carpeta de destino. No
repitan aquí la descripción completa de muestras, formatos, variables o criterios
de selección; enlacen el reporte. Si los datos no pueden publicarse, expliquen los
requisitos y el procedimiento autorizado para acceder a ellos.

EJEMPLO:
Los datos proceden de NCBI RefSeq. Para descargarlos, ejecuta:

    python scripts/download_data.py

Los archivos se guardarán en data/raw/. Los identificadores, versiones y criterios
de selección se documentan en docs/reporte-proyecto.md.
-->

Los datos proceden de *NCBI Bioproject*, asociado al paper de Bonizzoni M, Afrane Y, Dunn WA, Atieli FK, Zhou G, Zhong D, Li J, Githeko A, Yan G. Comparative transcriptome analyses of deltamethrin-resistant and -susceptible Anopheles gambiae mosquitoes from Kenya by RNA-Seq. PLoS One. 2012;7(9):e44607. doi: 10.1371/journal.pone.0044607. Epub 2012 Sep 7. PMID: 22970263; PMCID: PMC3436877. Con el identificador: *PRJNA170440*

Para su descraga este es el uso sugerido en la terminal:

```bash
    python src/download_data.py -p PJRNA170440 -e correo@ejemplo.com -a api_key -v
```
Los identificadores, versiones y criterios de selección se documentan en el
[reporte del proyecto](docs/reporte-proyecto.md).

Estos datos serán situados en las carpetas correspondientes:

```text
proyecto/
├── data/            # Datos de ejemplo o archivos pequeños
├──── raw/           # Datos crudos
├────── resistant/   # Datos de organismo resistente crudos
├────── susceptible/ # Datos de organismo susceptible crudos
```

## Uso

Para ejecutar el pipeline completo, se deben correr los módulos en el siguiente orden:

### 1. Descarga de datos

```bash
    python src/download_data.py -p PRJNA170440 -e correo@ejemplo.com -a api_key -v
```

Argumentos: -p (BioProject), -e (Correo NCBI), -a (API Key), -v (Verbose). Los metadatos se guardan en `docs/`, las lecturas en `data/raw/` y las gráficas en `results/`.

### 2. Reproducción de resultados completa 

Ejecuten los siguientes pasos en el orden indicado (desde la raíz del proyecto):

```bash
# 1. Descarga de datos
    python src/download_data.py -p PRJNA170440 -e correo@ejemplo.com -v

# 2. Control de calidad y Trimming
    python src/qc_trimming.py -w 8 -v

# 3. Alineamiento con HISAT2
    python src/align_reads.py -r data/genome/genome_Anopheles_gambiae.fasta -i data/index/hisat2/genome -gtf data/genome/genome_Anopheles_gambiae.gtf -w 8 -v

# 4. Conteo con featureCounts
    python src/count_features.py -gtf data/genome/genome_Anopheles_gambiae.gtf -v

# 5. Expresión diferencial (DESeq2 / Log2FC)
    python src/differential_expression.py -v

# 6. Anotación funcional (GO)
    python src/functional_annotation.py -v
```

## Documentación 

- El [reporte del proyecto](docs/reporte-proyecto.md). contiene el problema, las preguntas, la metodología y los resultados.
- La documentación técnica de cada módulo se encuentra en: [docs/docs_modules](docs/docs_modules).
- El formato para citar este repositorio se encuentra en: [CITATION](CITATION.cff).

## Equipo

- García Martínez Mónica — Colaborador
- Ledesma Gallegos Diego — Colaborador
- Muñoz Muñoz Luis Angel — Colaborador

## Citación 

Si utilizas este software, consulta [CITATION.cff](CITATION.cff) o la opción
**Cite this repository** de GitHub.

## Licencia

### **MIT License**

Consulta MIT License [LICENSE](LICENSE) para conocer los términos de uso.

## Agradecimientos

Proyecto desarrollado en la Licenciatura en Ciencias Genómicas, UNAM, como parte de la asignatura de Biooython.

Agradecemos, por su apoyo a las docentes: 
- Dra. Verónica Jiménez
- Dra. Leticia Vega
- M en C. Shirley Alquicira

Por su asesoría durante la realización del proyecto a: 
- Pablo Salazar 

Así como a: 
- Mi gata Daysi por su odio hacia mi persona. Me motiva el hecho de acabar con su paciencia para seguir adelante.

<!--
LISTA DE COMPROBACIÓN ANTES DE ENTREGAR:
- [ ] Se eliminaron o sustituyeron todos los textos entre corchetes.
- [ ] Los enlaces funcionan desde GitHub.
- [ ] Los comandos fueron probados en un ambiente limpio.
- [ ] El README permite instalar y ejecutar un ejemplo.
- [ ] No se incluyeron credenciales ni información sensible.
- [ ] El reporte contiene la explicación científica y evita duplicar el README.
- [ ] La versión, la licencia y la información de citación están actualizadas.
-->
