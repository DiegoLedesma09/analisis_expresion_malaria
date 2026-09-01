# Expresión diferencial entre organismos de Anopheles gambiae, mosquito vector de la malaria, resistentes y susceptibles a insecticidas

Se analizarán los datos disponibles en la base de datos del NCBI para comparar los perfiles de expresión génica (transcriptomas) de mosquitos *Anopheles gambiae* resistentes a insecticidas frente a los susceptibles. El objetivo es identificar qué genes se expresan de manera diferencial y que, por lo tanto, podrían estar implicados en la resistencia.

Para conocer el planteamiento, las preguntas de investigación, la metodología y
los resultados, consulta el [reporte de proyecto](docs/reporte-proyecto.md).

## Funcionalidades

- Descarga de los datos
- Control de calidad de las lecturas crudas.
- Trimming y limpieza
- Descarga del genoma de referencia y anotación
- Alineamiento del genoma 
- Cuantificación de expresión
- Control de calidad post-alineamiento
- Análisis de expresión diferencial
- Anotación funcional y enriquecimiento

<!--
AYUDA:
Describan solamente capacidades que ya funcionan. Comiencen cada elemento con un
verbo. Las funciones planeadas deben registrarse como issues y/o en GitHub Projects.

EJEMPLO:
- Descarga genomas mediante identificadores de NCBI.
- Lee archivos FASTA y GFF3.
- Genera una tabla de presencia y ausencia de genes.
-->

## Estructura del repositorio

```text
proyecto/
├── data/            # Datos de ejemplo o archivos pequeños
├──── curated/       # Datos curados
├────── resistent/   # Datos de organismo resistente curados
├────── susceptible/ # Datos de organismo susceptible curados
├──── raw/           # Datos crudos
├────── resistent/   # Datos de organismo resistente crudos
├────── susceptible/ # Datos de organismo susceptible crudos
├── docs/            # Reporte y documentación
├── results/         # Tablas y figuras generadas
├── src/             # Código fuente
├── tests/           # Pruebas
├── CITATION.cff     # Información para citar el software
├── LICENSE          # Licencia
└── README.md        # Introducción y guía rápida de uso
```
>

## Requisitos

- Python 3.14.3
- Git
- Conexión a internet para la descarga de datos

<!--
AYUDA:
Indiquen lo necesario antes de instalar o ejecutar el proyecto: versión de Python,
sistema operativo si es relevante, memoria, almacenamiento o acceso a servicios.

EJEMPLO:
- Python 3.12 o posterior.
- Git.
- 4 GB de memoria RAM.
- Conexión a internet para obtener datos.

No registren usuarios, contraseñas, tokens ni llaves privadas-->


## Datos # PENDIENTE

[Expliquen cómo obtener los datos y dónde colocarlos.]

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

Los datos proceden de *NCBI Sequence Read Archive (SRA)*, asociado al paper de Bonizzoni et al. 2015 (Parasites & Vectors, DOI: 10.1186/s13071-015-1083-z). Con el identificador: *SRP052073*

# Pendiente: Se ejecutará un módulo del script que busca los datos usando Bio.Entrez, pero aún no se ha determinado el script, por lo que, aún no tengo el uso de éste.

Los archivos se guardarán en data/raw/. Los identificadores, versiones y criterios de selección se documentan en docs/reporte-proyecto.md.

La procedencia y características detalladas se describen en el
[reporte del proyecto](docs/reporte-proyecto.md).

## Uso #PENDIENTE

[Expliquen qué hace el siguiente comando y qué entradas necesita.]

```bash
python src/main.py [argumentos]
```

[Indiquen dónde se guardan los resultados.]

<!--
AYUDA:
Incluyan al menos un ejemplo mínimo que pueda copiarse y ejecutarse. Sustituyan
los corchetes por valores reales. Expliquen entradas, opciones importantes y
archivos de salida sin describir toda la metodología.

EJEMPLO:
Para analizar los identificadores incluidos en data/accessions.txt:

    python src/main.py --input data/accessions.txt --output results/

El comando generará results/gene_matrix.csv y results/heatmap.png.
-->

## Reproducción de resultados # PENDIENTE

Ejecuten los siguientes pasos en el orden indicado:

```bash
python scripts/download_data.py
python src/run_analysis.py
python src/create_figures.py
```

Los resultados esperados se generarán en `results/`.

<!--
AYUDA:
Proporcionen la ruta más corta para regenerar el resultado principal desde los
datos originales. Los comandos deben indicar el orden correcto. Si el proceso
requiere parámetros o archivos de configuración, indíquenlos.
Las explicaciones científicas e interpretación de resultados pertenecen al reporte.

EJEMPLO:
Después de ejecutar los tres comandos se crearán la tabla comparativa y las
figuras utilizadas en el reporte. Sus nombres esperados deben indicarse aquí.
-->


## Documentación #PENDIENTE

- El reporte contiene el problema, las preguntas, la metodología y los resultados se encuentra en: (docs/reporte-proyecto.md)
- El reporte acerca de la manera correcta de citar este artículo se encuentra en: (CITATION.cff)

<!--
AYUDA:
Incluyan únicamente documentos que existan y comprueben sus enlaces. Agreguen
otros documentos sólo si evitan que el README sea demasiado extenso.

EJEMPLO:
- El reporte contiene el problema, las preguntas, la metodología y los resultados.
-->


## Equipo

- García Martínez Mónica — Colaborador
- Ledesma Gallegos Diego — Colaborador
- Muñoz Muñoz Luis Angel — Colaborador

<!--
AYUDA:
Identifiquen a las tres personas. Describan brevemente sus contribuciones generales.
Las actividades específicas y revisiones se consultan en los issues, Pull Requests,
GitHub Projects y el historial de Git.

EJEMPLO:
- Ana Pérez — procesamiento de datos y pruebas.
- Luis López — análisis y visualización.
- María García — documentación e integración.
-->

## Citación #PENDIENTE

Si utilizas este software, consulta [CITATION.cff](CITATION.cff) o la opción
**Cite this repository** de GitHub.

<!--
AYUDA:
Mantengan CITATION.cff actualizado con autores, título y versión. Si el proyecto
obtiene un DOI, añadan aquí la referencia recomendada.

EJEMPLO:
La forma recomendada de citar la versión v1.0.0 se encuentra en CITATION.cff.
-->

## Licencia

### **MIT License**

Consulta MIT License (LICENSE) para conocer los términos de uso

## Agradecimientos #PENDIENTE

Agradezco a mi gata Daysi por su odio hacia mi persona. Me motiva el hecho de acabar con su paciencia para seguir adelante

<!--
AYUDA:
Esta sección es opcional. Reconozcan laboratorios, docentes, instituciones,
proyectos o financiamientos que apoyaron el trabajo. No incluyan como autores a
personas que sólo deban aparecer en los agradecimientos.

EJEMPLO:
Proyecto desarrollado en la Licenciatura en Ciencias Genómicas, UNAM, como parte
de la asignatura [nombre].
-->

---

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
