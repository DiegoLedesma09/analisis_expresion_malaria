# Expresión diferencial entre poblaciones de *Anopheles gambiae*, mosquito vector de la malaria, resistentes y susceptibles a insecticidas

## Información general

| Dato | Nombre |
|:--|:--|
| Integrante 1 | García Martínez Mónica  monicagm@lcg.unam.mx |
| Integrante 2 | Ledesma Gallegos Diego  diegolg@lcg.unam.mx |
| Integrante 3 | Muñoz Muñoz Luis Angel  luismm@lcg.unam.mx |
| Fecha de creación | 25/08/2026 |
| Última actualización | 21/11/2026 |
| Repositorio | https://github.com/DiegoLedesma09/analisis_expresion_malaria.git |

## Resumen del proyecto

Se analizarán los datos disponibles en la base de datos del NCBI para comparar los perfiles de expresión génica (transcriptomas) de mosquitos *Anopheles gambiae* resistentes a insecticidas (deltametrina) frente a los susceptibles, con el objetivo es identificar qué genes se expresan de manera diferencial y que, por lo tanto, podrían estar implicados en la resistencia a insecticidas. El proyecto no ahondará en el análisis comparativo poblacional de metagenomas.

## 1. Contexto y antecedentes

La malaria (o paludismo) es una enfermedad causada por un parásito *Plasmodium*, el cual es trasmitido por la picadura de un mosquito infectado. Sólo el género *Anopheles* del mosquito transmite la malaria. Los síntomas de esta enfermedad pueden incluir fiebre, vómito y/o dolor de cabeza. La forma clásica de manifestación en el organismo es fiebre, sudoración y escalofríos, que aparecen 10 a 15 días después de la picadura del mosquito. Representa unproblema de salud publica en África subsahariana, ya que su morbilidad y mortalidad asociada afecta a esta población.

El control de la malaria se basa en el uso de intervenciones basadas en productos químicos en forma de mosquiteros tratados con insecticida de larga duración y pulverización residual en interiores.
En muchas zonas se están reforzando las medidas de control de los mosquitos, pero hay obstáculos importantes, tales como:

- El desarrollo de resistencia a insecticidas (IR) por parte del mosquito hacia los más clave como el DDT y los piretroides, sobre todo en África. 
- La inexistencia de insecticidas eficaces alternativos. 

No hay insecticidas alternativos que sean igual de eficaces y eficientes que el DDT y los piretroides, y el desarrollo de nuevos plaguicidas es una empresa cara y a largo plazo. Son imprescindibles prácticas de control de los vectores que hagan un uso juicioso de los insecticidas.
La detección de la resistencia a los insecticidas debe ser parte integral de las medidas nacionales de control, con el fin de garantizar que se están utilizando los métodos más eficaces de control de los vectores.

## 2. Planteamiento del problema

A pesar de que la resistencia a insecticidas en *Anopheles gambiae* es un problema crítico para el control de la malaria, los mecanismos moleculares no están completamente caracterizados, en parte porque los mecanismos varían entre poblaciones geográficas. Además, *Anopheles gambiae* está generando nuevos mecanismos para inhibir la acción del insecticidas; estos generan un mayor conflicto al caracterizar la vía que usan.

## 3. Justificación

La malaria es una enfermedad contagiosa en vectores cómunes en zonas del sudeste Africano,
la prevalencia del vector de transmición ha sido una problemática que produjo la creación
de insecticidas hacia ellos, pero debido a fuerzas de evolución, dichos organismos han evolucionado,
de tal manera que son resistentes a insecticidas cómunes. En este trabajo se busca encontrar los 
genes que promueven esta resistencia.

El análsis comparativo de datasets transcriptómicos requiere de un tratamiento estricto
de errores, por lo que, repetir el análisis transcriptómico con un tratamiento de *batch
error* podemos obtener posteriormente un perfil descriptivo basado en los genes de alta
expresividad en el trabajo.

## 4. Objetivo general

Tratar datos provenientes de secuenciación masiva (RNA-seq), identificar métricas de expresión transcriptómica y posteriormente identificar los genes que confieren la resistencia a insecticidas en *Anopheles gambiae*.

## 5. Preguntas de investigación

### Pregunta 1
**Pregunta:** ¿Qué genes se sobre expresan en mosquitos resistentes en presencia de insecticidas a comparación a los susceptibles?  
**Evidencia necesaria:** Genes con un Log2FoldChange positivo y significativo en los archivos de conteo de los organismos resistentes en comparación a organismos susceptibles.

### Pregunta 2
**Pregunta:** ¿Es posible relacionar los genes diferenciales con algún rasgo que les confiera la resistencia?  
**Evidencia necesaria:** Análisis de enriquecimiento funcional (Gene Ontology) para identificar procesos biológicos sobrerrepresentados en los genes candidatos.

## 6. Alcance y limitaciones

### Incluye
- Incluye un análisis de datasets de datos transcriptómicos de *Anopheles gambiae*, con organismos susceptibles y resistentes obtenidos de campo.
- Incluye un estudio de datos de expresión de los datos transcriptómicos.
- Incluye un enfoque con genes con el fin de una caracterización.
- Los resultados esperados son una lista de genes que promuevan la resistencia a insecticidas.

### No incluye
- No incluye un análisis comparativo poblacional de metagenomas.
- No incluye la secuenciación de los datos ni la obtención de los organismos.

### Limitaciones conocidas
- Podríamos encontrar genes que no estén asociados con una función pero tengan una sobreexpresión.
- El dataset público solo cuenta con una réplica biológica por condición (n=1), lo que limita el uso de pruebas estadísticas paramétricas tradicionales.
- El análisis se realizará en un plazo de 3 meses.

## 7. Propuesta de solución

Un programa modular (pipeline) que ejecuta las siguientes funciones:
- Descarga de los datos y genoma de referencia desde NCBI.
- Control de calidad de las lecturas crudas y limpieza (Trimming).
- Alineamiento de las lecturas contra el genoma de referencia.
- Cuantificación de expresión (matriz de conteos).
- Análisis de expresión diferencial y visualización de datos.
- Anotación funcional y enriquecimiento (Gene Ontology).

### 7.1 Resultado o producto esperado
Repositorio ejecutable con scripts, casos prueba, documentación, figuras regenerables, tests y datos accesibles.

## 8. Datos

### 8.1 Fuentes de datos

| Fuente | Identificador o versión | URL | Fecha de consulta | Licencia o condiciones |
|:--|:--|:--|:--|:--|
| NCBI Sequence Read Archive (SRA), asociado al BioProject PRJNA170440  | SRP052073 / PRJNA170440 | https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=170440 | 01/09/2026 | Público. Datos de RNA-seq de Anopheles gambiae resistentes y susceptibles a deltametrina. | 
| Publicación original (Bonizzoni et al., 2012) | PMID: 22970263 | https://pubmed.ncbi.nlm.nih.gov/22970263/ | 01/09/2026 | Artículo de acceso abierto (PLoS ONE). Describe el diseño experimental y los datos depositados en SRA. | 

### 8.2 Características de los datos

- **Organismo y taxonomía**:
  *Anopheles gambiae (forma S)*.
- **Diseño experimental / grupos**:
  Presentamos los análisis basados ​​en RNA-seq de los cambios transcripcionales en mosquitos   *Anopheles gambiae* de África Oriental, clasificados como resistentes o susceptibles a la deltametrina según la prueba de la OMS.
- **Origen geográfico**:
  Larvas colectadas en la Provincia Occidental de Kenia: localidades de Bungoma, Busia y Emutete.
  Año de colecta: 2012.
  Criadas hasta adultez en insectario del KEMRI, Kisumu.
- **Fenotipado**:
  WHO tube test estándar con deltametrina (0.05%).
  Resistentes: vivos 24h post-exposición.
  Susceptibles: derribados tempranamente (knock-down), sin signos de recuperación.
  Genotipado individual adicional de la mutación kdr (L1014S).
- **Datos moleculares**:
  Tipo de ensayo: RNA-seq.
  Plataforma: Illumina HiSeq2500.
  Tipo de lectura: pareada (paired-end).
  Longitud de lectura: 100 pb por extremo.
  Extracción de ARN: TRIzol.
  Lugar de secuenciación: DNA Technologies and Expression Analysis Core, UC Davis Genome Center.
- **Ensamblaje:** VectorBase AgamP3
- **Anotación:** AgamP3.7
- **Profundidad de secuenciación**: Entre 124,918,668 y 154,463,158 lecturas alineadas por librería, sin diferencias significativas entre muestras. 

### 8.3 Organización de los datos

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
```

### 8.4 Diccionario o formato de los datos

| Archivo o conjunto | Campo/columna | Tipo | Descripción | Valores o unidades |
|:--|:--|:--|:--|:--|
| SRR520428.fasta | - | FASTA | Lecturas transcriptómicas crudas | Pares de bases (pb) |
| SRR520427.fasta | - | FASTA | Lecturas transcriptómicas crudas | Pares de bases (pb) |

## 9. Metodología

### 9.1 Etapas del análisis o desarrollo

1. Descarga de datos crudos y genoma de referencia.
2. Control de calidad (FastQC) y limpieza de adaptadores (Trim Galore).
3. Alineamiento de lecturas al genoma (HISAT2) y procesamiento (samtools).
4. Cuantificación de expresión génica (featureCounts).
5. Análisis de expresión diferencial (InMoose / Log2FC).
6. Anotación funcional y enriquecimiento (goatools).

### 9.2 Herramientas y tecnologías

Para correr el programa, es necesaria la activación del ambiente **malaria_env**. [`malaria_env_linux.yml`](../envs/malaria_env_linux.yml)

El siguiente bloque de código explica cómo activar el ambiente.

```bash
conda env create -f envs/malaria_env.yml
conda activate malaria_env
```

Las herramientas contenidas en el ambiente son las siguientes:

| Herramienta | Versión | Propósito |
|:--|:--|:--|
| FASTQC (MultiQC) | v0.12.1. | Revisa calidad por base, contenido de adaptadores, duplicación y contenido de GC antes de seguir. |
| Trimmomatic | v0.41.| Eliminar adaptadores Illumina y bases de baja calidad en los extremos de las lecturas. |
| HISAT2 | v2.2.3 | Alineamiento de lecturas splice-aware contra el genoma. |
| samtools | v1.25 | Ordenamiento e indexación de archivos BAM. |
| featureCounts | v2.1.1 | Cuantificación de lecturas por gen. |
| InMoose (DESeq2) | v0.9.1 | Análisis de expresión diferencial (si hay réplicas). |
| goatools | v1.6.9 | Enriquecimiento funcional de Gene Ontology. | 

### 9.3 Estrategia de validación

Se realizaron pruebas con los datos de prueba descargados, verificando que el pipeline generara los archivos esperados en cada etapa (BAM indexados, matrices de conteos, gráficos). Además, se implementó un modo tolerante a fallos para el análisis de expresión diferencial, adaptándose al caso n=1.

## 10. Plan de trabajo

### 10.1 Distribución de responsabilidades

El desarrollo de este proyecto se llevó a cabo mediante un flujo de trabajo colaborativo y revisado por pares. Todos los integrantes participaron activamente en la creación de las funciones de cada módulo del pipeline. Para garantizar la calidad y la reproducibilidad del código, se implementó un flujo de integración basado en *Pull Requests* (PR): cada módulo fue desarrollado en una rama independiente, y ningún cambio se fusionó a la rama principal (`main`) sin antes haber sido revisado y aprobado por al menos otro miembro del equipo. 

*Nota: Todos los integrantes contribuyeron a la discusión del diseño de los módulos y a la validación de los resultados finales.*

| Integrante | Responsabilidad principal | Responsabilidad de revisión |
|:--|:--|:--|
| Mónica García | Codificación de las funciones en cada módulo, así como la redacción de la documentación.  | Revisión de los *Pull requests*.|
| Diego Ledesma | Gestión del repositorio: Creación de *Issues*, Ramas y *Pull requests*. Así como codificación. | Revisión del funcionamiento de cada módulo |
| Luis Ángel Muñoz | Codificación de las funciones en cada módulo, así como la redacción de la documentación. | Revisión de los *Pull requests*.  |

### 10.2 Riesgos y alternativas

| Riesgo | Probabilidad | Impacto | Prevención o alternativa |
|:--|:--|:--|:--|
| Falta de réplicas biológicas (n=1) | Alta | Medio | Se implementó un flujo alternativo (Log2FC directo y Volcano sin p-valor). |
| Archivos GTF no coinciden con el índice | Media | Alto | Validación manual de los ensamblajes antes del conteo. |

## 11. Resultados

A partir del análisis de enriquecimiento funcional de Gene Ontology (GO) sobre los genes diferencialmente expresados entre las condiciones resistente y susceptible de *Anopheles gambiae*, se obtuvieron perfiles funcionales diferenciados para ambos grupos.

### Enriquecimiento funcional en `up_resistant` (Sobreexpresados en resistentes)

En el grupo `up_resistant`, los términos GO enriquecidos se concentran en plegamiento, replegamiento y estabilización de proteínas (*protein folding*, *protein folding chaperone*, *protein refolding*, *ATP-dependent protein folding chaperone*, *protein stabilization*, *regulation of protein stability*) así como en respuesta a estímulos y estrés (*response to abiotic stimulus*, *response to temperature stimulus*, *response to heat*). 
Los conteos de genes (`k_estudio`) van de 3 a 15 genes, y los valores de `padj` oscilan entre 0.002 y 0.011 ver [`go_dotplot_up_resistant.png`](../results/go/go_dotplot_up_resistant.png) y [`go_barplot_up_resistant.png`](../results/go/go_barplot_up_resistant.png). 

### Enriquecimiento funcional en `up_susceptible`

En el grupo `up_susceptible`, destacan términos asociados a actividad proteolítica e hidrolítica (*serine-type endopeptidase activity*, *serine-type peptidase activity*, *serine hydrolase activity*, *endopeptidase activity*, *proteolysis*, *peptidase activity*, *hydrolase activity*), junto con regulación/inhibición de peptidasas y región extracelular (*peptidase inhibitor activity*, *peptidase regulator activity*, *extracellular region*). 
Los conteos de genes son mucho mayores, entre 13 y 65 genes, con valores de `padj` extremadamente bajos (del orden de 1e-09), indicando una alta significancia estadística ver [`go_dotplot_up_susceptible.png`](../results/go/go_dotplot_up_susceptible.png) y [`go_barplot_up_susceptible.png`](../results/go/go_barplot_up_susceptible.png). 

## 12. Discusión

Los resultados obtenidos muestran una convergencia biológica muy interesante. En el grupo de genes sobreexpresados en los mosquitos resistentes (`up_resistant`), el enriquecimiento en vías de **plegamiento de proteínas y respuesta al calor/estrés** sugiere que los mosquitos que sobreviven al insecticida activan mecanismos de protección celular y reparación de proteínas dañadas. Esto concuerda directamente con lo reportado por Bonizzoni et al. (2012), quienes observaron una sobrerrepresentación de términos asociados a respuesta a estrés y estímulos abióticos en mosquitos resistentes a deltametrina.

Por el contrario, los mosquitos susceptibles (`up_susceptible`) mostraron un enriquecimiento masivo en vías de **proteólisis y actividad peptidasa**. Biológicamente, esto sugiere un colapso celular o una respuesta de degradación de proteínas provocada por el efecto tóxico letal del insecticida. El paper original también reporta que los genes sobreexpresados en la cepa susceptible estaban enriquecidos en funciones de proteólisis, lo que refuerza la validez de nuestro pipeline.

**Limitaciones y aspectos visuales a considerar:** El enriquecimiento GO resume categorías funcionales y depende de la cobertura y calidad de las anotaciones disponibles para el organismo y de los umbrales de filtrado previos.

## 13. Conclusiones

El flujo de enriquecimiento funcional GO permitió contrastar los perfiles biológicos entre los genes sobreexpresados en mosquitos resistentes y susceptibles, generando visualizaciones (dotplot y barplot) que apoyan el análisis transcriptómico. 
Los resultados demuestran que es posible extraer conclusiones biológicas de un dataset público con limitaciones (n=1) mediante la implementación de un flujo de trabajo robusto. El pipeline modular desarrollado en Python permite automatizar el análisis transcriptómico desde la descarga hasta la anotación funcional. A pesar de la falta de réplicas biológicas, el uso de herramientas como Log2FoldChange y el análisis GO permitió identificar que la resistencia a insecticidas en *Anopheles gambiae* está asociada a respuestas de estrés y plegamiento proteico, mientras que la susceptibilidad se asocia a proteólisis. Esto aporta un marco reproducible para futuros análisis con datos más robustos.

## 14. Disponibilidad, licencia y citación

**Código:** [GitHub Repo](https://github.com/DiegoLedesma09/analisis_expresion_malaria)

**Datos:** [NCBI BioProject PRJNA170440](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=170440)

**Licencia del código:** [LICENSE](../LICENSE)

**Cómo citar:** [CITATION.cff](../CITATION.cff)

**Metadatos:** Consultar [codemeta.json](../codemeta.json)

**Versión o release:** [version ](https://github.com/DiegoLedesma09/analisis_expresion_malaria/tags)

## 15. Referencias

1. Organización Panamericana de la Salud. (s.f.). Malaria. Recuperado el 29 de agosto de 2026, de https://www.paho.org/es/temas/malaria
2. Bonizzoni M, Ochomo E, Dunn WA, Britton M, Afrane Y, Zhou G, Hartsel J, Lee MC, Xu J, Githeko A, Fass J, Yan G. RNA-seq analyses of changes in the Anopheles gambiae transcriptome associated with resistance to pyrethroids in Kenya: identification of candidate-resistance genes and candidate-resistance SNPs. Parasit Vectors. 2015 Sep 17;8:474. doi: 10.1186/s13071-015-1083-z. PMID: 26381877; PMCID: PMC4574070.
3. Zoh, M.G., Bonneville, JM., Laporte, F. et al. Deltamethrin and transfluthrin select for distinct transcriptomic responses in the malaria vector Anopheles gambiae. Malar J 22, 256 (2023). https://doi.org/10.1186/s12936-023-04673-5
4. Saizonou H, Impoinvil LM, Derilus D, Omoke D, Okeyo S, Dada N, Corredor C, Mulder N, Lenhart A, Ochomo E, Djogbénou LS. Transcriptomic analysis of Anopheles gambiae from Benin reveals overexpression of salivary and cuticular proteins associated with cross-resistance to pyrethroids and organophosphates. BMC Genomics. 2024 Apr 6;25(1):348. doi: 10.1186/s12864-024-10261-x. PMID: 38582836; PMCID: PMC10998338.
