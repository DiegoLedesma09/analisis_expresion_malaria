# Expresión diferencial entre poblaciones de *Anopheles gambiae*, mosquito vector de la malaria, resistentes y susceptibles a insecticidas

## Información general

| Dato | Nombre |
|:--|:--|
| Integrante 1 | García Martínez Mónica  monicagm@lcg.unam.mx |
| Integrante 2 | Ledesma Gallegos Diego  diegolg@lcg.unam.mx |
| Integrante 3 | Muñoz Muñoz Luis Angel  luismm@lcg.unam.mx |
| Fecha de creación | 25/08/2026 |
| Última actualización | 25/08/2026 |
| Repositorio | https://github.com/DiegoLedesma09/analisis_expresion_malaria.git |

## Resumen del proyecto

Se analizarán los datos disponibles en la base de datos del NCBI para comparar los perfiles de expresión génica (transcriptomas) de mosquitos *Anopheles gambiae* resistentes a insecticidas frente a los susceptibles. El objetivo es identificar qué genes se expresan de manera diferencial y que, por lo tanto, podrían estar implicados en la resistencia.

## 1. Contexto y antecedentes

<!-- AYUDA: Presenten la información necesaria para comprender el proyecto.
Definan conceptos biológicos y computacionales, describan qué se conoce y citen
trabajos, datos o herramientas relacionados.
PREGUNTAS: ¿Cuál es el fenómeno de interés? ¿Qué debe conocer quien lea el
reporte? ¿Qué métodos o herramientas se han utilizado anteriormente?
EJEMPLO: La resistencia antimicrobiana es un problema de salud pública. Aunque
existen bases especializadas, comparar varios genomas requiere integrar datos
procedentes de distintos archivos. -->

La malaria (o paludismo) es una enfermedad causada por un parásito Plasmodium, el cual es trasmitido por la picadura de un mosquito infectado. Sólo el género anófeles del mosquito transmite la malaria. Los síntomas de esta enfermedad pueden incluir fiebre, vómito y/o dolor de cabeza. La forma clásica de manifestación en el organismo es fiebre, sudoración y escalofríos, que aparecen 10 a 15 días después de la picadura del mosquito. Dado a esto es un problema de salud publica en África subsahariana con su morbilidad y mortalidad asociada afecta esta población.

El control de la malaria se basa en el uso de intervenciones basadas en productos químicos en forma de mosquiteros tratados con insecticida de larga duración y pulverización residual en interiores.

En muchas zonas se están reforzando las medidas de control de los mosquitos, pero hay obstáculos importantes, tales como:
- El desarrollo de resistencia a insecticidas (IR) por parte del mosquito hacia los más clave como el DDT y los piretroides, sobre todo en África. 
- La inexistencia de insecticidas eficaces alternativos. 

No hay insecticidas alternativos que sean igual de eficaces y eficientes que el DDT y los piretroides, y el desarrollo de nuevos plaguicidas es una empresa cara y a largo plazo. Son imprescindibles prácticas de control de los vectores que hagan un uso juicioso de los insecticidas.
La detección de la resistencia a los insecticidas debe ser parte integral de las medidas nacionales de control, con el fin de garantizar que se están utilizando los métodos más eficaces de control de los vectores.

## 2. Planteamiento del problema

<!-- AYUDA: Describan la dificultad, necesidad o vacío de conocimiento que
desean atender, a quién afecta y qué sucede si no se resuelve. No confundan el
problema con la herramienta ni con la solución.
EJEMPLO: La identificación manual de genes de resistencia en varios genomas es
lenta, propensa a errores y difícil de reproducir. -->

A pesar de que la resistencia a insecticidas en Anopheles gambiae es un problema crítico para el control de la malaria, los mecanismos moleculares no están completamente caracterizados, en parte porque los mecanismos varían entre poblaciones geográficas 

Ademas de que Anopheles gambiae esta generando nuevos mecanismos para inhibir la acción del insecticida estos generan un mayor conflicto al caracterizar la vía que usan


## 3. Justificación

<!-- AYUDA: Expliquen por qué vale la pena realizar el proyecto, cuál es su
relevancia biológica, científica, técnica o social y quién podría beneficiarse.
EJEMPLO: Un flujo automatizado reducirá errores y permitirá repetir el análisis
con los mismos datos, parámetros y versiones del software. -->

La malaria es una enfermedad contagiosa en vectores cómunes en zonas del sudeste Africano,
la prevalencia del vector de transmición ha sido una problemática que produjo la creación
de insecticidas hacia ellos, pero debido a fuerzas de evolución, dichos organismos han evolucionado,
de tal manera que son resistentes a insecticidas cómunes. En este trabajo se busca encontrar los 
genes que promueven esta resistencia.

El análsis comparativo de datasets transcriptómicos requiere de un tratamiento estricto
de errores, por lo que, repetir el análisis transcriptómico con un tratamiento de batch
error podemos obtener posteriormente un perfil descriptivo basado en los genes de alta
expresividad en el trabajo.

## 4. Objetivo general

<!-- AYUDA: Expresen el resultado global mediante un verbo en infinitivo. Debe
ser alcanzable durante el semestre.
EJEMPLO: Desarrollar un flujo reproducible en Python para identificar y comparar
genes de resistencia en un conjunto de genomas de E. coli. -->

Tratar datos provenientes de archivos FASTA, identificar datos de expresión transcriptómicos y
posteriormente identificar los genes que proveen la resistencia a insecticidas.

## 5. Preguntas de investigación

<!-- AYUDA: Formulen preguntas biológicas o computacionales que puedan
responderse con los datos y métodos disponibles. Indiquen qué evidencia sería
necesaria.
EJEMPLO: ¿Qué genes de resistencia aparecen en cada genoma? Evidencia: una tabla
de presencia y ausencia obtenida de las anotaciones. -->

### Pregunta 1

**Pregunta:** ¿Qué genes se sobre expresan en mosquitos resistentes en presencia de insecticidas a comparación a los susceptibles?  
**Evidencia necesaria:** Genes sobreexpresados en las lecturas de los archivos FASTA en los organismos resistentes en presencia de insecticidas en comparación a organismos susceptibles

### Pregunta 2

**Pregunta:** Relacionar los genes que encontremos con algún rasgo que les confiera la resistencia
**Evidencia necesaria:** Relacionar con deSeq2 para identificar los genes según los datos transcriptómicos

## 6. Alcance y limitaciones

<!-- AYUDA: Delimiten organismos, muestras, datos, análisis y resultado esperado. Indiquen
qué no se abordará y las restricciones de tiempo, cómputo, acceso o calidad.
EJEMPLO: Se analizarán como máximo 20 genomas completos de RefSeq. No se
utilizarán datos clínicos ni se realizará validación experimental. -->

### Incluye

- Incluye un análisis de 9 datasets de datos transcriptómicos de *Anopheles gambiae*, con organismos susceptibles y resistentes obtenidos de campo, así como organimos resistentes obtenidos de laboratorio.
- Incluye un estudio de datos de expresión de los datos transcriptómicos
- Incluye un enfoque con genes con el fin de una caracterización
- Los resultados esperados son una lista de genes que promuevan la resistencia a insecticidas.

### No incluye

- No incluye un análisis comparativo poblacional de metagenomas
- No incluye la secuenciación de los datos ni la obtención de los organismos

### Limitaciones conocidas

- Podriamos encontrar genes que no esten asociados con una función pero tengan una sobreexpresión.
- Podríamos tener datasets de mala calidad
- El análisis se realizará en un plazo de 3 meses

## 7. Propuesta de solución

<!-- AYUDA: Describan el producto o estrategia que podría resolver el problema.
Es una propuesta inicial y puede cambiar. Expliquen sus componentes, no sólo las
tecnologías.
EJEMPLO: Un programa modular recibirá identificadores, descargará archivos,
extraerá genes, almacenará resultados y generará visualizaciones. -->

Un programa modular que hará las siguientes funciones:

- Descarga de los datos
- Control de calidad de las lecturas crudas.
- Trimming y limpieza
- Descarga del genoma de referencia y anotación
- Alineamiento del genoma 
- Cuantificación de expresión
- Control de calidad post-alineamiento
- Análisis de expresión diferencial
- Anotación funcional y enriquecimiento

### 7.1 Resultado o producto esperado

<!-- AYUDA: Indiquen el entregable concreto: programa, paquete, flujo de análisis,
base de datos, visualizaciones u otro producto.
EJEMPLO: Repositorio ejecutable con scripts, datos de prueba, documentación,
tabla comparativa y figuras regenerables. -->

Repositorio ejecutable con scripts, casos prueba, documentación, figuras regenerables, tests y datos accesibles.

## 8. Datos

### 8.1 Fuentes de datos

<!-- AYUDA: Incluyan institución, base de datos, URL, identificador, versión o
fecha de consulta y condiciones de uso. No todos los proyectos usarán NCBI.
EJEMPLO: NCBI RefSeq, GCF_000005845.2, consultado el dd/mm/aaaa. -->

| Fuente | Identificador o versión | URL | Fecha de consulta | Licencia o condiciones |
|:--|:--|:--|:--|:--|
| NCBI Sequence Read Archive (SRA), asociado al paper de Bonizzoni et al. 2015 (Parasites & Vectors, DOI: 10.1186/s13071-015-1083-z) | SRP052073 | URL: https://www.ncbi.nlm.nih.gov/sra/?term=SRP052073 | - Consultado el 01/09/2026 | Mosquitos de campo resistentes a deltametrina (Provincia Occidental de Kenia) - Mosquitos de campo susceptibles a deltametrina (misma región) - Cepa de laboratorio Kisumu (susceptible de referencia, altamente endogámica) |

### 8.2 Características de los datos

<!-- AYUDA: Describan organismos, muestras, variables, formatos, versiones, tamaño y otros
atributos necesarios para interpretar los datos.
EJEMPLO: Archivos FASTA y GFF3 de 20 genomas completos de E. coli. -->

*Organismo y taxonomía*

Anopheles gambiae (forma S)

*Diseño experimental / grupos*

3 condiciones: resistente a deltametrina (campo), susceptible a deltametrina (campo), cepa Kisumu (susceptible de laboratorio, referencia)
9 librerías de RNA-seq en total
Cada librería = pool de ARN de 12 mosquitos individuales (no son mosquitos individuales secuenciados por separado)

*Origen geográfico*

Larvas colectadas en la Provincia Occidental de Kenia: localidades de Bungoma, Busia y Emutete
Año de colecta: 2012
Criadas hasta adultez en insectario del KEMRI, Kisumu

*Fenotipado*

WHO tube test estándar con deltametrina (0.05%)
Resistentes: vivos 24h post-exposición
Susceptibles: derribados tempranamente (knock-down), sin signos de recuperación
Genotipado individual adicional de la mutación kdr (L1014S) en el gen para

*Datos moleculares*

Tipo de ensayo: RNA-seq
Plataforma: Illumina HiSeq2500
Tipo de lectura: pareada (paired-end)
Longitud de lectura: 100 pb por extremo
Extracción de ARN: TRIzol
Lugar de secuenciación: DNA Technologies and Expression Analysis Core, UC Davis Genome Center

*Referencia usada en el análisis original*

Ensamblaje: VectorBase AgamP3
Anotación: AgamP3.7

Profundidad de secuenciación

Entre 124,918,668 y 154,463,158 lecturas alineadas por librería, sin diferencias significativas entre muestras


### 8.3 Organización de los datos

<!-- AYUDA: Muestren la estructura prevista. No suban datos sensibles, tokens,
contraseñas ni archivos grandes. Usen .gitignore y documenten cómo obtener lo
que no se guarde en Git.
EJEMPLO: data/raw conserva originales y data/processed los derivados. -->

```text
proyecto/
├── data/
│   ├── raw/
|      ├── resistent/
|      ├── susceptible/
│   ├── curated/
|      ├── resistent/
|      ├── susceptible/
├── docs/
├── notebooks/
├── results/
├── src/
└── tests/
```

### 8.4 Diccionario o formato de los datos

<!-- AYUDA: Describan campos o columnas relevantes. Incluyan fragmentos pequeños
cuando ayuden a comprender el formato, pero no archivos completos.
EJEMPLO: En GFF3, seqid identifica la secuencia; type indica gene, CDS, etc. -->

| Archivo o conjunto | Campo/columna | Tipo | Descripción | Valores o unidades |
|:--|:--|:--|:--|:--|
| SRR1763908.fasta | - | FASTA | Lecturas transcriptómicas | Lecturas de transcriptoma |

## 9. Metodología

<!-- AYUDA: Esta sección evolucionará. Primero describan el plan y después
actualícenla con lo que realmente ejecutaron, incluidos parámetros y decisiones. -->

### 9.1 Etapas del análisis o desarrollo

<!-- AYUDA: Describan la secuencia desde la obtención de datos hasta la validación
de resultados. Relacionen cada etapa con una pregunta u objetivo.
EJEMPLO: descarga, validación, transformación, análisis, visualización y pruebas. -->

1. Descarga de los datos
2. Control de calidad
3. Alineamiento del genoma
4. Cuantificación de expresión
5. Análisis de expresión diferencial
6. Anotación funcional

### 9.2 Herramientas y tecnologías

<!-- AYUDA: Registren lenguajes, bibliotecas y programas con sus versiones y
propósito. No incluyan credenciales.
EJEMPLO: Python 3.x; Biopython para leer formatos biológicos; Seaborn para
visualización. -->

| Herramienta | Versión | Propósito |
|:--|:--|:--|
| FASTQC (MultiQC) | v0.12.1. | Revisa calidad por base, contenido de adaptadores, duplicación y contenido de GC antes de seguir. |
| Trimmomatic | v0.41.| Eliminar adaptadores Illumina y bases de baja calidad en los extremos de las lecturas. |
| HISAT2 o STAR | v2.2.3 / v2.7.11b | Alinea las lecturas limpias splice-aware |
| featureCounts | v2.1.1 | Cuenta lecturas por gen |
| DESeq2 | v1.52.0  | Comparar resistentes vs. susceptibles, obteniendo log2FoldChange y p-valores ajustados (FDR) por gen |
| Profiler | v0.2.4. | Identificar categorías GO/KEGG sobrerrepresentadas | 

### 9.3 Estrategia de validación

<!-- AYUDA: Expliquen cómo comprobarán código y resultados: pruebas unitarias,
datos conocidos, comparación con otra herramienta o revisión manual.
EJEMPLO: Se compararán cinco anotaciones conocidas y se probarán entradas
válidas, identificadores inexistentes y archivos incompletos. -->

Se realizarán pruebas con datos conocidos, con el fin de reconocer la validez de nuestros pasos.

## 10. Plan de trabajo


### 10.1 Distribución de responsabilidades

<!-- AYUDA: Definan responsabilidades iniciales sin aislar a cada integrante.
Toda contribución importante debe ser revisada mediante Pull Request por otra
persona.
EJEMPLO: Ana desarrolla la descarga y revisa el módulo de visualización. -->

| Integrante | Responsabilidad principal | Responsabilidad de revisión |
|:--|:--|:--|
| [Nombre] | [Responsabilidad] | [Qué o a quién revisará] |

### 10.2 Riesgos y alternativas

<!-- AYUDA: Identifiquen situaciones que podrían impedir o retrasar el proyecto
y definan una alternativa.
EJEMPLO: Los datos requieren demasiado almacenamiento; alternativa: reducir el
número de genomas usando criterios documentados. -->

| Riesgo | Probabilidad | Impacto | Prevención o alternativa |
|:--|:--|:--|:--|
| [Riesgo] | Baja/Media/Alta | Bajo/Medio/Alto | [Acción] |

## 11. Resultados

<!-- AYUDA: Presenten resultados vinculados con preguntas y objetivos. Incluyan
tablas o figuras con títulos, leyendas y archivos de origen. Describan aquí lo
obtenido; interprétenlo en Discusión.
EJEMPLO: Tabla de presencia y ausencia generada por src/compare_genes.py.
PRIMERA SESIÓN: dejen esta sección como pendiente. -->

> Estado: pendiente. Se completará durante el desarrollo.



## 12. Discusión

<!-- AYUDA: Interpreten los resultados, expliquen si responden las preguntas,
compárenlos con los antecedentes y señalen limitaciones. No repitan únicamente
los valores.
EJEMPLO: La distribución observada sugiere..., aunque la interpretación está
limitada por la calidad de las anotaciones. -->

> Estado: pendiente. Se completará después de obtener resultados.

## 13. Conclusiones

<!-- AYUDA: Sinteticen qué se aprendió, qué preguntas se respondieron y si se
alcanzaron los objetivos. Incluyan aportes, limitaciones y trabajo futuro. No
introduzcan resultados nuevos.
EJEMPLO: El flujo permitió identificar..., pero será necesario incorporar... -->

> Estado: pendiente. Se completará al finalizar el proyecto.


## 14. Disponibilidad, licencia y citación

<!-- AYUDA: Indiquen dónde está el código, bajo qué licencia puede reutilizarse
y cómo citarlo. Relacionen esta sección con LICENSE, CITATION.cff, codemeta.json,
release final y, cuando corresponda, un DOI.
EJEMPLO: Código en GitHub bajo MIT; cita disponible en CITATION.cff. -->

**Código:** [URL]  
**Datos:** [URL, identificador o instrucciones]  
**Licencia del código:** [Licencia]  
**Cómo citar:** [Referencia o enlace a CITATION.cff]  
**Versión o release:** [URL]

## 15. Referencias

<!-- AYUDA: Registren publicaciones, datos, software y documentos consultados en
un formato uniforme. Incluyan DOI, URL o identificadores persistentes. Toda cita
del texto debe aparecer aquí.
EJEMPLO: Blattner, F. R. et al. (1997). The complete genome sequence of
Escherichia coli K-12. Science, 277(5331), 1453–1462.
https://doi.org/10.1126/science.277.5331.1453 -->

1. Organización Panamericana de la Salud. (s.f.). Malaria. Recuperado el 29 de agosto de 2026, de https://www.paho.org/es/temas/malaria
2. Bonizzoni M, Ochomo E, Dunn WA, Britton M, Afrane Y, Zhou G, Hartsel J, Lee MC, Xu J, Githeko A, Fass J, Yan G. RNA-seq analyses of changes in the Anopheles gambiae transcriptome associated with resistance to pyrethroids in Kenya: identification of candidate-resistance genes and candidate-resistance SNPs. Parasit Vectors. 2015 Sep 17;8:474. doi: 10.1186/s13071-015-1083-z. PMID: 26381877; PMCID: PMC4574070.
3. Zoh, M.G., Bonneville, JM., Laporte, F. et al. Deltamethrin and transfluthrin select for distinct transcriptomic responses in the malaria vector Anopheles gambiae. Malar J 22, 256 (2023). https://doi.org/10.1186/s12936-023-04673-5
4. Saizonou H, Impoinvil LM, Derilus D, Omoke D, Okeyo S, Dada N, Corredor C, Mulder N, Lenhart A, Ochomo E, Djogbénou LS. Transcriptomic analysis of Anopheles gambiae from Benin reveals overexpression of salivary and cuticular proteins associated with cross-resistance to pyrethroids and organophosphates. BMC Genomics. 2024 Apr 6;25(1):348. doi: 10.1186/s12864-024-10261-x. PMID: 38582836; PMCID: PMC10998338.


---

<!-- ORIENTACIÓN PARA LAS DOS PRIMERAS SESIONES:
Completen Información general, Resumen provisional, secciones 1 a 7, fuentes de
datos preliminares y plan de trabajo. Metodología, Resultados, Discusión,
Conclusiones Y Disponibilidad evolucionarán durante el
semestre. Sustituyan las indicaciones entre corchetes por contenido del equipo. -->