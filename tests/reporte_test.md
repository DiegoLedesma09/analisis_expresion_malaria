# Reporte de la prueba de testeo

## Objetivo

Se realiza la prueba de testeo para demostrar el objeto inicial de este proyecto el cual es crear un herramienta para realizar un análisis de expresión diferencial de 1 organismo comparando 2 condiciones una en respuesta a la otra.

## Contexto de la prueba

- Contamos con 6 lecturas del organismo *Aedes aegypti* provenientes de una tabla con los SRR a analizar y un taxid 7159 y un id Bioproject PRJNA801226 (estos últimos dos no presentes en la tabla).

## Diferencias entre el testeo y el proyecto original

- Se coloco el argumento de `--runs_file` para la función `read_runs_file` para indicarle el archivo a leer con los SRR que se buscaran y descargaran más adelante.
- Descarga las lecturas y sigue normalmente a la funcion de `qc_trimming`y `count_features`.
- Al contar con 3 o más muestras en el módulo `differential_expresion` nos permite generar el PCA y el volcanoplot con padj que en nuestro proyecto original por tener menos de 3 muestras nos impedia generarlo. (La mayor diferencia entre el testeo y el proyecto original)

## Resultados de cada módulo

### Descarga de datos

- `download_data`: Realizo la descarga del GTF, las lecturas SRR y el genoma de referencia

### Limpieza de lecturas

- `qc_trimming`: Realiza la limpieza de las 6 lecturas por medio de `trim_galore` y al momento de revisar las gráficas que genera tanto de las `susceptible`y `resistant` tienen una excelente calidad de secuencia por base el cual nos indica que estan en buenas condiciones y podemos seguir con el análisis.

### Alineacion de lecturas

- `align_reads`: Posterior a tener las lecturas limpias procedemos a alinearlas y ubicarlas respecto al genoma de referencia y al proporcionar un archivo GTF ubica exones y sitios de empalme conocidos de ese organismo, regresandonos un archivo con posiciones genómicas ubicadas en .bam.

### Expresión diferencial

- `differential_expresion`: A partir del archivo con las ubicaciones en .bam se crea una tabla respecto a la condición, elimina genes con baja expresión y comienza con la expresión diferencial y los divide por up de cada condición. Al tenerlos clasificados nos genera un PCA (donde podemos apreciar como se agrupan las lecturas de resistant y susceptible se agrupan por su condición), un volcano plot (que nos confirmara con valor estadistico los genes que se expresan) y un heatmap (el cual nos mostrara los 25 genes de up_resistant y up_susceptible).

### Anotación funcional

- `functional_annotation`: Usando la tabla de expresión diferencial los separa por condicion y mediante Gene Ontology y construye un universo estadistico usando los genes expresados durante la expresión diferencial y los compara contra el resto del genoma expresado para ver que procesos biológicos se repiten más de lo esperado por azar, agrupa funciones emparentadas y filtra falsos positivos y con las gráficas barplot y dotplot permite ver que mecanismos usa para sobrevivir a insecticida

## Problemas

- Si se cuenta con el Bioproject y taxid no se tendra ninguno.

## Conclusiones

Con este test queda demostrado que nuestra herramienta funciona correctamente contando con 2 o más de 3 lecturas y respecto a este parámetro decidira que ruta elegir y como ejemplificarlo gráficamente.