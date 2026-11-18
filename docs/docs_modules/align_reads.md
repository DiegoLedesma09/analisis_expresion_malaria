# Módulo de alineamiento de lecturas de RNA-seq (`align_reads.py`)

## Resumen 

Este módulo toma las lecturas limpias generadas por `qc_trimming.py` y las alinea contra el genoma de referencia de *Anopheles gambiae* usando **HISAT2**.
Posteriormente, convierte los alineamientos a formato BAM ordenado e indexado con **samtools**, y evalúa la calidad del alineamiento por muestra.

El módulo realiza lo siguiente:

1. Construye el índice del genoma de referencia con `hisat2-build` (omitiendo este paso si el índice ya existe).
2. Si se proporciona un archivo de anotación GTF, extrae automáticamente los sitios de empalme (*splice sites*) y exones para optimizar el mapeo de transcritos.
3. Itera sobre los grupos experimentales (`resistant` y `susceptible`) y alinea cada par de lecturas curadas (`_val_1.fq`, `_val_2.fq`) con `hisat2`, generando el archivo SAM y el reporte de métricas por muestra.
4. Convierte cada archivo SAM a formato BAM ordenado por coordenadas e indexado (`.bai`) mediante `samtools`, eliminando el archivo SAM temporal para optimizar el almacenamiento en disco.
5. Evalúa los resúmenes de alineamiento de cada grupo y emite una advertencia en el registro si alguna muestra queda por debajo del umbral mínimo de mapeo establecido.

## Requisitos 
- Python 3.12 
- **HISAT2** (`hisat2`, `hisat2-build`, `hisat2_extract_splice_sites.py` y `hisat2_extract_exons.py`)
- **samtools**
- Lecturas pareadas previamente curadas (`*_val_1.fq` y `*_val_2.fq`) organizadas por fenotipo en `data/curated/`

**Se recomienda ejecutar desde la raíz del proyecto, ya que las rutas por defecto son relativas (`data/`, `logs/`).**

## Estructura de archivos generados 

```text
analisis_expresion_malaria/
├── data/
│   ├── index/
│   │   └── hisat2/
│   │       ├── genome.exon          # Coordenadas de exones (si se usa -gtf)
│   │       ├── genome.ss            # Sitios de empalme (si se usa -gtf)
│   │       └── genome.[1-8].ht2     # Archivos binarios del índice HISAT2
│   └── aligned/
│       ├── resistant/
│       │   ├── SRR*.bam             # Alineamientos ordenados (cepa resistente)
│       │   ├── SRR*.bam.bai         # Índices binarios de los BAM
│       │   └── SRR*.summary.txt     # Resumen de tasas de alineamiento
│       └── susceptible/
│           ├── SRR*.bam             # Alineamientos ordenados (cepa susceptible)
│           ├── SRR*.bam.bai         # Índices binarios de los BAM
│           └── SRR*.summary.txt     # Resumen de tasas de alineamiento
└── logs/
    └── align_reads.log

## Pipeline 

```mermaid
flowchart TD
    A([Lecturas _val_1.fq y _val_2.fq]) --> B[build_index]
    B -.->|Si no existe| B1[(Índice .ht2)]
    B --> C[run_hisat2]
    C --> C1[(Archivo .sam)]
    C --> C2[(Archivo .summary.txt)]
    C1 --> D[convert_sam_to_bam]
    D --> D1[(BAM ordenado e indexado)]
    C2 --> E[get_alignment_rate]
```

## Componentes clave 

### Indexación del genoma

- build_index(fasta_path, index_prefix, threads, gtf_path): Construye el índice de HISAT2. Si el archivo .1.ht2 ya existe, omite el proceso para ahorrar tiempo. 
  Si se le pasa un GTF, extrae los sitios de empalme y exones.

### Alineamiento con HISAT2

- run_hisat2(fwd, rev, index_prefix, output_sam, summary_file, threads): Ejecuta HISAT2 en modo paired-end (-1 y -2). Guarda el alineamiento en formato SAM y el resumen de métricas en un .txt.

### Procesamiento y limpieza de BAM

- convert_sam_to_bam(sam_path, bam_path, threads=6) — transforma la salida cruda de HISAT2 mediante tres pasos secuenciales:

1. Ordena las lecturas por coordenadas genómicas con samtools sort -@ <threads> -o <bam> <sam>.

2. Construye el índice de acceso aleatorio (.bai) con samtools index -@ <threads> <bam>.

3. Elimina el archivo .sam original (os.remove) inmediatamente después de confirmar la indexación exitosa, evitando saturar el almacenamiento del servidor.

### Control de calidad de mapeo

- get_alignment_rate(summary_dir, suffix, umbral): Lee los resúmenes de HISAT2 con una expresión regular, extrae la tasa de alineamiento global y lanza un logger.warning si está por debajo del umbral (default: 70%).

## Uso por línea de comandos 

### Argumentos

| Argumento | Requerido | Default | Descripción |
| --- | --- | --- | --- |
| -r, --reference | Sí | — | Ruta al genoma de referencia (FASTA) |
| -i, --index | Sí | — | Ruta base y prefijo para el índice de HISAT2 |
| -gtf, --annotation | No | None | Ruta al archivo de anotación GTF |
| -d, --data_dir | No | data/curated/ | Directorio raíz de las lecturas curadas |
| -o, --outdir | No | data/aligned/ | Directorio de salida para los archivos BAM |
| -w, --workers | No | 6 | Hilos de procesamiento para hisat2 y samtools |
| -u, --umbral | No | 70.0 | Tasa mínima de alineamiento esperada (%) |
| -v, --verbose | No | False | Muestra mensajes DEBUG en consola |
| --log_file | No | logs/align_reads.log | Ruta del archivo de log |

## Ejemplos 

1. Ejecución con GTF para construir el índice desde cero

```bash
python src/align_reads.py -r data/genome/genome_Anopheles_gambiae.fasta -i data/index/hisat2/genome -gtf data/genome/genome_Anopheles_gambiae.gtf -v
```
