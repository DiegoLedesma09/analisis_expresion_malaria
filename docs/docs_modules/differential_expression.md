# Módulo de Expresión Diferencial (`differential_expression.py`)

## Resumen
Este módulo realiza el análisis de expresión diferencial entre los grupos `resistant` y `susceptible` a partir de la matriz de conteos generada por `featureCounts`. 

El módulo es **flexible**: detecta automáticamente cuántas réplicas biológicas hay. Si cuenta con 3 o más muestras por grupo, utiliza **InMoose (DESeq2)** para calcular p-valores ajustados. Si solo cuenta con 2 muestras (n=1 por grupo), calcula el `log2FoldChange` directamente a partir de datos normalizados (log2(CPM+1)). Finalmente, clasifica los genes, guarda los resultados y genera gráficos de alta calidad (PCA, Volcano Plot y Heatmap) adaptados a la cantidad de datos disponibles.

## Requisitos
- Python 3.12
- Librerías de Python: `inmoose`, `pandas`, `numpy`, `matplotlib`, `scikit-learn`, `seaborn`
- Matriz cruda generada por `featureCounts` (idealmente `data/counts/raw_counts.txt`).

## Estructura de archivos generados

```text
analisis_expresion_malaria/
├── results/
│   ├── exp_dif/
│   │   ├── resultados_completos.csv        # Tabla con todos los genes y su categoría
│   │   ├── up_resistant.txt                # Lista de genes sobreexpresados en resistentes
│   │   └── up_susceptible.txt              # Lista de genes sobreexpresados en susceptibles
│   └── images/
│       ├── pca.png                         # (Solo si hay >= 3 muestras totales)
│       ├── volcano.png                     # (Si hubo p-valores / DESeq2)
│       ├── volcano_sin_p.png               # (Si no hubo réplicas / solo log2FC)
│       └── heatmap.png                     # Top 25 genes up y 25 down (Z-score)
└── logs/
    └── differential_expression.log         # Archivo de registro
```

## Pipeline

```mermaid
flowchart TD
    A([raw_counts.txt]) --> B[load_counts]
    B --> C[filter_genes]
    A --> D[load_metadata]
    D -.->|Infiere grupos de la ruta| D1(Metadata)
    
    C --> E{¿Hay >= 3 réplicas\npor grupo?}
    D1 --> E
    
    E -->|Sí| F[DESeq2: run_differential_expression]
    E -->|No| G[Log2FC Directo: run_differential_expression]
    
    F --> H[classify_genes]
    G --> H
    
    H --> I[save_results]
    
    H --> J{¿padj es NaN?}
    J -->|Sí| K[plot_volcano_sin_p]
    J -->|No| L[plot_volcano]
    
    C --> M[plot_heatmap]
    D1 --> M
    
    D1 -.->|Si > 2 muestras| N[plot_pca]
```

## Componentes clave
load_counts: Lee y limpia la matriz de featureCounts, eliminando columnas de coordenadas para quedarse con una matriz de genes × muestras.

infer_conditions / load_metadata: Deduce el grupo experimental (resistant/susceptible) leyendo la ruta del archivo BAM original en el encabezado de featureCounts o mediante un archivo externo (--metadata).

run_differential_expression: Función central. Evalúa el número de muestras. Si hay réplicas suficientes, ejecuta DESeq2. Si no, aplica una normalización CPM y un logaritmo para calcular el log2FoldChange directo.

classify_genes: Etiqueta cada gen como up_resistant, up_susceptible o ns. Si existen p-valores, filtra por padj < 0.05 y |log2FC| >= 1.0. Si no existen, filtra solo por |log2FC| >= 1.0.

plot_pca: Genera un Análisis de Componentes Principales (si hay >2 muestras) para visualizar la varianza entre grupos.

plot_volcano / plot_volcano_sin_p: Genera gráficos de dispersión (Volcano plots). El primero usa los p-valores de DESeq2; el segundo adapta el eje Y a la expresión promedio cuando no hay estadística inferencial.

plot_heatmap: Genera un Mapa de Calor (con Seaborn) de los Top N genes diferenciales. Aplica un Z-score por gen para comparar visualmente patrones de expresión (rojo = alto, azul = bajo).

## Uso por línea de comandos

**Argumentos**
Argumento	Requerido	Default	Descripción
-m, --matrix	No	data/counts/raw_counts.txt	Ruta a la matriz de conteos cruda
-md, --metadata	No	None	Archivo con columnas sample y condition (opcional)
-o, --outdir	No	results/exp_dif	Directorio de salida de las tablas
-c1, --condition1	No	resistant	Condición de prueba (grupo de interés)
-c2, --condition2	No	susceptible	Condición de referencia (grupo base)
-pca, --pcaplot	No	True	Genera el PCA (usa --no-pcaplot para omitir)
-vol, --volcanoplot	No	True	Genera el Volcano Plot (usa --no-volcanoplot para omitir)
-hm, --heatmap	No	True	Genera el Heatmap (usa --no-heatmap para omitir)
-tn, --top_genes	No	25	Número de genes Top (up/down) a graficar en el Heatmap
-v, --verbose	No	False	Muestra mensajes DEBUG
--log_file	No	logs/differential_expression.log	Ruta del archivo de log

## Ejemplos

Ejecución básica (detecta automáticamente si usar DESeq2 o log2FC directo):
```bash
python src/differential_expression.py -v
```

Ejecución personalizando el número de genes para el Heatmap (ej. Top 50):
```bash
python src/differential_expression.py -tn 50 -v
```
