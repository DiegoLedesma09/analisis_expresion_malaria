import argparse as ap

import mysql.connector
import pandas as pd


def crear_base(host, user, password) -> None:
    """
    Creación de la base de datos del archivo de enriquecimiento funcional creado por
    functional_annotation.py

    Args:
        host (str) -> ID del local host de MySQL
        user (str) -> Usuario de MySQL
        password (str) -> Contraseña del usuario de MySQL

    """
    # 1. Leer el archivo TSV
    archivo_tsv = "results/go/go_enrichment_results.tsv"
    df = pd.read_csv(archivo_tsv, sep="\t")

    # Asegurar que 'significativo' sea booleano/entero (0 o 1) para MySQL
    df["significativo"] = df["significativo"].astype(bool).astype(int)

    # 2. Conectar al servidor MySQL (sin especificar BD aún)
    conexion = mysql.connector.connect(
        host= host,
        user= user,         
        password= password 
    )
    cursor = conexion.cursor()

    # 3. Crear la base de datos y seleccionarla
    cursor.execute("CREATE DATABASE IF NOT EXISTS go_enrichment_db CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;")
    cursor.execute("USE go_enrichment_db;")

    # 4. Crear la tabla
    sql_crear_tabla = """
    CREATE TABLE IF NOT EXISTS enriquecimiento_go (
        id INT AUTO_INCREMENT PRIMARY KEY,
        grupo VARCHAR(50) NOT NULL,
        go_id VARCHAR(20) NOT NULL,
        aspecto VARCHAR(10) NOT NULL,
        termino VARCHAR(255) NOT NULL,
        k_estudio INT NOT NULL,
        n_estudio INT NOT NULL,
        K_universo INT NOT NULL,
        N_universo INT NOT NULL,
        gene_ratio DOUBLE NOT NULL,
        p_value DOUBLE NOT NULL,
        padj DOUBLE NOT NULL,
        significativo BOOLEAN NOT NULL,
        profundidad INT NOT NULL,
        genes TEXT NOT NULL,
        INDEX idx_grupo (grupo),
        INDEX idx_go_id (go_id),
        INDEX idx_significativo (significativo)
    );
    """
    cursor.execute(sql_crear_tabla)

    # 5. Insertar los datos fila por fila
    sql_insertar = """
    INSERT INTO enriquecimiento_go (
        grupo, go_id, aspecto, termino, k_estudio, n_estudio,
        K_universo, N_universo, gene_ratio, p_value, padj,
        significativo, profundidad, genes
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    columnas = [
        "grupo", "go_id", "aspecto", "termino", "k_estudio", "n_estudio",
        "K_universo", "N_universo", "gene_ratio", "p_value", "padj",
        "significativo", "profundidad", "genes"
    ]

    # Reemplazar posibles NaN por None para evitar errores en MySQL
    datos = [
        tuple(None if pd.isna(val) else val for val in fila)
        for fila in df[columnas].itertuples(index=False, name=None)
    ]

    cursor.executemany(sql_insertar, datos)
    conexion.commit()

    print(f"Listo: se insertaron {cursor.rowcount} registros en 'go_enrichment_db.enriquecimiento_go'.")

    # 6. Cerrar conexión
    cursor.close()
    conexion.close()

def my_parser() -> ap.Namespace:
    """
    Retiene los argumentos deseados para el módulo desde línea de comandos.

    Args:
        None
    
    Returns:
        ap.Namespace: Argumentos indicados por el usuario
    """    
    parser = ap.ArgumentParser(
    description="Creación y retención de la base de datos obtenida del pipeline MOSQUITO DGE (Opcional)")
    parser.add_argument("-H", "--host", required=True, help="ID del local host de MySQL")
    parser.add_argument("-u", "--user", required=True, help="Nombre del usuario")
    parser.add_argument("-p", "--password", required=True, help="Contraseña del usuario")
    
    return parser.parse_args()

def main():
    """
    Ejecuta la petición de argumentos y el procedimiento para la creación de la base de datos
    asociada a los resultados obtenidos por el módulo functional_annotation.py
    """
    args = my_parser()

    crear_base(args.host, args.user, args.password)    

if __name__ == "__main__":
    main()
