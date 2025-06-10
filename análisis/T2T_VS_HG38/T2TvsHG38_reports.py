"""
Este script recibe un .csv procedente del script T2Tvshg38Compare.

Itera y filtra el dataframe con el objetivo de detectar variantes según diversos criterios de filtrado

Como resultado genera reportes con información general del estado de las variantes T2T y HG38 comparadas.

"""



import pandas as pd
import os
df = pd.read_csv("resultados/variantes_comparadas_twin3.csv", sep=";")
print(f"Número de variantes comparadas en el informe {len(df)}")

result_t2t = pd.read_csv(
    "T2T_converted/twin3/twin3_intervals.bed",
    sep='\t',
    header=None,
    names=['chrom', 'start', 'end', 'seqId', 'score', 'strand'],
    dtype={'chrom': str, 'start': int, 'end': int, 'seqId': str}
)





# --- Posición igual (Diff valor abs == 0) ---
df_pos0 = df[df["Diff valor abs"] == 0]

grupo1_mask = (
    (df_pos0["Clase-T2T"] == df_pos0["Clase-HG38"]) &
    (df_pos0["Subclase-T2T"] == df_pos0["Subclase-HG38"])
)
grupo1 = df_pos0[grupo1_mask]

grupo2_mask = (
    (df_pos0["Clase-T2T"] == df_pos0["Clase-HG38"]) &
    ~grupo1_mask
)
grupo2 = df_pos0[grupo2_mask]

grupo3_mask = (
    (df_pos0["Clase-T2T"] != df_pos0["Clase-HG38"]) &
    ~grupo1_mask &
    ~grupo2_mask
)
grupo3 = df_pos0[grupo3_mask]

# Imprimir resumen de variantes con posición igual
print("\n Variantes con posición igual (Diff valor abs == 0):")
print(f"Grupo 1: Clase y subclase iguales → {len(grupo1)} variantes")
print(f"Grupo 2: Solo clase igual (subclase distinta) → {len(grupo2)} variantes")
print(f"Grupo 3: Clase distinta → {len(grupo3)} variantes")

# --- Posición distinta (Diff valor abs ≠ 0 pero ≤10,000) ---
grupo4 = df[
    (df["Diff valor abs"] != 0) &
    (df["Diff valor abs"].abs() <= 10000)
]


# Subgrupo 4A: misma clase + subclase
subgrupo4A_mask = (
    (grupo4["Clase-T2T"] == grupo4["Clase-HG38"]) &
    (grupo4["Subclase-T2T"] == grupo4["Subclase-HG38"])
)
subgrupo4A = grupo4[subgrupo4A_mask]

# Subgrupo 4B: misma clase, subclase distinta
subgrupo4B_mask = (
    (grupo4["Clase-T2T"] == grupo4["Clase-HG38"]) &
    ~subgrupo4A_mask
)
subgrupo4B = grupo4[subgrupo4B_mask]

# Subgrupo 4C: clase distinta
subgrupo4C_mask = (
    (grupo4["Clase-T2T"] != grupo4["Clase-HG38"]) &
    ~subgrupo4A_mask &
    ~subgrupo4B_mask
)
subgrupo4C = grupo4[subgrupo4C_mask]

print("\n Variantes con posición distinta (|Diff valor abs| ≤ 10,000):")
print(f"Subgrupo 4A: Clase y subclase iguales → {len(subgrupo4A)} variantes")
print(f"Subgrupo 4B: Solo clase igual (subclase distinta) → {len(subgrupo4B)} variantes")
print(f"Subgrupo 4C: Clase distinta → {len(subgrupo4C)} variantes")

# --- Métricas estadísticas de desplazamiento ---
print("\nMétricas de desplazamiento para variantes con posición distinta:")
for nombre, grupo in [
    ("Subgrupo 4A (misma clase + subclase)", subgrupo4A),
    ("Subgrupo 4B (misma clase, subclase distinta)", subgrupo4B),
    ("Subgrupo 4C (clase distinta)", subgrupo4C),
]:
    media = grupo["Diff valor abs"].mean()
    std = grupo["Diff valor abs"].std()
    mediana = grupo["Diff valor abs"].median()
    print(f"{nombre}:")
    print(f"  - Variantes: {len(grupo)}")
    print(f"  - Media: {media:,.2f} pb")
    print(f"  - Desviación estándar: {std:,.2f} pb")
    print(f"  - Mediana: {mediana:,.2f} pb\n")


subclase_cambios = df[
    (df["Diff valor abs"] == 0) &
    (df["Subclase-T2T"] != df["Subclase-HG38"])
]
conteo_transiciones = subclase_cambios.groupby(
    ["Subclase-T2T", "Subclase-HG38"]
).size().reset_index(name="Frecuencia").sort_values("Frecuencia", ascending=False)

print(" Transiciones entre subclases (posición igual):")
print(conteo_transiciones)


clase_cambios = df[
    (df["Diff valor abs"] == 0) &
    (df["Clase-T2T"] != df["Clase-HG38"])
]
conteo_clase = clase_cambios.groupby(
    ["Clase-T2T", "Clase-HG38"]
).size().reset_index(name="Frecuencia").sort_values("Frecuencia", ascending=False)

print("\nTransiciones entre clases (posición igual):")
print(conteo_clase)



df_diferencias = df.drop(index=grupo1.index
                         .union(grupo2.index)
                         .union(grupo3.index)
                         .union(subgrupo4A.index)
                         .union(subgrupo4B.index)
                         .union(subgrupo4C.index))


print(f"Variantes NO clasificadas en ningún grupo: {len(df_diferencias)}")

# Mostrar las primeras filas (para verificar)
media = df_diferencias["Diff valor abs"].mean()
std = df_diferencias["Diff valor abs"].std()
mediana = df_diferencias["Diff valor abs"].median()

print(f"\nMétricas para variantes NO clasificadas:")
print(f"  - Media de desplazamiento: {media:,.2f} pb")
print(f"  - Mediana de desplazamiento: {mediana:,.2f} pb")
print(f"  - Desviación estándar: {std:,.2f} pb")
dif_clase = df_diferencias["Clase-T2T"] != df_diferencias["Clase-HG38"]
porc_dif_clase = dif_clase.mean() * 100
print(f" Variantes con CLASE distinta: {porc_dif_clase:.2f}%")
dif_subclase = df_diferencias["Subclase-T2T"] != df_diferencias["Subclase-HG38"]
porc_dif_subclase = dif_subclase.mean() * 100
print(f"Variantes con SUBCLASE distinta: {porc_dif_subclase:.2f}%")


bed_path = "T2T_converted/twin3/twin3_intervals.bed"

result_t2t = pd.read_csv(
    bed_path,
    sep='\t',
    header=None,
    names=['chrom', 'start', 'end', 'seqId'],  # Solo las columnas presentes
    dtype={'chrom': str, 'start': int, 'end': int, 'seqId': str}
)
print(len(result_t2t))

unique_insertions = result_t2t.drop_duplicates(subset=['chrom', 'start', 'end'])


total_unique = len(unique_insertions)


chromosomes_clasicos = {f"chr{i}" for i in range(1, 23)}.union({"chrX", "chrY"})


alt_contigs = unique_insertions[~unique_insertions["chrom"].isin(chromosomes_clasicos)]
num_alt_contigs = len(alt_contigs)


classic_chroms = unique_insertions[unique_insertions["chrom"].isin(chromosomes_clasicos)]
num_classic = len(classic_chroms)

# Mostrar resultados
print(f"Total de inserciones únicas: {total_unique}")
print(f"Contigs alternativos detectados: {num_alt_contigs}")
print(f"Inserciones en cromosomas clásicos: {num_classic}")
print("\nContigs alternativos únicos encontrados:")
print(alt_contigs['chrom'].unique())