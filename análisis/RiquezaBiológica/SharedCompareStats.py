import pandas as pd
"""
Este script analiza métricas sobre transposición a partir de ficheros csv creados a partir de los informes de retroinspector.
Elimina duplicados, creando una base de datos con toda la transposición compartida entre hermanos, lo que permite sacar métricas globales
sobre eventos representados en la transposición. 
"""
def chrT2TtoHG38(cromosoma):
    cp_to_chr = {
        "CP068277.2": "chr1",
        "CP068276.2": "chr2",
        "CP068275.2": "chr3",
        "CP068274.2": "chr4",
        "CP068273.2": "chr5",
        "CP068272.2": "chr6",
        "CP068271.2": "chr7",
        "CP068255.2": "chrX",
        "CP068269.2": "chr9",
        "CP068270.2": "chr8",
        "CP068267.2": "chr11",
        "CP068268.2": "chr10",
        "CP068266.2": "chr12",
        "CP068265.2": "chr13",
        "CP068264.2": "chr14",
        "CP068263.2": "chr15",
        "CP068262.2": "chr16",
        "CP068261.2": "chr17",
        "CP068260.2": "chr18",
        "CP068258.2": "chr20",
        "CP086569.2": "chrY",
        "CP068259.2": "chr19",
        "CP068256.2": "chr22",
        "CP068257.2": "chr21",
        "CP068254.1": "chrM"
    }
    return cp_to_chr[cromosoma]
ruta1= "twinsSHARE/HG38/twin1.csv"
ruta2= "twinsSHARE/HG38/twin2.csv"
ruta3= "twinsSHARE/HG38/twin3.csv"
dataframe1 = pd.read_csv(ruta1, sep=";")
dataframe2 = pd.read_csv(ruta2, sep=";")
dataframe3 = pd.read_csv(ruta3, sep=";")

dataframe1["Position"] = dataframe1["Position"].str.replace(r'\s+', '', regex=True)
dataframe1["Length(bp)"] = dataframe1["Length (bp)"].str.replace(r'\s+', '', regex=True)

dataframe2["Position"] = dataframe2["Position"].str.replace(r'\s+', '', regex=True)
dataframe2["Length(bp)"] = dataframe2["Length (bp)"].str.replace(r'\s+', '', regex=True)

dataframe3["Position"] = dataframe3["Position"].str.replace(r'\s+', '', regex=True)
dataframe3["Length(bp)"] = dataframe3["Length (bp)"].str.replace(r'\s+', '', regex=True)



dataframe1[['chr', 'pos_str']] = dataframe1['Position'].str.split(':', expand=True)
dataframe1['pos'] = dataframe1['pos_str'].astype(int)
dataframe1[['Class', 'Subclass']] =dataframe1['Type'].str.extract(r'^(?P<Class>.*?) \((?P<Subclass>.*?)\)$')
dataframe2[['chr', 'pos_str']] = dataframe2['Position'].str.split(':', expand=True)
dataframe2['pos'] = dataframe2['pos_str'].astype(int)
dataframe2[['Class', 'Subclass']] =dataframe2['Type'].str.extract(r'^(?P<Class>.*?) \((?P<Subclass>.*?)\)$')
dataframe3[['chr', 'pos_str']] = dataframe3['Position'].str.split(':', expand=True)
dataframe3['pos'] = dataframe3['pos_str'].astype(int)
dataframe3[['Class', 'Subclass']] =dataframe3['Type'].str.extract(r'^(?P<Class>.*?) \((?P<Subclass>.*?)\)$')

"""
Estas lineas de codigo desactivadas con # se utilizan para cuando las muestras procedan de T2T, cambiar el nombre de los cromosomas, de esta forma cuando se comparen en los reportes
hay consistencia en la nomenclatura.
"""
#dataframe1['chr'] = dataframe1['chr'].map(chrT2TtoHG38)
#dataframe2['chr'] = dataframe2['chr'].map(chrT2TtoHG38)
#dataframe3['chr'] = dataframe3['chr'].map(chrT2TtoHG38)
dataframe1['Pareja_ID'] = 'Pareja1'
dataframe2['Pareja_ID'] = 'Pareja2'
dataframe3['Pareja_ID'] = 'Pareja3'

dataframe1['locus_key'] = (
    dataframe1['chr'].astype(str) + ':' +
    dataframe1['pos'].astype(str) + '|' +
    dataframe1['Class'] + '|' +
    dataframe1['Subclass'] + '|' +
    dataframe1['Length (bp)'].astype(str)
)


grupo = (
    dataframe1[['locus_key', 'Present in sample']]
      .drop_duplicates()
      .groupby('locus_key')['Present in sample']
      .agg(lambda x: list(x))
      .reset_index()
      .rename(columns={'Present in sample': 'Gemelos_con_locus'})
)


grupo['Estado'] = grupo['Gemelos_con_locus'].apply(
    lambda x: 'Compartido' if len(x) > 1 else 'Único'
)

total_loci = grupo.shape[0]
loci_compartidos = grupo[grupo['Estado'] == 'Compartido'].shape[0]
loci_unicos = grupo[grupo['Estado'] == 'Único'].shape[0]


nombres_gemelos = dataframe1['Present in sample'].unique()
gemelo1, gemelo2 = nombres_gemelos[0], nombres_gemelos[1]  # se asume orden consistente

unicos_gemelo1 = grupo[
    (grupo['Estado'] == 'Único') &
    (grupo['Gemelos_con_locus'].apply(lambda x: x == [gemelo1]))
].shape[0]

unicos_gemelo2 = grupo[
    (grupo['Estado'] == 'Único') &
    (grupo['Gemelos_con_locus'].apply(lambda x: x == [gemelo2]))
].shape[0]


conteos_raw = dataframe1['locus_key'].value_counts().reset_index()
conteos_raw.columns = ['locus_key', 'veces_en_raw']
duplicados_raw = conteos_raw[conteos_raw['veces_en_raw'] > 1]


print(f"Total de loci distintos en la pareja: {total_loci}")
print(f"  · De ellos, {loci_compartidos} loci están compartidos por ambos gemelos.")
print(f"  · Y {loci_unicos} loci aparecen solo en uno de los gemelos:")
print(f"      - Únicos de {gemelo1}: {unicos_gemelo1}")
print(f"      - Únicos de {gemelo2}: {unicos_gemelo2}")

print()
print("Detalle de duplicados en el dataframe original (raw):")
print(duplicados_raw.head(10))


for df in [dataframe1, dataframe2, dataframe3]:
    df['locus_key'] = (
        df['chr'].astype(str) + ':' +
        df['pos'].astype(str) + '|' +
        df['Class'] + '|' +
        df['Subclass'] + '|' +
        df['Length (bp)'].astype(str)
    )

# Dado que los informes de retroinspector reportan variantes comunes como duplicadas es necesario filtrar
dataframe1 = dataframe1.drop_duplicates(subset='locus_key')
print(len(dataframe1))
dataframe2 = dataframe2.drop_duplicates(subset='locus_key')
dataframe3 = dataframe3.drop_duplicates(subset='locus_key')

df_all = pd.concat([dataframe1, dataframe2, dataframe3], axis=0, ignore_index=True)
conteo_por_chr = df_all.groupby(['Pareja_ID','chr']).size().reset_index(name='count')
matriz_chr = conteo_por_chr.pivot(index='Pareja_ID', columns='chr', values='count').fillna(0)

print(matriz_chr)

stats_chr = pd.DataFrame({
    'Media_por_chrom':   matriz_chr.mean(axis=0),
    'DesvEst_por_chrom': matriz_chr.std(axis=0),
    'NumParejas_observ': matriz_chr.count(axis=0)  # cuántas parejas aportaron valor >0 en cada crom
})

print(stats_chr.sort_values(by='Media_por_chrom', ascending=False))




conteo_por_class = df_all.groupby(['Pareja_ID','Class']).size().reset_index(name='count')


matriz_class = conteo_por_class.pivot(index='Pareja_ID', columns='Class', values='count').fillna(0)

print(matriz_class)

stats_class = pd.DataFrame({
    'Media_por_class':   matriz_class.mean(axis=0),
    'DesvEst_por_class': matriz_class.std(axis=0),
    'NumParejas_observ': matriz_class.count(axis=0)
})

print(stats_class.sort_values(by='Media_por_class', ascending=False))




conteo_por_sub = df_all.groupby(['Pareja_ID','Subclass']).size().reset_index(name='count')


matriz_sub = conteo_por_sub.pivot(index='Pareja_ID', columns='Subclass', values='count').fillna(0)

print(matriz_sub.head(10))


stats_sub = pd.DataFrame({
    'Media_por_sub':   matriz_sub.mean(axis=0),
    'DesvEst_por_sub': matriz_sub.std(axis=0),
    'NumParejas_observ': matriz_sub.count(axis=0)
})


print(stats_sub.sort_values(by='Media_por_sub', ascending=False))




df_all['Locus_ID'] = df_all['chr'].astype(str) + ':' + df_all['pos'].astype(str)

locus_counts = df_all.groupby('Locus_ID')['Pareja_ID'].nunique().reset_index(name='Parejas_con_locus')


recurrentes = locus_counts[locus_counts['Parejas_con_locus'] > 1]
print(f"Total de loci diferenciales únicos: {locus_counts.shape[0]}")
print(f"Loci que aparecen en ≥2 parejas: {recurrentes.shape[0]}")


replicas = df_all[['Locus_ID','Pareja_ID']].drop_duplicates()

locus_parejas = replicas.groupby('Locus_ID')['Pareja_ID'].apply(lambda x: ','.join(x.unique())).reset_index()
locus_parejas = locus_parejas.merge(locus_counts, on='Locus_ID')


locus_recurrentes = locus_parejas[locus_parejas['Parejas_con_locus'] > 1]
print(locus_recurrentes.head(10))




reps = df_all[['Locus_ID','Pareja_ID']].drop_duplicates()


mat_pareja_locus = pd.crosstab(reps['Locus_ID'], reps['Pareja_ID'])


print(mat_pareja_locus.head(10))