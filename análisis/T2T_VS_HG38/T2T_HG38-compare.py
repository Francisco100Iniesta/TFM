"""
Este script recibe como argumentos 3 ficheros:
 -[Variantes de retroinspector Genoma version hg38]
 -[Variantes de retroinspector Genoma T2t]
 -[Bed de intervalos con posiciones cromosómicas mapeadas del Genoma T2T al hg38]

 Itera por id buscando en Variantes T2T variantes del BED que se pudieron mapear a hg38. Las variantes coincidentes por Id, heredan la nueva posicion de Start y End.
 Una vez generado la lista de variantes T2T mapeadas, se evalua cada una de forma individual sobre la lista de variantes Hg38 bajo dos criterios:
 *Coincidencia de cromosoma
 *Localizar una inserción coincidente o más cercana en hg38

 Como resultado vuelve un .csv donde se compara cada posición del T2T respecto hg38, indicando valores clave como:
 -distancia de solapamiento
 -clase
 -subclase
 -tamaño del elemento
 -Muestra en la que esta presente
 -Id de los informes de retroinspector
TODOS ESTOS CAMPOS SE REPORTAN PARA HG38 Y T2T


"""
import os
# Crear la carpeta 'resultados' si no existe
os.makedirs("resultados", exist_ok=True)
import pandas as pd
ruta_hg38= "twinsShare/HG38/twin3.csv"
ruta_t2t="twinsShare/T2T/twin3.csv"

result_t2t = pd.read_csv(
    "T2T_converted/twin3/twin3_intervals.bed",
    sep='\t',
    header=None,
    names=['chrom', 'start', 'end', 'seqId', 'score', 'strand'],
    dtype={'chrom': str, 'start': int, 'end': int, 'seqId': str}
)
#nombre= os.path.splitext(os.path.basename(ruta))[0]
dataframe_hg38 = pd.read_csv(ruta_hg38, sep=";")
dataframe_t2t = pd.read_csv(ruta_t2t, sep=";")
dataframe_hg38["Position"] = dataframe_hg38["Position"].str.replace(r'\s+', '', regex=True)
dataframe_t2t["Position"] = dataframe_t2t["Position"].str.replace(r'\s+', '', regex=True)
dataframe_t2t["Length(bp)"] = dataframe_t2t["Length (bp)"].str.replace(r'\s+', '', regex=True)
dataframe_hg38["Length(bp)"] = dataframe_hg38["Length (bp)"].str.replace(r'\s+', '', regex=True)
cols_variantes = ['Position', 'seqId', 'Length (bp)', 'Type']
dataframe_t2t = dataframe_t2t.drop_duplicates(subset=cols_variantes)
dataframe_hg38 = dataframe_hg38.drop_duplicates(subset=cols_variantes)

dataframe_t2t[['Class', 'Subclass']] = dataframe_t2t['Type'] \
        .str.extract(r'^(?P<Class>.*?) \((?P<Subclass>.*?)\)$')
dataframe_hg38[['Class', 'Subclass']] = dataframe_hg38['Type'] \
        .str.extract(r'^(?P<Class>.*?) \((?P<Subclass>.*?)\)$')
print(len(dataframe_hg38))
print(len(dataframe_t2t))
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


def comparador_posiciones(end_t2t, end_hg38):
    diferencia=end_t2t-end_hg38
    aguas=""
    if diferencia <0:
        aguas="-"
        verificacion=""
    elif diferencia==0:
        verificacion="Insercción Identica"
    else:
        aguas="+"
        verificacion = ""
    comparacion=(f"{verificacion} Diferencia respecto hg38: {aguas}{abs(diferencia)}")
    return comparacion,abs(diferencia),diferencia

def construir_tablas(dataframe1, dataframe2, resultados):
    seleccion = []
    for i, fila1 in dataframe1.iterrows():
        id1 = fila1.iloc[1]
        datos = None
        for j, res in resultados.iterrows():
            if res.iloc[3] == id1:
                datos = {"Cromosoma": res.iloc[0], "Posicion": res.iloc[2]}
                break
        if datos is None:
            continue
        cromosoma1 = datos["Cromosoma"]
        start1 = int(datos["Posicion"])
        longitud1 = int(str(fila1.iloc[2]).replace('\u202f', '').replace(" ", "").replace(",", ""))
        tipo1 = fila1.iloc[3]
        paciente1 = fila1.iloc[5]
        clase1=fila1.iloc[7]
        subclase1=fila1.iloc[8]
        diferencia_min = float('inf')
        mejor_candidato = None
        for j, fila2 in dataframe2.iterrows():
            crom2, start2_str = fila2.iloc[0].split(':')
            start2 = int(start2_str)
            if crom2 != cromosoma1:
                continue
            id2 = fila2.iloc[1]
            longitud2 = int(str(fila2.iloc[2]).replace('\u202f', '').replace(" ", "").replace(",", ""))
            tipo2 = fila2.iloc[3]
            paciente2 = fila2.iloc[5]
            clase2 = fila2.iloc[7]
            subclase2 = fila2.iloc[8]
            etiqueta, diff_pos,posicion = comparador_posiciones(start2, start1)
            if diff_pos < diferencia_min:
                diferencia_min = diff_pos
                mejor_candidato = {
                    "Cromosoma": cromosoma1,
                    "Posicion-T2T": start1,
                    "Posicion-HG38": start2,
                    "Longitud-T2T": longitud1,
                    "Longitud2": longitud2,
                    "Diferencia_tamaño": abs(longitud2 - longitud1),
                    "Clase-T2T": clase1,
                    "Clase-HG38": clase2,
                    "Subclase-T2T":subclase1,
                    "Subclase-HG38": subclase2,
                    "id-T2T": id1,
                    "id-HG38": id2,
                    "Variación": etiqueta,
                    "Diff valor abs":posicion,
                    "paciente-T2T": paciente1,
                    "paciente-HG38": paciente2
                }
        if mejor_candidato:
            seleccion.append(mejor_candidato)
    return seleccion


a=construir_tablas(dataframe_t2t,dataframe_hg38,result_t2t)
# Convertir a DataFrame
df_a = pd.DataFrame(a)
# Guardar a CSV
df_a.to_csv("resultados/variantes_comparadas_twin3.csv", sep=";", index=False)
