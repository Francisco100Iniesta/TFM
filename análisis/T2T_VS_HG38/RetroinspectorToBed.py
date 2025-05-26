import pandas as pd
import os
ruta = "RetroinspectorAnnotation/t2t.csv"
nombre= os.path.splitext(os.path.basename(ruta))[0]
dataframe = pd.read_csv(ruta, sep=";")
"""
Este script recibe un csv con separador; construido a partir de copiar y pegar anotaciones de retroinspector 
 construyendo un bed de posiciones[start=inserccion-1][end=inserccion], id del propio informe y sustituye los cromosomas 
 del ensamblaje t2t por los clásicos de los contigs tradicionales para el hg38.
"""
dataframe["Position"] = dataframe["Position"].str.replace(r'\s+', '', regex=True)
with open(f"{nombre}.bed","w")as archivo:
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

    for i in dataframe.iterrows():
        cromosoma=(i[1][0])
        cromosoma,start=cromosoma.split(':')
        end = int(start)
        start=int(start)-1
        id=(i[1][1])
        archivo.write(f"{cp_to_chr[cromosoma]}\t{start}\t{end}\t{id}\n")

