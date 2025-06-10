"""
Este script recibe archivos .csv con información generada por el script ShareCompareStats.
Específicamente recibe: media de cromosomas,clases y subclases. PARA LAS VERSIONES T2T y HG38

Realiza estadistica comparativa entre las medias T2T/ HG38, y reporta información comparativa y de significancia.
"""
from scipy.stats import wilcoxon, spearmanr, pearsonr
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)
ruta1= "shared/Cromosomas_-_Shared_h38.csv"
ruta2= "shared/Cromosomas_-_Shared_T2T.csv"
ruta3= "shared/Clases_-_Shared_h38.csv"
ruta4= "shared/Clases_-_Shared_T2T.csv"
ruta5= "shared/Subclases_-_Shared_h38.csv"
ruta6= "shared/Subclases_-_Shared_T2T.csv"

dataframeCHR_hg38 = pd.read_csv(ruta1, sep=",")
dataframeCHR_T2T= pd.read_csv(ruta2, sep=",")
dataframeClass_hg38= pd.read_csv(ruta3, sep=",")
dataframeClass_T2T= pd.read_csv(ruta4, sep=",")
dataframeSUBCLASS_hg38= pd.read_csv(ruta5, sep=",")
dataframeSUBCLASS_T2T= pd.read_csv(ruta6, sep=",")

print(dataframeCHR_hg38)

merged_chr = dataframeCHR_hg38.merge(
    dataframeCHR_T2T,
    on="chr",
    suffixes=('_hg38', '_T2T')
)
print(len(merged_chr))
# Test de Wilcoxon pareado
stat, p = wilcoxon(merged_chr["Media_por_chrom_hg38"], merged_chr["Media_por_chrom_T2T"])
print(f"Wilcoxon pareado - Media por cromosoma\nEstadístico: {stat}, p-valor: {p:.4f}")

# Correlación
pearson_corr, _ = pearsonr(merged_chr["Media_por_chrom_hg38"], merged_chr["Media_por_chrom_T2T"])
spearman_corr, _ = spearmanr(merged_chr["Media_por_chrom_hg38"], merged_chr["Media_por_chrom_T2T"])
print(f"Pearson: {pearson_corr:.3f} | Spearman: {spearman_corr:.3f}")

# Gráfico scatter hg38 vs T2T
plt.figure(figsize=(8,6))
plt.scatter(merged_chr["Media_por_chrom_hg38"], merged_chr["Media_por_chrom_T2T"])
plt.plot([merged_chr["Media_por_chrom_hg38"].min(), merged_chr["Media_por_chrom_hg38"].max()],
         [merged_chr["Media_por_chrom_hg38"].min(), merged_chr["Media_por_chrom_hg38"].max()],
         linestyle='--', color='gray')
plt.xlabel("hg38 - Media por cromosoma")
plt.ylabel("T2T - Media por cromosoma")
plt.title("Comparación de medias por cromosoma (hg38 vs T2T)")
plt.grid(True)
plt.tight_layout()
plt.show()


deltas = merged_chr["Media_por_chrom_hg38"] - merged_chr["Media_por_chrom_T2T"]
from scipy.stats import shapiro
stat, p = shapiro(deltas)
print(f"Shapiro-Wilk para diferencias (hg38 - T2T): estadístico={stat:.3f}, p-valor={p:.4f}")

from scipy.stats import ttest_rel

stat, p = ttest_rel(
    merged_chr["Media_por_chrom_hg38"],
    merged_chr["Media_por_chrom_T2T"]
)

print(f"t-test pareado: estadístico={stat:.3f}, p-valor={p:.4f}")


# Diferencias entre hg38 y T2T
merged_chr["diferencia"] = merged_chr["Media_por_chrom_hg38"] - merged_chr["Media_por_chrom_T2T"]

# 1. Media y desviación de las diferencias
media_diff = merged_chr["diferencia"].mean()
std_diff = merged_chr["diferencia"].std()

# 2. Cohen's d
cohen_d = media_diff / std_diff

print(f"Diferencia media: {media_diff:.2f}")
print(f"Desviación estándar de diferencias: {std_diff:.2f}")
print(f"Cohen's d: {cohen_d:.2f}")


"""
Diferencias entre clases
"""

merged_class = dataframeClass_hg38.merge(
    dataframeClass_T2T,
    on="Class",
    suffixes=('_hg38', '_T2T')
)

# Test de Wilcoxon
from scipy.stats import wilcoxon, ttest_rel, shapiro, pearsonr, spearmanr

stat, p = wilcoxon(merged_class["Media_por_class_hg38"], merged_class["Media_por_class_T2T"])
print(f"Wilcoxon pareado - Media por clase\nEstadístico: {stat}, p-valor: {p:.4f}")

# Correlaciones
pearson_corr, _ = pearsonr(merged_class["Media_por_class_hg38"], merged_class["Media_por_class_T2T"])
spearman_corr, _ = spearmanr(merged_class["Media_por_class_hg38"], merged_class["Media_por_class_T2T"])
print(f"Pearson: {pearson_corr:.3f} | Spearman: {spearman_corr:.3f}")

# Normalidad de las diferencias
deltas = merged_class["Media_por_class_hg38"] - merged_class["Media_por_class_T2T"]
stat, p = shapiro(deltas)
print(f"Shapiro-Wilk para diferencias (hg38 - T2T): estadístico={stat:.3f}, p-valor={p:.4f}")

# t-test pareado
stat, p = ttest_rel(
    merged_class["Media_por_class_hg38"],
    merged_class["Media_por_class_T2T"]
)
print(f"t-test pareado: estadístico={stat:.3f}, p-valor={p:.4f}")

# Cohen's d
merged_class["diferencia"] = deltas
media_diff = deltas.mean()
std_diff = deltas.std()
cohen_d = media_diff / std_diff
print(f"Diferencia media: {media_diff:.2f}")
print(f"Desviación estándar de diferencias: {std_diff:.2f}")
print(f"Cohen's d: {cohen_d:.2f}")



"""
Diferencias entre SUBCLASES
"""

merged_subclass = dataframeSUBCLASS_hg38.merge(
    dataframeSUBCLASS_T2T,
    on="Subclass",
    suffixes=('_hg38', '_T2T')
)
stat, p = wilcoxon(merged_subclass["Media_por_sub_hg38"], merged_subclass["Media_por_sub_T2T"])
print(f"Wilcoxon pareado - Media por subclase\nEstadístico: {stat}, p-valor: {p:.4f}")
from scipy.stats import pearsonr, spearmanr

pearson_corr, _ = pearsonr(merged_subclass["Media_por_sub_hg38"], merged_subclass["Media_por_sub_T2T"])
spearman_corr, _ = spearmanr(merged_subclass["Media_por_sub_hg38"], merged_subclass["Media_por_sub_T2T"])
print(f"Pearson: {pearson_corr:.3f} | Spearman: {spearman_corr:.3f}")
deltas = merged_subclass["Media_por_sub_hg38"] - merged_subclass["Media_por_sub_T2T"]
stat, p = shapiro(deltas)
print(f"Shapiro-Wilk para diferencias (hg38 - T2T): estadístico={stat:.3f}, p-valor={p:.4f}")
tat, p = ttest_rel(
    merged_subclass["Media_por_sub_hg38"],
    merged_subclass["Media_por_sub_T2T"]
)
print(f"t-test pareado: estadístico={stat:.3f}, p-valor={p:.4f}")

merged_subclass["diferencia"] = merged_subclass["Media_por_sub_hg38"] - merged_subclass["Media_por_sub_T2T"]
media_diff = merged_subclass["diferencia"].mean()
std_diff = merged_subclass["diferencia"].std()
cohen_d = media_diff / std_diff

print(f"Diferencia media: {media_diff:.2f}")
print(f"Desviación estándar de diferencias: {std_diff:.2f}")
print(f"Cohen's d: {cohen_d:.2f}")


subclasses_hg38 = set(dataframeSUBCLASS_hg38["Subclass"])
subclasses_T2T = set(dataframeSUBCLASS_T2T["Subclass"])

hg38_only = subclasses_hg38 - subclasses_T2T
T2T_only = subclasses_T2T - subclasses_hg38
common_subclasses = subclasses_hg38 & subclasses_T2T

# Cálculos
num_hg38_only = len(hg38_only)
num_T2T_only = len(T2T_only)
num_common = len(common_subclasses)
total_hg38 = len(dataframeSUBCLASS_hg38)
total_T2T = len(dataframeSUBCLASS_T2T)
total_unique = len(subclasses_hg38 | subclasses_T2T)

# Mostrar resultados
print("Subclasses sólo en hg38:", num_hg38_only, hg38_only)
print("Subclasses sólo en T2T:", num_T2T_only, T2T_only)
print("Subclasses en común:", num_common, common_subclasses)
print("Total filas hg38:", total_hg38)
print("Total filas T2T:", total_T2T)
print("Total únicas (hg38 ∪ T2T):", total_unique)


merged_subclass = dataframeSUBCLASS_hg38.merge(
    dataframeSUBCLASS_T2T,
    on="Subclass",
    suffixes=('_hg38', '_T2T')
)
print("Total filas mergeadas (elementos a comparar):", len(merged_subclass))