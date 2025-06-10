
library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

#Leo directamente del .tsv que genera methylartist
data <- read.table("segmeth_twin2.tsv",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#  Busco que mínimo que mis islas tengan 10 reads de soporte
cov1 <- data$alignments_LSA.2024.0039_GS1375.I.Sangre_m_readcount
cov2 <- data$alignments_LSA.2024.0040_GS1375.H.Sangre_m_readcount
keep <- cov1 >= 10 & cov2 >= 10
data_filt <- data[keep, ]

#Además puedo extraer un bed de coordenadas de islas cpg para intentar realizar enriquecimiento funcional.
data_filt$score <- "."
bed_df <- data_filt[, c("seg_chrom", "seg_start", "seg_end", "seg_name", "score", "seg_strand")]

# Guardar como BED
write.table(
  bed_df,
  file = "cpg_diff.bed",
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

# 3) Busco la columna que me dice la fracción de metilación para la marca (5mc o 5hmc)
twin1 <- data_filt$alignments_LSA.2024.0039_GS1375.I.Sangre_m_methfrac
twin2 <- data_filt$alignments_LSA.2024.0040_GS1375.H.Sangre_m_methfrac
data_filt$diff  <- twin1 - twin2
diff  <- twin1 - twin2
data_filt$twin1<-data_filt$alignments_LSA.2024.0039_GS1375.I.Sangre_m_methfrac
data_filt$twin2<-data_filt$alignments_LSA.2024.0040_GS1375.H.Sangre_m_methfrac

#creo un dataset filtrado por cromosoma para realizar la comparación enter muestras
wilcox_por_chr <- data_filt %>%
  group_by(seg_chrom) %>%
  group_split() %>%
  map_df(function(df_chr) {
    chr <- unique(df_chr$seg_chrom)
    
  
# Wilcoxon por parejas
    test <- wilcox.test(df_chr$twin1, df_chr$twin2, paired = TRUE)
    
    data.frame(
      seg_chrom = chr,
      p_value = test$p.value
    )
  })

# Ajuste por múltiples pruebas (FDR Benjamini-Hochberg)
wilcox_por_chr$p_adj <- p.adjust(wilcox_por_chr$p_value, method = "BH")

wilcox_por_chr <- wilcox_por_chr %>% arrange(p_adj)
signif <- wilcox_por_chr %>% filter(p_adj < 0.05)

print(wilcox_por_chr)
#Representamos visualemente la significancia
library(ggplot2)

ggplot(wilcox_por_chr, aes(x = reorder(seg_chrom, p_adj), y = -log10(p_adj))) +
  geom_bar(stat = "identity", fill = ifelse(wilcox_por_chr$p_adj < 0.05, "red", "grey")) +
  coord_flip() +
  labs(
    title = "Significancia por cromosoma (Wilcoxon, FDR ajustado)",
    x = "Cromosoma",
    y = expression(-log[10](p[adj]))
  ) +
  theme_minimal()

cromosomas_significativos <- wilcox_por_chr %>%
  filter(p_adj < 0.05) %>%
  pull(seg_chrom)

data_signif <- data_filt %>%
  filter(seg_chrom %in% cromosomas_significativos)

#Filtramos un nuevo dataset con los cromosomas diferenciales par plotear de forma comparativa las diferencias en sus medias 
medias_por_chr <- data_signif %>%
  pivot_longer(cols = c(twin1, twin2), names_to = "twin", values_to = "methfrac") %>%
  group_by(seg_chrom, twin) %>%
  summarise(mean_meth = mean(methfrac), .groups = "drop") %>%
  pivot_wider(names_from = twin, values_from = mean_meth) %>%
  mutate(
    perc_diff = 100 * (twin1 - twin2) / ((twin1 + twin2)/2)
  )
medias_long <- medias_por_chr %>%
  pivot_longer(cols = c(twin1, twin2), names_to = "twin", values_to = "mean_meth")
ggplot(medias_long, aes(x = seg_chrom, y = mean_meth, fill = twin)) +
  geom_col(position = "dodge") +
  geom_text(
    data = medias_por_chr,
    aes(x = seg_chrom, y = pmax(twin1, twin2) + 0.01, label = paste0(round(perc_diff, 2), "%")),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0
  ) +
  labs(
    title = "Fracción de metilación por cromosoma (media por gemelo)",
    subtitle = "Con diferencia porcentual entre twin1:GEMELOS_21.061 y twin2:GEMELOS_21.065",
    x = "Cromosoma", y = "Media de metilación", fill = "Gemelo"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# 4) Para ver diferencias reales selecciono un umbral del 20% como mínimo, en diferencias para su fraccion de metilación
umbral <- 0.2
sel    <- abs(diff) >= umbral
data_sig <- data_filt[sel, ]

# 5) Construyo un data frame para guardar la información de interés y poder representarla
resultado <- data.frame(
  seg_id   = data_sig$seg_id,
  chrom    = data_sig$seg_chrom,
  start    = data_sig$seg_start,
  end      = data_sig$seg_end,
  name     = data_sig$seg_name,
  meth_t1  = data_sig$alignments_LSA.2024.0039_GS1375.I.Sangre_m_methfrac,
  meth_t2  = data_sig$alignments_LSA.2024.0040_GS1375.H.Sangre_m_methfrac,
  diff     = diff[sel]
)
#Dado que quiero analizar si las islas diferenciales presentan inserciones diferenciales de elementos transponibles guardo el dataframe.
write.csv(resultado, file = "islacpg-Retro/CpGDiferenciales/CpG_diferencialesTwin1.csv", row.names = FALSE)

#Voy a  plotear con un barplot la distribución de islas CpG diferenciales por cromosoma
counts <- table(resultado$chrom)
counts <- sort(counts, decreasing = TRUE)

pct <- counts / sum(counts) * 100
cols <- heat.colors(length(counts))
bp <- barplot(
  counts,
  col         = cols,
  main        = "Distribución de metilaciones diferenciales por cromosoma",
  xlab        = "Cromosomas",
  ylab        = "Número de diferencias filtradas",
  las         = 1,
  cex.names   = 0.9,
  cex.axis    = 0.9,

)
legend(
  "topright", 
  legend = paste0(round(pct, 1), "%"),   # ej. "12.3%", "8.7%", ...
  fill   = cols,
  cex    = 0.8,
  bty    = "n"   # sin cuadro alrededor
)


# Tambien verifico como se reparten mis diferencias mediante boxplot
muestra1 <- resultado$meth_t1
muestra2 <- resultado$meth_t2


# Hacer un boxplot comparativo
pdf("boxplot_twins.pdf", width = 6, height = 6)
boxplot(
  muestra1, muestra2,
  names = c("Twin1", "Twin2"),
  main  = "Comparación de metilación 5mC",
  ylab  = "Fracción metilada",
  col   = c("darkgreen", "darkred"),
  border= "darkgrey",
  notch = TRUE
)
dev.off()
median1=median(muestra1)
median2=median(muestra2)
mean1=mean(muestra1)
mean2=mean(muestra2)
diffmean=(abs(median1-median2)*100)
print(diffmean)
diffmedian=(abs(mean1-mean2)*100)
print(diffmedian)



pdf("density_diff_twins_zoom.pdf", width=6, height=4)
d <- density(diff)

plot(d, 
     main = "Densidad de Δ metilación (5mC)", 
     xlab = "Twin1 – Twin2", 
     ylab = "Densidad",
     xlim = c(-0.3, 0.3),     # recorta ±0.3
     ylim = c(0, 10),         # recorta altura para ver la cola
     lwd  = 2)

abline(v = median(diff), col = "red", lty = 2, lwd = 2)  
dev.off()




library(ggplot2)
library(reshape2)

seleccion <- 106
top_regions <- resultado[order(-abs(resultado$diff)), ][1:seleccion, ]

heat_data <- data.frame(
  seg_id = top_regions$seg_id,      
  Twin1  = top_regions$meth_t1,
  Twin2  = top_regions$meth_t2,
  stringsAsFactors = FALSE
)

heat_data$seg_id <- factor(
  heat_data$seg_id,
  levels = rev(heat_data$seg_id)
)
library(reshape2)
heat_long <- melt(
  heat_data,
  id.vars       = "seg_id",
  variable.name = "Gemelo",
  value.name    = "MethFrac"
)

library(ggplot2)
pdf("heatmap_top50_diff.pdf", width = 5, height = 8)
ggplot(heat_long, aes(x = Gemelo, y = seg_id, fill = MethFrac)) +
  geom_tile() +
  scale_fill_gradient2(
    low      = "lightblue",
    mid      = "white",
    high     = "red",
    midpoint = median(heat_long$MethFrac)
  ) +
  theme_minimal(base_size = 12) +
  labs(
    x     = NULL,
    y     = "Isla CpG",
    fill  = "5mC frac.",
    title = paste0("Top ", seleccion, " islas CpG diferenciales (Diferencias del ", umbral*100, "%)")
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 8),
    plot.title   = element_text(hjust = 0.5)
  )
dev.off()

