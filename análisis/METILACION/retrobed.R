library(data.table)


dt <- fread("twin2.csv", sep = ";")

# Extraer chrom y posición
dt[, c("chrom","pos_str") := tstrsplit(Position, ":", fixed=TRUE)]
dt[, pos_str := gsub("[^0-9]", "", pos_str)]
dt[, pos := as.integer(pos_str)]

# creo start=end=pos
dt[, `:=`(start = pos, end = pos)]

# renombrar las muestras
dt[, Sample := fifelse(
  Sample == "reads_LSA-2024-0039_GS1375-I-Sangre", "twin1",
  fifelse(Sample == "reads_LSA-2024-0040_GS1375-H-Sangre", "twin2", NA_character_)
)]

#esto se comento en la metodología dado que se reeportan duplicados para inserciones compartidas en los informes de retroinspector.

dt <- dt[!duplicated(dt, by = c("chrom","pos"))]


# se calcula el centro de la insercion
dt[, center := floor((start + end) / 2)]

# Se crea una ventana de 300 bases aguas arriba y abajo
dt[, `:=`(
  start = pmax(0L, center - 300L),
  end   = center + 300L
)]

# guardamos el elemento para no perderlo en metylartist y añadimos que se busque en ambos sentidos con *
dt[, Type := paste0(Type, ":", Sample)]
dt[, strand := "*"]

# generamos bed
dt[, c("Position","seqId","Length","pos","pos_str","Sample") := NULL]

fwrite(
  dt[, .(chrom, start, end, Type, strand)],
  "unshared_insertions_twin2.bed",
  sep       = "\t",
  col.names = FALSE,
  quote     = FALSE
)


