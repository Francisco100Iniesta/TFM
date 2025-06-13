suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(data.table)
  library(org.Hs.eg.db)
})

gff_path <- "genomic.gff"      #archivo gff descargado desde la web del ncbi                 
out_rds  <- "chm13v2.0_annotation.rds"          

# Construir TxDb con la libreria
txdb <- makeTxDbFromGFF(gff_path, format="gff3")

contig_map <- c(
  "CP068277.2" = "NC_060925.1",  # chr1
  "CP068276.2" = "NC_060926.1",  # chr2
  "CP068275.2" = "NC_060927.1",  # chr3
  "CP068274.2" = "NC_060928.1",  # chr4
  "CP068273.2" = "NC_060929.1",  # chr5
  "CP068272.2" = "NC_060930.1",  # chr6
  "CP068271.2" = "NC_060931.1",  # chr7
  "CP068270.2" = "NC_060932.1",  # chr8
  "CP068269.2" = "NC_060933.1",  # chr9
  "CP068268.2" = "NC_060934.1",  # chr10
  "CP068267.2" = "NC_060935.1",  # chr11
  "CP068266.2" = "NC_060936.1",  # chr12
  "CP068265.2" = "NC_060937.1",  # chr13
  "CP068264.2" = "NC_060938.1",  # chr14
  "CP068263.2" = "NC_060939.1",  # chr15
  "CP068262.2" = "NC_060940.1",  # chr16
  "CP068261.2" = "NC_060941.1",  # chr17
  "CP068260.2" = "NC_060942.1",  # chr18
  "CP068259.2" = "NC_060943.1",  # chr19
  "CP068258.2" = "NC_060944.1",  # chr20
  "CP068257.2" = "NC_060945.1",  # chr21
  "CP068256.2" = "NC_060946.1",  # chr22
  "CP068255.2" = "NC_060947.1",  # chrX
  "CP086569.2" = "NC_060948.1",  #chrY
  "CP068254.1" = NA            # chrM (sin RefSeq asignado)
)
rev_map <- setNames(names(contig_map), contig_map)
old_ns <- intersect(seqlevels(txdb), names(rev_map))
new_ns <- rev_map[old_ns]
names(new_ns) <- old_ns

# Renombrar cromosoma de refseq a genBank
txdb <- renameSeqlevels(txdb, new_ns)

# A partir de este punto es buscar y aplicar funciones de genomic features y  Genomic cranges para extraer información

# Promotores
tx_gr     <- transcripts(txdb, columns=c("TXID","GENEID","TXNAME"))
prom_gr   <- promoters(tx_gr, upstream=1000, downstream=100)
prom_gr$symbol <- mcols(tx_gr)$GENEID
prom_gr$tx_id   <- mcols(tx_gr)$TXNAME
prom_gr$class    <- "genes"
prom_gr$subclass    <- "promoters"
# Exones
ex_by_gene    <- exonsBy(txdb, by="gene")
exons_gr      <- unlist(ex_by_gene)
exons_gr$symbol <- names(exons_gr)
exons_gr$class    <- "genes"
exons_gr$subclass    <- "exons"

#  Intrones
intr_by_tx       <- intronsByTranscript(txdb, use.names=TRUE)
introns_gr       <- unlist(intr_by_tx)
names(introns_gr) <- NULL
intr_counts       <- elementNROWS(intr_by_tx)
introns_gr$tx_id <- rep(names(intr_by_tx), intr_counts)
introns_gr$class  <- "genes"
introns_gr$subclass  <- "introns"
tx2gene          <- select(txdb,
                           keys=keys(txdb, keytype="TXID"),
                           keytype="TXID",
                           columns=c("GENEID","TXNAME"))
introns_gr$symbol <- tx2gene$GENEID[ match(introns_gr$tx_id, tx2gene$TXNAME) ]

#  UTRs 5'
utr5_by_tx      <- fiveUTRsByTranscript(txdb, use.names=TRUE)
utr5_gr         <- unlist(utr5_by_tx)
names(utr5_gr)  <- NULL
utr5_counts     <- elementNROWS(utr5_by_tx)
utr5_gr$tx_id   <- rep(names(utr5_by_tx), utr5_counts)
utr5_gr$class    <- "genes"
utr5_gr$subclass    <- "5UTRs"
utr5_gr$symbol <- tx2gene$GENEID[ match(utr5_gr$tx_id, tx2gene$TXNAME) ]

# UTRs 3'
utr3_by_tx      <- threeUTRsByTranscript(txdb, use.names=TRUE)
utr3_gr         <- unlist(utr3_by_tx)
names(utr3_gr)  <- NULL
utr3_counts     <- elementNROWS(utr3_by_tx)
utr3_gr$tx_id   <- rep(names(utr3_by_tx), utr3_counts)
utr3_gr$class    <- "genes"
utr3_gr$subclass    <- "3UTRs"
utr3_gr$symbol <- tx2gene$GENEID[ match(utr3_gr$tx_id, tx2gene$TXNAME) ]



#CDS
cds_by_tx <- cdsBy(txdb, by="tx", use.names=TRUE)

# 2) Filtrar fuera los elementos con name = NA
valid_tx   <- !is.na(names(cds_by_tx))
cds_by_tx  <- cds_by_tx[valid_tx]
cds_gr     <- unlist(cds_by_tx, use.names=FALSE)
tx_names   <- rep(names(cds_by_tx), elementNROWS(cds_by_tx))
cds_gr$tx_name <- tx_names
tx2gene <- select(
  txdb,
  keys    = unique(tx_names),
  keytype = "TXNAME",
  columns = "GENEID"
)

cds_gr$gene_id <- tx2gene$GENEID[ match(cds_gr$tx_name, tx2gene$TXNAME) ]
cds_gr$class    <- "genes"
cds_gr$subclass <- "cds"
cds_gr$symbol <- cds_gr$gene_id


#LAS ANOTACIONES ESTAN EN UN FORMATE granges (parecido a una lista de texto), entonces se fusionan todas en una sola 
all_annots <- c(prom_gr, exons_gr, introns_gr, utr5_gr, utr3_gr, cds_gr)
names(all_annots) <- NULL

# Se necesita para retroinspector un id de subclase para que se indique en los informes, entonces se crea una columna con la subclase del elemento para que se reporte en dichos informes
all_annots$id <- paste0(all_annots$subclass, "_", seq_along(all_annots))

# Eliminamos columnas sin id dado que nos reportaría como nada en el pipeline y la colapsaría
if (is.null(all_annots$tx_id)) {
  all_annots$tx_id <- NA
}
# pasamos a datafreme
dt=(as.data.frame(all_annots))
# guarsdamos como .rds standard este luego hay que subirlo a la carpeta data del pipeline
saveRDS(all_annots, file = out_rds)
