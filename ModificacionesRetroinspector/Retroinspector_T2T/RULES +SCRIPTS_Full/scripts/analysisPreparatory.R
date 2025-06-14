options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

library(data.table)
data.table::setDTthreads(threads = 2)
library(stringi)
library(ggplot2)
library(magrittr)
library(qqman) %>% suppressMessages()

# Patients and folders
# ----
samples = snakemake@params[["samples"]]
# ----

if (F) {
  
  saveRDS(snakemake, paste0(snakemake@config[["outputPath"]], "/snk.rds"))
  q()
  snakemake = readRDS("~/path/snk2.rds")
  setwd("~/path/")
}
if (F) {
  snakemake = readRDS("~/ontology/proj/public5/snk.rds")
}

annotate <- function(tableToAnnotate) {
  library(annotatr)
  library(GenomicRanges)
  library(data.table)

  annot_path <- snakemake@input[["annotationRds"]]
  annotation <- readRDS(annot_path)
  saveRDS(annotation, snakemake@output[["annotation"]])

  message("Built annotation")
  message("Annotating...")

  ann_gr <- annotatr::annotate_regions(
    regions = makeGRangesFromDataFrame(
      tableToAnnotate,
      seqnames.field     = "seqnames",
      keep.extra.columns = TRUE
    ),
    annotations   = annotation,
    ignore.strand = FALSE,
    quiet         = TRUE
  )

  annotatedInsertions <- data.table(as.data.frame(ann_gr, row.names = NULL, check.names = FALSE))

  for (col in c("annot.strand", "annot.seqnames", "annot.width", "annot.id")) {
    if (col %in% names(annotatedInsertions)) {
      annotatedInsertions[, (col) := NULL]
    }
  }
  
   if ("annot.symbol" %in% names(annotatedInsertions)) {
    annotatedInsertions[, annot.symbol := as.character(annot.symbol)]
  }
  if ("annot.gene_id" %in% names(annotatedInsertions)) {
    annotatedInsertions[, annot.gene_id := as.character(annot.gene_id)]
  }
  if ("symbol" %in% names(annotatedInsertions)) {
    annotatedInsertions[, symbol := as.character(symbol)]
  }
  if ("gene_id" %in% names(annotatedInsertions)) {
    annotatedInsertions[, gene_id := as.character(gene_id)]
  }


  message("Finished annotation")
  return(annotatedInsertions)
}


repeatMaskerTable = fread(
    snakemake@input[["repeatMaskerBed"]],
    col.names = c(
    "seqnames", "start", "end", "name", "score", "strand", "seqId", 
    "sw_score", 
    # "rename_score", 
    "repeat.class",
    "repeat.subclass",
    "repeat.start",
    "repeat.end",
    "repeat.left",
    "repeat.strains",
    "repeat.divergence_percentage",
    "repeat.deletion_percentage",
    "repeat.insertion_percentage",
    "vcf_alt",
    "vcf_info",
    paste0("id", samples %>% stri_replace_all_fixed(., pattern = "-", replacement = "."))
    )
  )
# Generating the usable datasets: min3 and trusty, lax and strict critera, respectively
genoFields = c("GT", "PSV", "LN", "DR", "ST", "QV", "TY", "ID", "RAL", "AAL", "CO", "SC")
genoFieldsOfInt = c("GT", "PSV", "DR") # Used for counting occurrence and genotyping
pos = which(genoFields %in% genoFieldsOfInt) # The position of "PSV" and "DR" in the genotype columns

for (col in repeatMaskerTable %>% colnames() %>%
      .[. %in% 
      paste0("id", samples %>% stri_replace_all_fixed(., pattern = "-", replacement = "."))
      ]) {
  # create columns for GT, PSV and DR for each patient
  repeatMaskerTable[, paste0(col, "_", genoFieldsOfInt) := tstrsplit(.SD[[col]], ":", fixed = T, keep = pos)]
  # keep second value of DR, which is the one supporting the INS
  # Variant callers dont consistently report the support for the reference allele
  repeatMaskerTable[, paste0(col, "_DR") := tstrsplit(.SD[[paste0(col, "_DR")]], ",", fixed = T, keep = 2)]
  # Fix NAs
  repeatMaskerTable[,
      paste0(col, "_DR") :=
        fifelse(.SD[[paste0(col, "_DR")]] == ".", 0, as.numeric(.SD[[paste0(col, "_DR")]]) )
      ] # This will warn about NA by coercion, which are not included in the data.table
}

# Take the DR columns, create a bool vector with DR >=3, make it numeric (1 or 0), assign as column
# this creates a column with the lax criterion
# stored as numeric for addition
repeatMaskerTable[, SUPP_VEC_min3 := apply(.SD, 2, `>=`, y = snakemake@params[["re"]]) %>% apply(., 2, as.integer) %>% apply(., 1, paste0, collapse = ""), .SDcols = grep(pattern = "_DR$", colnames(repeatMaskerTable))]
repeatMaskerTable[, SUPP_min3 := stri_count_fixed(SUPP_VEC_min3, "1")] # Count SUPPORT with the min3 criteria

# Initialize the column with only the two callers criteria
# Since a data.table is a list of columns (like a data.frame)
# apply with margin = 2 iterates rows instead of columns, and the opposite for margin = 1
# because it is "transposed"
repeatMaskerTable[
  , SUPP_VEC_trusty := apply(.SD, 2, `==`, y = "11") %>%
    apply(., 2, as.integer) %>%
    apply(., 1, paste0, collapse = ""),
  .SDcols = patterns("_PSV$")]
# The value for the column is provisional, and missing the RE >= 3 part

# ("00", "11") => "00"
# ("01", "11") => "01"
# ("11", "11") => "11"
bitwiseAndStr = function (a, b) {
  lapply (1:length(a), function (i) {
    a2 = unlist( strsplit(a[[i]], split = "", fixed = T) )
    b2 = unlist( strsplit(b[[i]], split = "", fixed = T) )
    paste0(ifelse(a2 == b2 & a2 == "1", "1", "0"), collapse = "")
  }
  )
}
# Now we do both callers AND min3 -> strict criterion
correctTrusty = unlist( bitwiseAndStr(repeatMaskerTable$SUPP_VEC_min3, repeatMaskerTable$SUPP_VEC_trusty) )
repeatMaskerTable[, SUPP_VEC_trusty := correctTrusty]
rm(correctTrusty)
repeatMaskerTable[, SUPP_trusty := stri_count_fixed(SUPP_VEC_trusty, "1")]

# We will also need the INS' length
repeatMaskerTable[, SVLEN := as.numeric(gsub(".*SVLEN=([^;]+);.*", "\\1", vcf_info))]

# We extract the genotype from the columns and store it as numeric
# 0 is hom. ref., 1 is het., 2 is hom. alt.
for (col in repeatMaskerTable %>% colnames() %>% 
      grep(pattern = "_GT$", x = ., value = T)) {
  repeatMaskerTable[, (col) := get(col) %>%
                                  lapply(function (x) {
                                    x %>%
                                      strsplit(split = "/", fixed = T) %>%
                                      unlist() %>%
                                      as.numeric() %>%
                                      sum()
                                  }) %>% unlist()]
}
repeatMaskerTable[
  , genotype := .SD %>% 
    simplify2array() %>% 
    t() %>% 
    split(., rep(1:ncol(.), each = nrow(.))) %>% 
    unname(), 
  .SDcols = grep(pattern = "_GT$", colnames(repeatMaskerTable))]
repeatMaskerTable[, SUPP_mask_geno := genotype %>% lapply(`>=`, y = 1)]
repeatMaskerTable[, SUPP_VEC_geno := SUPP_mask_geno %>% lapply(as.numeric) %>% lapply(paste0, collapse = "") %>% unlist()]
repeatMaskerTable[, SUPP_geno := SUPP_mask_geno %>% lapply(sum) %>% unlist()]
repeatMaskerTable[, maf := genotype %>% lapply(function(x) sum(x) / (length(samples) * 2 ) ) %>% unlist() ] # humans are diploid, thus the length() * 2
# Create vector bool columns
repeatMaskerTable[, SUPP_mask_min3 := SUPP_VEC_min3 %>%
            lapply(strsplit, split = "") %>%
            lapply(function(x)as.logical(as.numeric(unlist(x))))]
repeatMaskerTable[, SUPP_mask_trusty := SUPP_VEC_trusty %>%
            lapply(strsplit, split = "") %>%
            lapply(function(x)as.logical(as.numeric(unlist(x))))]

repeatMaskerTable[, # split level2-level3 and replace NA for level 3 with "-"
  c("repeat.subclass", "repeat.family") := stri_match_all_regex(repeat.subclass, "([^-]+)(?:-([^-]+))?$") %>% rapply(., f = `[`, ... = 2:3, how = "r") %>%
      lapply(function (x) {x[is.na(x)] = "-"; x}) %>% data.table::transpose() ]

# Remove unnecessary columns
repeatMaskerTable[, grep("^id", colnames(repeatMaskerTable)) := NULL]
repeatMaskerTable[, vcf_info := NULL]




#
# Filtering
#
# ----
chrs = c(
  "CP068277.2", "CP068276.2", "CP068275.2", "CP068274.2",
  "CP068273.2", "CP068272.2", "CP068271.2", "CP068270.2",
  "CP068269.2", "CP068268.2", "CP068267.2", "CP068266.2",
  "CP068265.2", "CP068264.2", "CP068263.2", "CP068262.2",
  "CP068261.2", "CP068260.2", "CP068259.2", "CP068258.2",
  "CP068257.2", "CP068256.2", "CP068255.2", "CP086569.2"
)


class1 = c("Retroposon", "SINE", "LINE", "LTR") # Class I
class2 = c("DNA")   # Class II elements

repeatMaskerTable[, repeat.percentage := (repeat.end - repeat.start + 1) / SVLEN] # Percentage of INS covered by RepeatMasker annotation
repeatMaskerTableMin3 = repeatMaskerTable[SUPP_min3 > 0]
allIns = repeatMaskerTable # An unfiltered version, with only one row per INS (so no annotation)
# rm(annotatedInsertions); gc()
# We filter out the ones that are not >=85% repetitive elements (includes TEs and non TEs),
# the ones that are not TEs and the ones not in chr1-22, X or Y
repeatMaskerTableMin3 = repeatMaskerTableMin3[repeat.percentage >= 0.85 & repeat.class %in% c(class1, class2) & seqnames %in% chrs]
print("Filtered repeatMaskerTableMin3")
annotatedInsertionsMin3 = annotate(repeatMaskerTableMin3)
# insertionsTable = annotatedInsertionsMin3[, .SD[1], by = seqId] # Only one annotation per INS, a "unique" version of the above
# print("Created <insertionsTable> table")
# insertionsTable[, grep("^annot\\.", colnames(insertionsTable)) := NULL] # Remove the annotation data (incomplete in this table)
print("Removed annotation columns from <insertionsTable>")
saveRDS(annotatedInsertionsMin3, snakemake@output[["annotatedInsertionsMin3"]])
saveRDS(repeatMaskerTableMin3, snakemake@output[["insertionsTable"]])
saveRDS(allIns, snakemake@output[["allIns"]])
message("Annotating...delection")

meDeletionsMin3 = fread(snakemake@input[["repeatMaskerVcfDel"]], skip = "#CHROM")
meDeletionsMin3 = meDeletionsMin3[
  !INFO %like% "ME=[^;]*\\?[^;]*" &
  !INFO %like% "MEFAM=[^;]*\\?[^;]*" &
  `#CHROM` %in% chrs
]

meDeletionsMin3[, `:=`(
  survivorId = ID,
  ID         = paste0("surv", seq_len(.N)),
  svlen      = abs(as.numeric(sub(".*SVLEN=-([0-9]+).*", "\\1", INFO)))
)]

genoFields = meDeletionsMin3[1, 9] %>% unname() %>% stri_split_fixed(pattern = ":") %>% unlist()
genoFieldsOfInt = c("GT", "PSV")
pos = which(genoFields %in% genoFieldsOfInt)
for (col in meDeletionsMin3 %>% colnames() %>% 
     grep(pattern = "^Sample", x = ., ignore.case = T, value = T) ) {
  meDeletionsMin3[, paste0(col, "_", genoFieldsOfInt) := tstrsplit(.SD[[col]], ":", fixed = T, keep = pos)]
}
meDeletionsMin3[, SUPP_VEC_min3 := apply(.SD, 2, `!=`, y = "NaN") %>% apply(., 2, as.integer) %>% apply(., 1, paste0, collapse = ""), .SDcols = grep(pattern = "_PSV$", colnames(meDeletionsMin3))]
meDeletionsMin3[, SUPP_min3 := stri_count_fixed(SUPP_VEC_min3, "1")]
meDeletionsMin3[, SUPP_VEC_trusty := apply(.SD, 2, `==`, y = "11") %>% apply(., 2, as.integer) %>% apply(., 1, paste0, collapse = ""), .SDcols = patterns("_PSV$")]
meDeletionsMin3[, SUPP_trusty := stri_count_fixed(SUPP_VEC_trusty, "1")]
meDeletionsMin3[, SVLEN := as.numeric(gsub(".*SVLEN=([^;]+);.*", "\\1", INFO)) %>% abs()]
# extract genotype for deletions
for (col in meDeletionsMin3 %>% colnames() %>% 
     grep(pattern = "_GT$", x = ., value = T)) {
  meDeletionsMin3[, (col) := get(col) %>% ifelse(. == "./.", "0/0", .) %>%
                    lapply(function (x) {
                      x %>%
                        strsplit(split = "/", fixed = T) %>%
                        unlist() %>%
                        as.numeric() %>%
                        sum()
                    }) %>% unlist()]
}
meDeletionsMin3[
  , genotype := .SD %>% 
    simplify2array() %>% 
    t() %>% 
    split(., rep(1:ncol(.), each = nrow(.))) %>% 
    unname(), 
  .SDcols = grep(pattern = "_GT$", colnames(meDeletionsMin3))]
meDeletionsMin3[, SUPP_VEC_geno := genotype %>% lapply(`>=`, y = 1)]
meDeletionsMin3[, SUPP_geno := SUPP_VEC_geno %>% lapply(sum) %>% unlist()]
meDeletionsMin3[, maf := genotype %>% lapply(function(x) sum(x) / (length(samples) * 2) ) %>% unlist() ] # humans are diploid, thus the length() * 2

meDeletionsMin3[, grep("^(Sample|FORMAT)", colnames(meDeletionsMin3), ignore.case = T) := NULL]
meDeletionsMin3 = meDeletionsMin3[
  , .(
    seqnames = `#CHROM`,
    start = POS,
    end = sub(pattern = ".*;END=([^;]+).*", replacement = "\\1", x = INFO, perl = T) %>% as.numeric(),
    id = ID,
    SUPP_min3,
    SUPP_VEC_min3,
    SUPP_trusty,
    SUPP_VEC_trusty,
    genotype,
    maf,
    SUPP_geno,
    SUPP_VEC_geno,
    svlen,
    # info = INFO,
    name = sub(pattern = ".*ME=([^;]+).*", replacement = "\\1", x = INFO, perl = T),
    repeat.class = sub(pattern = ".*MEFAM=([^;]+).*", replacement = "\\1", x = INFO, perl = T),
    repeat.coords = sub(pattern = ".*MECOORDS=([^;]+).*", replacement = "\\1", x = INFO, perl = T)
  )
][, c("repeat.class", "repeat.subclass") := repeat.class %>% 
    stri_match_all_regex("([^/]+)(?:/(.+))?$") %>% 
    rapply(., f = `[`, ... = 2:3, how = "r") %>% transpose()][
      , repeat.subclass := fifelse(is.na(repeat.subclass), "-", repeat.subclass)
    ][
      , c("repeat.subclass", "repeat.family") := 
        stri_match_all_regex(repeat.subclass, "([^-]+|-)(?:-(.+))?$") %>% 
        rapply(., f = `[`, ... = 2:3, how = "r") %>% transpose()
    ][, repeat.family := fifelse(is.na(repeat.family), "-", repeat.family)]

saveRDS(meDeletionsMin3, snakemake@output[["meDeletionsMin3"]])