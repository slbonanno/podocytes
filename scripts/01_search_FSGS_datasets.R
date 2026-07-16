# the first step in this project is to find data that is FSGS and podocyte relevant
# FSGS is a medically unambiguous acronym; metadata with this classifier is a good filter


###########################################################################################
# 1. find scRNA-seq datasets that contain podocytes on entrez_search
###########################################################################################

install.packages("rentrez")
library(rentrez)
library(dplyr)
library(GEOquery)

# entrez_search is going to pull max 30 records with these filters
# this does not contain much metadata, we will run this out through summary to get that.
# tried at first without Organism filter - best to keep to Human, since want to look at Human dysregulated
res_kidney <- entrez_search(
  db = "gds",
  term = 'kidney[All Fields] AND podocyte[All Fields] AND single cell[All Fields] AND Homo sapiens[Organism] AND gse[Filter]',
  retmax = 30
)

# this outputs a complex df with records.  cmds below tell you how many records pulled (should be <= 30) and their ids
# res_kidney$count
# res_kidney$ids
# in this query, using Human as a filter, we get 29 records returned.

# explore one record
one <- entrez_summary(db = "gds", id = res_kidney$ids[1])
one$title
one$taxon
one$suppfile

###########################################################################################
# 2. loop to run entrez_summary on all records (queries NCBI)
###########################################################################################

# in R appending to a vector repeatedly inside a loop is slow; pre-allocate empty storage is faster
titles    <- character(length(res_kidney$ids))
taxons    <- character(length(res_kidney$ids))
suppfiles <- character(length(res_kidney$ids))
n_samps   <- character(length(res_kidney$ids))
gses      <- character(length(res_kidney$ids))

# loop through entrez_summary for each entry
for (i in seq_along(res_kidney$ids)) {
  rec <- entrez_summary(db = "gds", id = res_kidney$ids[i])
  titles[i]    <- rec$title
  taxons[i]    <- rec$taxon
  suppfiles[i] <- rec$suppfile
  n_samps[i]   <- rec$n_samples
  gses[i]      <- rec$gse
}

df_kidney <- data.frame(gse = gses, title = titles, taxon = taxons,
                        n_samples = n_samps, suppfile = suppfiles)
View(df_kidney)

# this provides 29 records, some of which are a better match for the task at hand:
# robust scRNA-seq dataset with enough podocytes and other cell types
# to define marker genes for podocytes
# looked at the list of 29 and chose 5 that look most likely (by descr) to help here
# but when I looked more carefully, only ended up with 1 record that looked worth following
# aiming for at least 2 or 3, so inspect all 29 bioinformatically instead

# this will output description of the 29 records
for (i in seq_len(nrow(df_kidney))) {
  cat(sprintf("[%d] GSE%s | n=%s | %s\n    suppfile: %s\n\n",
              i, df_kidney$gse[i], df_kidney$n_samples[i],
              df_kidney$title[i], df_kidney$suppfile[i]))
}

# I tried to filter these by the below, but it dropped records that are useful
# exclude_terms <- c("organoid", "iPSC", "embryonic", "mouse", "murine", "E18.5", "DNase", "chromatin accessibility", "methylation", "GWAS")
# it is better to inspect manually.  Did this and decided to keep records [7] and [22]

# here instead, searching for records with "single cell" or "single-cell"
df_kidney_single_cell <- filter(
  df_kidney,
  grepl("single[- ]cell", title, ignore.case = TRUE)
)

# look at the retained records - there are 13 out of 29
View(df_kidney_single_cell)
nrow(df_kidney_single_cell)

# look at details for the retained 13 records
for (i in seq_len(nrow(df_kidney_single_cell))) {
  cat(sprintf("GSE%s | n=%s | %s\n    suppfile: %s\n\n",
              df_kidney_single_cell$gse[i], df_kidney_single_cell$n_samples[i],
              df_kidney_single_cell$title[i], df_kidney_single_cell$suppfile[i]))
}

###########################################################################################
# 3. extend download and unzip process to a few other datasets.
###########################################################################################

# 1. GSE209781: 3 healthy, 3 diabetic kidney disease, 10x scRNA-seq — downloaded, extracted, confirmed standard barcodes/features/matrix format
# 2. GSE131882: 3 healthy, 3 early diabetic nephropathy, snRNA-seq (10x) — known-good from literature review, not yet verified via GEO search (missed by "single cell" keyword since it's single-nucleus); file structure TBC
# 3. GSE195797: human, n=2, parietal epithelial cell (PEC) crescent formation study, 10x scRNA-seq — glomerular tissue, likely contains podocytes; small sample size, file structure TBC
# 4. GSE270701: human, n=1, pediatric SRNS (coenzyme Q10 nephropathy) podocytopathy, 10x scRNA-seq — disease-relevant but single patient; file structure TBC

library(GEOquery)

# use GEO supp files from GEOquery to confirm what data is provided (processed or raw only?)
# the below confirms .tar is provided, can't look inside without decompressing

# ---- 1: list files for all 3, don't download yet ----
# files_209781 <- getGEOSuppFiles("GSE209781", fetch_files = FALSE) # already done
files_131882 <- getGEOSuppFiles("GSE131882", fetch_files = FALSE)
files_195797 <- getGEOSuppFiles("GSE195797", fetch_files = FALSE)
files_270701 <- getGEOSuppFiles("GSE270701", fetch_files = FALSE)

# files_209781 # already done in figuring out how to download and unzip
files_131882
files_195797
files_270701

# ---- 2: download all three (run only after reviewing Chunk 1 output) ----
getGEOSuppFiles("GSE131882", fetch_files = TRUE, baseDir = "data/raw")
getGEOSuppFiles("GSE195797", fetch_files = TRUE, baseDir = "data/raw")
getGEOSuppFiles("GSE270701", fetch_files = TRUE, baseDir = "data/raw")

# these were then unzipped using terminal, tar -xvf etc
# note that getGEOSuppFiles sometimes errors - keep track and delete and redownload when nec

###########################################################################################
# 4. continue extracting data - not all filteypes are the same (not all std 10x scRNA-seq)
# i.e. GSE...882 is snRNA-seq and not providing 10x triplet files as usual, but maybe high-value biologically
###########################################################################################


################################ GSE195797 ###############################################
library(Matrix)

mat_195797 <- readMM("data/raw/GSE195797/GSE195797_matrix.mtx.gz")
dim(mat_195797)

barcodes_195797 <- readLines(gzfile("data/raw/GSE195797/GSE195797_barcodes.tsv.gz"))
length(barcodes_195797)
head(barcodes_195797)

# Extract just the suffix after the second dash, check how many unique sample tags exist
sample_tags <- sub("^[ACGT]+-1-", "", barcodes_195797)
table(sample_tags)

# RPC_CTL (control) and RPC_TR (treated/disease, presumably)
# 3,283 and 2,825 cells respectively, matching the n=2 from the series metadata.

################################ GSE131882 ###############################################

## note that this was supplied as RDS but gz twice, so normal readRDS was failing
# gunzip first in bash, then readRDS works

one_sample <- readRDS("data/raw/GSE131882/GSM3823939_control.s1.dgecounts_final.rds")
str(one_sample, max.level = 2)

# classic zUMIs output structure: umicount vs readcount, exon / inex (pre-mRNA) / intron
# So the dataset-appropriate pick is likely umicount$inex, not the naive umicount$exon
# worth flagging since defaulting to matching 10x's convention would actually be the wrong call for nucleus data specifically.

str(one_sample$umicount$inex, max.level = 2)

# one_sample$umicount$inex$all is the sparse gene × cell matrix (dgCMatrix, the exact same class Read10X() returns for the other datasets)
# and downsampling$downsampled_ is a QC-normalized-depth variant you don't need

mat_131882 <- one_sample$umicount$inex$all
dim(mat_131882)
mat_131882[1:5, 1:5]


# =========================================================================
# SUMMARY — 4 datasets confirmed downloadable + loadable, ready
# for notebook 02 (Seurat clustering & podocyte identification)
# =========================================================================
#
# GSE209781 — 3 healthy (NM), 3 DKD, 10x scRNA-seq
#   Load: Read10X() per sample folder (NM01-03, DKD01-03)
#   Gene IDs: symbols
#
# GSE131882 — 3 control, 3 diabetic, snRNA-seq (zUMIs pipeline)
#   Load: readRDS() -- files are DOUBLE-gzipped, need manual
#         double gunzip first (GEO packaging quirk)
#         then extract $umicount$inex$all (UMI counts,
#         exon+intron combined -- appropriate for snRNA-seq)
#   Gene IDs: Ensembl (ENSG...) -- WILL NEED ID CONVERSION
#             before searching by marker gene symbol
#
# GSE195797 — n=2 (RPC_CTL, RPC_TR), 10x scRNA-seq, PEC/crescent study
#   Load: readMM() on single combined matrix, split samples via
#         barcode suffix ("-RPC_CTL" / "-RPC_TR")
#   Gene IDs: symbols
#
# GSE270701 — n=1 (MCD1), 10x scRNA-seq, pediatric podocytopathy
#   Load: Read10X() on single sample folder
#   Gene IDs: symbols
#   NOTE: sample labeled "MCD1" (minimal change disease) --
#         verify this matches expected phenotype before use,
#         does not obviously match "CoQ10 nephropathy" from
#         series title
# =========================================================================







