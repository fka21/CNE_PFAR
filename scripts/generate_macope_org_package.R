#!/usr/bin/env Rscript

### SCRIPT METADATA ###
# Title: Generate Org Package for Macropodus opercularis
# Author: Ferenc Kagan
# Date: 2024.08.24
# Description: This script creates an organism package for Macropodus opercularis using `AnnotationForge`,
#              incorporating gene annotations, chromosome information, and Gene Ontology (GO) terms.

### LIBRARY LOADING ###
library(AnnotationForge)
library(GenomicFeatures)
library(dplyr)
library(readr)

### INPUT FILES ###
# Paths to input files
gff_file <- "input/macOpe2_annotations/macOpe2_annotated_p2u.gff"
sym_file <- "input/macOpe2_annotations/pfish_macOpe2_geneID_UCSC_ann.gff"
go_file <- "input/macOpe2_annotations/emapper/gene2go.tsv"

### READ IN DATA ###
# Create TxDb object from GFF file
macope_txdb <- makeTxDbFromGFF(gff_file, organism = "Macropodus opercularis")

# Read gene symbol and ID mapping
macope_sym <- read_tsv(sym_file, col_names = FALSE) %>%
  dplyr::rename(GID = X1, ALIAS = X2) %>%
  filter(!ALIAS %in% c("-", "NA")) %>%
  filter(!is.na(ALIAS))

# Read Gene Ontology (GO) annotations
macope_go <- read_tsv(go_file, col_names = FALSE) %>%
  mutate(EVIDENCE = "ISO") %>%
  dplyr::rename(GID = X1, GO = X2)

### DATA PROCESSING ###
# Extract chromosome information
macope_chr <- select(macope_txdb, 
                     keys = keys(macope_txdb, "GENEID"), 
                     keytype = "GENEID", 
                     columns = c("GENEID", "CDSCHROM")) %>%
  dplyr::rename(GID = GENEID, CHROM = CDSCHROM)

### GENERATE ORGANISM PACKAGE ###
makeOrgPackage(
  gene_info = macope_sym, 
  chromosome = macope_chr, 
  go = macope_go,
  version = "0.1",
  author = "Ferenc Kagan <ferenc.kagan@mpinat.mpg.de>",
  maintainer = "Ferenc Kagan <ferenc.kagan@mpinat.mpg.de>",
  outputDir = "~/Documents/Projects/CNE_PFARs/annotforge/",
  tax_id = "158451",
  genus = "Macropodus",
  species = "opercularis",
  goTable = "go"
)
