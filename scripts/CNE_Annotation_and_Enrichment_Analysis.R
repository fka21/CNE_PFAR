### SCRIPT METADATA ###
# Title: Genomic Range Analysis and Annotation
# Author: Ferenc Kagan
# Date: 2024.08.24
# Description: This script performs a detailed analysis of conserved non-coding elements (CNEs) and associated genomic annotations
#              in the Macropodus opercularis genome, including data cleaning, exploratory data analysis, and enrichment analysis.

### LIBRARIES ###
# Install and load necessary libraries
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

# Load all required libraries using pacman
pacman::p_load(
  GenomicRanges, tidyverse, rtracklayer, gprofiler2, rGREAT, regioneR,
  patchwork, RColorBrewer, Biostrings, ChIPseeker, GenomicFeatures
)
### SET WORKING DIRECTORY ###
setwd("~/Documents/Projects/CNE_PFARs/")

### FUNCTIONS ###
# Define a function to read and clean data files
read_and_clean <- function(file_path, col_names = TRUE, remove_suffix = NULL) {
  data <- read_tsv(file_path, col_names = col_names)
  if (!is.null(remove_suffix)) {
    data <- data %>%
      mutate(across(everything(), ~str_remove_all(., remove_suffix)))
  }
  return(data)
}

### READ IN DATA ###
mope_cne <- read_and_clean("input/mope_most-cons.bed", col_names = FALSE)
mope_cne_supp <- read_tsv("input/mope_phyloP", col_names = TRUE) %>%
  mutate(pval = p.adjust(pval))

mope_sizes <- read_and_clean("input/mope_sizes.txt", col_names = FALSE, remove_suffix = "\\.[0-9]$")
colnames(mope_sizes) <- c("chrom", "length")

mope_anno <- GenomicFeatures::makeTxDbFromGFF("input/macOpe2_annotations/macOpe2_annotated_p2u.gff",
                                              organism = "Macropodus opercularis")

mope_mask <- read_and_clean("input/mope_masks.bed", col_names = FALSE)

### EXPLORATORY DATA ANALYSIS (EDA) ###
# Convert CNEs and mask data to GRanges objects
mope_cne_gr <- GRanges(seqnames = mope_cne$X1,
                       ranges = IRanges(start = mope_cne$X2, end = mope_cne$X3))

mope_mask_gr <- GRanges(seqnames = mope_mask$X1,
                        ranges = IRanges(start = mope_mask$X2, end = mope_mask$X3))

# Filter significant annotations and convert to GRanges
anar <- mope_cne_supp %>%
  filter(pval <= 0.05)
anar_gr <- GRanges(seqnames = anar$chr,
                   ranges = IRanges(start = anar$start,
                                    end = anar$end))

# Ensure consistent sequence levels between objects
seqlevelsStyle(mope_cne_gr) <- seqlevelsStyle(mope_anno)
common_seqlevels <- intersect(seqlevels(mope_cne_gr), seqlevels(mope_anno))
mope_cne_gr <- keepSeqlevels(mope_cne_gr, common_seqlevels, pruning.mode = "coarse")
mope_anno <- keepSeqlevels(mope_anno, common_seqlevels, pruning.mode = "coarse")

# Remove conserved sequences overlapping with exons and mask regions
exons <- GenomicFeatures::exons(mope_anno)
mope_cne_gr <- mope_cne_gr %>%
  subsetByOverlaps(exons, invert = TRUE) %>%
  subsetByOverlaps(mope_mask_gr, invert = TRUE)

anar_gr <- anar_gr %>%
  subsetByOverlaps(exons, invert = TRUE) %>%
  subsetByOverlaps(mope_mask_gr, invert = TRUE) %>%
  subsetByOverlaps(mope_cne_gr, invert = FALSE, maxgap = 100)

# Annotate peaks and create visualizations
peak_anno_list <- lapply(list(mope_cne_gr, anar_gr), annotatePeak, TxDb = mope_anno,
                         tssRegion = c(-3000, 3000),
                         addFlankGeneInfo = TRUE, 
                         genomicAnnotationPriority = c("Intergenic", "Downstream", "Promoter", "5UTR", "3UTR", "Intron", "Exon"))
names(peak_anno_list) <- c("Teleost_CNE", "ANAR")

# Plot and save peak annotation distributions
p1 <- plotAnnoBar(peak_anno_list) + ggtitle(NULL)
p2 <- plotDistToTSS(peak_anno_list) + ggtitle(NULL) + ylab("CNE (%) (5' -> 3')")
combined_plot <- p1 / p2 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 18))
ggsave("output/peak_annotation_distribution.png", combined_plot, units = 'in', width = 10, height = 10, dpi = 320)

### OVERLAP WITH DIFFERENTIALLY EXPRESSED GENES (DGE) ###
dge <- read_tsv("input/20240621-edgeR_LabOrg-Gill.tabular") %>%
  filter(FDR <= 0.05)
crossref <- read_tsv("input/pfish_macOpe2_geneID_UCSC_ann.gff", col_names = FALSE)

# Merge ANAR annotation with DGE and cross-reference data
anar_tibble <- peak_anno_list$ANAR@anno %>%
  as_tibble() %>%
  separate_rows(flank_geneIds, sep = ";") 

merged_data <- anar_tibble %>%
  left_join(dge, by = c("flank_geneIds" = "GeneID")) %>%
  left_join(crossref, by = c("flank_geneIds" = "X1")) %>%
  dplyr::select(-c("transcriptId", "flank_txIds")) %>%
  distinct()

# Write output data to file
write.table(merged_data, file = "output/anar_annotations.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

# Create BED file of unique annotations
asd_bed <- merged_data %>% distinct() %>% dplyr::select(1:3)
write.table(asd_bed, "output/ANAR.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Calculate genomic coverage percentages
cne_coverage <- sum(width(mope_cne_gr)) / sum(as.integer(mope_sizes$length)) * 100
anar_coverage <- sum(width(anar_gr)) / sum(as.integer(mope_sizes$length)) * 100

# Set sequence lengths for ANAR GRanges object
seqlengths(anar_gr) <- as.integer(mope_sizes$length[match(names(seqlengths(anar_gr)), mope_sizes$chrom)])

# Plot ANAR chromosomal distribution
cne_distribution_plot <- CNEr::plotCNEDistribution(anar_gr, chrScale = "Mb", chrs = paste0("chr", 1:24)) +
  geom_point(size = 0.01, alpha = 0.5) +
  theme(legend.key = element_blank(), 
        text = element_text(size = 7),
        strip.background = element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("output/ANAR_chromosomal_distribution.png", cne_distribution_plot, width = 10, height = 10, units = 'in', dpi = 300)

### ANNOTATION WIDTH DISTRIBUTION PLOT ###
# Prepare data for plotting annotation widths
plot_data1 <- peak_anno_list$ANAR@anno %>%
  as_tibble() %>%
  mutate(Annotation = str_replace_all(annotation, "Intron .*", "Intron"),
         Annotation = str_replace_all(annotation, "Promoter .*", "Promoter"),
         Category = "ANAR") 

plot_data2 <- peak_anno_list$Teleost_CNE@anno %>%
  as_tibble() %>%
  mutate(Annotation = str_replace_all(annotation, "Intron .*", "Intron"),
         Annotation = str_replace_all(annotation, "Promoter .*", "Promoter"), 
         Category = "CNE")

# Create and save density plot
density_plot <- bind_rows(plot_data1, plot_data2) %>%
  ggplot(aes(x = width, color = Category)) +
  geom_density(size = 2, alpha = 0.5) +
  facet_wrap(~Annotation) +
  scale_x_log10() +
  theme_bw() +
  ylab("Density") +
  xlab("log10(nt)") +
  theme(text = element_text(size = 15),
        legend.position = c(0.065, 0.85),
        legend.background = element_rect(fill = "white", color = "black"))
ggsave("output/cne_widths.png", density_plot, height = 10, width = 10, units = 'in', dpi = 300)

### ENRICHMENT ANALYSIS ###
cnes <- peak_anno_list$Teleost_CNE@anno %>%
  as_tibble() %>%
  group_by(geneId) %>%
  summarise(CNE = n())

anars <- peak_anno_list$ANAR@anno %>%
  as_tibble() %>%
  group_by(geneId) %>%
  summarise(ANAR = n())

enrichment_results <- full_join(cnes, anars) %>%
  mutate(CNE_genome = length(mope_cne_gr),
         ANAR_genome = length(anar_gr)) %>%
  mutate(enrich = phyper(ANAR - 1, ANAR_genome, CNE_genome - ANAR_genome, CNE, lower.tail = FALSE)) %>%
  mutate(enrich = p.adjust(enrich)) %>%
  filter(enrich <= 0.05) %>%
  left_join(crossref, by = c("geneId" = "X1")) %>%
  left_join(dge, by = c("geneId" = "GeneID"))

write.table(enrichment_results, "output/genes_enriched_in_ANARs.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
