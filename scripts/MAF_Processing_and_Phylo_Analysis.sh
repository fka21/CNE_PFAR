#!/bin/bash

### SCRIPT METADATA ###
# Title: MAF Processing and Phylogenetic Analysis
# Author: Ferenc Kagan
# Date: 2024.08.24
# Description: This script processes MAF files generated from HAL, performs filtering and sorting, extracts 4D sites,
#              generates a background model, and runs phylogenetic analyses using phyloFit and phastCons.

### VARIABLES ###
# Directory containing MAF files
MAF_DIR="/projects/ferenc.kagan/CNE_PFARs/03_whole_genome_alignment/cactus_new-mope"  # Update to the correct directory

# Features file (GFF)
FEATURES_FILE="/projects/ferenc.kagan/CNE_PFARs/03_whole_genome_alignment/cactus_new-mope/macOpe2_annotated_p2u.gff"  # Update to the correct path

# Number of cores to use for processing
CORES=50

# Tree structure for phylogenetic analysis
TREE_STRUCTURE="((amex,drer),(trub,(((olat,nfur),onil),(gacu,(ates,(mope,bspl))labirynth))))"

### PREPROCESSING ###
# Ensure the necessary output directories exist
mkdir -p "${MAF_DIR}/filtered_maf" "${MAF_DIR}/filtered_sorted_maf" "${MAF_DIR}/features" "${MAF_DIR}/4d" "${MAF_DIR}/4d_sites"

### STEP 1: Generate and Process MAF Files ###
cactus-hal2maf ./js_hal2maf ./teleost.hal ./teleost_mope-ref.maf \
  --maxCores $CORES --noAncestors --filterGapCausingDupes --refGenome mope \
  --dupeMode single --chunkSize 10000 --logFile hal2maf.log --coverage

# Remove duplicates and sort MAF files
mafDuplicateFilter -k -m teleost_mope-ref.maf > teleost_mope-ref.dedup.maf
maf_sort teleost_mope-ref.dedup.maf mope > teleost_mope-ref.dedup.sorted.maf
mafSplit -useFullSequenceName ./ teleost_mope-ref.dedup.sorted.maf

# Clean up MAF filenames
rename 's/mope//g' ./mope*
mv *maf ../split_maf

### STEP 2: Process Each MAF File Individually ###
for maf_file in "$MAF_DIR"/split_maf/*.maf; do
  base_name=$(basename "$maf_file" .maf)
  echo "Processing $maf_file"

  # Filter by length, sort, and extract features
  mafFilter -l 3 --maf "$maf_file" > "${MAF_DIR}/filtered_maf/${base_name}_filtered.maf"
  maf_sort "${MAF_DIR}/filtered_maf/${base_name}_filtered.maf" mope > "${MAF_DIR}/filtered_sorted_maf/${base_name}_filtered_sorted.maf"

  grep -P "${base_name}\t" "$FEATURES_FILE" > "${MAF_DIR}/features/${base_name}_features.gff"
  sed -i 's/^/mope./g' "${MAF_DIR}/features/${base_name}_features.gff"

  # Convert GFF to SS format and extract 4D sites
  msa_view "${MAF_DIR}/filtered_sorted_maf/${base_name}_filtered_sorted.maf" --in-format MAF --4d --features "${MAF_DIR}/features/${base_name}_features.gff" > "${MAF_DIR}/4d/${base_name}_4d.ss"
  msa_view "${MAF_DIR}/4d/${base_name}_4d.ss" --in-format SS --out-format SS --tuple-size 1 > "${MAF_DIR}/4d_sites/${base_name}_4d_sites.ss"

  echo "Finished processing $maf_file"
done

### STEP 3: Clean Up Empty Files and Aggregate 4D Sites ###
find "${MAF_DIR}/4d_sites/" -type f -empty -print -delete

msa_view --unordered-ss --in-format SS --out-format SS --aggregate mope,gacu,bspl,ates,onil,nfur,olat,trub,drer,amex "${MAF_DIR}/4d_sites/*_4d_sites.ss" > aggregated_4d_sites.ss

### STEP 4: Generate Background Model and Run Phylogenetic Analysis ###
phyloFit --tree "$TREE_STRUCTURE" --subst-mod REV --out-root general --msa-format SS aggregated_4d_sites.ss

for maf_file in "$MAF_DIR"/split_maf/*.maf; do
  base_name=$(basename "$maf_file" .maf)
  echo "Processing $maf_file"

  # Run phyloFit for each MAF file
  phyloFit --tree "$TREE_STRUCTURE" --subst-mod HKY85 --features "$FEATURES_FILE" --catmap "NCATS = 0" --out-root "$base_name" "$maf_file"

  # Run phastCons
  phastCons --target-coverage 0.25 --expected-length 50 --rho 0.4 --most-conserved "${base_name}_most-cons.bed" --score --msa-format MAF "$maf_file" "${base_name}.background.mod" > "${base_name}_scores.wig"

  # Run phyloP
  phyloP --features "${base_name}_most-cons.bed" --msa-format MAF --subtree labirynth --mode ACC "${base_name}.background.mod" "${base_name}.maf" > "${base_name}_phyloP"

  echo "Finished processing $maf_file"
done

### STEP 5: Concatenate Results ###
cat *_most-cons.bed > concat_most-cons.bed
