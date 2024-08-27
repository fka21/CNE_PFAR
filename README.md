# Conserved non-coding element analysis in _M.opercularis_
Repository containing scripts used for analysing conserved non-coding elements in the _M.opercularis_ genome.

## Aim

The aim of the project is to find conserved non-coding elements (CNE) in the _M.opercularis_ genome which can be associated with the evolution of the labyrinth organ. This is done by whole genome alignments of multiple teleost species with and without the labyrinth organ, calling the CNEs and testing acceleration of the CNE sequences in the lineage with labyrinth organ. Some subsequent analyses for mostly visualization are performed with custom R scripts.

## Repo structure

`scripts/` - contains the scripts used during the analysis.

`misc/` - contains ancilliary files for the analysis.

`annoforge/` - contains the AnnotationForge package constructed from the _M.opercularis_ genome.


### Order of scripts

1. OrthoFinder2 was used to generate the species tree:
  ```
  orthofinder -t 20 -a 20 -f /projects/ferenc.kagan/CNE_PFARs/00_raw_data/orthofinder
  ```
2. `scripts/MAF_Processing_and_Phylo_Analysis.sh`
3. `scripts/CNE_Annotation_and_Enrichment_Analysis.R`
