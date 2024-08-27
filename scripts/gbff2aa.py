#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 14:47:48 2024

@author: ferenc.kagan
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Define the input and output file names
input_file = "~/Downloads/GCA_030770545.1_macOpe2_genomic.gbff.gz"
output_file = "~/Documents/Projects/CNE_PFARs/input/macOpe2_annotations/macOpe2_proteins.fasta"

# Open the output file for writing
with open(output_file, "w") as out_handle:
    # Parse the GBFF file
    for record in SeqIO.parse(input_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if 'translation' in feature.qualifiers:
                    protein_seq = feature.qualifiers['translation'][0]
                    protein_id = feature.qualifiers.get('protein_id', ['unknown'])[0]
                    description = feature.qualifiers.get('product', ['unknown'])[0]
                    seq_record = SeqRecord(Seq(protein_seq), id=protein_id, description=description)
                    SeqIO.write(seq_record, out_handle, "fasta")
                else:
                    print(f"No translation available for {feature.qualifiers.get('locus_tag', ['unknown'])[0]}")
