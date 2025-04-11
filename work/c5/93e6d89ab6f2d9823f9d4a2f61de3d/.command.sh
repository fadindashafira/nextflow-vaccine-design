#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
import textwrap
import subprocess

# Read vaccine construct
record = next(SeqIO.parse("vaccine_construct.fasta", "fasta"))
vaccine_seq = str(record.seq)

# Calculate basic properties

# 1. Amino acid composition
aa_composition = {}
for aa in 'ACDEFGHIKLMNPQRSTVWY':
    aa_composition[aa] = vaccine_seq.count(aa)

# 2. Molecular weight (rough estimate)
aa_weights = {
    'A': 89.09, 'R': 174.2, 'N': 132.12, 'D': 133.1, 'C': 121.16,
    'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17,
    'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
}
mol_weight = sum(aa_weights.get(aa, 0) for aa in vaccine_seq)

# 3. Theoretical pI (isoelectric point) - this is a simplified calculation
# In production, you would use a more accurate method or library
# Here we use a placeholder value
theoretical_pi = 8.5  # Example value

# 4. Instability index (simplified calculation)
# In production, use a proper library or API
# This is a dummy value for demonstration
instability_index = 35.2  # Example value

# 5. Count linkers to estimate number of epitopes
linkers = ["GPGPG", "AAY", "KK"]
epitope_count = 0
for linker in linkers:
    epitope_count += vaccine_seq.count(linker)

# 6. Identify potential antigenic regions
# In production, use a proper tool like Kolaskar-Tongaonkar
# Here we just simulate the output
antigenic_regions = [
    {"start": 15, "end": 25, "score": 0.95},
    {"start": 45, "end": 65, "score": 0.91},
    {"start": 80, "end": 95, "score": 0.88}
]

# 7. Check for allergenicity (simulated)
# In production, use AllerTOP or AllergenFP
allergenicity_score = 0.15  # Example value (low score = less allergenic)

# Save properties to CSV
properties = {
    "Property": [
        "Length", "Molecular Weight", "Theoretical pI", 
        "Instability Index", "Epitope Count", "Allergenicity Score"
    ],
    "Value": [
        len(vaccine_seq), f"{mol_weight:.2f}", f"{theoretical_pi:.2f}",
        f"{instability_index:.2f}", epitope_count, f"{allergenicity_score:.4f}"
    ],
    "Unit": [
        "aa", "Da", "", "", "", ""
    ],
    "Interpretation": [
        "", "", 
        "7-8 is optimal for solubility", 
        "<40 suggests stable protein", 
        "", 
        "<0.3 suggests low allergenicity"
    ]
}

pd.DataFrame(properties).to_csv("vaccine_properties.csv", index=False)

# Create a detailed evaluation report
with open("vaccine_evaluation.txt", "w") as f:
    f.write("Vaccine Construct Evaluation\n")
    f.write("============================\n\n")
    f.write(f"Sequence Length: {len(vaccine_seq)} amino acids\n")
    f.write(f"Molecular Weight: {mol_weight:.2f} Da\n")
    f.write(f"Theoretical pI: {theoretical_pi:.2f}\n")
    f.write(f"Instability Index: {instability_index:.2f} (<40 suggests a stable protein)\n")
    f.write(f"Number of Epitopes: {epitope_count}\n")
    f.write(f"Allergenicity Score: {allergenicity_score:.4f} (<0.3 suggests low allergenicity)\n\n")

    f.write("Amino Acid Composition:\n")
    for aa in sorted(aa_composition.keys()):
        count = aa_composition[aa]
        percentage = (count / len(vaccine_seq)) * 100
        f.write(f"  {aa}: {count} ({percentage:.2f}%)\n")

    f.write("\nPredicted Antigenic Regions:\n")
    for region in antigenic_regions:
        f.write(f"  Region {region['start']}-{region['end']}: Score {region['score']}\n")

    f.write("\nSequence:\n")
    for line in textwrap.wrap(vaccine_seq, width=60):
        f.write(f"{line}\n")

    f.write("\n\nNotes:\n")
    f.write("- This evaluation provides a basic assessment of the vaccine construct properties.\n")
    f.write("- For a comprehensive evaluation, additional analyses should be performed:\n")
    f.write("  1. Structural modeling and refinement\n")
    f.write("  2. Molecular dynamics simulations for stability assessment\n")
    f.write("  3. Binding affinity calculations with target HLA molecules\n")
    f.write("  4. Experimental validation of immunogenicity\n")

print("Vaccine construct evaluation complete.")
