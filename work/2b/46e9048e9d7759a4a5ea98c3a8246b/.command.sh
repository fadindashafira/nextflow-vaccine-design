#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv("combined_epitopes.csv")
top_epitopes = df.sort_values("consensus_score", ascending=False).head(10)
construct_seq = "GPGPG".join(top_epitopes["sequence"].dropna().tolist())

with open("vaccine_construct.fasta", "w") as f:
    f.write(">Vaccine_Construct\n" + construct_seq + "\n")
