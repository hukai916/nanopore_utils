"""
Given two input folders: fasta and ref, return the fasta-ref matching info. This is for get_eventalign.py step.

Usage:
python match_fasta_ref.py fasta_dir ref_dir out_dir
"""

import sys
import os
import pathlib
import glob

fasta_dir = sys.argv[1]
ref_dir   = sys.argv[2]
outdir    = sys.argv[3]

match_fasta_ref = {}
fasta_files = glob.glob(fasta_dir + "/read_seqs_*.fasta")
ref_files   = glob.glob(ref_dir + "/ref_*.fasta")

for fasta in fasta_files:
    _5mer = fasta.split("_")[-1].split(".fasta")[0]
    if not _5mer in match_fasta_ref:
        match_fasta_ref[_5mer] = [os.path.realpath(fasta)]

for ref in ref_files:
    _5mer = ref.split("_")[-1].split(".fasta")[0]
    if _5mer in match_fasta_ref:
        match_fasta_ref[_5mer].append(os.path.realpath(ref))

pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
for key in match_fasta_ref:
    if len(match_fasta_ref[key]) == 2:
        print("Outputing for ", key, " ...")
        outname = os.path.join(outdir, "match_fasta_ref_" + key + ".txt")

        with open(outname, "w") as f:
            f.write("\t".join([key, match_fasta_ref[key][0], match_fasta_ref[key][1]]))
