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
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

match_fasta_ref = {}
fasta_files = glob.glob(fasta_dir + "/reads_fasta_*.fasta")
ref_files   = glob.glob(ref_dir + "/ref_*.fasta")
print(fasta_files)
print(fastq_dir)

for fasta in fasta_files:
    _5mer = fasta.split("_")[-1].split(".fasta")[0]
    if not _5mer in match_fasta_ref:
        match_fasta_ref[_5mer] = [fasta]

for ref in ref_files:
    _5mer = ref.split("_")[-1].split(".fasta")[0]
    if _5mer in match_fasta_ref:
        match_fasta_ref[_5mer].append(ref)

for key in match_fasta_ref:
    print(key, match_fasta_ref[key])
