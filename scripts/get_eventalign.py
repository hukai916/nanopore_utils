"""
Output the eventalign values given the following input:
Raw Fast5 folder, match_fasta_ref file, outdir

Usage:

python get_eventalign.py fast5_dir match_fasta_ref.txt outdir hpc/local
"""

import sys
import os
import subprocess

fast5_dir = sys.argv[1]
match_fasta_ref_file = sys.argv[2]
outdir    = sys.argv[3]
executor  = sys.argv[4]

if not executor == "hpc":
    executor = "local"

_5mer, fasta, ref = match_fasta_ref_file.split()

nanopolish_index = "nanopolish index -d " + fast5_dir + " " + fasta
minimap2         = "minimap2 -ax map-ont " + ref + " " + fasta + " | samtools sort -o " + outdir + "/" + _5mer + ".sorted.bam -T " + _5mer + ".tmp"
samtools_index   = "samtools index " + outdir + "/" + _5mer + ".sorted.bam"
nanopolish_align = "nanopolish eventalign --reads " + fasta + " --bam " + outdir + "/" + _5mer + ".sorted.bam --genome " + ref + " --scale-events > " + outdir + "/" +_5mer + ".eventalign.txt"

commands = ["module load minimap2/2.17", "module load samtools/1.4.1", "module load nanopolish/v0.13.2", nanopolish_index, minimap2, samtools_index, nanopolish_align]

if executor == "local":
    command = " && ".join(commands)
    print(command)
