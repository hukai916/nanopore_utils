"""
Output the eventalign values given the following input:
Raw Fast5 folder, match_fasta_ref file, outdir

Usage:

python get_eventalign.py fast5_dir match_fasta_ref.txt outdir hpc/local
"""

import sys
import os
import subprocess
import pathlib

fast5_dir = sys.argv[1]
match_fasta_ref_file = sys.argv[2]
outdir    = sys.argv[3]
executor  = sys.argv[4]

if not executor == "hpc":
    executor = "local"

with open(match_fasta_ref_file) as f:
    for line in f:
        _5mer, fasta, ref = line.split()

nanopolish_index = "nanopolish index -d " + fast5_dir + " " + fasta
minimap2         = "minimap2 -ax map-ont " + ref + " " + fasta + " | samtools sort -o " + outdir + "/" + _5mer + ".sorted.bam -T " + _5mer + ".tmp"
samtools_index   = "samtools index " + outdir + "/" + _5mer + ".sorted.bam"
nanopolish_align = "nanopolish eventalign --reads " + fasta + " --bam " + outdir + "/" + _5mer + ".sorted.bam --genome " + ref + " --scale-events > " + outdir + "/" +_5mer + ".eventalign.txt"

commands = ["module load minimap2/2.17", "module load samtools/1.4.1", "module load nanopolish/v0.13.2", nanopolish_index, minimap2, samtools_index, nanopolish_align]
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
command  = " && ".join(commands)

if executor == "local":
    print(_5mer, " running ...")
    subprocess.call(command, shell = True)
    print(_5mer, " done")
else if executor == "hpc":
    command = 'bsub -q short -n 1 -W 4:00 -R select[rh=6] -R rusage[mem=4000] -o ' + outdir + "/" + _5mer + ".log " + '"' + command  + '"'
    print(command) 
