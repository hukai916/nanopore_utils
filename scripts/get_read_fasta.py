"""
Given read_ids folder, fastq file, and output dir, for each read_ids.txt, extract the fasta and save to output dir.

Usage:
python get_read_fasta.py read_ids_dir fastq outdir
"""

import sys
import pathlib
import os
import glob
from Bio import SeqIO

read_ids_dir = sys.argv[1]
fastq  = sys.argv[2]
outdir = sys.argv[3]

pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
read_ids_list = glob.glob(read_ids_dir + '/read_ids_*.txt')

read_id_fasta = {}

for read_ids in read_ids_list:
    _5mer = read_ids.split(".txt")[0].split("read_ids_")[1]
    outfile = os.path.join(outdir, "read_seqs_" + _5mer + ".fasta")
    with open(outfile, "w") as f:
        pass
    with open(read_ids) as f:
        for line in f:
            id = line.strip()
            if not id in read_id_fasta:
                read_id_fasta[id] = outfile

for record in SeqIO.parse(fastq, "fastq"):
    if record.name in read_id_fasta:
        with open(read_id_fasta[record.name][1], "a") as f:
            f.write(">" + record.name + "\n")
            f.write(record.seq + "\n")
