# nanopore_utils

## scripts/
### get_refs.py
Return reference fasta by concatenating combinations of the 5 nucleotides.

Usage:
python get_refs.py outdir

### bam_parser.py
Parse aligned BAM file and return read info as csv file.

Usage:
python bam_parser.py input.bam ref.fasta # use reads as reference and barcodes as queries.

### get_read_ids.py
Parse csv from bam_parser.py and output reference fasta to designated folder.

Usage:
python get_read_ids.py input.csv outdir

### get_read_fasta.py
Given read_ids folder, fastq file, and output dir, for each read_ids.txt, extract the fasta and save to output dir.

Usage:
python get_read_fasta.py read_ids_dir fastq

Dependencies:
- Bio

### match_fasta_ref.py
Given two input folders: fasta and ref, return the fasta-ref matching info. This is for get_eventalign.py step.

Usage:
python match_fasta_ref.py fasta_dir ref_dir out_dir


### get_eventalign.py

TODO
