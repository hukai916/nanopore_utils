# nanopore_utils

## scripts/
### get_refs.py
Return reference fasta by combining the 5 nucleotides.

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

###


bash/get_current.sh
  Extract the current value.
