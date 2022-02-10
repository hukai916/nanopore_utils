"""
Get read ids for match_5 reads.

Usage:
python get_read_ids.py input.csv outdir

Ref:
1. Accepted csv format:
col1: read_id
col2: barcoded vs called
        match_5: 5xBC matches with called 5mer
        match_x: XxBC matches with called Xmer (other than 5)
        mis_5:   called 5mer don't match with 5xBC
col3: oligo_pattern
col4: barcoded_nucleotides before acceptor
col5: reverse of col4 and also replace U with T
col6: called nucleotides
col7: rightmost mapped coordinate of acceptor on reference
col8: leftmost mapped coordinate of donor/donor_rem on reference

"""

import os
import pathlib

csv        = sys.argv[1]
outdir     = sys.argv[2]
