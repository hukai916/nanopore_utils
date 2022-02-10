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
import sys
import pathlib

csv        = sys.argv[1]
outdir     = sys.argv[2]

match5 = {}

pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

with open(csv) as f:
    for line in f:
        cols  = line.split()
        if cols[1] == "match_5":
            _5mer = cols[4][::-1]
            outfile = "read_ids_" + _5mer + ".txt"
            if not _5mer in match5:
                match5[_5mer] = [cols[0]]
                with open(os.path.join(outdir, outfile), "w") as f:
                    f.write(_5mer + "\n")
            else:
                match5[_5mer].append(cols[0])
                with open(os.path.join(outdir, outfile), "a") as f:
                    f.write(_5mer + "\n")

for key in match5:
    print(key, len(match5[key]))
