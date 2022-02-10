"""
Output all combination of references (5mer).

Usage:
python get_refs.py out_dir

"""
import sys
import os
import pathlib

outdir = sys.argv[1]

acceptor    = "GCTATTGTCTGCCCATGTGGCGCCCCAATTAGTGACCGCACAAGAACAGTAAGGA"
_5t         = "TGGTCTAGAG"
_3t         = "CACAGTTGCCATTCCATAGTCTCTCCTATGGAATGGCAACTGTG"
donor       = "TGTATGGATACTAGGGTACGAAGAGCTAGGTCGACTAGC"
donor_rem   = "TGTATGGA"
donor_3t    = "TGTATGGATACTAGGGTACACAGTTGCCATTCCATAG" # why Will used this?
l_odd       = "AGAG" # linker, nubmering based on insersion order
l_even      = "ACCA"

barcodes    = {
    # barcode: [barcode_seq, encoded, sequenced]
    "bc1": ["GGAGGATCTCAGGTAGGACTAACCGCTAGATCTTGG", "A", "A"],
    "bc2": ["GGAGGAAGCTAGCACCATTGTGCCAAGGAATCTTGG", "C", "C"],
    "bc3": ["GGAGGACCTATCCACTAGACGAGGTAATTGTCTTGG", "G", "G"],
    "bc4": ["GGAGGAGGCACAGAATGTGACCACAATGTGTCTTGG", "T", "T"],
    "bc5": ["GGAGGAGTAGGTAGTGCCATTCTTCGTGGATCTTGG", "U", "T"]
}

pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

for insert1 in barcodes.keys():
    for insert2 in barcodes.keys():
        for insert3 in barcodes.keys():
            for insert4 in barcodes.keys():
                for insert5 in barcodes.keys():
                    ref_name = "ref_" + barcodes[insert1][1] + barcodes[insert2][1] + barcodes[insert3][1] + barcodes[insert4][1] + barcodes[insert5][1]
                    outfile = ref_name + ".fasta"
                    with open(os.path.join(outdir, outfile), "w") as f:
                        f.write(">" + ref_name + "\n")
                        seq = ''.join([_5t, l_even, barcodes[insert5][0], l_odd, barcodes[insert4][0], l_even, barcodes[insert3][0], l_odd, barcodes[insert2][0], l_even, barcodes[insert1][0], l_odd, acceptor, barcodes[insert1][2], barcodes[insert2][2], barcodes[insert3][2], barcodes[insert4][2], barcodes[insert5][2], donor_3t])
                        f.write(seq)
