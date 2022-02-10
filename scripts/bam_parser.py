"""
Parse BAM (sorted) to output a csv file containing reads (reference) info.
Usage:

python bam_parser.py xxx.bam ref.fasta

Result:
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

query_name:
bc1-bc5, acceptor, donor, donor_rem
cigartuples: [(0, 28), (1,1)] # Consume reference: 0, 2, 3, 7, 8 (better use get_reference_positions[-1])

Valid reads: 5x(bc)-acceptor-5nt(must match 5xbc)-donor

References:
  CIGAR:
    M 0 alignment match (can be a sequence match or mismatch) yes yes
    I 1 insertion to the reference yes no
    D 2 deletion from the reference no yes
    N 3 skipped region from the reference no yes
    S 4 soft clipping (clipped sequences present in SEQ) yes no
    H 5 hard clipping (clipped sequences NOT present in SEQ) no no
    P 6 padding (silent deletion from padded reference) no no
    = 7 sequence match yes yes
    X 8 sequence mismatch yes yes
  Read attributes:
  'aend', 'alen', 'aligned_pairs', 'bin', 'blocks', 'cigar', 'cigarstring', 'cigartuples', 'compare', 'flag', 'from_dict', 'fromstring', 'get_aligned_pairs', 'get_blocks', 'get_cigar_stats', 'get_forward_qualities', 'get_forward_sequence', 'get_overlap', 'get_reference_positions', 'get_reference_sequence', 'get_tag', 'get_tags', 'has_tag', 'header', 'infer_query_length', 'infer_read_length', 'inferred_length', 'is_duplicate', 'is_paired', 'is_proper_pair', 'is_qcfail', 'is_read1', 'is_read2', 'is_reverse', 'is_secondary', 'is_supplementary', 'is_unmapped', 'isize', 'mapping_quality', 'mapq', 'mate_is_reverse', 'mate_is_unmapped', 'mpos', 'mrnm', 'next_reference_id', 'next_reference_name', 'next_reference_start', 'opt', 'overlap', 'pnext', 'pos', 'positions', 'qend', 'qlen', 'qname', 'qqual', 'qstart', 'qual', 'query', 'query_alignment_end', 'query_alignment_length', 'query_alignment_qualities', 'query_alignment_sequence', 'query_alignment_start', 'query_length', 'query_name', 'query_qualities', 'query_sequence', 'reference_end', 'reference_id', 'reference_length', 'reference_name', 'reference_start', 'rlen', 'rname', 'rnext', 'seq', 'setTag', 'set_tag', 'set_tags', 'tags', 'template_length', 'tid', 'tlen', 'to_dict', 'to_string', 'tostring'
"""

import sys
import pysam
from Bio import SeqIO

bam = sys.argv[1]
ref = sys.argv[2]

# print("Reading reads into memory ...")
ref_dict = {}
with open(ref) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        ref_dict[record.id] = record.seq
# print("Done.")

samfile = pysam.AlignmentFile(bam, "rb")
barcode_dict = {"bc1": "A", "bc2": "C", "bc3": "G", "bc4": "T", "bc5": "U", "": ""}

current_ref = {}
done = 0
for read in samfile.fetch():
    query_name  = read.query_name
    if not read.reference_name in current_ref:
        if not (len(current_ref) == 0):
            last_read       = current_ref[list(current_ref.keys())[0]]
            last_read_ref   = list(current_ref.keys())[0]
            oligos          = "_".join(last_read[0])
            if "acceptor" in oligos:
                barcodes        = oligos.split("acceptor")[0].split("_") # only count barcodes before acceptor
            else:
                barcodes = "NA"
            # print(oligos)
            # print(barcodes)
            # print(oligos)
            # print(last_read_ref)
            barcoded_nt     = "".join([barcode_dict[x] for x in barcodes if x in barcode_dict])
            if barcoded_nt == "":
                barcoded_nt = "NA"
                barcoded_nt_rev = "NA"
                barcoded_nt_rev_no_U = "NA"
            else:
                barcoded_nt_rev = barcoded_nt[::-1]
                barcoded_nt_rev_no_U = barcoded_nt_rev.replace("U", "T")

            adapter         = last_read[1]
            adapter_rightmost = last_read[2]
            donor           = last_read[3]
            donor_leftmost  = last_read[4]
            assert  last_read_ref in ref_dict, "Record not in reference!"
            called_nt       = last_read[5]
            valid           = "NA"
            if (called_nt == barcoded_nt_rev_no_U) and (len(called_nt) == 5):
                valid       = "match_5" # perfectly valid
            elif called_nt == barcoded_nt_rev_no_U and not called_nt == "NA":
                valid       = "match_x" # matched but not 5 nt
            elif len(called_nt) == 5:
                valid       = "mis_5" # called nt is 5, but not match with barcodes

            # print(last_read)
            print(last_read_ref, valid, oligos, barcoded_nt, barcoded_nt_rev_no_U, called_nt, adapter_rightmost, donor_leftmost)
            # print(last_read_ref, oligos, adapter, adapter_rightmost, donor_leftmost, donor, called_nt)
            # exit()
            # print(current_ref)
        current_ref = {}
        current_ref[read.reference_name] = [[], "NA", "NA", "NA", "NA", "NA"]
        # [[oligo list], "exist | not_exist [acceptor]", rightmost_pos_acceptor, "exist | not_exist [donor|donor_exist]", leftmost_pos_donor, called_nt]


    if query_name in ["bc1", "bc2", "bc3", "bc4", "bc5"]:
        current_ref[read.reference_name][0].append(query_name)
    elif query_name == "acceptor" and not current_ref[read.reference_name][1] == "exist":
        # only track the first acceptor
        current_ref[read.reference_name][0].append(query_name)
        current_ref[read.reference_name][1] = "exist"
        current_ref[read.reference_name][2] = read.get_reference_positions()[-1]
    elif query_name in ["donor", "donor_rem"] and not current_ref[read.reference_name][3] == "exist": # only track the first donor
        current_ref[read.reference_name][0].append(query_name)
        current_ref[read.reference_name][3] = "exist"
        current_ref[read.reference_name][4] = read.reference_start
    if current_ref[read.reference_name][1] == "exist" and current_ref[read.reference_name][3] == "exist":
        if current_ref[read.reference_name][4] >= current_ref[read.reference_name][2]:
            called_nt = ref_dict[read.reference_name][current_ref[read.reference_name][2] + 1: current_ref[read.reference_name][4]]
        else:
            called_nt = "INVAID" # donor in front of acceptor
        current_ref[read.reference_name][5] = called_nt

samfile.close()
