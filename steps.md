# Preprocessing for nanopore sequencing data
For detailed usage of each .py code, check out the respective docstring comments on top.

## Overview: two major parts:
### Part1: obtaining match_5 reads
Strategy:  
Map oligos to all reads with BWA/bowie2 and extract valid match_5 reads by scanning the oligo mapping pattern.  
In practice, the mapping is really slow given millions of reads as reference, therefore, split the reads into 100k chunks and perform mapping for each chunk separately.  
Sort the BAM files and parse using custom Python script, which output a text file containing read_id, oligos_pattern, called 5mer vs barcoded_5mer, etc..;
A read is labeled as “match_5” only if:  
i. It has an oligo mapping pattern as “5xBC_acceptor_5nt_donor” (column3)  
ii. 5xBC encoded nucleotides match with called 5nt.

### Part2: extract eventalign values
Strategy:  
Use nanopolish eventalign command. Need to prepare reads_fasta and ref_fasta for each 5mer per eventalign command requirement.

## Part1: step-by-step
**Step1:** split reference reads into chunks, each with 100k reads.
```
cd /home/kh45w/project/umw_paul_kaufman/kai/for_will/kai
mkdir split_100k
split --lines=200000 5R_BrdU split_100k/
```

**Step2:** map oligos to reference chunks with bowtie2.
```
cd /home/kh45w/project/umw_paul_kaufman/kai/for_will/kai/split_100k; bsub < bowtie2.bsub
```
A copy of the bsub code is as below. Note that the resulting BAM must be position sorted per requirement of bam_parser.py.
```
#!/bin/bash
#BSUB -q short
#BSUB -W 4:00
#BSUB -n 4
#BSUB -R "rusage[mem=10000]"
#BSUB -R "span[hosts=1]"
#BSUB -J "bowtie2[1-15]"
#BSUB -o bowtie2_log/out.%J.%I

module load samtools/1.9
module load bowtie2/2.4.1

samples=(aa  ab	 ac ad ae af ag ah ai aj ak al am an ao)
sample=${samples[$LSB_JOBINDEX - 1]}

cd /home/kh45w/project/umw_paul_kaufman/kai/for_will/kai/split_100k
bowtie2-build -f ${sample} ${sample}
bowtie2 -a -N 1 -L 20 -x  ${sample} -f ../oligos.fa -S bowtie2_${sample}.sam

samtools sort bowtie2_${sample}.sam -o bowtie2_${sample}.bam
samtools index bowtie2_${sample}.bam

python bam_parser.py bowtie2_${sample}.bam ${sample} > bowtie2_${sample}.txt

```
Instruction on how to read output from `bam_parser.py`:
```
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
```
A piece of sample output of `bam_parser.py`.
```
eb0f752b-7c9d-46e8-ac3d-6708514c745d NA bc5 NA NA NA NA NA
ac68ab14-b5bd-4ff0-9639-e59e3c415bf8 mis_5 bc5_bc4_bc3_bc5_bc3_acceptor_donor UTGUG GTGTT ATGTT 290 296
dbb667b3-c150-40ef-a7df-719cf43886e2 NA bc3_bc1_bc2 NA NA NA NA NA
a60cf84a-76e1-4232-b9f5-32e28cd46ce4 match_5 bc1_bc1_bc1_bc2_bc2_acceptor_donor_rem_bc1 AAACC CCAAA CCAAA 287 293
5cf0c283-b17e-412f-9017-d14b61299b6e NA bc3_bc4_bc1_bc5_bc3 NA NA NA NA NA
65025da6-0b02-44da-b076-10259563533d NA bc4_bc2_bc3_bc2_bc1_acceptor_donor_rem TCGCA ACGCT AC 292 295
```

**Step3:** merge all results
```
cd /home/kh45w/project/umw_paul_kaufman/kai/for_will/kai/split_100k
cat bowtie2_*.txt > /project/umw_william_flavahan/Kai/nanopore/res/bowtie2_5r_BrdU.txt

# 5r_dU results are here: /home/kh45w/project/umw_paul_kaufman/kai/for_will/kai/5r_dU, merged result also saved to /project/umw_william_flavahan/Kai/nanopore/res
```

## Part2: step-by-step
**Step1:** prepare reference fasta for each 5mer
```
cd /home/kh45w/project/umw_willian_flavahan/nanopore/nanopore_utils/scripts
python get_refs.py refs
```
This creates a refs folder containing all possible 5mer references.  
The reference has naming format like `ref_AGUTA.fasta`.

**Step2:** prepare read fasta for each 5mer  
**Step2-1:** prepare read id list for each 5mer
```
cd /home/kh45w/project/umw_willian_flavahan/nanopore/nanopore_utils/scripts
python get_read_ids.py /project/umw_william_flavahan/Kai/nanopore/res/bowtie2_5r_BrdU.txt reads_id_5r_BrdU
```
This takes as input the output of bam_parser.py and output read_ids for each 5mer into reads_id_5r_BrdU folder.  
The reads_id files have naming format like `read_ids_CGTGU.txt`.  

**Step2-2:** prepare read fasta for each 5mer
```
cd /home/kh45w/project/umw_willian_flavahan/nanopore/nanopore_utils/scripts
cat /project/umw_william_flavahan/Nanopore/analysis/mm2/SPBS/Runs_for_basecaller/5r_BrdU/fastq/*.fastq > 5r_BrdU.fastq

python get_read_fasta.py reads_id_5r_BrdU 5r_BrdU.fastq reads_fasta_5r_BrdU
```
This takes as input the reads_id folder from Step2-2, the fastq file, and output reads_fasta into reads_fasta_5r_BrdU for each 5mer.   
The reads_fasta files have naming format like `read_seqs_AGCAA.fasta`.

**Step3:** pair up read fasta with ref for Step4  
This is to prepare for the input to Step5.
```
cd /home/kh45w/project/umw_willian_flavahan/nanopore/nanopore_utils/scripts
python match_fasta_ref.py reads_fasta_5r_BrdU refs match_fasta_ref_5r_BrdU
```
This creates paired reads_fasta with its ref (both with full path) into individual file.  
The output files have structure like:
```
more  match_fasta_ref_5r_BrdU/match_fasta_ref_CAUTU.txt

# below is command output
CAUTU	/project/umw_william_flavahan/Kai/nanopore/nanopore_utils/scripts/reads_fasta_5r_BrdU/read_seqs_CAUTU.fasta	/project/umw_william_flavahan/Kai/nanopore/nanopore_utils/
scripts/refs/ref_CAUTU.fasta
```
**Step4:**: get subset fast5 of match_5 reads  
This is to prepare for the input to Step5.
To subset fast5, the tools are installed in docker: hukai916/nanopore:0.1.
Note that fast5 are compressed in VBZ format, which is not compatible with downstream eventaligner unless you install the VBZ plugin. To solve it, choose to convert fast5 to gzip format. But a better solution would be to use a container with the plugin installed (pending task).
```
cat /project/umw_william_flavahan/Kai/nanopore/res/bowtie2_5r_BrdU.txt | grep "match_5" | cut -f 1 -d " " > read_ids_match5_5r_BrdU.txt
singularity exec ../../nanopore_0.1.sif fast5_subset -i /project/umw_william_flavahan/Nanopore/analysis/mm2/SPBS/Runs_for_basecaller/5r_BrdU/fast5 -s fast5_match5_5r_BrdU --read_id_list read_ids_match5_5r_BrdU.txt

singularity exec ../../nanopore_0.1.sif compress_fast5 -i fast5_match5_5r_BrdU -s fast5_match5_5r_BrdU_gzip -c gzip

```

**Step5**: get eventalign files:
```
cd /home/kh45w/project/umw_willian_flavahan/nanopore/nanopore_utils/scripts
for x in match_fasta_ref_5r_BrdU/*; do python get_eventalign.py fast5_match5_5r_BrdU_gzip $x eventalign_5r_BrdU hpc; done

```
Basically, get_eventalign.py takes as input of the fast5_match5 folder, the full paths to paired read_fasta and ref that are stored in `match_fasta_ref_xxxxx.txt`, and which executor to leverage (either local or hpc).
The eventalign results will be stored at `eventalign_5r_BrdU` folder.
