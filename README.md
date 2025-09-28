# Project Title: SNP Analysis Pipeline
## Course: Biological Computing in Python (MSc Bioinformatics) 

## Overview
This project contains a Python script to analyze Single Nucleotide Polymorphisms of Plasmodium falciparum from a VCF file and classify them into different mutation categories. It uses information from genome feature format (GFF) and sequence (FASTA) files to determine the effect of SNPs at the RNA and protein level. Plasmodium falciparum is the parasite responsible for the most severe form of malaria in humans. Understanding genetic variation, such as SNPs, within its genome is important for studying drug resistance, parasite biology, and potential treatment strategies.

Getting Started

Clone this repository and navigate into the project directory.

Install dependencies from the requirements.txt file:

pip install -r requirements.txt


Run the script as shown below.

Usage

Run the script from the project root directory:

python3 scripts/AssignmentBCPYComplete.py \
  --vcfFile data/assessmentData.vcf.gz \
  --gffFile data/PlasmoDB-54_Pfalciparum3D7.gff \
  --fastaFile data/PlasmoDB-54_Pfalciparum3D7_Genome.fasta

Workflow

Logging setup
A log file is created to capture errors, warnings, statistics, and file locations.

Input validation
The script verifies that all input files exist and are in the expected formats using the custom function isExistingInRightFormat().

VCF record processing

Only high-quality VCF records are considered.

Each SNP is checked to determine if it lies in a coding region (CDS).

The transcript containing the SNP is reconstructed, and both the reference and alternate transcripts are generated.

Transcripts are translated into amino acids, and the SNP’s amino acid position is extracted.

Based on whether the reference and alternate amino acids are identical or different, the SNP is categorized as synonymous or non-synonymous.

Output generation
Results are written to multiple output files (see below).

Output Files

The script produces three key outputs:

Log File

Confirms validity of input files.

Displays statistics about VCF record quality.

Summarizes counts of each mutation type.

Lists file locations and any warnings or errors.

Bar Graph

Visualizes high-quality variants.

Categories: Non-coding, Synonymous, Non-synonymous.

Tab-Separated Table

Contains detailed mutation information with the following 9 columns:

Column	Description
CHROM	Sequence where the variant is located
POS	SNP position on the sequence
REF	Reference allele base (from VCF)
ALT	Alternate (mutated) allele base (from VCF)
Type	Mutation type: Non-coding, Synonymous, or Non-synonymous
Transcript	Transcript (mRNA) ID where the mutation occurs
Protein Location	Amino acid coordinate of the SNP in the protein
Ref AA	Reference amino acid at the SNP location
Alt AA	Alternate amino acid at the SNP location
License

This project is distributed for educational and research purposes.

Would you like me to also add an “Example Outputs” section at the end (e.g., results.log, mutation_stats.tsv, mutation_plot.png) so users know what filenames to expect?
