# Project Title: SNP Analysis Pipeline
## Course: Biological Computing in Python (MSc Bioinformatics) 

## Overview
This project contains a Python script to analyze Single Nucleotide Polymorphisms of Plasmodium falciparum from a VCF file and classify them into different mutation categories. It uses information from genome feature format (GFF) and sequence (FASTA) files to determine the effect of SNPs at the RNA and protein level. Plasmodium falciparum is the parasite responsible for the most severe form of malaria in humans. Understanding genetic variation, such as SNPs, within its genome is important for studying drug resistance, parasite biology, and potential treatment strategies.

## Getting Started

1. Clone this repository and navigate into the project directory.
2. ```pip install -r requirements.txt```
3. Run the script as ```python3 Script/Variant Sorter.py --vcfFile Data/assessmentData.vcf.gz --gffFile Data/PlasmoDB-54_Pfalciparum3D7.gff --fastaFile Data/PlasmoDB-54_Pfalciparum3D7_Genome.fasta```

## Output Files

1. **Log File**

   * Confirms validity of input files.
   * Displays statistics about VCF record quality.
   * Summarizes counts of each mutation type.
   * Lists file locations and any warnings or errors.

2. **Bar Graph**

   * Visualizes high-quality variants.
   * Categories: *Non-coding*, *Synonymous*, *Non-synonymous*.

3. **Tab-Separated Table**

   * Contains detailed mutation information such as type of mutation, reference and alt variants and their translated amino acids
