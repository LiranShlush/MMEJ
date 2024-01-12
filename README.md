# Genomic Sequence Analysis Toolkit

This repository contains a set of Python and MATLAB scripts for genomic sequence analysis. The toolkit provides functionalities for identifying homologies in genomic sequences, processing MMEJ search sequences, and analyzing BAM files to calculate the number of mutant reads and coverage for each MMEJ deletion event.

## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Usage](#usage)
- [File Descriptions](#file-descriptions)


## Introduction

The toolkit includes three main scripts:

1. **Main_Exome_100_Shay.m**: A MATLAB script for identifying homologies in genomic sequences within specified windows.

2. **find_homolog_in_str_Shay.m**: A supporting MATLAB function for finding homologies in a given sequence using dot plot analysis.

3. **Del_Read.py**: A Python script for processing MMEJ search sequences, analyzing BAM files, and outputting the number of mutant reads and coverage for each MMEJ deletion event.

   

## Requirements

- **MATLAB Environment:** 
- **Python Environment:** 
      pandas==0.24.2
      numpy==1.16.3
      regex==2021.9.24
      biopython==1.73  
      pysam==0.11.2.2
  
## Usage
Run the Del_Read.py script by providing the "homolog file" generated from the MATLAB scripts, BAM file per chromosome, output_file_name,  and chromosome  as command-line arguments.

**Command-line:**
python Del_Read.py homolog_file.tsv input_chr1.bam.bam output_file.tsv chromosome

## File Descriptions

Main_Exome_100.m: The main MATLAB script for identifying homologies in genomic sequences.

find_homolog_in_str.m: A supporting MATLAB function for finding homologies in a given sequence.

Del_Read.py: The Python script for processing MMEJ search sequences, and searching them in the BAM files, and outputting results.

/Hg19Exome/: Directory containing exome coordinates.

/Hg19Genome/: Directory containing genome coordinates.








