# Genomic Sequence Analysis Toolkit

This repository contains a set of Python and MATLAB scripts for genomic sequence analysis. The toolkit provides functionalities for identifying homologies in genomic sequences, processing MMEJ search sequences, and analyzing BAM files to calculate the number of mutant reads and coverage for each MMEJ deletion event.

## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Usage](#usage)
- [File Descriptions](#file-descriptions)


## Introduction

The toolkit includes three main scripts:

1. **Main_Exome_100.m**: A MATLAB script for identifying homologies in genomic sequences with the maximum distance in the homologies being 100 bases.
The script uses a sliding window to examine exonic sequences, pinpointing possible homologies that suggest Microhomology-Mediated End-joining (MMEJ) events. After finding these homologies, the script adjusts the coordinates to match the reference genome and checks them against it. If needed, an optional sequence validation step compares the genomic sequence, highlighting any differences. The results are organized, removing duplicates, fixing coordinates, and filtering out events that are too large or overlapping. The final refined MMEJ events are then outputted to a designated file. 

2. **find_homolog_in_str.m**: A supporting MATLAB function for finding homologies in a given sequence using dot plot analysis.
   The script "Main_Exome_100.m" uses the "find_homolog_in_str.m" function to find homologies. It utilizes a dot plot approach, calling the my_seqdotplot function with specific parameters to generate a similarity matrix. Subsequently, the script processes the matrix to isolate homologies based on a diagonal-match in a specified window size. It addresses potential artifacts, such as truncated homologs, by identifying and removing them. The final results are presented as a matrix containing information about the start and end positions and the length of the identified homologs. Users can customize the script by adjusting parameters, providing a versatile tool for homology detection in biological sequences.

3. **Del_Read.py**: A Python script for processing MMEJ search sequences, analyzing BAM files, and outputting the number of mutant reads and coverage for each MMEJ deletion event.
   
 The input file should be a homology file with specific columns detailing the coordinates and sequences related to the homologs. The script utilizes the PySAM library to work with BAM files and extract coverage information. The script iterates through each line in the input file, calculating the start and end positions of the deletion event. It then retrieves coverage information from the BAM file for regions around the deletion event. The script counts the occurrences of canonical, imperfect type A, and imperfect type B MMEJ search sequences within the specified regions.


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

1)Main_Exome_100_Shay(chromosome)

2)python {Del_Read.py} {homolog_file.tsv} {input_chr1.bam} {output_file.tsv} {chromosome}
Replace homolog_file.tsv, input_chr1.bam, output_file.tsv, and chromosome with the actual file paths and chromosome identifier.

## File Descriptions

Main_Exome_100.m: The main MATLAB script for identifying homologies in genomic sequences.

find_homolog_in_str.m: A supporting MATLAB function for finding homologies in a given sequence.

Del_Read.py: The Python script for processing MMEJ search sequences, searching them in the BAM files, and outputting results.

/Hg19Exome/: Directory containing exome coordinates.

/Hg19Genome/: Directory containing genome coordinates.








