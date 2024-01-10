# Genomic Sequence Analysis Toolkit

This repository contains a set of Python and MATLAB scripts for genomic sequence analysis. The toolkit provides functionalities for identifying homologies in genomic sequences, processing MMEJ search sequences, and analyzing BAM files to calculate the number of mutant reads and coverage for each MMEJ deletion event.

## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Usage](#usage)
- [File Descriptions](#file-descriptions)
- [License](#license)
- [Contact](#contact)

## Introduction

The toolkit includes three main scripts:

1. **Main_Exome_100_Shay.m**: A MATLAB script for identifying homologies in genomic sequences within specified windows.

2. **find_homolog_in_str_Shay.m**: A supporting MATLAB function used for finding homologies in a given sequence using dot plot analysis.

3. **MMEJ_Deletion_Analysis.py**: A Python script for processing MMEJ search sequences, analyzing BAM files, and outputting the number of mutant reads and coverage for each MMEJ deletion event.

   Running the file:  python Del_Read.pu {f} {b} {o} {c}

## Requirements

- **MATLAB Environment:** The MATLAB scripts require a MATLAB environment.
- **Python Environment:** The Python script requires a Python environment with the `pysam` and `numpy` libraries.

## Usage

### Genomic Sequence Analysis (MATLAB)

1. Clone the repository to your local machine.

```bash
git clone https://github.com/your-username/genomic-sequence-analysis.git
