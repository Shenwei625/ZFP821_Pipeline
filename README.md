# ZFP821_Pipeline

This repository contains analysis pipelines for RNA-seq, ChIP-seq, and ATAC-seq data from the ZFP821 project, covering preprocessing, alignment, quality control, and downstream analyses.

## Requirements

Operating system:
- Linux (tested on Ubuntu 20.04)

Software dependencies:
- Fastp 0.23.4
- Bowtie2 2.5.4
- STAR 2.7.11b
- Samtools 1.21
- Deeptools 3.5.6
- featureCounts 2.0.7
- TEtranscripts 2.2.3
- R 4.4.1
- DESeq2 1.44.0

No non-standard hardware is required.
The pipeline can be run on a standard workstation.

## Input

- FASTQ files in ./data/
- Expression matrix files in ./data/

## Usage

Just follow the process in the script fold

## Output

The pipeline generates:
- BAM files
- bigwig files
- Read count tables
- Differential expression results
