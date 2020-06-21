---
layout: post
title: Single cell
subtitle: Single-cell analysis
gh-repo: leequn/leequn.github.io
gh-badge: [star, fork, follow]
tags: [Bioinformatics, Single Cell]
comments: true
---

# Single-cell analysis

## 1. Introduction
### 1.1 Why scRNA-seq
- explore which cell types are present in a tissue
- identify unknown/rare cell types or states
- elucidate the changes in gene expression during differentiation processes or across time or states
- identify genes that are differentially expressed in particular cell types between conditions (e.g. treatment or disease)
- explore changes in expression among a cell type while incorporating spatial, regulatory, and/or protein information
<br/>
![pic1](/assets/img/scrnaseq/1.jpg)

### 1.2 Challenges of scRNA-seq analysis
- Large volume of data
- Low depth of sequencing per cell
- Biological variability across cells/samples
	- Transcriptional bursting
	- Varying rates of RNA processing
	- Continuous or discrete cell identities (e.g. the pro-inflammatory potential of each individual T cell)
	- Environmental stimuli
	- Temporal changes
- Technical variability across cells/samples
	- Cell-specific capture efficiency
	- Library quality
	- Amplification bias
	- Batch effects

## 2. scRNA-seq methods
![pic2](/assets/img/scrnaseq/2.jpg)

## 3. scRNA-seq flow (BIOLOGICAL REPLICATES ARE STILL NEEDED!)
- **Generation of the count matrix (method-specific steps)**: formating reads, demultiplexing samples, mapping and quantification
- **Quality control of the raw counts**: filtering of poor quality cells
- **Clustering of filtered counts**: clustering cells based on similarities in transcriptional activity (cell types = different clusters)
- **Marker identification**: identifying gene markers for each cluster
- **Optional downstream steps**

![pic3](/assets/img/scrnaseq/3.jpg)

### 3.1 Raw data to count matrix
RNA sequences (also referred to as reads or tags), will be derived either from the 3’ ends (or 5’ ends) of the transcripts (10X Genomics, CEL-seq2, Drop-seq, inDrops) or from full-length transcripts (Smart-seq)
- 3’ (or 5’)-end sequencing:
	- More accurate quantification through use of unique molecular identifiers distinguishing biological duplicates from amplification (PCR) duplicates
	- Larger number of cells sequenced allows better identity of cell type populations
	- Cheaper per cell cost
	- Best results with > 10,000 cells
- Full length sequencing:
	- Detection of isoform-level differences in expression
	- Identification of allele-specific differences in expression
	- Deeper sequencing of a smaller number of cells
	- Best for samples with low number of cells

Many of the same analysis steps need to occur for 3’-end sequencing as for full-length, but 3’ protocols have been increasing in popularity and consist of a few more steps in the analysis. Therefore, our materials are going to detail the analysis of data from these 3’ protocols with a focus on the droplet-based methods (inDrops, Drop-seq, 10X Genomics).

#### 3.1.1 What information is present in each of the reads (3’-end reads (includes all droplet-based methods)
![pic4](/assets/img/scrnaseq/4.jpg)
- Sample index: determines which sample the read originated from
	- Added during library preparation - needs to be documented
- Cellular barcode: determines which cell the read originated from
	- Each library preparation method has a stock of cellular barcodes used during the library preparation
- Unique molecular identifier (UMI): determines which transcript molecule the read originated from
	- The UMI will be used to collapse PCR duplicates
	- Reads with different UMIs mapping to the same transcript were derived from different molecules and are biological duplicates - each read should be counted.
	- Reads with the same UMI originated from the same molecule and are technical duplicates - the UMIs should be collapsed to be counted as a single read.
	- In image below, the reads for ACTB should be collapsed and counted as a single read, while the reads for ARL1 should each be counted.

![pic6](/assets/img/scrnaseq/6.jpg)

- Sequencing read1: the Read1 sequence
- Sequencing read2: the Read2 sequence

The analysis workflow for scRNA-seq is similar for the different droplet-based scRNA-seq methods, but the parsing of the UMIs, cell IDs, and sample indices, will differ between them. For example, below is a schematic of the 10X sequence reads, where the indices, UMIs and barcodes are placed differently:

![pic5](/assets/img/scrnaseq/5.jpg)

#### 3.1.2 Generation of count matrix
After sequencing, the sequencing facility will either output the raw sequencing data as BCL or FASTQ format or will generate the count matrix. If the reads are in BCL format, then we will need to convert to FASTQ format. There is a useful command-line tool called bcl2fastq that can easily perform this conversion.

![pic7](/assets/img/scrnaseq/7.jpg)

- Formatting reads and filtering noisy cellular barcodes
- Demultiplexing sample reads
- Mapping/pseudo-mapping to cDNAs
- Collapsing UMIs and quantification of reads

If using 10X Genomics library preparation method, then the [Cell Ranger pipeline](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) would be used for all of the above steps.

### 3.2 Quality control set-up
![pic8](/assets/img/scrnaseq/8.jpg)

After quantifying gene expression we need to bring this data into R to generate metrics for performing QC. In this lesson we will talk about the format(s) count data can be expected in, and how to read it into R so we can move on to the QC step in the workflow. We will also discuss the dataset we will be using and the associated metadata.

#### 3.2.1 Exploring the example dataset
In this paper, the authors present a a computational algorithm that harnesses genetic variation (eQTL) to determine the genetic identity of each droplet containing a single cell (singlet) and identify droplets containing two cells from different individuals (doublets) [Kang et al, 2017](https://www.nature.com/articles/nbt.4042).
The data used to test their algorithm is comprised of pooled Peripheral Blood Mononuclear Cells (PBMCs) taken from eight lupus patients, split into control and interferon beta-treated (stimulated) conditions.
![pic9](/assets/img/scrnaseq/9.jpg)

##### 3.2.1.1 Raw data
This dataset is available on GEO (GSE96583), however the available counts matrix lacked mitochondrial reads, so we downloaded the BAM files from the SRA (SRP102802). These BAM files were converted back to FASTQ files, then run through Cell Ranger to obtain the count data that we will be using.

**NOTE**: The counts for this dataset is also freely available from 10X Genomics and is used as part of the [Seurat tutorial](https://satijalab.org/seurat/v3.0/immune_alignment.html).

##### 3.2.1.2 Metadata
In addition to the raw data, we also need to collect information about the data; this is known as metadata. There is often a temptation to just start exploring the data, but it is not very meaningful if we know nothing about the samples that this data originated from.

Some relevant metadata for our dataset is provided below:

- The libraries were prepared using 10X Genomics version 2 chemistry
- The samples were sequenced on the Illumina NextSeq 500
- PBMC samples from eight individual lupus patients were separated into two aliquots each.
	- One aliquot of PBMCs was activated by 100 U/mL of recombinant IFN-β for 6 hours.
	- The second aliquot was left untreated.
	- After 6 hours, the eight samples for each condition were pooled together in two final pools (stimulated cells and control cells). We will be working with these two, pooled samples.
- 12,138 and 12,167 cells were identified (after removing doublets) for control and stimulated pooled samples, respectively.

- Since the samples are PBMCs, we will expect immune cells, such as:
	- B cells
	- T cells
      	- NK cells
	- monocytes
	- macrophages
	- possibly megakaryocytes

**It is recommended that you have some expectation regarding the cell types you expect to see in a dataset prior to performing the QC. This will inform you if you have any cell types with low complexity (lots of transcripts from a few genes) or cells with higher levels of mitochondrial expression. This will enable us to account for these biological factors during the analysis workflow.**

- **Control sample** LINK: https://pan.baidu.com/s/1ChLI7QuTcXIXz9bXSdqbtQ   PASSWORD: 7ezu
- **Stimulated sample** LINK: https://pan.baidu.com/s/12SRTbCiWyetSFI4DOjKbSA  PASSWORD: oaw2 

#### 3.3 Setting up the R environment
##### 3.3.1 Your working directory

```
single_cell_rnaseq/
├── data
│   ├──ctrl_raw_feature_bc_matrix
│   │   ├──barcodes.tsv
│   │   ├──features.tsv
│   │   ├──matrix.mtx
│   ├──stim_raw_feature_bc_matrix
│   │   ├──barcodes.tsv
│   │   ├──features.tsv
│   │   ├──matrix.mtx
├── results
└── figures
```
##### 3.3.2 Loading libraries
Save the Rscript as `quality_control.R`

```R
#load packages
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
```
##### 3.3.3 Loading single-cell RNA-seq count data
Regardless of the technology or pipeline used to process your single-cell RNA-seq sequence data, the output will generally be the same. That is, for each individual sample you will have the following three files:

1. a file with the cell IDs, representing all cells quantified
2. a file with the gene IDs, representing all genes quantified
3. a matrix of counts per gene for every cell

data/ctrl_raw_feature_bc_matrix folder:


## References

{: .box-note}
- [Introduction to Single-cell RNA-seq](https://hbctraining.github.io/scRNA-seq/schedule/)
- [single-cell-tutorial](https://github.com/theislab/single-cell-tutorial)
