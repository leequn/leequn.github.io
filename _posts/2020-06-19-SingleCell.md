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

## Introduction
### Why scRNA-seq
- explore which cell types are present in a tissue
- identify unknown/rare cell types or states
- elucidate the changes in gene expression during differentiation processes or across time or states
- identify genes that are differentially expressed in particular cell types between conditions (e.g. treatment or disease)
- explore changes in expression among a cell type while incorporating spatial, regulatory, and/or protein information
<br/>
![pic1](/assets/img/scrnaseq/1.jpg)

### Challenges of scRNA-seq analysis
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

## scRNA-seq methods
![pic2](/assets/img/scrnaseq/2.jpg)

## scRNA-seq flow
![pic3](/assets/img/scrnaseq/3.jpg)
### raw data to count matrix
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

## References

{: .box-note}
- [Introduction to Single-cell RNA-seq](https://hbctraining.github.io/scRNA-seq/schedule/)
- [single-cell-tutorial](https://github.com/theislab/single-cell-tutorial)
