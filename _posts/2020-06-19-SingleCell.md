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

- **Control sample** LINK: `https://pan.baidu.com/s/1ChLI7QuTcXIXz9bXSdqbtQ`   PASSWORD: 7ezu
- **Stimulated sample** LINK: `https://pan.baidu.com/s/12SRTbCiWyetSFI4DOjKbSA`  PASSWORD: oaw2 

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
library(stringr)
library(ggplot2)
```
##### 3.3.3 Loading single-cell RNA-seq count data
Regardless of the technology or pipeline used to process your single-cell RNA-seq sequence data, the output will generally be the same. That is, for each individual sample you will have the following three files:

1. a file with the cell IDs, representing all cells quantified
2. a file with the gene IDs, representing all genes quantified
3. a matrix of counts per gene for every cell
**DONOT UNZIP THESE FILES WHEN YOU DO scRNA-seq PROCESSING, YOU CAN UNZIP THESE FILES TO SEE INFORMATION OF THESE FILES**
For example: data/ctrl_raw_feature_bc_matrix folder:
- `barcodes.tsv:` This is a text file which contains all cellular barcodes present for that sample. Barcodes are listed in the order of data presented in the matrix file (i.e. these are the column names).

![pic10](/assets/img/scrnaseq/10.jpg)

- `features.tsv:` This is a text file which contains the identifiers of the quantified genes. The source of the identifier can vary depending on what reference (i.e. Ensembl, NCBI, UCSC) you use in the quantification methods, but most often these are official gene symbols. The order of these genes corresponds to the order of the rows in the matrix file (i.e. these are the row names).

![pic11](/assets/img/scrnaseq/11.jpg)

- `matrix.mtx:` This is a text file which contains a matrix of count values. The rows are associated with the gene IDs above and columns correspond to the cellular barcodes. Note that there are many zero values in this matrix.

![pic12](/assets/img/scrnaseq/12.jpg)

Loading this data into R requires us to **use functions that allow us to efficiently combine these three files into a single count matrix**. However, instead of creating a regular matrix data structure, the functions we will use create a **sparse matrix** to improve the amount of space, memory and CPU required to work with our huge count matrix.

Different methods for reading in data include:
- `readMM()`: This function is from the Matrix package and will turn our standard matrix into a sparse matrix. The features.tsv file and barcodes.tsv must first be individually loaded into R and then they are combined. For specific code and instructions on how to do this please see our additional material.
- `Read10X()`: This function is from the Seurat package and will use the Cell Ranger output directory as input. In this way individual files do not need to be loaded in, instead the function will load and combine them into a sparse matrix for you. *We will be using this function to load in our data!*

##### 3.3.3.1 Reading in a single sample (`read10X()`)
When working with 10X data and its proprietary software Cell Ranger, you will always have an `outs` directory. Within this directory you will find a number of different files including:
- **web_summary.html**: report that explores different QC metrics, including the mapping metrics, filtering thresholds, estimated number of cells after filtering, and information on the number of reads and genes per cell after filtering.
- **BAM alignment files**: files used for visualization of the mapped reads and for re-creation of FASTQ files, if needed
- **filtered_feature_bc_matrix**: folder containing all files needed to construct the count matrix using data filtered by Cell Ranger
- **raw_feature_bc_matrix**:  folder containing all files needed to construct the count matrix using the raw unfiltered data

We are mainly interested in the `raw_feature_bc_matrix` as we wish to perform our own QC and filtering while accounting for the biology of our experiment/biological system.

If we had a single sample, we could generate the count matrix and then subsequently create a [Seurat object](https://github.com/satijalab/seurat/wiki/Seurat):

```R
#How to read in 10X data for a single sample (output is a sparse matrix)
ctrl_counts <- Read10X(data.dir="data/ctrl_raw_feature_bc_matrix")
#Turn count matrix into a Seurat object (output is a Seurat object)
ctrl <- CreateSeuratObject(counts = ctrl_counts, min.features = 100)
```
**NOTE**: The `min.features` argument specifies the minimum number of genes that need to be detected per cell. This argument will filter out poor quality cells that likely just have random barcodes encapsulated without any cell present. Usually, cells with less than 100 genes detected are not considered for analysis.

**Seurat automatically creates some metadata** for each of the cells when you use the Read10X() function to read in data. This information is stored in the meta.data slot within the Seurat object (see more in the note below).

The Seurat object is a custom list-like object that has well-defined spaces to store specific information/data. You can find more information about the slots in the Seurat object at this **[link](https://github.com/satijalab/seurat/wiki/Seurat).**

```
#Explore the metadata
head(ctrl@meta.data)
```

metadata:

```
                    orig.ident nCount_RNA nFeature_RNA
AAACATACAATGCC-1 SeuratProject       2344          874
AAACATACATTTCC-1 SeuratProject       3125          896
AAACATACCAGAAA-1 SeuratProject       2578          725
AAACATACCAGCTA-1 SeuratProject       3261          979
AAACATACCATGCA-1 SeuratProject        746          362
AAACATACCTCGCT-1 SeuratProject       3519          866
```

What do the columns of metadata mean?
- **orig.ident**: this often contains the sample identity if known, but will default to “SeuratProject”
- **nCount_RNA**: number of UMIs per cell
- **nFeature_RNA**: number of genes detected per cell


##### 3.3.3.2 Reading in multiple samples with a `for loop`
In practice, you will likely have several samples that you will need to read in data for using one of the 2 functions we discussed earlier (`Read10X()` or `readMM()`). So, to make the data import into R more efficient we can use a `for loop`, that will interate over a series of commands for each of the inputs given.
In R, it has the following structure/syntax:
```
## DO NOT RUN
for (variable in input){
	command1
	command2
	command3
}
```
The for loop we will be using today will iterate over the two sample “files” and execute two commands for each sample - **(1)** read in the count data (`Read10X()`) and **(2)** create the Seurat objects from the read in data (`CreateSeuratObject()`):
```
# Create each individual Seurat object for every sample
for (file in c("ctrl_raw_feature_bc_matrix", "stim_raw_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}
```

Now that we have created both of these objects, let’s take a quick look at the metadata to see how it looks:

```
# Check the metadata in the new Seurat objects
head(ctrl_raw_feature_bc_matrix@meta.data)
head(stim_raw_feature_bc_matrix@meta.data)
```

metadata:

```
> head(ctrl_raw_feature_bc_matrix@meta.data)
                                 orig.ident nCount_RNA nFeature_RNA
AAACATACAATGCC-1 ctrl_raw_feature_bc_matrix       2344          874
AAACATACATTTCC-1 ctrl_raw_feature_bc_matrix       3125          896
AAACATACCAGAAA-1 ctrl_raw_feature_bc_matrix       2578          725
AAACATACCAGCTA-1 ctrl_raw_feature_bc_matrix       3261          979
AAACATACCATGCA-1 ctrl_raw_feature_bc_matrix        746          362
AAACATACCTCGCT-1 ctrl_raw_feature_bc_matrix       3519          866
> head(stim_raw_feature_bc_matrix@meta.data)
                                 orig.ident nCount_RNA nFeature_RNA
AAACATACCAAGCT-1 stim_raw_feature_bc_matrix       1221          606
AAACATACCCCTAC-1 stim_raw_feature_bc_matrix       1782          807
AAACATACCCGTAA-1 stim_raw_feature_bc_matrix       1451          605
AAACATACCCTCGT-1 stim_raw_feature_bc_matrix       1549          747
AAACATACGAGGTG-1 stim_raw_feature_bc_matrix       1303          558
AAACATACGCGAAG-1 stim_raw_feature_bc_matrix       5445         1330
```

Next, we need to merge these objects together into a single Seurat object. This will make it easier to run the QC steps for both sample groups together and enable us to easily compare the data quality for all the samples.

We can use the `merge()` function from the Seurat package to do this:

```
# Create a merged Seurat object
merged_seurat <- merge(x = ctrl_raw_feature_bc_matrix, 
                       y = stim_raw_feature_bc_matrix, 
                       add.cell.id = c("ctrl", "stim"))
```

Because the same cell IDs can be used for different samples, we add a **sample-specific prefix** to each of our cell IDs using the `add.cell.id` argument. If we look at the metadata of the merged object we should be able to see the prefixes in the rownames:
```
# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)
```
#### 3.3.4 Generating quality metrics
Remember that Seurat automatically creates some metadata for each of the cells:

![pic13](/assets/img/scrnaseq/13.jpg)

We need to calculate some additional metrics for plotting:
- **number of genes detected per UMI**: this metric with give us an idea of the complexity of our dataset (more genes detected per UMI, more complex our data)
- **mitochondrial ratio**: this metric will give us a percentage of cell reads originating from the mitochondrial genes

The number of genes per UMI for each cell is quite easy to calculate, and we will log10 transform the result for better comparison between samples.

```
# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
```

Seurat has a convenient function that allows us to calculate the **proportion of transcripts mapping to mitochondrial genes**. The `PercentageFeatureSet()` will take a pattern and search the gene identifiers. For each column (cell) it will take the sum of the counts slot for features belonging to the set, divide by the column sum for all features and multiply by 100. `Since we want the ratio value for plotting, we will reverse that step by then dividing by 100.`

`**NOTE:**` The pattern provided (“^MT-“) works for human gene names. You may need to adjust depending on your organism of interest. If you weren’t using gene names as the gene ID, then this function wouldn’t work. We have [code available](https://github.com/hbctraining/scRNA-seq/blob/master/lessons/mitoRatio.md) to compute this metric on your own.

`We need to add additional information to our metadata for our QC metrics, such as the cell IDs, condition information, and various metrics.` While it is quite easy to **add information directly to the metadata slot in the Seurat object using the $ operator**, we will extract the dataframe into a separate variable instead. In this way we can continue to insert additional metrics that we need for our QC analysis without the risk of affecting our `merged_seurat` object.

We will create the metadata dataframe by extracting the meta.data slot from the Seurat object:

```
#Create metadata dataframe
metadata <- merged_seurat@meta.data
```

You should see each cell ID has a `ctrl_` or `stim_` prefix as we had specified when we merged the Seurat objects. These prefixes should match the sample listed in `orig.ident`. Let’s begin by **adding a column with our cell IDs** and **changing the current column names** to be more intuitive:

```
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- dplyr::rename(metadata, seq_folder = orig.ident,nUMI = nCount_RNA,nGene = nFeature_RNA)
```

Now, let’s get **sample names for each of the cells** based on the cell prefix:

```
# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"
````

Now you are **all setup with the metrics you need to assess the quality of your data!** Your final metadata table will have rows that correspond to each cell, and columns with information about those cells:

![pic14](https://raw.githubusercontent.com/leequn/leequn.github.io/master/assets/img/scrnaseq/14.jpg)

**Saving the updated metadata to our Seurat object**

Before we assess our metrics we are going to save all of the work we have done thus far back into our Seurat object. We can do this by simply assigning the dataframe into the meta.data slot:

```
# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata
# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat.RData")
```

#### 3.3.5 Assessing the quality metrics
Now that we have generated the various metrics to assess, we can explore them with visualizations. We will assess various metrics and then decide on which cells are low quality and should be removed from the analysis:
- Cell counts
- UMI counts per cell
- Genes detected per cell
- UMIs vs. genes detected
- Mitochondrial counts ratio
- Novelty

**What about doublets?** In single-cell RNA sequencing experiments, doublets are generated from two cells. They typically arise due to errors in cell sorting or capture, especially in droplet-based protocols involving thousands of cells. Doublets are obviously undesirable when the aim is to characterize populations at the single-cell level. In particular, they can incorrectly suggest the existence of intermediate populations or transitory states that do not actually exist. Thus, it is desirable to remove doublet libraries so that they do not compromise interpretation of the results.

**Why aren’t we checking for doublets?** Many workflows use maximum thresholds for UMIs or genes, with the idea that a much higher number of reads or genes detected indicate multiple cells. While this rationale seems to be intuitive, it is not accurate. Also, many of the tools used to detect doublets tend to get rid of cells with intermediate or continuous phenotypes, although they may work well on datasets with very discrete cell types. [Scrublet](https://github.com/AllonKleinLab/scrublet) is a popular tool for doublet detection, but we haven’t adequately benchmarked it yet. Currently, we recommend not including any thresholds at this point in time. When we have identified markers for each of the clusters, we suggest exploring the markers to determine whether the markers apply to more than one cell type.

##### 3.3.5.1 Cell counts
The cell counts are determined by the number of unique cellular barcodes detected. For this experiment, between 12,000 -13,000 cells are expected.

In an ideal world, you would expect the number of unique cellular barcodes to correpsond to the number of cells you loaded. However, this is not the case as capture rates of cells are only a proportion of what is loaded. For example, the inDrops cell **capture efficiency** is higher (70-80%) compared to 10X which can drop is between 50-60%.

**NOTE**: The capture efficiency could appear much lower if the cell concentration used for library preparation was not accurate. Cell concentration should NOT be determined by FACS machine or Bioanalyzer (these tools are not accurate for concentration determination), instead use a hemocytometer or automated cell counter for calculation of cell concentration.

The cell numbers can also vary by protocol, **producing cell numbers that are much higher than what we loaded**. For example, during the inDrops protocol, the cellular barcodes are present in the hydrogels, which are encapsulated in the droplets with a single cell and lysis/reaction mixture. While each hydrogel should have a single cellular barcode associated with it, occasionally a hydrogel can have more than one cellular barcode. Similarly, with the 10X protocol there is a chance of obtaining only a barcoded bead in the emulsion droplet (GEM) and no actual cell. Both of these, in addition to the presence of dying cells can lead to a higher number of cellular barcodes than cells.

```
# Visualize the number of cell counts per sample
ggplot(metadata,aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
```
![pic15](/assets/img/scrnaseq/15.jpg)







## References

{: .box-note}
- [Introduction to Single-cell RNA-seq](https://hbctraining.github.io/scRNA-seq/schedule/)
- [single-cell-tutorial](https://github.com/theislab/single-cell-tutorial)
