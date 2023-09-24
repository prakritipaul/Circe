# Circe: Suite of Computational Methods to Predict Reprogramming Cocktails  

**Circe is a suite of computational methods to predict reprogramming Cocktails for neuronal cell types in *Ciona intestinalis***, a marine invertebrate that serves as a simple model organism for chordate development. **A reprogramming Cocktail is a combination of transcription factors** that when overexpressed, can induce a particular cellular identity in the cell type in which it was introduced. This process of transforming one cell type into another is referred to as **direct reprogramming**. 

![Circe](https://cocktails-r-objects.s3.amazonaws.com/method_figure.png)

This repository contains **all relevant code and analyses** pertaining to the direct reprogramming of three neuronal cell types in Ciona intestinalis- the Bipolar Tail Neurons (BTNs), Palp Sensory Cells (PSCs), and Decussating Neurons (ddNs)- **using Circe**. 

**Circe's predictions and its experimental validations** are present in Chapter 3 of Prakriti Paul Chacha's PhD Thesis "Computational Approaches to Reprogram Neuronal Cell Identities in *Ciona intestinalis*" [1]. 

**Readers are encouraged to clone this repository and run the raw Rmarkdown files in RStudio.** Each step of Circe is documented in detail in these files, so the reader will better understand the various parts of the pipeline and interact with the data structures created such as lists, data frames, etc. The contents of the Directories below will also make more sense. 



## Python Code
The **Materials and Methods Section is especially** helpful to understand what each script implements. 

### "get_pseudostates.py"
- Implements "Identification of Pseudostates" method.

###  "GRNBoost2_pipeline.py"
-  Makes 50 GRNs for a given cell state/terminal cell type, which are used for subsequent analyses, namely getting "sorted_total_importance_dfs" and "sorted_percent_top_list".
-  Refer to Computational Supplementary Figures and Results. 
*GRN is short for Gene Regulatory Network*.

### "subsampling_prediction_confidence.py"
- Gets a percentage/confidence for the overall prediction by counting the number of times that a TF in the predicted Cocktail appears in a subsample's inferred Cocktail.

## R Code
> "(BTN/PSC/ddN)_cocktail_pipeline.(Rmd/md)" consists of entire analysis pipeline using Circe:
- Identification of Pseudostates
- Calculation of Differentially Expressed Genes
- Group Lasso
- Comparison to other Models
	-- (Linear) Gaussian Multivariate Multi-Response Linear Regression
	-- (Nonlinear) GRNBoost2 

These files are also referred to as "cocktail_pipeline files".

## Directories
### DEG_lists
Contains the DEG lists generated for all three cell types: BTNs, PSCs, and ddNs.
*Note: DEG is short for Differentially Expressed Gene*.

The Subsection **"Biological Concepts Underlying the Circe Framework" of the Results Section is especially helpful** in understanding what each of these lists pertain to. 

#### These DEG lists include:
- Candidate Cocktail TFs
- DEG TFs in Mother compared to its Siblings
- Precursor DEGs compared to its Sibling (Precursor Identity Genes)
- DEG TFs in Precursor compared to Mother 
- DEG TFs in Precursor compared to its Siblings
- Terminal Candidate Cocktail TFs
- DEG TFs in Terminal Cell Type compared to its Siblings 
- Terminal DEGs compared to its Siblings


### GRNs 
- Subdirectories consist of the 50 GRNs built for precursor cell states and terminal cell types by GRNBoost2 [2].
- This directory also consists of inputs to GRNBoost2 (cocktail khids and GRN_df).
- These inputs were generated from the respective cocktail_pipeline files. 


### helper_files
- "Ciona_khid_TF.tsv" contains khid and human ortholog name. It was generated from Ghost database (ghost.zool.kyoto-u.ac.jp).
- "final_helper_functions.R" and "final_main_pipeline_start_script.R" are starter scripts for the cocktail_pipeline files.


### PC_pseudotime_dfs
- Consists of "PC_pseudotime_dfs", which contain the PC coordinate values and pseudotimes for cells belonging to a given cell population segment. 
- They were generated in the cocktail_pipeline files. 

### subsampling
- Contains "group_lassoed_tf_df" and "combined_subsampled_tf_df" csv files that are used as inputs to "subsampling_prediction_confidence.py".

## Notes
1. R and python code are used together in the cocktail_pipeline files.

2. Relevant objects such as the Seurat and URD objects are provided in and are automatically downloaded from Amazon S3 bucket "cocktails-r-objects".

3. All environments generated from running the various pipelines are provided in the same Amazon S3 bucket. 

## AWS Links to R environments
These are the R Environments generated from running the "cocktail_pipeline.Rmd" files on the various neuronal cell types. 

https://cocktails-r-objects.s3.amazonaws.com/Final_BTN_cocktail_pipeline_environment.RData.zip

https://cocktails-r-objects.s3.amazonaws.com/Final_ddN_cocktail_pipeline_environment.RData.zip

https://cocktails-r-objects.s3.amazonaws.com/Final_PSC_cocktail_pipeline_environment.RData.zip 

## Computational Supplementary Results
In addition to the Chapter 3 Supplementary Appendix present in the Thesis (which focuses solely on the BTNs), a "Computational Supplementary Results" document is available (https://cocktails-r-objects.s3.amazonaws.com/Computational+Supplementary+Results.pdf). Among its contents are comparisons of Circe to other linear and non-linear methods, TF groups selected at different values of lambda, and analyses of the terminally differentiated states of the other cell types (PSCs and ddNs).

## Citations
[1] Chacha, Prakriti Paul. Computational Approaches to Reprogram Neuronal Cell Identities in Ciona intestinalis. Diss. Princeton University, 2023.

[2] Moerman, Thomas, et al. "GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks." Bioinformatics 35.12 (2019): 2159-2161.
