# MMR-D Detection Using Albayrak *et al.* (2020) Algorithm

We implemented an algorithm to determine the mismatch repair status of a patient based on NGS data. The algorithm was described by [Albayrak *et al.* (2020)](https://doi.org/10.1200/PO.20.00185) and we provide an implementation written in R. The implementation was kindly provided by Vanessa Klauenberg. 

Basically, our implementation can be applied to any genome version and the analysis has to be restricted to covered regions in the genome by providing a bed file. 

## Installation and dependencies

To use this code, you need to have R installed on your system (tested version: 4.4 or later). Additionally, you’ll need several R/Bioconductor packages, which you can install by running the following commands in your R console.

**Step 1: Install Required CRAN Packages**

Install the necessary CRAN packages by running the following command:

```
install.packages(c("writexl")) # Add more packages as needed
```

**Step 2: Install Required Bioconductor Packages**

Some functions in this code use Bioconductor packages. To install these, first install the Bioconductor Manager if you haven’t already:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

Then use `BiocManager` to install the Bioconductor packages:

```
bmpackages <- c("Biostrings", "GenomicRanges", "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install(bmpackages)
```

**Step 3: Download the Code**

Clone the repository to your local machine:

```
git clone https://github.com/your-username/mmr-detection.git
```

## Usage

To run the MMR-D detection algorithm on your dataset, load the R scripts and pass your data to the function as shown below. The code is provided as a RStudio project can be used independently of RStudio as well.

`01_homopplymers.R`: This script can be used to extract homopolymers of a user-defined length (`N > 1`) from a given reference genome. The implementation is not optimized for memory usage. Therefore, we provided precomputed regions for hg19 (`res_hg19.Rdata`) and hg38 (`res_hg38.Rdata`) for homopolymers of length 5 or longer.

`02_bedHomopolymers.R`: This script can be used to subset the homopolymers to the covered regions (e.g. panel targets, whole exome). The script needs a bed file and an R object containing the homopolymers as input.

`03_MMRstatus.R`: This script takes the output of `02_bedHomopolymers.R` and a MAF file with one or multiple samples (Tumor_Sample_Barcodes) to perform the calculation of the MMR status. The ouptut is an Excel file with MMR status (proficient, P; deficient, D) for each sample. Additionally, the number of variants, indels, and single base indels are given. Furthermore, a likelihood score is given for the MMR status, which depends on the number of the indel density in the data.

## References

1. Albayrak, A. *et al.* Clinical Pan-Cancer Assessment of Mismatch Repair Deficiency Using Tumor-Only, Targeted Next-Generation Sequencing. *JCO Precision Oncology* 1084–1097 (2020) doi:10.1200/PO.20.00185.

