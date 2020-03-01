# RSEQNORM
 MIXnorm and SMIXnorm for RNA-seq data normalization
 
 ## Introduction
MIXnorm and SMIXnorm are normalization methods designed for Formalin-Fixed Paraffin-Embedded (FFPE) RNA-sequencing (RNA-seq) data. MIXnorm relies on a two-component mixture model, which models non-expressed genes by zero-inflated Poisson distributions and models expressed genes by truncated normal distributions. SMIXnorm is a simplified version of MIXnorm, which uses an simplified mixture model and requires less computation. We recommend using SMIXnorm for FFPE RNA-seq data normalization for faster computation when the number of samples is larger than 25. Though designed specifically for FFPE RNA-seq data, MIXnorm and SMIXnorm are directly applicable to normalize fresh-frozen (FF) RNA-seq data as a special case of FFPE RNA-seq normalization. To obtain the maximum likelihood estimates, we developed a nested EM algorithm, in which closed-form updates are available in each iteration.

## Dependencies

R version: 3.5.3 (2019-3-11)

Platform: x86_64-apple-darwin15.6.0 (64-bit)

Running under: macOS Mojave 10.14.5

R Packages: truncnorm_1.0-8

## Input files
Input file should be a raw read count matrix in gene level with (J genes)*(I samples).

## Example
### Install
You can install RSEQNORM from `github` using the `devtool`. 

```{r}
install.packages("devtools")
library(devtools)
install_github("S-YIN/RSEQNORM")
library(RSEQNORM)
```

### Example data
The ccRCC.RData is the clear cell renal cell carcinoma (ccRCC) data from  Eikrem et al. (2016), which contains 18,458 protein coding genes and 32 FFPE RNA-seq data.

```{r}
data(ccRCC)
head(ccRCC)
```

### Usage
RSEQNORM is implemented in R. The scripts are under folder [R](https://github.com/S-YIN/RSEQNORM/tree/master/R).  `MIXnorm` and `SMIXnorm` are the core functions that produce the normalized expression matrix, proportion of expressed genes and the probability of being expressed for each gene. 

```{r}
smix <- SMIXnorm(dat = ccRCC, max_iter = 20, tol = 0.01, appr = TRUE)
mix <- MIXnorm(dat = ccRCC, max_iter = 20, tol = 0.01, appr = TRUE)
normalized.by.smix <- smix$SMIX_normalized_log
normalized.by.mix <- mix$MIX_normalized_log
express.gene.smix <- rownames(ccRCC)[smix$D > 0.5]
express.gene.mix <- rownames(ccRCC)[mix$D > 0.5]
```

### Details
* The input data must be raw read count matrix of dimension `J genes * I samples`.
* The default maximum number of nested EM iteration (max_iter) is `20`, recommend range `(10, 50)`.
* The default convergency criteria (tol) is `0.01`, recommend range `(1e-5, 1)`.
* The default setting uses an approximate version of normalization (`appr=TRUE`). The exact normalization uses the posterior probabilities of genes being expressed or not to produce the normalized data. The approximate version uses a cut-off value of `0.5` on those probabilities to classify genes as expressed or not, then non-expressed genes will be normalized to exact `0`.  







