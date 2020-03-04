# RSEQNORM
RSEQNORM implements MIXnorm and SMIXnorm for RNA-seq data normalization.
 
 ## Introduction
MIXnorm and SMIXnorm are normalization methods designed for Formalin-Fixed Paraffin-Embedded (FFPE) RNA-sequencing (RNA-seq) data. MIXnorm relies on a two-component mixture model, which models non-expressed genes by zero-inflated Poisson distributions and models expressed genes by truncated normal distributions. SMIXnorm is a simplified version of MIXnorm, which uses a simplified mixture model and requires less computation. We recommend using SMIXnorm for FFPE RNA-seq data normalization for faster computation when the number of samples is larger than 25. Though designed specifically for FFPE RNA-seq data, MIXnorm and SMIXnorm are directly applicable to normalize fresh-frozen (FF) RNA-seq data as well. To obtain the maximum likelihood estimates, we developed a nested EM algorithm, in which closed-form updates are available in each iteration.

## Dependencies

R version: > 3.5.3 (2019-3-11)

Platform: x86_64-apple-darwin15.6.0 (64-bit)

Running under: macOS Mojave 10.14.5

R Packages: truncnorm_1.0-8

## Input files
Input file should be a J-by-I raw read count matrix in gene level with (J genes)*(I samples). We encourage to supply gene names as row names and sample names as column names. An example of the input data in R is shown below.

```{r}
> dim(example.data)[1] 18458    32> example.data[1:6,1:6]          NF9   TF9  NF10  TF10  NF11  TF11A1BG        0     0     1     2     0     0A1CF      595    67   292    52   259   112A2M     45347 56829 15779 39418 33654 32715A2ML1       2     0     1     3     0     0A3GALT2     3     5     2     4     1     1A4GALT    497   382   248   429   351   312
```

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
RSEQNORM is implemented in R. The scripts are under folder [R](https://github.com/S-YIN/RSEQNORM/tree/master/R).  `MIXnorm` and `SMIXnorm` are the core functions that produce the normalized expression matrix (`$MIX_normalized_log` and `$SMIX_normalized_log`), proportion of expressed genes (`$phi`) and the probability of being expressed for each gene (`$D`). 

#### Run SMIXnorm and MIXnorm
```{r}
smix <- SMIXnorm(dat = ccRCC, max_iter = 20, tol = 0.01, appr = TRUE)
mix <- MIXnorm(dat = ccRCC, max_iter = 20, tol = 0.01, appr = TRUE)
```
See ?SMIXnorm and ?MIXnorm for additional documentation.

#### Obtain normalized expression
```{r}
normalized.by.smix <- smix$SMIX_normalized_log
normalized.by.mix <- mix$MIX_normalized_log
```

The normalized expression is a matrix in the same dimension as the input matrix.
```{r}
> dim(normalized.by.smix)[1] 18458    32> normalized.by.smix[1:6,1:6]              NF9        TF9      NF10       TF10      NF11       TF11A1BG    0.0000000  0.0000000 0.0000000  0.0000000 0.0000000  0.0000000A1CF    0.8885373 -0.8772782 0.9724398 -1.1300976 0.5749542 -0.2046440A2M     5.2204180  5.8510337 4.9587658  5.4816137 5.4381895  5.4635877A2ML1   0.0000000  0.0000000 0.0000000  0.0000000 0.0000000  0.0000000A3GALT2 0.0000000  0.0000000 0.0000000  0.0000000 0.0000000  0.0000000A4GALT  0.7088968  0.8512491 0.8097201  0.9633957 0.8779038  0.8141714
```

#### Identify expressed genes
```{r}
express.gene.smix <- rownames(ccRCC)[smix$D > 0.5]
express.gene.mix <- rownames(ccRCC)[mix$D > 0.5]
```
`$D` is the probability of being expressed for each gene. `$phi$ is the proportion of expressed genes among the data. In the example data, both MIXnorm and SMIXnorm identify 78.7\% genes as being expressed.
```{r}
> smix$D[1:6]        A1BG         A1CF          A2M        A2ML1      A3GALT2       A4GALT 6.562856e-47 1.000000e+00 1.000000e+00 1.410672e-48 1.538995e-27 1.000000e+00

> round(smix$phi,3)[1] 0.787> round(mix$phi,3)[1] 0.787
```


### Details
* The input data must be raw read count matrix of dimension `J genes * I samples`.
* The default maximum number of nested EM iteration (max_iter) is `20`, recommend range `(10, 50)`.
* The default convergency criteria (tol) is `0.01`, recommend range `(1e-5, 1)`.
* The default setting uses an approximate version of normalization (`appr=TRUE`). The exact normalization uses the posterior probabilities of genes being expressed or not to produce the normalized data. The approximate version uses a cut-off value of `0.5` on those probabilities to classify genes as expressed or not, then non-expressed genes will be normalized to exact `0`.  

### Citation
Yin, S., Wang, X., Jia, G., and Xie, Y. (2020). MIXnorm: Normalizing Gene Expression Data from RNA Sequencing of Formalin-Fixed Paraffin-Embedded Samples. Bioinformatics. In Press. DOI: 10.1093/bioinformatics/btaa153.





