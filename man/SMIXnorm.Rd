% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMIXnorm.R
\name{SMIXnorm}
\alias{SMIXnorm}
\title{Produce SMIXnorm normalized expression matrix}
\usage{
SMIXnorm(dat, max_iter = 20, tol = 0.01, appr = T)
}
\arguments{
\item{dat}{input raw read count matrix. dim of dat = J genes * I samples.}

\item{max_iter}{maximum number of iterations for the nested EM algorithm default is 20, recommend range (10, 50).}

\item{tol}{convergency criteria, default is 1e-2, recommend range (1e-5,1).}

\item{appr}{binary True of False, indicates if the approximate version of normalization should be used.}
}
\value{
A list contains the normalized expression matrix (MIX_normalized_log), proportion of expressed genes (phi) and probabilities of being expressed for all genes (D).
}
\description{
Calls the SMIXnorm.mle function to obtain MLE of the simplified mixture model,
then produces the normalized expression matrix.
}
