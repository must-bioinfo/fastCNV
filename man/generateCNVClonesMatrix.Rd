% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateCNVClonesMatrix.R
\name{generateCNVClonesMatrix}
\alias{generateCNVClonesMatrix}
\title{Generate CNV Matrix for CNV Clusters by Chromosome Arm}
\usage{
generateCNVClonesMatrix(
  seuratObj,
  healthyClusters = NULL,
  values = "scores",
  cnv_thresh = 0.15
)
}
\arguments{
\item{seuratObj}{A Seurat object containing CNV data and metadata.}

\item{healthyClusters}{A numeric vector or \code{NULL}. If provided, clusters specified in this vector
will be labeled as "Benign" instead of "Clone". Default is \code{NULL}.}

\item{values}{one of 'scores' or 'calls'. 'scores' returns the mean CNV score per cluster,
while 'calls' uses \code{cnv_thresh} to establish a cut-off for gains and losses, returning a matrix
of CNV calls (0=none, 1=gain, -1=loss).}
}
\value{
A matrix of CNVs with row names corresponding to the clone or benign labels and columns representing
the chromosome arms, with values corresponding to CNV scores or CNV calls.
}
\description{
This function generates a matrix of metacells where each metacell corresponds to a CNV cluster.
The CNV matrix is calculated by chromosome arm. If specified, certain clusters will be labeled as "Benign"
rather than "Clone".
}
