% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNVanalysis.R
\name{CNVanalysis}
\alias{CNVanalysis}
\title{CNVanalysis
Runs the CNV functions given a seurat object or a list of seurat objects}
\usage{
CNVanalysis(
  object,
  referenceVar = NULL,
  referenceLabel = NULL,
  doRecapPlot = TRUE,
  pooledReference = TRUE,
  scaleOnReferenceLabel = TRUE,
  assay = NULL,
  thresholdPercentile = 0.01,
  genes = getGenes(),
  windowSize = 100,
  windowStep = 20,
  topNGenes = 7000
)
}
\arguments{
\item{object}{The Seurat object or list of seurat objects containing the data to analyze}

\item{referenceVar}{The variable name of the annotations in the Seurat metadata}

\item{referenceLabel}{The label given to the observations wanted as reference (can be any type of annotation)}

\item{doRecapPlot}{Default \code{TRUE}. Will output the CNV heatmaps by annotation if \code{TRUE}.}

\item{pooledReference}{Default \code{TRUE}. Will build a pooled reference across all samples if \code{TRUE}.}

\item{scaleOnReferenceLabel}{If you want to scale the results depending on the normal observations}

\item{assay}{Name of the assay to run the CNV on. Takes the results of prepareCountsForCNVAnalysis by default if available}

\item{thresholdPercentile}{Which quantiles to take (if 0.01 it will take 0.01-0.99). Background noise appears with higher numbers.}

\item{genes}{List of genes from ensembl}

\item{windowSize}{Size of the genomic windows}

\item{windowStep}{Step between the genomic windows}

\item{topNGenes}{Number of top expressed genes to keep}
}
\value{
This function returns the genomic scores per genomic window per seurat object. If given a list with more than one seurat object, and their annotations, it will output a heatmap per cell type given
}
\description{
CNVanalysis
Runs the CNV functions given a seurat object or a list of seurat objects
}
