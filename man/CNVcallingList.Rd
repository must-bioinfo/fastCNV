% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNVCallingList.R
\name{CNVcallingList}
\alias{CNVcallingList}
\title{CNVcalling but for a list of seurat objects}
\usage{
CNVcallingList(
  seuratList,
  assay = NULL,
  referenceVar = NULL,
  referenceLabel = NULL,
  scaleOnReferenceLabel = TRUE,
  thresholdPercentile = 0.01,
  genes = getGenes(),
  windowSize = 100,
  windowStep = 20,
  topNGenes = 7000
)
}
\arguments{
\item{seuratList}{List of seurat objects}

\item{assay}{Name of the assay to run the CNV on. Takes the results of \code{prepareCountsForCNVAnalysis} by default if available}

\item{referenceVar}{The variable name of the annotations in the Seurat metadata}

\item{referenceLabel}{The label given to the observations wanted as reference (can be any type of annotation)}

\item{scaleOnReferenceLabel}{If you want to scale the results depending on the reference observations}

\item{thresholdPercentile}{Which quantiles to take (if 0.01 it will take 0.01-0.99). Background noise appears with higher numbers.}

\item{genes}{List of genes from ensembl}

\item{windowSize}{Size of the genomic windows}

\item{windowStep}{Step between the genomic windows}

\item{topNGenes}{Number of top expressed genes}
}
\value{
This function returns a list of the genomic scores per genomic window
}
\description{
CNVcalling but for a list of seurat objects
}
