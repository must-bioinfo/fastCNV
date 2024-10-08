% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCNVResults.R
\name{plotCNVResults}
\alias{plotCNVResults}
\title{A function to plot the CNV results into a heatmap}
\usage{
plotCNVResults(
  seuratObj,
  referenceVar = NULL,
  splitPlotOnVar = referenceVar,
  savePath = ".",
  downsizePlot = FALSE
)
}
\arguments{
\item{seuratObj}{Seurat object containing the data used to get the genomicScores.}

\item{referenceVar}{The variable name of the annotations in the Seurat metadata}

\item{splitPlotOnVar}{The variable name on which to split the heatmap rows.}

\item{savePath}{Path to save the pdf heatmap. If \code{NULL}, plot won't be saved (default = \code{.}).}

\item{downsizePlot}{Subset the observations to speed up the plotting process (default = \code{FALSE}).}
}
\value{
This function builds a heatmap and saves it as a .pdf file in your working directory
}
\description{
A function to plot the CNV results into a heatmap
}
