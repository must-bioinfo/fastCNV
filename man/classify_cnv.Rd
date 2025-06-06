% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNVclassify.R
\name{classify_cnv}
\alias{classify_cnv}
\title{CNV Classification for a Given Vector
Classifies each CNV value in a vector as "loss", "gain", or "no alteration" based on specified thresholds.}
\usage{
classify_cnv(cnv_vector, peaks)
}
\arguments{
\item{cnv_vector}{A numeric vector of CNV scores for a given chromosome arm.}

\item{peaks}{A numeric vector containing the thresholds for classification.
The default is \code{c(-0.1, 0, 0.1)}:
\itemize{
\item Loss: CNV values below the first threshold (\code{peaks[1]}).
\item Gain: CNV values above the third threshold (\code{peaks[3]}).
\item No alteration: CNV values between the first and third thresholds.
}}
}
\value{
A character vector of classifications for each CNV value in the input vector,
with possible values of \code{"loss"}, \code{"gain"}, or \code{"no_alteration"}.
}
\description{
CNV Classification for a Given Vector
Classifies each CNV value in a vector as "loss", "gain", or "no alteration" based on specified thresholds.
}
\keyword{internal}
