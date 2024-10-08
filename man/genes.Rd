% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genesData.R
\docType{data}
\name{genes}
\alias{genes}
\title{Genes data from ensemble version 111}
\format{
An object of class list.
}
\source{
https://www.ensembl.org/index.html
}
\usage{
data(genes)
}
\description{
Data downloaded from the ensembl website version 111
Data containing ensembl gene id, HUGO nomenclature and entrez gene id,
location on chromosome, biotype and length for each gene (about 76 000 genes)
}
\examples{
data(genes)
hgnc <- genes$hgnc_symbol
entrez <- genes$entrezgene_id
}
\keyword{datasets}
