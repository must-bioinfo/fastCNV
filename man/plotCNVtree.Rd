% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCNVtree.R
\name{plotCNVtree}
\alias{plotCNVtree}
\title{Plot an Annotated Phylogenetic Tree with CNV Events}
\usage{
plotCNVtree(tree_data, clone_cols = NULL)
}
\arguments{
\item{tree_data}{A data frame containing tree structure and annotations,
typically produced by \code{annotateCNVtree}.}

\item{clone_cols}{a color palette to color the clones. If NULL, points are
not colored. If TRUE, clones are colored using default color palette. If a
palette is given, clones are colored following the palette, with
values passed to \href{ggplot::scale_color_manual}{\code{scale_color_manual}}.}
}
\value{
A \code{ggplot} object representing the annotated phylogenetic tree.
}
\description{
This function generates a plot of an annotated phylogenetic tree using \code{ggtree}.
It displays tip labels, tip points, and labels for CNV events associated with
each node.
}
\examples{
cnv_matrix <- structure(c(0.2, 0.4, 0, 0, 0.1, 0, 0.1, 0.2, 0.2), dim = c(3L,
3L), dimnames = list(c("Clone 1", "Clone 2", "Clone 3"), c("Region 1",
                                                          "Region 2", "Region 3")))
tree <- CNVtree(cnv_matrix)
tree_data <- annotateCNVtree(tree, cnv_matrix)
plotCNVtree(tree_data)

}
