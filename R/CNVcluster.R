#' Perform CNV Clustering with Seurat Object
#'
#' The `CNVcluster` function performs hierarchical clustering on a genomic score matrix extracted from a Seurat object.
#' It provides options for plotting a dendrogram, an elbow plot for optimal cluster determination,
#' and cluster visualization on the dendrogram. The resulting cluster assignments are stored in the Seurat object.
#'
#' @param seuratObj A Seurat object containing a "genomicScores" assay with a matrix of genomic scores for clustering.
#' @param k Optional. The number of clusters to cut the dendrogram into. If `NULL`, the optimal number of clusters is determined automatically using the elbow method.
#' @param h Optional. The height at which to cut the dendrogram for clustering. If both `k` and `h` are provided, `k` takes precedence.
#' @param referenceVar Optional. The name of the metadata column in the Seurat object that contains reference annotations. If `NULL`, all observations will be clustered.
#' @param cellTypesToCluster Optional. Cell type annotations to use for clustering. If `NULL` (default), all observations will be clustered.
#'
#' @details
#' The function computes a Manhattan distance matrix and performs hierarchical clustering using the Ward.D2 method.
#' If `k` is not provided, the elbow method is applied to determine the optimal number of clusters based on the within-cluster sum of squares (WSS).
#'
#' The clusters are assigned to the Seurat object under the metadata column `cnv_clusters`.
#'
#' @return A Seurat object with an additional metadata column, `cnv_clusters`, containing the cluster assignments.
#'
#' @importFrom proxy dist as.dist
#' @importFrom utils tail
#' @importFrom graphics abline text
#' @importFrom stats hclust cutree rect.hclust
#'
#' @export


CNVCluster <- function(seuratObj,
                       k = NULL,
                       h = NULL,
                       referenceVar = NULL,
                       cellTypesToCluster = NULL) {

  if (is.null(k)){kDetection = "automatic"}
  if (!is.null(k)){kDetection = "manual"}

  if (!is.null(referenceVar) && !is.null(cellTypesToCluster)) {
    if (length(which(seuratObj@meta.data[[referenceVar]] %in% cellTypesToCluster)) == 0) {
      message("No observations found corresponding to the cellTypesToCluster. Clustering all observations instead.")
    } else {
      seuratObj_orig <- seuratObj
      seuratObj <- suppressWarnings(suppressMessages(subset(seuratObj, cells = Seurat::Cells(seuratObj)[which(seuratObj@meta.data[[referenceVar]] %in% cellTypesToCluster)])))
    }
  } else {
    seuratObj_orig <- seuratObj
  }

  if (length(Seurat::Cells(seuratObj,assay = "genomicScores")) < 30000) {
    genomicMatrix <- t(as.matrix(Seurat::GetAssay(seuratObj, assay = "genomicScores")$data))

    dist_cos <- proxy::dist(genomicMatrix, method = "Manhattan")

    dist_matrix <- proxy::as.dist(dist_cos)
    hc <- stats::hclust(dist_matrix, method = "ward.D2")

    dist_matrix_full <- as.matrix(dist_matrix)

    find_elbow <- function(x, y, sensitivity = 0.05) {
      x_norm <- (x - min(x)) / (max(x) - min(x))
      y_norm <- (y - min(y)) / (max(y) - min(y))

      slopes <- diff(y_norm) / diff(x_norm)
      slope_changes <- diff(slopes)

      significant_changes <- which(abs(slope_changes) > sensitivity * max(abs(slope_changes)))

      if (length(significant_changes) == 0) {
        return(x[ceiling(length(x)/2)])
      }

      first_quarter <- ceiling(length(x) / 4)
      elbow_point <- significant_changes[significant_changes > first_quarter][1]

      if (is.na(elbow_point)) {
        elbow_point <- utils::tail(significant_changes, 1)
      }

      return(x[elbow_point + 1])
    }

    k_values <- 1:15
    wss <- sapply(k_values, function(k) {
      clusters <- cutree(hc, k = k)
      sum(sapply(unique(clusters), function(cl) {
        sum(dist_matrix_full[which(clusters == cl), which(clusters == cl)]^2) / (2 * sum(clusters == cl))
      }))
    })

    if(kDetection == "automatic") {k <- find_elbow(k_values, wss)}

    clusters <- cutree(hc, k = k, h = h)

    seuratObj$cnv_clusters = as.factor(clusters)
    seuratObj_orig@meta.data$cnv_clusters <- 0
    seuratObj_orig@meta.data[Seurat::Cells(seuratObj), "cnv_clusters"] = seuratObj$cnv_clusters
    seuratObj_orig$cnv_clusters <- as.factor(seuratObj_orig$cnv_clusters)

  } else {

    mat <- as.matrix(Seurat::GetAssayData(seuratObj, assay = "genomicScores", layer = "data"))
    pca_res <- prcomp(t(mat), center = TRUE, scale. = TRUE)
    X <- pca_res$x[, 1:30]

    set.seed(42)

    if (kDetection == "automatic") {
      wss <- numeric(15)
      for (ktest in 1:15) {
        km_temp <- kmeans(X, centers = ktest, iter.max = 100)
        wss[ktest] <- km_temp$tot.withinss
      }
      diff1 <- diff(wss)
      diff2 <- diff(diff1)
      k <- which.min(diff2) + 1
    }

    km <- kmeans(X, centers = k, iter.max = 100)

    seuratObj_orig$cnv_clusters <- NA
    seuratObj_orig$cnv_clusters[colnames(seuratObj_orig@assays$genomicScores)] <- 0
    seuratObj_orig$cnv_clusters[colnames(mat)] <- km$cluster
    seuratObj_orig$cnv_clusters <- as.factor(seuratObj_orig$cnv_clusters)

  }

  return(seuratObj_orig)
}




