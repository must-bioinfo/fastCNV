#' fastCNV_10XHD calls all of the internal functions needed to compute the putative CNV on a Seurat Visium HD object or a list of Seurat Visium HD objects
#'
#' This function orchestrates the CNV analysis on a Seurat Visium HD object (or multiple objects). It calls internal functions such as
#' `CNVAnalysis` and `PlotCNVResults` to compute the CNVs and generate heatmaps. The results are saved in the metadata of the Seurat object(s), with options for
#' generating and saving plots.
#'
#' @param seuratObjHD Seurat object or list of Seurat objects to perform the CNV analysis on.
#' @param sampleName Name of the sample or a list of names corresponding to the samples in the `seuratObj`.
#' @param referenceVar The variable name of the annotations in the Seurat metadata to be used as reference.
#' @param referenceLabel The label given to the observations you want as reference (can be any type of annotation).
#' @param assay Name of the assay to run the CNV on. Takes the results of `prepareCountsForCNVAnalysis` by default if available.
#' @param pooledReference Default is `TRUE`. Will build a pooled reference across all samples if `TRUE`.
#' @param denoise If `TRUE`, the data will be denoised (default = `TRUE`).
#' @param scaleOnReferenceLabel If `TRUE`, scales the results depending on the normal observations (default = `TRUE`).
#' @param thresholdPercentile Which quantiles to take (default 0.01). For example, `0.01` will take quantiles between 0.01-0.99. Background noise appears with higher numbers.
#' @param geneMetadata List of genes and their metadata (default uses genes from Ensembl version 113).
#' @param chrArmsToForce A chromosome arm (e.g., `"8p"`, `"3q"`) or a list of chromosome arms (e.g., `c("3q", "8p", "17p")`) to force into the analysis.
#' @param genesToForce A list of genes to force into the analysis (e.g. `c("FOXP3","MUC16","SAMD15")`).
#' @param regionToForce Chromosome region to force into the analysis (vector containing chr, start, end).
#' @param windowSize Size of the genomic windows for CNV analysis (default = 150).
#' @param windowStep Step between the genomic windows (default = 10).
#' @param saveGenomicWindows If `TRUE`, saves the information of the genomic windows in the current directory (default = `FALSE`).
#' @param topNGenes Number of top expressed genes to keep (default = 7000).
#' @param getCNVPerChromosomeArm If `TRUE`, will save the CNV per chromosome arm into the metadata.
#' @param getCNVClusters If `TRUE`, will perform clustering on the CNV scores and save them in the metadata of the Seurat object as `cnv_clusters`.
#' @param tumorLabel The label within `referenceVar` that specifies the tumor/malignant population (can be any type of annotation).
#' @param k_clusters Optional. Number of clusters to cut the dendrogram into. If `NULL`, the optimal number of clusters is determined automatically using the elbow method.
#' @param h_clusters Optional. The height at which to cut the dendrogram for clustering. If both `k` and `h` are provided, `k` takes precedence.
#' @param plotDendrogram Logical. Whether to plot the dendrogram (default = `FALSE`).
#' @param plotClustersOnDendrogram Logical. Whether to highlight clusters on the dendrogram (default = `FALSE`).
#' @param plotElbowPlot Logical. Whether to plot the elbow plot used for determining the optimal number of clusters (default = `FALSE`).
#' @param doPlot If `TRUE`, will build a heatmap for each of the samples (default = `TRUE`).
#' @param clustersVar The name of the metadata column containing cluster information (default = `"cnv_clusters"`).
#' @param clusters_palette A color palette for `clustersVar`.
#' You can provide a custom palette as a vector of color codes (e.g., `c("#F8766D", "#A3A500", "#00BF7D")`).
#' @param printPlot If `TRUE`, the heatmap will be printed in the console (default = `FALSE`, the plot will only be saved in a PDF).
#' @param savePath Path to save the heatmap plot. If `NULL`, the plot won't be saved (default = `.`).
#' @param outputType Specifies the file format for saving the plot, either `"png"` or `"pdf"` (default = `"png"`).
#' @param splitPlotOnVar The name of the metadata column to split the observations during the `plotCNVResults` step, if different from `referenceVar`.
#' @param referencePalette A color palette for `referenceVar`.
#' You can provide a custom palette as a vector of color codes (e.g., `c("#FF0000", "#00FF00")`).
#'
#' @import Seurat
#'
#' @return A Seurat object or a list of Seurat objects after all the analysis is complete. Heatmaps of the CNVs for every object in `seuratObj` are generated and saved in the specified path (default = current working directory).
#'
#' @export

fastCNV_10XHD <- function(seuratObjHD,
                          sampleName,
                          referenceVar = NULL,
                          referenceLabel = NULL,
                          assay = "Spatial.016um",

                          pooledReference = TRUE,
                          denoise = TRUE,
                          scaleOnReferenceLabel = TRUE,
                          thresholdPercentile = 0.01,
                          geneMetadata = getGenes(),
                          windowSize = 150,
                          windowStep = 10,
                          saveGenomicWindows = FALSE,
                          topNGenes = 7000,
                          chrArmsToForce = NULL,
                          genesToForce = NULL,
                          regionToForce = NULL,

                          getCNVPerChromosomeArm = TRUE,

                          getCNVClusters = FALSE,
                          tumorLabel = NULL,
                          k_clusters = NULL,
                          h_clusters = NULL,
                          plotDendrogram = FALSE,
                          plotClustersOnDendrogram = FALSE,
                          plotElbowPlot = FALSE,

                          doPlot = TRUE,
                          printPlot = FALSE,
                          savePath = ".",
                          outputType = "png",
                          clustersVar = "cnv_clusters",
                          splitPlotOnVar = clustersVar,
                          referencePalette = "default"){

  if(!length(seuratObjHD)==length(sampleName)) stop(crayon::red("seuratObjHD & sampleName should have the same length"))

  options(future.globals.maxSize = 8000*1024^2)

  message(crayon::yellow(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Preparing HD object...")))
  if (length(seuratObjHD) == 1) {
    Seurat::DefaultAssay(seuratObjHD) <- assay
    assaysCells <- Seurat::Cells(seuratObjHD)
    newHDobj <- suppressWarnings(suppressMessages(subset(seuratObjHD, cells = assaysCells)))
    newHDobj@project.name = sampleName
  } else if (length(seuratObjHD) > 1) {
    newHDobj <- list()
    for (i in 1:length(seuratObjHD)) {
      Seurat::DefaultAssay(seuratObjHD[[i]]) <- assay
      assaysCells <- Seurat::Cells(seuratObjHD[[i]])
      newHDobj[[i]] <- suppressWarnings(suppressMessages(subset(seuratObjHD[[i]], cells = assaysCells)))
      newHDobj[[i]]@project.name = sampleName[i]
    }
  }
  message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done !")))
  invisible(gc())

  ## Do CNV Analysis on the seurat / list of seurat Visium HD Objects
  newHDobj <- CNVanalysis(object = newHDobj,
                          referenceVar = referenceVar,
                          referenceLabel = referenceLabel,
                          assay = assay,
                          pooledReference = pooledReference,
                          denoise = denoise,
                          scaleOnReferenceLabel = scaleOnReferenceLabel,
                          thresholdPercentile = thresholdPercentile,
                          geneMetadata = geneMetadata,
                          windowSize = windowSize,
                          windowStep = windowStep,
                          saveGenomicWindows = saveGenomicWindows,
                          topNGenes = topNGenes,
                          chrArmsToForce = chrArmsToForce,
                          genesToForce = genesToForce,
                          regionToForce = regionToForce)
  invisible(gc())

  seuratObjHD[["genomicScores"]] = Seurat::GetAssay(newHDobj, assay = "genomicScores")
  seuratObjHD$cnv_fraction = NA
  seuratObjHD$cnv_fraction[rownames(newHDobj@meta.data)] = newHDobj$cnv_fraction
  Seurat::DefaultAssay(seuratObjHD) = "Spatial.016um"
  rm(newHDobj) ; invisible(gc())

  if (getCNVPerChromosomeArm == TRUE){
    message(crayon::yellow(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Computing CNV per chromosome arm...")))
    seuratObjHD <- CNVPerChromosomeArm(seuratObjHD)
    invisible(gc())
    message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done!")))
  }

  if (getCNVClusters == TRUE){
    message(crayon::yellow(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," CNV clustering may crash with large Visium HD samples. Turn `getCNVClusters` to `FALSE` to avoid this crashing.")))
    message(crayon::yellow(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Running CNV clustering...")))
    Seurat::DefaultAssay(seuratObjHD) = assay
    seuratObjHD <- CNVcluster(seuratObj = seuratObjHD,
                           referenceVar = referenceVar,
                           tumorLabel = tumorLabel,
                           k = k_clusters,
                           h = h_clusters,
                           plotDendrogram = plotDendrogram,
                           plotClustersOnDendrogram = plotClustersOnDendrogram,
                           plotElbowPlot = plotElbowPlot)
    invisible(gc())
    message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done!")))
  }

  if (length(seuratObjHD) == 1) {
    if (doPlot == TRUE) {
      plotCNVResultsHD(seuratObjHD = seuratObjHD,
                       printPlot = printPlot,
                       savePath = savePath,
                       outputType = outputType,
                       referenceVar = referenceVar,
                       clustersVar = clustersVar,
                       splitPlotOnVar = splitPlotOnVar,
                       referencePalette = referencePalette)
      invisible(gc())
    }

  } else if (length(seuratObjHD) > 1) {
    if (doPlot == TRUE) {
      for (i in 1:length(seuratObjHD)) {
        plotCNVResultsHD(seuratObjHD = seuratObjHD[[i]],
                         printPlot = printPlot,
                         savePath = savePath,
                         outputType = outputType,
                         referenceVar = referenceVar,
                         clustersVar = clustersVar,
                         splitPlotOnVar = splitPlotOnVar,
                         referencePalette = referencePalette)
        invisible(gc())
      }
    }
  }

  return(seuratObjHD)
}
