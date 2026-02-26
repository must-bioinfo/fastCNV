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
#' @param denoise If `TRUE`, the denoised data will be used in the heatmap (default = `TRUE`).
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
#' @param k_clusters Optional. Number of clusters to cut the dendrogram into. If `NULL`, the optimal number of clusters is determined automatically using the elbow method.
#' @param h_clusters Optional. The height at which to cut the dendrogram for clustering. If both `k` and `h` are provided, `k` takes precedence.
#' @param cellTypesToCluster Optional. Cell type annotations to use for clustering. If `NULL` (default), all observations will be clustered.
#' @param mergeCNV Logical. Whether to merge the highly correlated CNV clusters.
#' @param mergeThreshold A numeric value between 0 and 1. Clusters with correlation greater than this threshold will be merged. Default is 0.98.
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
#' @param raster_resize_mat Whether resize the matrix to let the dimension of the matrix the same as the dimension of the raster image.
#' Default is TRUE.
#'
#' @import Seurat
#' @importFrom crayon red yellow green black
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

                          getCNVClusters = TRUE,
                          k_clusters = NULL,
                          h_clusters = NULL,
                          cellTypesToCluster = NULL,

                          mergeCNV = TRUE,
                          mergeThreshold = 0.98,

                          doPlot = TRUE,
                          denoise = TRUE,
                          printPlot = FALSE,
                          savePath = ".",
                          outputType = "png",
                          clustersVar = "cnv_clusters",
                          clusters_palette = "default",
                          splitPlotOnVar = clustersVar,
                          referencePalette = "default",
                          raster_resize_mat = TRUE){

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
  newHDobj <- CNVAnalysis(object = newHDobj,
                          referenceVar = referenceVar,
                          referenceLabel = referenceLabel,
                          assay = assay,
                          pooledReference = pooledReference,
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

  if (length(newHDobj) == 1) {
    seuratObjHD[["rawGenomicScores"]] = Seurat::GetAssay(newHDobj, assay = "rawGenomicScores")
    seuratObjHD[["genomicScores"]] = Seurat::GetAssay(newHDobj, assay = "genomicScores")
    seuratObjHD$cnv_fraction = NA
    seuratObjHD$cnv_fraction[rownames(newHDobj@meta.data)] = newHDobj$cnv_fraction
    Seurat::DefaultAssay(seuratObjHD) = assay
    seuratObjHD@project.name = sampleName
    rm(newHDobj) ; invisible(gc())
  } else if (length(seuratObjHD) > 1) {
    for (i in 1:length(newHDobj)){
      mid <- newHDobj[[i]]
      seuratObjHD[[i]][["rawGenomicScores"]] = Seurat::GetAssay(mid, assay = "rawGenomicScores")
      seuratObjHD[[i]]$cnv_fraction = NA
      seuratObjHD[[i]]$cnv_fraction[rownames(mid@meta.data)] = mid$cnv_fraction
      Seurat::DefaultAssay(seuratObjHD[[i]]) = assay
      seuratObjHD[[i]]@project.name = sampleName[[i]]
      rm(mid) ; invisible(gc())
    }
    rm(newHDobj) ; invisible(gc())
  }


  if (getCNVPerChromosomeArm == TRUE){
    message(crayon::yellow(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Computing CNV per chromosome arm...")))
    if (length(seuratObjHD) == 1) {
      seuratObjHD <- CNVPerChromosomeArm(seuratObjHD)
      invisible(gc())
    } else if (length(seuratObjHD) > 1) {
      for (i in 1:length(seuratObjHD)){
        seuratObjHD[[i]] <- CNVPerChromosomeArm(seuratObjHD[[i]])
        invisible(gc())
      }
    }
    message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done!")))
  }

  if (getCNVClusters == TRUE){
    message(crayon::yellow(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Running CNV clustering...")))
    if (length(seuratObjHD) == 1) {
      Seurat::DefaultAssay(seuratObjHD) = assay
      seuratObjHD <- CNVCluster(seuratObj = seuratObjHD,
                                k = k_clusters,
                                h = h_clusters,
                                referenceVar = referenceVar,
                                cellTypesToCluster = cellTypesToCluster)
      invisible(gc())
    } else if (length(seuratObjHD) > 1) {
      for (i in 1:length(seuratObjHD)){
        Seurat::DefaultAssay(seuratObjHD[[i]]) = assay
        seuratObjHD[[i]] <- CNVCluster(seuratObj = seuratObjHD[[i]],
                                  k = k_clusters,
                                  h = h_clusters,
                                  referenceVar = referenceVar,
                                  cellTypesToCluster = cellTypesToCluster)
        invisible(gc())
      }
    }

    if (mergeCNV == TRUE) {
      if(length(seuratObjHD) == 1) {
        seuratObjHD <- mergeCNVClusters(seuratObj = seuratObjHD, mergeThreshold = mergeThreshold)
      } else if (length(seuratObjHD) > 1) {
        for (i in 1:length(seuratObjHD)){
          seuratObjHD[[i]] <- mergeCNVClusters(seuratObj = seuratObjHD[[i]], mergeThreshold = mergeThreshold)
        }
      }
    }
    message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done!")))
  }

  if (length(seuratObjHD) == 1) {
    if (doPlot == TRUE) {
      tryCatch(
        {
          plotCNVResultsHD(
            seuratObjHD = seuratObjHD,
            denoise = denoise,
            printPlot = printPlot,
            savePath = savePath,
            outputType = outputType,
            referenceVar = referenceVar,
            clustersVar = clustersVar,
            clusters_palette = clusters_palette,
            splitPlotOnVar = splitPlotOnVar,
            referencePalette = referencePalette,
            raster_resize_mat = raster_resize_mat
          )
          invisible(gc())
        },
        error = function(e) {
          # fallback if error
          plotCNVResultsHD(
            seuratObjHD = seuratObjHD,
            denoise = denoise,
            printPlot = printPlot,
            savePath = savePath,
            outputType = outputType,
            referenceVar = referenceVar,
            clustersVar = clustersVar,
            clusters_palette = clusters_palette,
            splitPlotOnVar = splitPlotOnVar,
            referencePalette = referencePalette,
            raster_resize_mat = FALSE
          )
          invisible(gc())
        }
      )
    }

  } else if (length(seuratObjHD) > 1) {
    if (doPlot == TRUE) {
      for (i in seq_along(seuratObjHD)) {
        tryCatch(
          {
            plotCNVResultsHD(
              seuratObjHD = seuratObjHD[[i]],
              denoise = denoise,
              printPlot = printPlot,
              savePath = savePath,
              outputType = outputType,
              referenceVar = referenceVar,
              clustersVar = clustersVar,
              clusters_palette = clusters_palette,
              splitPlotOnVar = splitPlotOnVar,
              referencePalette = referencePalette,
              raster_resize_mat = raster_resize_mat
            )
            invisible(gc())
          },
          error = function(e) {
            plotCNVResultsHD(
              seuratObjHD = seuratObjHD[[i]],
              denoise = denoise,
              printPlot = printPlot,
              savePath = savePath,
              outputType = outputType,
              referenceVar = referenceVar,
              clustersVar = clustersVar,
              clusters_palette = clusters_palette,
              splitPlotOnVar = splitPlotOnVar,
              referencePalette = referencePalette,
              raster_resize_mat = FALSE
            )
            invisible(gc())
          }
        )
      }
    }
  }

  return(seuratObjHD)
}
