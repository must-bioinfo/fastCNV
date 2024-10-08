#' fastCNV calls all of the internal functions needed to compute the putative CNV on a seurat object or a list of seurat objects
#'
#' @param seuratObj Seurat object or list of seurat objects
#' @param sampleName Name of the sample or list of names corresponding to the samples in the `seuratObj`
#' @param referenceVar The variable name of the annotations in the Seurat metadata
#' @param referenceLabel The label given to the observations wanted as reference (can be any type of annotation)
#' @param splitPlotOnVar The name given to the meta.data slot to split the observations on during the `plotCNVResults` step, if different from `referenceVar`.
#' @param downsizePlot Subset the observations to speed up the plotting process (default = `FALSE`).
#' @param doPlot If `TRUE` will plot a heatmap for each of the samples (default = `TRUE`)
#' @param doRecapPlot Default `TRUE`. Will output the CNV heatmaps by annotation if `TRUE`.
#' @param prepareCounts If `FALSE` will not run the `prepareCountsForCNVAnalysis` function. Default `TRUE`.
#' @param pooledReference Default `TRUE`. Will build a pooled reference across all samples if `TRUE`.
#' @param assay Name of the assay to run the CNV on. Takes the results of `prepareCountsForCNVAnalysis` by default if available
#' @param getCNVFractionPerChromosome If `TRUE`, will save into the metadata the CNV fraction per chromosome arm
#' @param savePath Default `.`. If `NULL` the heatmap won't be saved as a pdf.
#' @param aggregFactor The number of counts per spot desired (default `30 000`). If less than `1000`, will not run the `prepareCountsForCNVAnalysis` function.
#' @param seuratClusterResolution The resolution wanted for the seurat clusters (default 0.8)
#' @param aggregateByVar If `referenceVar` is given, whether to use it to pool the observations
#' @param reClusterSeurat Whether to re-cluster if the Seurat object given already has a `seurat_clusters` slot in its meta.data
#' @param scaleOnReferenceLabel If you want to scale the results depending on the normal observations
#' @param thresholdPercentile Which quantiles to take (if 0.01 it will take 0.01-0.99). Background noise appears with higher numbers.
#' @param genes List of genes (default genes from ensembl version 111)
#' @param windowSize Size of the genomic windows
#' @param windowStep Step between the genomic windows
#' @param topNGenes Number of top expressed genes to keep
#'
#' @return This function returns a list of the seurat objects after all the analysis and it creates the heatmaps of the CNVs of every object in seuratObj and saves them in a pdf in the current working directory
#'
#' @export
#'


fastCNV <- function (seuratObj, sampleName, referenceVar = NULL, referenceLabel = NULL,
                     splitPlotOnVar = referenceVar, downsizePlot = FALSE, doPlot = TRUE,
                     doRecapPlot = TRUE, prepareCounts = TRUE, pooledReference = TRUE, assay = NULL,
                     getCNVFractionPerChromosome = TRUE, savePath = ".",
                     aggregFactor=15000, seuratClusterResolution = 0.8, aggregateByVar = T,
                     reClusterSeurat = F, scaleOnReferenceLabel = TRUE, thresholdPercentile = 0.01,
                     genes=getGenes(), windowSize=100, windowStep=20, topNGenes=7000){

  if(!length(seuratObj)==length(sampleName)) stop("error-fastCNV : seuratObj & sampleName should have the same length")

  if(length(seuratObj) == 1){
    seuratObj <- list(seuratObj) ; names(seuratObj) <- sampleName
  }
  if (prepareCounts == TRUE & aggregFactor>=1000) {
    print("Aggregating counts matrix.")
    for (i in 1:length(seuratObj)) {
      seuratObj[[i]] <- prepareCountsForCNVAnalysis(seuratObj[[i]], sampleName = sampleName[[i]],
                                                   referenceVar = referenceVar,
                                                   aggregateByVar = aggregateByVar ,
                                                   aggregFactor=aggregFactor,
                                                   seuratClusterResolution = seuratClusterResolution,
                                                   reClusterSeurat = reClusterSeurat  )
    }
  } else {
    for (i in 1:length(seuratObj)) {
      seuratObj[[i]]@project.name = sampleName[[i]]
    }
  }


  if (length(seuratObj) > 1) {
    print("Running CNVAnalysis")
    seuratObj <- CNVanalysis(seuratObj, referenceVar = referenceVar, referenceLabel = referenceLabel,
                              doRecapPlot = doRecapPlot, pooledReference = pooledReference,
                              scaleOnReferenceLabel = scaleOnReferenceLabel, assay = assay,
                              thresholdPercentile = thresholdPercentile, genes = genes,
                              windowSize = windowSize, windowStep = windowStep, topNGenes = topNGenes)
    print("CNVAnalysis done!")
    if (doPlot == TRUE) {
      print("Plotting CNV results. This step may take some time.")
      for (i in 1:length(seuratObj)) {
        if (Seurat::Project(seuratObj[[i]]) == "SeuratProject") {Seurat::Project(seuratObj[[i]]) = paste0("Sample",i)}
          plotCNVResults(seuratObj[[i]], referenceVar = referenceVar, splitPlotOnVar = splitPlotOnVar, savePath = savePath, downsizePlot = downsizePlot)
        }
    }

  } else {
    print("Running CNVAnalysis")
    seuratObj <- CNVanalysis(seuratObj[[1]], referenceVar = referenceVar, referenceLabel = referenceLabel,
                         doRecapPlot = doRecapPlot, pooledReference = pooledReference,
                         scaleOnReferenceLabel = scaleOnReferenceLabel, assay = assay,
                         thresholdPercentile = thresholdPercentile, genes = genes,
                         windowSize = windowSize, windowStep = windowStep, topNGenes = topNGenes)
    print("CNVAnalysis done!")
    if (doPlot == TRUE) {
      print("Plotting CNV results. This step may take some time.")
      plotCNVResults(seuratObj = seuratObj, referenceVar = referenceVar, splitPlotOnVar = splitPlotOnVar, savePath = savePath, downsizePlot = downsizePlot)
    }
  }

  if (getCNVFractionPerChromosome == TRUE) {
    if (length(seuratObj) == 1) {
          seuratObj <- CNVfractionPerChromosome(seuratObj)
    } else {
      for (i in 1:length(seuratObj)) {
        seuratObj[[i]] <- CNVfractionPerChromosome(seuratObj[[i]])
      }
    }
  }

  return (seuratObj)
}
