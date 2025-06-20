
## Introduction <img src = "https://must-bioinfo.github.io/fastCNV/articles/logo.png" align = "right" width = "200">

`fastCNV` is a package that helps you detect, plot and analyse the putative Copy Number Variations (CNVs) in single cell (scRNA-seq) data or Spatial Transcriptomics (ST) data. Built on `SeuratObject`, it is easily integrated into scRNA-seq or ST pipelines. 

*WARNING:* Project is still under construction and function usage may change.

## Installation

To install fastCNV, run in R:

```         
remotes::install_github("must-bioinfo/fastCNV")
```

## Usage

Extensive tutorials to run `fastCNV`on scRNA-seq and spatial transcriptomics data are available to get started [here](https://must-bioinfo.github.io/fastCNV/articles/index.html).


## NEW : RUN fastCNV ON YOUR VISIUM HD DATA

You can now run `fastCNV` on visium HD data using `fastCNV_10XHD()`

![](https://must-bioinfo.github.io/fastCNV/articles/fastCNV_HD_files/figure-html/cnvheatmapt_clusters-1.png)

You can also see the cnv clusters on the slide : 

<img src=https://must-bioinfo.github.io/fastCNV/articles/fastCNV_HD_files/figure-html/plot_cnv_clusters-1.png width = 70%>

And see the clonality tree : 

<img src=https://must-bioinfo.github.io/fastCNV/articles/fastCNV_HD_files/figure-html/tree_cnv-1.png width = 70%>


## Capabilities

`fastCNV` can plot a heatmap of inferred CNVs: 

![](https://must-bioinfo.github.io/fastCNV/articles/fastCNV_sc_files/figure-html/heatmap_sc.png)

It also calculates a `cnv_fraction`, which can be plotted with `Seurat` standard plotting functions:

![](https://must-bioinfo.github.io/fastCNV/articles/fastCNV_sc_files/figure-html/plot_cnv_umap-1.png)

`cnv_fractions` can also be visualized spatially for Spatial Transcriptomics samples:

![](https://must-bioinfo.github.io/fastCNV/articles/fastCNV_ST_files/figure-html/plot_cnv_fraction-1.png )

And `cnv_fractions` can be used to obtain clonal clusters (`cnv_clusters`):

<img src="https://must-bioinfo.github.io/fastCNV/articles/fastCNV_sc_files/figure-html/plot_cnv_clusters-1.png" alt="Clonal clusters for sc data" width="50%">

The subclonality tree can be plotted too : 

<img src="https://must-bioinfo.github.io/fastCNV/articles/fastCNV_sc_files/figure-html/tree_cnv-1.png" alt="Subclonality tree" width="500">

