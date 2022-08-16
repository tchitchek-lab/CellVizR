# 1. Introduction

Cytometry data are now classically analyzed using non-linear
dimensionality reduction approaches, but it is still challenging to
easily handle the whole pipeline of computational analyses.

UMAPVizR allows the statistical analysis and visualization of
high-dimensional cytometry data using manifold algorithms and clustering
methods. Especially, several key analysis steps are available to perform
data importation, manifold generation, cell cluster identification,
statistical analyses, cluster visualization, and quality controls of
generated results.

UMAPVizR can import cell events from FCS or txt file formats using
different transformation, down-sampling, and normalization approaches.
Manifold representations can be generated using the UMAP, tSNE or
LargeVis algorithms to project cell events into a lower dimensionality
space. Cell clusters can be identified using multiple clustering
algorithms, depending on the user’s assumptions. The characteristics of
cell clusters can be visualized using scatter plots, categorical heatmap
of marker expressions, or using parallel coordinates representations.
Cell clusters having abundances differently expressed between biological
conditions can be identified using several statistical tests.
Statistical results can be visualized using volcano plots or heatmaps.

## 1.1 Workflow overview

In the `UMAPVizR` workflow, an S4 object is created to store data and
sample information is implemented for analysis. This stored information
will allow performing the statistics and visualization of the dataset.

<img src="README/figures/workflow.png" width="90%" style="display: block; margin: auto;" />

*Figure 1: Workflow of UMAPVizR*

*The analysis in UMAPVizR consists of 5 main steps: (1) importing the
data in FCS or txt format resulting in the creation of an S4 UMAPdata
object; (2) assigning the metadata (sample information) into the
UMAPdata object; and (3) generating the manifold and clustering. The
computed results can be (4) visualized in different manners and (5)
analyzed using statistical approaches.*

## 1.2 Input data

The following conditions must be respected to analyze data with
`UMAPVizR`:

-   **Type and format of data**: The cytometry data that can be analyzed
    and integrated with `UMAPVizR` are flow, mass or spectral cytometry
    data. The input files can be in standard cytometry format (FCS) or
    txt format.
-   **Compensation**: Before starting an analysis with `UMAPVizR`,
    performing the compensation steps for flow cytometry and spectral
    data with conventional software (FlowJo, Kaluza, etc) is necessary.
-   **Cleaning and gating**: It is recommended to remove debris, dead
    cells and doublets before the analysis. A pre-gating on a cell
    population of interest (e.g.lymphocytes, B cells, etc.) can be
    performed.

# 2. Quick start

In this section, the main analysis steps of `UMAPVizR` are presented.

These steps cover several aspects, such as:

-   Installing the package
-   Importing the data and creating an `UMAPdata` object
-   Creating the manifold and clustering
-   Generating basic visualization

## 2.1 Installation

To download `UMAPVizR` it is required `devtools`:

``` r
install.packages("devtools")
library("devtools")
install_github("tchitchek-lab/UMAPVizR")
```

The `UMAPVizR` package automatically downloads the necessary packages
for its operation such as: `coin`, `concaveman`, `dendextend`,
`flowCore`, `ggdendro`, `gglot2`, `gridExtra`, `MASS`, `plyr`,
`reshape`, `reshape2`, `rstatix`, `Rtsne`, `scales`, `stats`, `stringr`,
`uwot`. If not, the packages are available on the `CRAN`, except
`flowCore` which is available on `Bioconductor`.

Once installed, `UMAPVizR` can be loaded using the following command:

``` r
library("UMAPVizR")
```

## 2.2 Importing cell expression profiles (import)

The `import` function allows importing the expression matrix of the
cytometry files into a `UMAPdata` object.

The files to be loaded must be in FCS or txt format. The `import`
function is used as below:

``` r
# creation of a vector containing the names of the files 
files <- list.files("C:/Users/GWMA/Documents/Transreg/03_Kaluza_exports_renamed/Panel_03_NK/", 
                    pattern = "fcs", full.names = TRUE)

# import the FCS files into a UMAPdata object 
UMAPV <- import(files, 
                filetype = "fcs", 
                transform = "logicle", 
                exclude.markers = c("FS","FS.1","FS.2", "SS","SS.1","SS.2", "Time"), 
                downsampling = 5000)
```

The main arguments of the `import` function are:

-   the `filetype` argument, which allows defining the data file type
    (`fcs` or `txt`)
-   the `transform` argument, which allows choosing the type of
    transformation to apply to the data. Possible values are: `none`,
    `logicle`, `arcsinh` and `logarithmic`. It is advised to use a
    `logicle` transform for flow cytometry, and to use an `arcsinh`
    transform for mass cytometry,
-   the `exclude_markers` argument, which is used to remove the
    irrelevant channels

After importing the dataset, the `plotCellCounts` function allows you to
see the number of cells in each sample to be displayed as follows:

``` r
plot(plotCellCounts(UMAPV, 
                    stats = c("min","median","mean","q75","max"),
                    samples = NULL,
                    sort = TRUE))
```

![](README/figure-markdown_github/plotCellCounts-1.png)

## 2.3 Assigning meta-information of biological samples (assignMetadata)

The metadata (information about the biological samples) can be assigned
to each sample in the dataset. These metadata are then used by the
different visualization methods to properly represent biological
conditions, timepoints, and individuals. The metadata argument must be a
dataframe that contains exclusively the following column names:

-   individual: corresponds to the sample identifier,
-   condition: corresponds to the biological condition of the sample,
-   timepoint: corresponds to the timepoint of the sample (optional).

Here is an example of a metadata assignment:

``` r
# creation of the dataframe 
metadata <- data.frame("individual"= c("10105LA","10209HE","10306CG","10307BR","10503DC",
                                       "10707BL","11204CD","20208AA","20210RF",
                                       "10105LA","10209HE","10304KJ","10306CG","10309BR",
                                       "10503DC","11204CD","20208AA","20210RF",
                                       "21203AS",
                                       "10105LA","10207BL","10209HE","10304KJ","10306CG",
                                       "10503DC","10807BR","10904VP","11204CD","20208AA",
                                       "20210RF",
                                       "10105LA","10207BL","10209HE","10304KJ","10306CG",
                                       "10503DC","10807BR","10904VP","11204CD","20208AA",
                                       "20210RF","21203AS"),
                       "condition"= rep("HV", 42),
                       "timepoint"= c(rep("V1", 9), rep("V6", 10), rep("V7", 11), rep("V8", 12))
)

rownames(metadata) = paste0(metadata$timepoint,"_", metadata$individual)

# assign the dataframe 
UMAPV <- assignMetadata(UMAPV, 
                        metadata = metadata)
```

## 2.4 Manifold construction and clustering

This section consists in generating the manifold using different
algorithms combined with cell cluster identification.

Two methods are available, depending on the parameters selected:

-   The manifold is generated first, followed by cell cluster
    identification
-   Cell cluster identification is performed followed by the manifold

In the example below, the first method has been performed.

### 2.4.1 Generating a manifold of cell events (generateManifold)

The first step is to compute the manifold on the dataset by following
the instructions below:

``` r
# Perform Manifold from the "UMAPdata" object
UMAPV <- generateManifold(UMAPV, 
                          markers = c("TCRgd", "NKP44", "HLADR", "NKp30", "NKp46",
                                      "NKG2D", "CD3", "CD16", "CD56", "CD8"), 
                          type = "UMAP", 
                          n_neighbors = 15,
                          n_components = 2,
                          metric = "euclidean",
                          n_epochs = NULL,
                          n_threads = 40, 
                          n_sgd_threads = 1,
                          scale = FALSE)
```

    ## Manifold markers are: TCRgd, NKP44, HLADR, NKp30, NKp46, NKG2D, CD3, CD16, CD56, CD8

    ## 14:23:28 UMAP embedding parameters a = 1.896 b = 0.8006

    ## 14:23:28 Read 193322 rows and found 10 numeric columns

    ## 14:23:28 Using Annoy for neighbor search, n_neighbors = 15

    ## 14:23:28 Building Annoy index with metric = euclidean, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 14:23:43 Writing NN index file to temp file C:\Users\GWMA\AppData\Local\Temp\Rtmp6ZjvPa\file2a606008fe9
    ## 14:23:44 Searching Annoy index using 40 threads, search_k = 1500
    ## 14:23:53 Annoy recall = 100%
    ## 14:23:53 Commencing smooth kNN distance calibration using 40 threads
    ## 14:23:58 Initializing from normalized Laplacian + noise
    ## 14:24:03 Commencing optimization for 200 epochs, with 3859266 positive edges using 1 thread
    ## 14:26:22 Optimization finished

The main arguments of the `generateManifold` function are:

-   the `markers` argument, which specifies the markers to be used for
    the manifold generation
-   the `type` argument, which specifies the manifold method to use

### 2.4.2 Identifying cell clusters having similar marker expression (identifyClusters)

The second step is to identify cell clusters by following the
instructions below:

``` r
# Clustering computation from the manifold 
UMAPV <- identifyClusters(UMAPV, 
                          space = "manifold", 
                          method = "kmeans", 
                          centers = 120, 
                          nstart = 3)
```

    ## Identifying cell clusters...

    ## computing cell clusters boundaries...

    ## computing cell cluster count matrix...

    ## computing cell cluster abundance matrix...

The main arguments of the `identifyClusters` function are:

-   the `space` argument, which determines if the clustering is done on
    the markers or the manifold coordinates
-   the `method` argument, which specifies the clustering algorithm to
    use

After clustering, the `plotClustersCounts` function allows to visualize
the cells of each sample in the clusters as follows:

``` r
plot(plotClustersCounts(UMAPV, 
                        clusters = NULL,
                        sort = TRUE))
```

![](README/figure-markdown_github/plotClustersCounts-1.png)

## 2.5 Basic Visualization

Once the manifold has been generated and cell clusters have been
identified, it is possible to perform different types of visualization
which are detailed below.

### 2.5.1 Representation of a computed manifold (PlotManifold)

The `plotManifold` function displays a computed manifold representation
for a given analysis. Cell clusters are delimited by black lines on the
manifold.

The main argument of the `plotManifold` function is the `markers`
argument which is used to specify the colour of the cells. If the
`density` value is used, then a UMAP representation showing the
distribution of the cell density for all samples will be shown as below:

``` r
# Display manifold overlay by 'density' 
plotManifold(UMAPV, 
             markers = "density")
```

![](README/figure-markdown_github/PlotManifold-1.png)

If the name of the marker is used, then the intensity of marker
expression, overlaid on the manifold (e.g. CD8), will be shown as below:

``` r
# Display manifold overlay by 'markers'  
plotManifold(UMAPV, 
             markers = "CD8")
```

![](README/figure-markdown_github/PlotManifold2-1.png)

It is possible to specify the biological samples to be displayed in the
representation using the `samples` argument as below:

``` r
# Display manifold overlay by 'density' by sample 
plotManifold(UMAPV, 
             markers = "density", 
             samples = "V1_10105LA")
```

![](README/figure-markdown_github/PlotManifold3-1.png)

### 2.5.2 Heatmap of cell marker expressions (plotHmExpressions)

The `plotHmExpressions` function shows marker median relative
expressions for all clusters in the whole dataset.

The mean of the median expression of each marker is classified into 4
categories (the number of categories can be changed by users, `nb.cat`
parameters). Hierarchical clustering is performed at both the marker and
cluster levels and is represented using dendrograms (the hierarchical
clustering parameters can be changed by users `method.hclust`
parameters).

This function is used as below:

``` r
# Heatmap of expression markers 
hm.exp <- plotHmExpressions(UMAPV)
```

![](README/figure-markdown_github/PlotHMExpressions-1.png) It is
possible to customize the `plotHmExpressions` with these parameters:

-   the `markers` argument, which specifies the markers to be displayed
-   the `clusters` argument, which specifies the identifiers of the
    clusters to be displayed

These parameters can be used independently of each other as in the
following example:

``` r
# Heatmap of expression markers 
hm.exp <- plotHmExpressions(UMAPV, 
                            markers = c("NKP44", "NKp30", "NKp46", "NKG2D"), 
                            clusters = c(1:50))
```

![](README/figure-markdown_github/plotHmExpressions2-1.png)

# 3. Statistics and visualization

## 3.1 Compute differential abundance analyses

Once the cell clustering performed, it is possible to do a differential
analysis of cell cluster abundances to identify relevant cell clusters.

The `computeStatistics` function allows to perform the such operation
and several parameters must be taken into consideration:

-   the `condition` argument, which specifies the biological condition
    to be compared
-   the `ref.condition` argument, which specifies the reference
    biological condition
-   the `test.statistics` argument, which specifies the name of the
    statistical test to use
-   the `paired` argument, which specifies if samples are paired in the
    statistical comparison

This function is used as follows:

``` r
# Compute statistics 
baseline = "V1"
list.conditions <- c("V6", "V7", "V8")

for (condition in list.conditions) {
  UMAPV <- computeStatistics(UMAPV, 
                             condition = paste0(condition), 
                             ref.condition = paste0(baseline),
                             test.statistics = "wilcox.test",
                             paired = FALSE)
}
```

    ## Computing of the wilcox.test for: V6 vs. V1

    ## Computing of the wilcox.test for: V7 vs. V1

    ## Computing of the wilcox.test for: V8 vs. V1

## 3.2 Visualisation of statistical analysis

### 3.2.1 Volcano plot of statistical analysis (plotVolcanoPlot)

The `plotVolcanoPlot` function shows the clusters whose number of
associated cells is statistically different between two biological
conditions and/or timepoints.

For each cluster, the p-value (indicated by -log10(p-value)) is
represented on the Y-axis and the cell abundance fold-change (indicated
by log2(fold-change)) is represented on the X-axis. The thresholds for
the p-value (`th.pv` parameter) and the fold-change (`th.fc` parameter)
are shown as dotted lines. Cell clusters down-represented are shown in
green and cell clusters up-represented are shown in red.

Here is an example for generating such representation:

``` r
# Volcano plot for differential analysis 
plotVolcanoPlot(UMAPV,
                comparison = ("V7 vs. V1"),
                th.pv = 1.3,
                th.fc = 1.5)
```

![](README/figure-markdown_github/plotVolcanoPlot-1.png)

### 3.2.2 Heatmap of statistical analysis results (plotHmStatistics)

The `plotHmStatistics` function shows the differences in abundance
between different conditions for each cluster.

For each cluster, the p-value, the log2(fold-change) and the effect size
(`statistics` parameters) can be represented. Down-represented clusters
are represented in orange, and up-represented clusters are represented
in blue. Furthermore, it is possible to choose the clusters to be
represented with the `clusters` parameter.

Here is an example for generating such representation:

``` r
# Heatmap of statistics
hm.stats <- plotHmStatistics(UMAPV, 
                             clusters = NULL,
                             statistics = "pvalue")
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](README/figure-markdown_github/plotHmStatistics-1.png)

## 3.3 Visualisation of cell cluster abundances

### 3.3.1 Heatmap of cell cluster abundances (plotHmAbundances)

The `plotHmAbundances` function shows the cellular distribution of
samples within a given cluster.

The more the sample is represented within the cluster, the redder the
tile. If the sample is not represented in the cluster, then the tile
will be black. The `plotHmAbundances` function can be interesting to
visualize the abundance of statistically different clusters between two
conditions, as in the following example:

``` r
#Samples to study
V1 <- unique(UMAPV@samples)[grepl("V1", unique(UMAPV@samples))]
V6 <- unique(UMAPV@samples)[grepl("V6", unique(UMAPV@samples))]
samples = c(V1, V6)

#Statistically different clusters
stats <- UMAPV@statistic[UMAPV@statistic$comparison == "V6 vs. V1",]
clusters = stats[stats$pvalue<=0.01 & abs(stats$lfc)>log(1.5)/log(2),]$clusters

# Heatmap of abundances
plotHmAbundances(UMAPV, 
                 clusters = clusters,
                 samples = samples,
                 rescale = TRUE)
```

![](README/figure-markdown_github/plotHmAbundances-1.png)

### 3.3.2 Cell cluster abundances using a boxplot representation (plotBoxplot)

The `plotBoxPlot` function shows the cell distribution between several
biological conditions and/or timepoints for a single cluster or for a
combined set of clusters.

This display shows the abundances of the user-defined cell clusters
(`clusters` parameter). It is possible to observe the cell abundance as
a function of the biological condition or timepoint (`obervation`
parameter). In addition, statistical tests can be performed and
displayed directly on the boxplot.

Here is an example for generating such representation:

``` r
# Boxplot for differential analysis
plotBoxplot(UMAPV, 
            clusters = "31",
            samples = NULL,
            observation = "timepoint", 
            test.statistics = "wilcox.test")
```

![](README/figure-markdown_github/plotBoxplot-1.png)

Other possible parameters to customize the `plotBoxPlot` are:

-   the `samples` argument, which specifies the biological samples to be
    displayed
-   the `paired` argument, which specifies if samples are paired in the
    statistical comparison

### 3.3.3 MDS representation based on cell cluster abundances (plotMDS)

The `plotMDS` function shows similarities between samples or clusters
based on cell cluster abundances.

Each point represents a sample or a cluster (`levels` parameter) and the
distance between the points is proportional to the Euclidean distance
between these objects. It is possible to observe the cell abundance as a
function of the biological condition or timepoint (`condition.samples`
parameter)

Here is an example for generating such representation:

``` r
# MDS
plotMDS(UMAPV, 
        levels = "samples", 
        condition.samples = "timepoint", 
        clusters = NULL, 
        samples = NULL)
```

    ## Warning: ggrepel: 29 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](README/figure-markdown_github/plotMDS-1.png)

Other possible parameters to customize the `plotMDS` are:

-   the `clusters` argument, which specifies the identifiers of the
    clusters to be displayed
-   the `samples` argument, which specifies the biological samples to be
    displayed

### 3.3.4 PCA representation based on cell cluster abundances (plotPCA)

The `plotPCA` function shows similarities between samples or clusters
based on cell cluster abundances.

Each point represents a sample or a cluster (`levels` parameter). It is
possible to observe the cell abundance as a function of the biological
condition or timepoint (`condition.samples` parameter)

``` r
# PCA
plotPCA(UMAPV, 
        levels = "clusters", 
        clusters = NULL, 
        samples = NULL, 
        condition.samples = "condition")
```

![](README/figure-markdown_github/plotPCA-1.png)

Other possible parameters to customize the `plotPCA` are:

-   the `clusters` argument, which specifies the identifiers of the
    clusters to be displayed
-   the `samples` argument, which specifies the biological samples to be
    displayed

# 4. Quality control

The `UMAPVizR` package allows to perform quality control of generated
results.

The quality control can be performed:

-   on the input dataset to check the names and range expression of the
    markers of each sample
-   on the generated results, to check the quality of the cell
    clustering.

## 4.1 Quality control of the dataset

The input dataset can be checked in two ways. The first method checks
the concordance of the markers names between the different samples.

Here is an example of generating such quality control:

``` r
# Check for marker concordance
QCN <- QCMarkerNames(files)
```

    ##            nb_cells FS-H FS-A FS-W SS-H SS-A SS-W       FL1-A    FL2-A  FL3-A
    ## V1_10105LA     5768   FS   FS   FS   SS   SS   SS TCR gd-FITC NKP44-PE DR-ECD
    ## V1_10209HE     4944   FS   FS   FS   SS   SS   SS TCR gd-FITC NKP44-PE DR-ECD
    ## V1_10306CG     4746   FS   FS   FS   SS   SS   SS TCR gd-FITC NKP44-PE DR-ECD
    ## V1_10307BR     3615   FS   FS   FS   SS   SS   SS TCR gd-FITC NKP44-PE DR-ECD
    ## V1_10503DC     5877   FS   FS   FS   SS   SS   SS TCR gd-FITC NKP44-PE DR-ECD
    ## V1_10707BL     7823   FS   FS   FS   SS   SS   SS TCR gd-FITC NKP44-PE DR-ECD
    ##                 FL4-A      FL5-A     FL6-A    FL7-A     FL8-A      FL9-A FL10-A
    ## V1_10105LA NKp30-Pcy5 NKp46-Pcy7 NKG2D-APC CD3-A700 CD16-A750 CD56-BV421 CD8-KO
    ## V1_10209HE NKp30-Pcy5 NKp46-Pcy7 NKG2D-APC CD3-A700 CD16-A750 CD56-BV421 CD8-KO
    ## V1_10306CG NKp30-Pcy5 NKp46-Pcy7 NKG2D-APC CD3-A700 CD16-A750 CD56-BV421 CD8-KO
    ## V1_10307BR NKp30-Pcy5 NKp46-Pcy7 NKG2D-APC CD3-A700 CD16-A750 CD56-BV421 CD8-KO
    ## V1_10503DC NKp30-Pcy5 NKp46-Pcy7 NKG2D-APC CD3-A700 CD16-A750 CD56-BV421 CD8-KO
    ## V1_10707BL NKp30-Pcy5 NKp46-Pcy7 NKG2D-APC CD3-A700 CD16-A750 CD56-BV421 CD8-KO

If the marker names are not the same for each sample, they can be
corrected using the `renameMarkers` as below:

``` r
# Rename markers if necessary
UMAPV <- renameMarkers(UMAPV, marker.names = c("TCRgd", "NKP44", "HLADR", "NKp30", "NKp46",
                                               "NKG2D", "CD3", "CD16", "CD56", "CD8"))
```

The second method computes the 5 centiles and 95 centiles expression
values for each marker of each sample:

``` r
# Check the expression values for markers
QCR <- QCMarkerRanges(files)
```

    ##                  FS       FS       FS       SS       SS       SS TCR gd-FITC
    ## V1_10105LA 4.836585 4.989121 4.353023 4.246526 4.410427 4.347033   1.2229851
    ## V1_10209HE 4.830411 4.955269 4.347033 4.233577 4.376030 4.334798   1.4042910
    ## V1_10306CG 4.947511 5.075207 4.347033 4.306356 4.453486 4.340958   1.6249665
    ## V1_10307BR 4.884884 5.021064 4.347033 4.316584 4.488392 4.353023   1.6031246
    ## V1_10503DC 4.829002 4.962259 4.353023 4.289901 4.446210 4.347033   1.6368241
    ## V1_10707BL 4.815139 4.954592 4.364761 4.337542 4.477769 4.347033   0.9124178
    ##             NKP44-PE   DR-ECD NKp30-Pcy5 NKp46-Pcy7 NKG2D-APC  CD3-A700
    ## V1_10105LA 1.2894598 1.271010  1.1037258   1.989170  2.056610 1.9222306
    ## V1_10209HE 0.9668970 1.579406  1.1496302   1.846668  2.414580 0.9294985
    ## V1_10306CG 0.8354852 1.698004  0.8368412   2.761478  3.339898 1.2062602
    ## V1_10307BR 1.3685657 1.532980  1.1084658   2.161430  2.873332 0.6686698
    ## V1_10503DC 1.0512715 1.318435  0.8781855   2.253348  2.616961 0.8506702
    ## V1_10707BL 0.6770187 1.399496  1.4884913   2.546774  3.139840 2.0235424
    ##            CD16-A750 CD56-BV421   CD8-KO
    ## V1_10105LA  2.265913   2.708021 1.823717
    ## V1_10209HE  2.694609   2.973497 1.934102
    ## V1_10306CG  2.824189   2.967019 1.673800
    ## V1_10307BR  2.326850   3.037347 1.498806
    ## V1_10503DC  2.686252   2.916044 2.133573
    ## V1_10707BL  2.834552   3.253349 2.259912

    ##                  FS       FS       FS       SS       SS       SS TCR gd-FITC
    ## V1_10105LA 5.158144 5.305961 4.443943 4.621098 4.813464 4.506733    2.280538
    ## V1_10209HE 5.173436 5.307780 4.403514 4.584182 4.724468 4.413979    2.118809
    ## V1_10306CG 5.248328 5.374601 4.392790 4.614367 4.754747 4.413979    2.195158
    ## V1_10307BR 5.180815 5.317134 4.408778 4.653551 4.830843 4.448741    2.347433
    ## V1_10503DC 5.137767 5.276119 4.413979 4.626658 4.789096 4.434183    2.290403
    ## V1_10707BL 5.141039 5.285225 4.419119 4.696846 4.834221 4.413979    2.117568
    ##            NKP44-PE   DR-ECD NKp30-Pcy5 NKp46-Pcy7 NKG2D-APC CD3-A700 CD16-A750
    ## V1_10105LA 2.698542 3.751472   3.070920   4.164283  3.928043 3.178521  4.747111
    ## V1_10209HE 2.408619 3.910715   3.057855   3.956062  3.914080 2.842833  4.826450
    ## V1_10306CG 2.412678 3.838884   3.186744   4.040787  4.063245 3.043769  4.817522
    ## V1_10307BR 2.504910 4.329139   3.344018   3.950651  3.987541 2.864951  4.703335
    ## V1_10503DC 2.414765 3.107674   3.167953   3.917157  3.895880 2.806967  4.751921
    ## V1_10707BL 2.377488 3.506573   3.267490   3.982711  3.964351 3.146767  4.549099
    ##            CD56-BV421   CD8-KO
    ## V1_10105LA   4.359073 3.387718
    ## V1_10209HE   4.428926 3.528046
    ## V1_10306CG   4.336497 3.744981
    ## V1_10307BR   4.392716 3.813623
    ## V1_10503DC   4.264874 3.901291
    ## V1_10707BL   4.446657 3.906904

## 4.2 Control quality of the cell clustering result

The quality control of clustering can be checked in two ways.

The first method allows the identification of small clusters,
i.e.clusters whose number of cells is below a specific threshold. The
results can be represented as a heatmap. On the left are the
contributions of each sample and on the right are the contribution of
the whole dataset. If the tile is red then the cluster is less than the
specified number of cells, if the tile is green, the cluster is greater
than or equal to the specified number of cells. The percentage of
clusters with a small number of cells among all clusters is shown at the
top of the heatmap.

The function is as below:

``` r
# QC for small clusters 
QCS <- QCSmallClusters(UMAPV,
                       th.size = 50, 
                       plot.device = TRUE)
```

![](README/figure-markdown_github/QCSmallClusters-1.png)

    ##      V1_10105LA V1_10209HE V1_10306CG V1_10307BR V1_10503DC V1_10707BL
    ## [1,]      FALSE       TRUE      FALSE       TRUE       TRUE      FALSE
    ## [2,]       TRUE       TRUE       TRUE       TRUE      FALSE       TRUE
    ## [3,]       TRUE      FALSE      FALSE       TRUE       TRUE       TRUE
    ## [4,]       TRUE       TRUE       TRUE      FALSE       TRUE       TRUE
    ## [5,]       TRUE      FALSE       TRUE      FALSE      FALSE       TRUE
    ## [6,]      FALSE      FALSE       TRUE       TRUE      FALSE       TRUE
    ##      V1_11204CD V1_20208AA V1_20210RF V6_10105LA V6_10209HE V6_10304KJ
    ## [1,]       TRUE       TRUE       TRUE      FALSE      FALSE       TRUE
    ## [2,]      FALSE      FALSE      FALSE       TRUE      FALSE       TRUE
    ## [3,]       TRUE       TRUE      FALSE      FALSE      FALSE       TRUE
    ## [4,]       TRUE       TRUE       TRUE       TRUE      FALSE       TRUE
    ## [5,]       TRUE      FALSE      FALSE       TRUE      FALSE       TRUE
    ## [6,]      FALSE      FALSE      FALSE      FALSE      FALSE       TRUE
    ##      V6_10306CG V6_10309BR V6_10503DC V6_11204CD V6_20208AA V6_20210RF
    ## [1,]      FALSE       TRUE       TRUE       TRUE      FALSE      FALSE
    ## [2,]      FALSE       TRUE      FALSE      FALSE      FALSE      FALSE
    ## [3,]       TRUE       TRUE       TRUE       TRUE       TRUE      FALSE
    ## [4,]      FALSE       TRUE       TRUE      FALSE      FALSE      FALSE
    ## [5,]       TRUE       TRUE       TRUE       TRUE      FALSE      FALSE
    ## [6,]       TRUE       TRUE       TRUE       TRUE       TRUE       TRUE
    ##      V6_21203AS V7_10105LA V7_10207BL V7_10209HE V7_10304KJ V7_10306CG
    ## [1,]      FALSE       TRUE      FALSE       TRUE       TRUE      FALSE
    ## [2,]       TRUE      FALSE       TRUE      FALSE      FALSE      FALSE
    ## [3,]       TRUE      FALSE       TRUE      FALSE       TRUE      FALSE
    ## [4,]       TRUE       TRUE       TRUE      FALSE      FALSE      FALSE
    ## [5,]       TRUE      FALSE       TRUE      FALSE      FALSE       TRUE
    ## [6,]       TRUE      FALSE       TRUE      FALSE      FALSE      FALSE
    ##      V7_10503DC V7_10807BR V7_10904VP V7_11204CD V7_20208AA V7_20210RF
    ## [1,]       TRUE       TRUE       TRUE       TRUE       TRUE      FALSE
    ## [2,]      FALSE       TRUE       TRUE       TRUE      FALSE      FALSE
    ## [3,]      FALSE       TRUE       TRUE       TRUE      FALSE      FALSE
    ## [4,]       TRUE      FALSE       TRUE       TRUE      FALSE      FALSE
    ## [5,]      FALSE       TRUE       TRUE       TRUE      FALSE       TRUE
    ## [6,]      FALSE       TRUE       TRUE       TRUE       TRUE      FALSE
    ##      V8_10105LA V8_10207BL V8_10209HE V8_10304KJ V8_10306CG V8_10503DC
    ## [1,]       TRUE      FALSE       TRUE       TRUE      FALSE       TRUE
    ## [2,]      FALSE       TRUE       TRUE       TRUE       TRUE      FALSE
    ## [3,]      FALSE       TRUE      FALSE       TRUE      FALSE      FALSE
    ## [4,]       TRUE       TRUE       TRUE       TRUE       TRUE       TRUE
    ## [5,]      FALSE       TRUE      FALSE      FALSE       TRUE       TRUE
    ## [6,]      FALSE       TRUE      FALSE      FALSE       TRUE      FALSE
    ##      V8_10807BR V8_10904VP V8_11204CD V8_20208AA V8_20210RF V8_21203AS
    ## [1,]       TRUE       TRUE       TRUE       TRUE      FALSE      FALSE
    ## [2,]       TRUE       TRUE      FALSE      FALSE      FALSE      FALSE
    ## [3,]       TRUE       TRUE       TRUE       TRUE      FALSE      FALSE
    ## [4,]       TRUE       TRUE       TRUE       TRUE       TRUE      FALSE
    ## [5,]      FALSE       TRUE       TRUE      FALSE       TRUE      FALSE
    ## [6,]      FALSE       TRUE       TRUE      FALSE       TRUE      FALSE
    ##      total.cells
    ## [1,]       FALSE
    ## [2,]       FALSE
    ## [3,]       FALSE
    ## [4,]       FALSE
    ## [5,]       FALSE
    ## [6,]       FALSE

The second method allows to identify the uniform clusters, i.e.those
with unimodal expression and low dispersion of expression for all its
markers.

The most important parameter of the `QCUniformClusters` function is
`uniform.test`, three possibilities:

-   `uniform` corresponds to the verification of the unimodal
    distribution of markers with Hartigan’s test (`th.pvalue`
    parameter),
-   `IQR` corresponds to the verification of the distribution of markers
    so that they are not below the IQR threshold (`th.IQR` parameter)
-   `both` correspond to the combination of the two parameters: uniform
    and IQR

The results can be represented as a heatmap. If the tile is green then
the cell clusters have the uniform phenotype, if the tile is red, the
cell clusters have the phenotype that is not uniform. The percentage of
clusters having a uniform phenotype among all clusters is shown at the
top of the heatmap. If the score is high, it indicates that the
clustering is good.

The function is as below:

``` r
# QC for uniform clusters
QCU <- QCUniformClusters(UMAPV,
                         uniform.test = "both",
                         th.pvalue = 0.05,
                         th.IQR = 2,
                         plot.device = TRUE)
```

    ## Using clusters as id variables
    ## Using clusters as id variables

![](README/figure-markdown_github/QCUniformClusters-1.png)

    ##   clusters markers    pv_dip       IQR passed
    ## 1        1    CD16 0.9936791 0.3635300   TRUE
    ## 2        1     CD3 0.9988149 0.3031861   TRUE
    ## 3        1    CD56 0.9927818 0.3802556   TRUE
    ## 4        1     CD8 0.9924886 0.4070512   TRUE
    ## 5        1   HLADR 0.9976176 0.4741935   TRUE
    ## 6        1   NKG2D 0.9655686 0.2297570   TRUE

# 5. Advanced usage

## 5.1 Upsampling

The `performUpsampling` function allows the data set to be implemented
if downsampling has been performed.

This function is used after performing the manifold and clustering (Step
2.4). After calculating the centroids from the existing clusters, the
implemented cells will be associated according to their expression
similarity with the centroid.

The procedure is as follows:

``` r
UMAPV <- performUpsampling(UMAPV,
                           files = files)
```

## 5.2 Export

The `export` function allows extracting of the dataset in FCS or txt
format with some parameters such as UMAP coordinates and clusters.

Please note that if downsampling and upsampling have been performed,
only the downsampled cells will be extracted.

With the following method:

``` r
export(UMAPV,
       filename = "Analyses_NK_K100.fcs",
       clusters = NULL,
       samples = NULL)
```

    ## [1] "Analyses_NK_K100.fcs"
