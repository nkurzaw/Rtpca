# Rtpca: an R package for (differential) thermal proximity coaggregation analysis

[![Build Status](https://travis-ci.org/nkurzaw/Rtpca.svg?branch=master)](https://travis-ci.org/nkurzaw/Rtpca) [![codecov](https://codecov.io/gh/nkurzaw/Rtpca/branch/master/graph/badge.svg)](https://codecov.io/gh/nkurzaw/Rtpca)

# Installation

Install the development version of the package from Github.
```{r, eval = FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(“nkurzaw/Rtpca”)
```

2. Load the package into your R session.
```{r Load, message=FALSE}
library(Rtpca)
```

# Introduction

The `Rtpca` package contains functions to analyze thermal proteome profiling (TPP) experiments in the light of protein-protein interactions (PPIs). It makes use of a method termed thermal proximity coaggregation analysis (TPCA) and also allows to test for changes in PPIs across different conditions.

This vignette only represents a minimal example. To have a look at a more realistic example feel free to check out [this more realistic example](https://github.com/nkurzaw/Rtpca_analysis/blob/master/Hashimoto_et_al_analysis.pdf).

We also load the `TPP` package to illustrate how to import TPP data with the Bioconductor package and then input in into the `Rtpca` functions.
```{r, message=FALSE, warning=FALSE}
library(TPP)
```

## Import Thermal proteome profiling data using the TPP package
```{r}
data("hdacTR_smallExample")
```

Filter hdacTR_data to speed up computations
```{r}
set.seed(123)
random_proteins <- sample(hdacTR_data[[1]]$gene_name, 300)
```

```{r}
hdacTR_data_fil <- lapply(hdacTR_data, function(temp_df){
    filter(temp_df, gene_name %in% random_proteins)
})
```

We can now import our small example dataset using the import function from the `TPP` package:
```{r}
trData <- tpptrImport(configTable = hdacTR_config, data = hdacTR_data_fil)
```

```{r}
data("ori_et_al_complex_ppis")
```

## Run TPCA on data from a single condition

We can run TPCA for protein-protein interactions like this:
```{r}
vehTPCA <- runTPCA(
    objList = trData,
    ppiAnno = ori_et_al_complex_ppis
)
```
Note: it is not necessary that your data has the format of the TPP package (ExpressionSet), you can also supply the function with a list of matrices of data frames (in the case of data frames you need to additionally indicate with column contains the protein or gene names).

We can also run TPCA to test for coaggregation of protein complexes like this:
```{r}
data("ori_et_al_complexes_df")

vehComplexTPCA <- runTPCA(
    objList = trData,
    complexAnno = ori_et_al_complexes_df,
    minCount = 2
)
```

We can plot a ROC curve for how well our data captures protein-protein interactions:
```{r}
plotPPiRoc(vehTPCA, computeAUC = TRUE)
```

And we can also plot a ROC curve for how well our data captures protein complexes:

```{r}
plotComplexRoc(vehComplexTPCA, computeAUC = TRUE)

```


## Run differential TPCA on two conditions 

In order to test for protein-protein interactions that change significantly between both conditions, we can run the `runDiffTPCA` as illustrated below:

```{r}
diffTPCA <- 
    runDiffTPCA(
        objList = trData[1:2], 
        contrastList = trData[3:4],
        ppiAnno = ori_et_al_complex_ppis)
```

We can then plot a volcano plot to visualize the results:
```{r}
plotDiffTpcaVolcano(diffTPCA)
```

As mentioned above, this vignette includes only a very minimal example, have a look at a more extensive example [here](https://github.com/nkurzaw/Rtpca_analysis/blob/master/Hashimoto_et_al_analysis.pdf).
