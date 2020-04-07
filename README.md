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

Thermal proteome profiling (TPP) (Savitski et al.,
2014; Mateus et al., 2020) is a mass spectrometry-based, proteome-wide implemention of the cellular thermal shift assay (Martinez Molina et al., 2013). It was originally developed to study drug-(off-)target engagement. However, it
was realized that profiles of interacting protein pairs appeared more
similar than by chance which was coined as 'thermal proximity
co-aggregation' (TPCA) (Tan et al., 2018). The R package `Rtpca` enables
analysis of TPP datasets using the TPCA concept for studying
protein-protein interactions and protein complexes and also allows to
test for differential protein-protein interactions across different
conditions.

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

```{r}
sessionInfo()
```

## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Mojave 10.14.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] TPP_3.12.0          magrittr_1.5        Rtpca_0.0.99       
    ## [4] tidyr_1.0.0         dplyr_0.8.3         Biobase_2.44.0     
    ## [7] BiocGenerics_0.30.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.2           lattice_0.20-38      utf8_1.1.4          
    ##  [4] assertthat_0.2.1     zeallot_0.1.0        digest_0.6.22       
    ##  [7] foreach_1.4.7        sme_1.0.2            R6_2.4.0            
    ## [10] plyr_1.8.4           futile.options_1.0.1 backports_1.1.5     
    ## [13] stats4_3.6.1         evaluate_0.14        ggplot2_3.2.1       
    ## [16] pillar_1.4.2         rlang_0.4.1          lazyeval_0.2.2      
    ## [19] VennDiagram_1.6.20   data.table_1.12.6    rmarkdown_1.16      
    ## [22] nls2_0.2             labeling_0.3         splines_3.6.1       
    ## [25] stringr_1.4.0        RCurl_1.95-4.12      munsell_0.5.0       
    ## [28] compiler_3.6.1       xfun_0.10            pkgconfig_2.0.3     
    ## [31] htmltools_0.4.0      tidyselect_0.2.5     tibble_2.1.3        
    ## [34] gridExtra_2.3        codetools_0.2-16     fansi_0.4.0         
    ## [37] crayon_1.3.4         MASS_7.3-51.4        bitops_1.0-6        
    ## [40] grid_3.6.1           gtable_0.3.0         lifecycle_0.1.0     
    ## [43] formatR_1.7          pROC_1.15.3          scales_1.1.0        
    ## [46] zip_2.0.4            cli_1.1.0            stringi_1.4.3       
    ## [49] farver_2.0.3         reshape2_1.4.3       fdrtool_1.2.15      
    ## [52] doParallel_1.0.15    limma_3.40.6         futile.logger_1.4.3 
    ## [55] ellipsis_0.3.0       vctrs_0.2.0          openxlsx_4.1.0.1    
    ## [58] lambda.r_1.2.4       RColorBrewer_1.1-2   iterators_1.0.12    
    ## [61] tools_3.6.1          glue_1.3.1           purrr_0.3.3         
    ## [64] yaml_2.2.0           colorspace_1.4-1     VGAM_1.1-2          
    ## [67] knitr_1.25

References
==========

Martinez Molina, D., Jafari, R., Ignatushchenko, M., Seki, T., Larsson,
E.A., Dan, C., Sreekumar, L., Cao, Y., and Nordlund, P. (2013).
Monitoring drug target engagement in cells and tissues using the
cellular thermal shift assay. Science *341*, 84–87.

Mateus, A., Kurzawa, N., Becher, I., Sridharan, S., Helm, D., Stein, F.,
Typas, A., and Savitski, M.M. (2020). Thermal proteome profiling for
interrogating protein interactions. Molecular Systems Biology 16, e9232.

Savitski, M.M., Reinhard, F.B.M., Franken, H., Werner, T., Savitski,
M.F., Eberhard, D., Martinez Molina, D., Jafari, R., Dovega, R.B.,
Klaeger, S., et al. (2014). Tracking cancer drugs in living cells by
thermal profiling of the proteome. Science *346, 6205*, 1255784.

Tan, C.S.H., Go, K.D., Bisteau, X., Dai, L., Yong, C.H., Prabhu, N.,
Ozturk, M.B., Lim, Y.T., Sreekumar, L., Lengqvist, J., et al. (2018).
Thermal proximity coaggregation for system-wide profiling of protein
complex dynamics in cells. Science *359, 6380*, 1170–1177.