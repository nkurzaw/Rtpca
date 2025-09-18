## Tests for plotPPiRoc
## library(Rtpca); library(testthat)
context("plotPPiRoc")

#hdf5M <- writeHDF5Array(x = 
ppiRocDf <- tibble(
    TPR = c(0, 0.1, 0.2, 0.4, 0.5, 0.7, 0.9, 1),
    FPR = c(0, 0.05, 0.1, 0.2, 0.5, 0.7, 0.9, 1)
)#, chunkdim = c(8, 2))

#colnames(hdf5M) <- c("TPR", "FPR")

tpcaTest <- new(
    "tpcaResult",
    PPiRocTable = ppiRocDf #hdf5M
    )

p <- plotPPiRoc(tpcaTest)

expect_equal(
    class(p),
    c("ggplot2::ggplot", "ggplot", "ggplot2::gg", "S7_object", "gg")
)
