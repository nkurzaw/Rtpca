## Tests for plotPPiRoc
## library(Rtpca); library(testthat)
context("plotPPiRoc")

tpcaTest <- new(
    "tpcaResult",
    PPiRocTable = tibble(
        TPR = c(0, 0.1, 0.2, 0.4, 0.5, 0.7, 0.9, 1),
        FPR = c(0, 0.05, 0.1, 0.2, 0.5, 0.7, 0.9, 1)
    ))

p <- plotPPiRoc(tpcaTest)

expect_equal(
    class(p),
    c("gg", "ggplot")
)
