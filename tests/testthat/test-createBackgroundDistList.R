## Tests for .createBackgroundDistList
## library(Rtpca); library(testthat)
context(".createBackgroundDistList")

library(dplyr)
library(Biobase)

dist_mat <- as.matrix(dist(matrix(1:20, ncol = 4)))

set.seed(1)
backSamplList <- Rtpca:::.createBackgroundDistList(
    distMat = dist_mat, nMemVec = 3:4)

expect_equal(
    round(sapply(backSamplList, mean), 3),
    c("3" = 4.002, "4" = 3.997)
)
