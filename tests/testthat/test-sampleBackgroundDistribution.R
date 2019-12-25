## Tests for .sampleBackgroundDistribution
## library(Rtpca); library(testthat)
context(".sampleBackgroundDistribution")

library(dplyr)
library(Biobase)

dist_mat <- as.matrix(dist(matrix(1:20, ncol = 4)))

set.seed(1)
backSampl <- Rtpca:::.sampleBackgroundDistribution(
    distMat = dist_mat, nMem = 3)

expect_equal(
    mean(backSampl),
    4.002
)
