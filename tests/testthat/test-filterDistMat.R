## Tests for .filterDistMat
## library(Rtpca); library(testthat)
context(".filterDistMat")

library(dplyr)
library(Biobase)

m1 <- matrix(1:28, ncol = 4)
m2 <- matrix(2:25, ncol = 4)
m3 <- matrix(c(2:10, 1:19), ncol = 4)

rownames(m1) <- 1:7
rownames(m2) <- 3:8
rownames(m3) <- 2:8

mat_list <- list(
    m1, m2, m3
)

dist_mat <- Rtpca:::createDistMat(
    mat_list
)

ppi_anno <- tibble(
    x = "3",
    y = "4",
    combined_score = 700,
    pair = "3:4")

ref_mat <- matrix(2)
dimnames(ref_mat) <- list("3", "4")

expect_identical(
    Rtpca:::.filterDistMat(dist_mat, ppi_anno = ppi_anno),
    ref_mat
)
