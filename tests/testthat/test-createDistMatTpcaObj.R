## Tests for .createDistMatTpcaObj
## library(Rtpca); library(testthat)
context(".createDistMatTpcaObj")

library(dplyr)
library(Biobase)

m1 <- matrix(1:12, ncol = 4)
m2 <- matrix(2:13, ncol = 4)
m3 <- matrix(c(2:10, 1:7), ncol = 4)

rownames(m1) <- 1:3
rownames(m2) <- 2:4
rownames(m3) <- 2:5

colnames(m1) <- paste0("X", 1:4)
colnames(m2) <- paste0("X", 1:4)
colnames(m3) <- paste0("X", 1:4)

mat_list <- list(
    m1, m2, m3
)

ppi_anno <- tibble(
    gene_name1 = "2",
    gene_name2 = "3",
    combined_score = 700,
    pair = "2:3")

tpcaTest1 <- new(
    "tpcaResult",
    ObjList = mat_list,
    ComplexAnnotation = ppi_anno
)

tpcaTest2 <- Rtpca:::.createDistMatTpcaObj(
    tpcaObj = tpcaTest1)

#M2 <- writeHDF5Array(matrix(c(0, 2, 0, 2), ncol = 2))
M2 <- matrix(c(0, 2, 2, 0), byrow = TRUE, ncol = 2)

expect_equivalent(
    tpcaTest2@DistMat,
    M2
)
