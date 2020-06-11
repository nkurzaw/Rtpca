## Tests for tpcaResult class generic ObjList
## library(Rtpca); library(testthat)
context("ObjList")

m1 <- matrix(1:12, ncol = 4)
m2 <- matrix(2:13, ncol = 4)
m3 <- matrix(c(2:10, 1:7), ncol = 4)

rownames(m1) <- 1:3
rownames(m2) <- 2:4
rownames(m3) <- 2:5

mat_list <- list(
    m1, m2, m3
)
tpcaObj <- new("tpcaResult", ObjList = mat_list)

expect_equal(ObjList(tpcaObj),
             mat_list)