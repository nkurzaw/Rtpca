## Tests for tpcaResult class generic CommonFeatures
## library(Rtpca); library(testthat)
context("CommonFeatures")

m1 <- matrix(1:12, ncol = 4)
m2 <- matrix(2:13, ncol = 4)
m3 <- matrix(c(2:10, 1:7), ncol = 4)

rownames(m1) <- 1:3
rownames(m2) <- 2:4
rownames(m3) <- 2:5

mat_list <- list(
    m1, m2, m3
)

ppi_anno <- tibble(
    x = "2",
    y = "3",
    combined_score = 700,
    pair = "2:3")

tpcaObj <- runTPCA(
    objList = mat_list,
    complexAnno = NULL,
    ppiAnno = ppi_anno
)

expect_equal(
    CommonFeatures(tpcaObj),
    c("2", "3")
)