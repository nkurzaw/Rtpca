## Tests for runDiffTPCA
## library(Rtpca); library(testthat)
context("runDiffTPCA")

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

c1 <- matrix(29:2, ncol = 4)
c2 <- matrix(26:3, ncol = 4)
c3 <- matrix(c(11:3, 20:2), ncol = 4)

rownames(c1) <- 1:7
rownames(c2) <- 3:8
rownames(c3) <- 2:8

contrast_list <- list(
    c1, c2, c3
)

ppi_anno <- tibble(
    x = c("3", "3"),
    y = c("5", "7"),
    pair = c("3:5", "3:7"))

ref_df <- tibble(
    pair = c("3:5", "3:7"),
    valueC2 = c(4, 8)
)

diff_tpca <- runDiffTPCA(mat_list, contrast_list, 
                         ppiAnno = ppi_anno)

expect_identical(
    diff_tpca@diffTpcaResultTable[,c(1,3)],
    ref_df
)
