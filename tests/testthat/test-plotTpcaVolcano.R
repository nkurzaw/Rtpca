## Tests for plotTpcaVolcano
## library(Rtpca); library(testthat)
context("plotTpcaVolcano")

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

complex_anno <- tibble(
    protein = c("3", "4", "5",
       "4", "5", "6", "7"),
    id = c(rep("1", 3), rep("2", 4)),
    count = c(rep(3, 3), rep(4, 4)))

tpca_result <- runTPCA(
  mat_list, complexAnno = complex_anno)

p <- plotTpcaVolcano(tpca_result)

expect_equal(
    class(p),
    c("gg", "ggplot")
)
