## Tests for .distMat2AnnotatedPPiDf
## library(Rtpca); library(testthat)
context(".distMat2AnnotatedPPiDf")

library(dplyr)
library(Biobase)

m <- matrix(c(0, 2, 2, 0), ncol = 2)
dimnames(m) <- list(c("2", "3"),
                    c("2", "3"))

ppi_anno <- tibble(
    gene_name1 = "2",
    gene_name2 = "3",
    combined_score = 700,
    pair = "2:3")

tst_df <- Rtpca:::.distMat2AnnotatedPPiDf(
    dist_mat = m,
    ppi_anno = ppi_anno
)

ctrl_df <- tibble(
    rowname = "3",
    key = "2",
    value = 2,
    pair = "2:3"
)

expect_equal(
    tst_df,
    ctrl_df
)
