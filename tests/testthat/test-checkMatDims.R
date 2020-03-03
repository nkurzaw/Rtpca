## Tests for .checkMatDims
## library(Rtpca); library(testthat)
context(".checkMatDims")

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

expect_error(
    Rtpca:::.checkMatDims(
        mat_list 
    )
)

filt_list <- Rtpca:::.getMatList(
    mat_list, 
    Rtpca:::.getCommonRownames(mat_list)
)

expect_invisible(
    Rtpca:::.checkMatDims(
        filt_list 
    )
)
