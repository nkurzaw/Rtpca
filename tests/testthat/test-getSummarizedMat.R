## Tests for .getCommonRownames
## library(Rtpca); library(testthat)
context(".getSummarizedMat")

m1 <- matrix(1:12, ncol = 4)
m2 <- matrix(2:13, ncol = 4)
m3 <- matrix(c(2:10, 1:7), ncol = 4)

rownames(m1) <- 1:3
rownames(m2) <- 2:4
rownames(m3) <- 2:5

colnames(m1) <- paste0("X", 1:4)
colnames(m2) <- paste0("X", 1:4)
colnames(m3) <- paste0("X", 1:4)

m1_new <- m1[-1,]
m2_new <- m2[-3,]
m3_new <- m3[-c(3,4),]

mat_list <- list(
    m1, m2, m3
)

mat_list_new <- list(
    m1_new, m2_new, m3_new
)

ctrl_m <- matrix(
    c(2, 3, 5, 6, 8, 9, 11, 12),
    ncol = 4)
rownames(ctrl_m) <- c("2", "3")
colnames(ctrl_m) <- paste0("X", 1:4)

expect_equal(
    Rtpca:::.getSummarizedMat(mat_list_new),
    ctrl_m
)

expect_error(
    Rtpca:::.getSummarizedMat(mat_list),
    "checkMatDims: unequal matrix dimensions!"
)
