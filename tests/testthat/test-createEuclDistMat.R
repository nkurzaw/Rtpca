## Tests for .getCommonRownames
## library(Rtpca); library(testthat)
context("createEuclDistMat")

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

ctrl_dist_mat <- matrix(c(0, 2 ,2, 0),
                        ncol = 2)
rownames(ctrl_dist_mat) <- 
    colnames(ctrl_dist_mat) <- c("2", "3")

expect_identical(
    createEuclDistMat(mat_list),
    ctrl_dist_mat
)

expr1 <- ExpressionSet(m1)
expr2 <- ExpressionSet(m2)
expr3 <- ExpressionSet(m3)

exprSet_list <- list(
    expr1, expr2, expr3
)

expect_identical(
    createEuclDistMat(exprSet_list),
    ctrl_dist_mat
)

df_list <- list(
    data.frame(m1), 
    data.frame(m2), 
    data.frame(m3)
)

expect_identical(
    createEuclDistMat(exprSet_list), 
    ctrl_dist_mat
)

tbl_list <- list(
    data.frame(m1) %>% 
        tbl_df %>% 
        mutate(gene_name = 1:3), 
    data.frame(m2) %>% 
        tbl_df %>% 
        mutate(gene_name = 2:4), 
    data.frame(m3) %>% 
        tbl_df %>% 
        mutate(gene_name = 2:5)
)

expect_identical(
    createEuclDistMat(tbl_list, 
                      rownameCol = "gene_name"), 
    ctrl_dist_mat
)
