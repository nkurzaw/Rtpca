## Tests for .getCommonRownames
## library(Rtpca); library(testthat)
context(".getCommonRownames")

library(dplyr)
library(Biobase)

m1 <- matrix(1:12, ncol = 4)
m2 <- matrix(2:13, ncol = 4)
m3 <- matrix(c(2:10, 1:3), ncol = 4)

rownames(m1) <- 1:3
rownames(m2) <- 2:4
rownames(m3) <- 2:4

mat_list <- list(
    m1, m2, m3
)

expect_identical(
    Rtpca:::.getCommonRownames(mat_list), 
    c("2", "3")
)

df_list <- list(
    data.frame(m1), 
    data.frame(m2), 
    data.frame(m3)
)

expect_identical(
    Rtpca:::.getCommonRownames(df_list), 
    c("2", "3")
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
        mutate(gene_name = 2:4)
)

expect_identical(
    as.character(
        Rtpca:::.getCommonRownames(tbl_list, 
                           rownameCol = "gene_name")
    ), 
    c("2", "3")
)

expr1 <- ExpressionSet(m1)
expr2 <- ExpressionSet(m2)
expr3 <- ExpressionSet(m3)

exprSet_list <- list(
    expr1, expr2, expr3
)

expect_identical(
    Rtpca:::.getCommonRownames(exprSet_list), 
    c("2", "3")
)

expect_error(
    Rtpca:::.getCommonRownames(c("2", "3")), 
    paste("getCommonRownames: Supplied object is neither", 
          "ExpressionSet nor matrix or data.frame!")
)
