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

m1_new <- m1[-1,]
m2_new <- m2[-3,]
m3_new <- m3[-c(3,4),]

mat_list <- list(
    m1, m2, m3
)

mat_list_new <- list(
    m1_new, m2_new, m3_new
)

expect_identical(
    .getMatList(mat_list, commonRownames = c("2", "3")),
    mat_list_new
)

df_list <- list(
    data.frame(m1), 
    data.frame(m2), 
    data.frame(m3)
)

expect_equal(
    .getMatList(df_list, commonRownames = c("2", "3")), 
    mat_list_new
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
    .getMatList(tbl_list, commonRownames = c(2, 3),
                rownameCol = "gene_name"), 
    mat_list_new
)
