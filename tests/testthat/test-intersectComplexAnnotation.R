## Tests for .intersectComplexAnnotation
## library(Rtpca); library(testthat)
context(".intersectComplexAnnotation")

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

complexAnno <- tibble(
    id = "complex1",
    protein = as.character(c(1:3, 5)),
    count = 4
)

tpcaTest1 <- 
    new("tpcaResult", 
    ObjList = mat_list,
    ComplexAnnotation = complexAnno)

tpcaTest1 <- Rtpca:::.createDistMatTpcaObj(
    tpcaTest1
)

tpcaTest2 <- Rtpca:::.intersectComplexAnnotation(
    tpcaTest1,
    minCount = 2
)

expect_equal(
    tpcaTest2@ComplexAnnotation$protein,
    tpcaTest2@CommonFeatures
)

tpcaTest3 <- Rtpca:::.intersectComplexAnnotation(
    tpcaTest1,
    minCount = 3
)

expect_equal(
    tpcaTest3@ComplexAnnotation$protein,
    character(0)
)
