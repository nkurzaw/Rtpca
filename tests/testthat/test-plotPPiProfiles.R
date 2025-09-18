## Tests for plotPPiProfiles
## library(Rtpca); library(testthat)
library(Biobase)
context("plotPPiProfiles")

set.seed(12)
m1 <- matrix(rnorm(50), ncol = 10)
m2 <- matrix(rnorm(50), ncol = 10)

rownames(m1) <- letters[1:5]
rownames(m2) <- letters[1:5]

colnames(m1) <- paste("fc", 1:10, sep = "_")
colnames(m2) <- paste("fc", 1:10, sep = "_")

pheno <- data.frame(
    temperature = seq(37, 67, length.out = 10))
rownames(pheno) <- paste("fc", 1:10, sep = "_")

eset1 <- ExpressionSet(
    assayData = m1,
    phenoData = AnnotatedDataFrame(pheno)
)

eset2 <- ExpressionSet(
    assayData = m2,
    phenoData = AnnotatedDataFrame(pheno)
)

tpcaObj1 <- new("tpcaResult",
                ObjList = list(eset1),
                ContrastList = list(eset2),
                CtrlCondName = "control",
                ContrastCondName = "treatment")

p1 <- plotPPiProfiles(tpcaObj1, pair = c("b", "d"))

expect_equal(
    class(p1),
    c("ggplot2::ggplot", "ggplot", "ggplot2::gg", "S7_object", "gg")
)

# objList is list of matrices

attributes(m1)$temperature <- seq(37, 67, length.out = 10)
attributes(m2)$temperature <- seq(37, 67, length.out = 10)

tpcaObj2 <- new("tpcaResult",
               ObjList = list(m1),
               ContrastList = list(m2),
               CtrlCondName = "control",
               ContrastCondName = "treatment")

p2 <- plotPPiProfiles(tpcaObj2, pair = c("b", "d"))

expect_equal(
    class(p2),
    c("ggplot2::ggplot", "ggplot", "ggplot2::gg", "S7_object", "gg")
)
