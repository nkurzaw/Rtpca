% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{runDiffTPCA}
\alias{runDiffTPCA}
\title{Run differential TPCA analysis}
\usage{
runDiffTPCA(
  objList,
  contrastList,
  ctrlCondName = "control",
  contrastCondName = "treatment",
  ppiAnno = NULL,
  complexAnno = NULL,
  rownameCol = NULL,
  summaryMethodStr = "median",
  distMethodStr = "euclidean",
  n = 10000,
  p_adj_method = "BH"
)
}
\arguments{
\item{objList}{input list of objects, e.g.
ExpressionSets retrieved after TPP data import
or matrices or data frames}

\item{contrastList}{input list of objects for
comparison at e.g. different treatment condtion,
same file formats work as for objList}

\item{ctrlCondName}{character string indicating the name
of the control condition, default is "control"}

\item{contrastCondName}{character string indicating the name
of the contrast condition, default is "treatment"}

\item{ppiAnno}{data frame annotation known
protein-protein interactions (PPI) to test}

\item{complexAnno}{data frame annotating
known protein complexes of interest to test}

\item{rownameCol}{in case the input objects are tibbles
this parameter takes in the name (character) of the column 
specifying protein names or ids}

\item{summaryMethodStr}{character string indicating a method 
to use to summarize measurements across replicates, 
default is "median", other options are c("mean", "rbind")}

\item{distMethodStr}{method to use within dist function,
default is 'euclidean'}

\item{n}{number of random protein pair draws to obtain
empirical p-value, default is 10000}

\item{p_adj_method}{method to be used for multiple
testing adjusment, default is "BH"}
}
\value{
an object of class tpcaResult
with the following slots:
1) ObjList: containing the supplied list of
objects
}
\description{
Run differential TPCA analysis
}
\examples{
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

diff_tpca <- Rtpca:::runDiffTPCA(
  mat_list, contrast_list, ppiAnno = ppi_anno)

}
