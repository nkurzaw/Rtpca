#' S4 TPCA Result Class
#'
#' @slot ObjList list.
#' @slot ContrastList list.
#' @slot DistMat matrix.
#' @slot ContrastDistMat matrix.
#' @slot CommonFeatures vector.
#' @slot ComplexAnnotation data.frame.
#' @slot ComplexBackgroundDistributionList list.
#' @slot PPiAnnotation data.frame.
#' @slot PPiRocTable data.frame.
#' @slot summaryFUN function.
#' @slot distMethod character.
#'
#' @return an object of class tpcaResult
#' with the following slots:
#' 1) ObjList: containing the supplied list of
#' objects
#'
#' @export
#'
#' @examples
#' m1 <- matrix(1:12, ncol = 4)
#' m2 <- matrix(2:13, ncol = 4)
#' m3 <- matrix(c(2:10, 1:7), ncol = 4)
#' 
#'rownames(m1) <- 1:3
#' rownames(m2) <- 2:4
#' rownames(m3) <- 2:5
#' 
#' mat_list <- list(
#'     m1, m2, m3
#' )
#' tpcaObj <- new("tpcaResult", ObjList = mat_list)
tpcaResult <- setClass("tpcaResult",
                      slots = list(
                          ObjList = "list",
                          ContrastList = "list",
                          CommonFeatures = "vector",
                          DistMat = "matrix",
                          ContrastDistMat = "matrix",
                          ComplexAnnotation = "data.frame",
                          ComplexBackgroundDistributionList = "list",
                          PPiAnnotation = "data.frame",
                          PPiRocTable = "data.frame",
                          summaryFUN = "function",
                          distMethod = "character"))
