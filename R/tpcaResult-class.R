#' S4 TPCA Result Class
#'
#' @slot ObjList list.
#' @slot ContrastList list.
#' @slot DistMat matrix.
#' @slot ContrastDistMat matrix.
#' @slot ComplexAnnotation data.frame.
#' @slot PPiAnnotation data.frame.
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
#' tpcaObj <- tpcaResult(ObjList = mat_list)
tpcaResult <- setClass("tpcaResult",
                      slots = list(
                          ObjList = "list",
                          ContrastList = "list",
                          DistMat = "matrix",
                          ContrastDistMat = "matrix",
                          ComplexAnnotation = "data.frame",
                          PPiAnnotation = "data.frame",
                          PPiRocTable = "data.frame",
                          summaryFUN = "function",
                          distMethod = "character"))
