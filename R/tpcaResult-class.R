#' S4 TPCA Result Class
#'
#' @slot ObjList list.
#' @slot ContrastList list.
#' @slot CtrlCondName character.
#' @slot ContrastCondName character.
#' @slot DistMat matrix.
#' @slot ContrastDistMat matrix
#' @slot CommonFeatures vector.
#' @slot ComplexAnnotation data.frame.
#' @slot ComplexBackgroundDistributionList list.
#' @slot PPiAnnotation data.frame.
#' @slot PPiRocTable data.frame
#' @slot PPiRocTableAnno data.frame
#' @slot ComplexRocTable data.frame
#' @slot summaryMethod character.
#' @slot distMethod character.
#' @slot tpcaResultTable data.frame.
#' @slot diffTpcaResultTable data.frame.
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
                          CtrlCondName = "character",
                          ContrastCondName = "character",
                          CommonFeatures = "vector",
                          DistMat = "matrix", #"DelayedMatrix",
                          ContrastDistMat = "matrix",
                          ComplexAnnotation = "data.frame",
                          ComplexBackgroundDistributionList = "list",
                          PPiAnnotation = "data.frame",
                          PPiRocTable = "data.frame", #"DelayedMatrix",
                          PPiRocTableAnno = "data.frame", #"DelayedMatrix",
                          ComplexRocTable = "data.frame",
                          summaryMethod = "character",
                          distMethod = "character",
                          tpcaResultTable = "data.frame",
                          diffTpcaResultTable = "data.frame"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "tpcaResult",
          function(object)
          {
              cat("class:", class(object), "\n")
              
              cat('Slot "ObjList":', names(object@ObjList), 
                  "of class", class(object@ObjList), 
                  "and length", length(object@ObjList), "\n")
              
              cat('Slot "ContrastList":', names(object@ContrastList), 
                  "of class", class(object@ContrastList),
                  "and length", length(object@ContrastList),"\n")
              
              cat('Slot "DistMat" with dimension:', dim(object@DistMat), "\n") 
              # cat('head(DistMat):', "\n")
              # max_col <- ifelse(ncol(object@DistMat) < 5, ncol(object@DistMat), 5)
              # max_row <- ifelse(nrow(object@DistMat) < 5, nrow(object@DistMat), 5)
              # for(i in seq_len(max_row)){
              #     cat(object@DistMat[i, 1:max_col], "\n")
              # }
              
              cat('Slot "ContrastDistMat" with dimension:', dim(object@ContrastDistMat), "\n") 
              
              cat('Slot "ComplexAnnotation" of class:', class(object@ComplexAnnotation), 
                  'with dim:', dim(object@ComplexAnnotation), "\n")
              
              cat('Slot "ComplexBackgroundDistributionList" of class:', class(object@ComplexBackgroundDistributionList), 
                  'with length:', length(object@ComplexBackgroundDistributionList), "\n")
              
              cat('Slot "PPiAnnotation" of class:', class(object@PPiAnnotation), 
                  'with dim:', dim(object@PPiAnnotation), "\n")
              
              cat('Slot "PPiRocTable" of class:', class(object@PPiRocTable), 
                  'with dim:', dim(object@PPiRocTable), "\n")
              
              cat('Slot "PPiRocTableAnno" of class:', class(object@PPiRocTableAnno), 
                  'with dim:', dim(object@PPiRocTableAnno), "\n")
              
              cat('Slot "ComplexRocTable" of class:', class(object@ComplexRocTable), 
                  'with dim:', dim(object@ComplexRocTable), "\n")
              
              cat('Slot "summaryMethod":', object@summaryMethod, "\n")
              
              cat('Slot "distMethod":', ifelse(!is.na(object@distMethod), 
                                               object@distMethod, "not specified") , "\n")
              
              cat('Slot "tpcaResultTable" of class:', class(object@tpcaResultTable), 
                  'with dim:', dim(object@tpcaResultTable), "\n")
              
              cat('Slot "diffTpcaResultTable" of class:', class(object@diffTpcaResultTable), 
                  'with dim:', dim(object@diffTpcaResultTable), "\n")
          })
