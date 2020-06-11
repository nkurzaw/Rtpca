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
#' objects (e.g. a list of Expression Sets 
#' summarizing a TPP experiment)
#' 2) ContrastList: containing the supplied list 
#' of constrast objects (if supplied for
#' performance of a differential Rtpca analysis)
#' 3) CtrlCondName: character string indicating
#' the control condition, e.g. "control"
#' 4) ContrastCondName: character string indicating
#' the contrast condition, e.g. "drug treatment"
#' 5) DistMat: a matrix containing all pairwise
#' protein-protein distances obtained from comparing
#' their melting curves in the control condition
#' 6) ContrastDistMat: a matrix containing all pairwise
#' protein-protein distances obtained from comparing
#' their melting curves in the contrast condition
#' 7) CommonFeatures: a vector containing the features
#' (proteins) found in common between control and 
#' contrast condition
#' 8) ComplexAnnotation: a data frame supplied by the
#' user annotating protein to protein complexes
#' 9) ComplexBackgroundDistributionList: a list of 
#' distances drawn for random groups of proteins
#' with different number of members
#' 10) PPiAnnotation: a data frame supplied by the
#' user annotating protein-protein interactions
#' 11) PPiRocTable: data frame containing false positive
#' rate and true positive rate based on ranking the TPCA
#' analysis results by euclidean distance of melting curves
#' of protein pairs, annotated PPIs are considered true 
#' positives
#' 12) PPiRocTableAnno: annotation to PPiRocTable
#' 13) ComplexRocTable: data frame containing false positive
#' rate and true positive rate based on ranking the TPCA
#' analysis results by euclidean distance of melting curves
#' of proteins within annotated complexes, annotated complexes
#' are considered true positives, proteins in randomly permuted 
#' complex annotations are considered false positives
#' 14) summaryMethod: character string of summarization method
#' used to summarize data across replicates
#' 15) distMethod: character string of distance method
#' used to compare melting curves of proteins
#' 16) tpcaResultTable: data frame containing the results
#' from a tpca analysis
#' 17) diffTpcaResultTable: data frame containing the results
#' from a differential tpca analysis
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

              cat('Slot "ContrastDistMat" with dimension:', 
                  dim(object@ContrastDistMat), "\n") 
              
              cat('Slot "ComplexAnnotation" of class:', 
                  class(object@ComplexAnnotation), 
                  'with dim:', dim(object@ComplexAnnotation), "\n")
              
              cat('Slot "ComplexBackgroundDistributionList" of class:', 
                  class(object@ComplexBackgroundDistributionList), 
                  'with length:', length(
                      object@ComplexBackgroundDistributionList), "\n")
              
              cat('Slot "PPiAnnotation" of class:', 
                  class(object@PPiAnnotation), 
                  'with dim:', dim(object@PPiAnnotation), "\n")
              
              cat('Slot "PPiRocTable" of class:', 
                  class(object@PPiRocTable), 
                  'with dim:', dim(object@PPiRocTable), "\n")
              
              cat('Slot "PPiRocTableAnno" of class:', 
                  class(object@PPiRocTableAnno), 
                  'with dim:', dim(object@PPiRocTableAnno), "\n")
              
              cat('Slot "ComplexRocTable" of class:', 
                  class(object@ComplexRocTable), 
                  'with dim:', dim(object@ComplexRocTable), "\n")
              
              cat('Slot "summaryMethod":', object@summaryMethod, "\n")
              
              cat('Slot "distMethod":', 
                  ifelse(!is.na(object@distMethod), 
                         object@distMethod, "not specified") , "\n")
              
              cat('Slot "tpcaResultTable" of class:', 
                  class(object@tpcaResultTable), 
                  'with dim:', dim(object@tpcaResultTable), "\n")
              
              cat('Slot "diffTpcaResultTable" of class:', 
                  class(object@diffTpcaResultTable), 
                  'with dim:', dim(object@diffTpcaResultTable), "\n")
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generics
###

#### Get ObjList
setGeneric("ObjList", function(object)
    standardGeneric("ObjList"))

#' Extract ObjList
#' @aliases ObjList
#' @param object and object of class tpcaResult
#' @return an object list containing TPP data
#' @export
#' @examples
#' m1 <- matrix(1:12, ncol = 4)
#' m2 <- matrix(2:13, ncol = 4)
#' m3 <- matrix(c(2:10, 1:7), ncol = 4)
#' 
#' rownames(m1) <- 1:3
#' rownames(m2) <- 2:4
#' rownames(m3) <- 2:5
#' 
#' mat_list <- list(
#'     m1, m2, m3
#' )
#' tpcaObj <- new("tpcaResult", ObjList = mat_list)
#' ObjList(tpcaObj)
setMethod("ObjList", "tpcaResult", function(object){
    object@ObjList
})


#### Get CommonFeatures
setGeneric("CommonFeatures", function(object)
    standardGeneric("CommonFeatures"))

#' Extract CommonFeatures
#' @aliases CommonFeatures
#' @param object and object of class tpcaResult
#' @return a vector of common features across
#' replicates
#' @export
#' @examples
#' m1 <- matrix(1:12, ncol = 4)
#' m2 <- matrix(2:13, ncol = 4)
#' m3 <- matrix(c(2:10, 1:7), ncol = 4)
#' 
#' rownames(m1) <- 1:3
#' rownames(m2) <- 2:4
#' rownames(m3) <- 2:5
#' 
#' mat_list <- list(
#'     m1, m2, m3
#' )
#' 
#' ppi_anno <- tibble(
#'     x = "2",
#'     y = "3",
#'     combined_score = 700,
#'     pair = "2:3")
#' 
#' tpcaObj <- runTPCA(
#'     objList = mat_list,
#'     complexAnno = NULL,
#'     ppiAnno = ppi_anno
#' )
#' 
#' CommonFeatures(tpcaObj)
setMethod("CommonFeatures", "tpcaResult", function(object){
    object@CommonFeatures
})

#### Set CommonFeatures
setGeneric("SetCommonFeatures", 
           function(object, commonFeatures)
    standardGeneric("SetCommonFeatures"))

#' Set CommonFeatures
#' @aliases SetCommonFeatures
#' @param object and object of class tpcaResult
#' @param commonFeatures a vector of characters
#' indicating the common features across replicates
#' @return a vector of common features across
#' replicates
#' @export
#' @examples
#' m1 <- matrix(1:12, ncol = 4)
#' m2 <- matrix(2:13, ncol = 4)
#' m3 <- matrix(c(2:10, 1:7), ncol = 4)
#' 
#' rownames(m1) <- 1:3
#' rownames(m2) <- 2:4
#' rownames(m3) <- 2:5
#' 
#' mat_list <- list(
#'     m1, m2, m3
#' )
#' 
#' tpcaObj <- new("tpcaResult", ObjList = mat_list)
#' SetCommonFeatures(tpcaObj, c("2", "3"))
setMethod("SetCommonFeatures", "tpcaResult", 
          function(object, commonFeatures){
    object@CommonFeatures <- commonFeatures
    return(object)
})

#### Get tpcaResultTable
setGeneric("tpcaResultTable", function(object)
    standardGeneric("tpcaResultTable"))

#' Extract tpcaResultTable
#' @aliases tpcaResultTable
#' @param object and object of class tpcaResult
#' @return a data frame containing the results
#' from a tpca analysis
#' @export
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
#' tpcaResultTable(tpcaObj)
setMethod("tpcaResultTable", "tpcaResult", function(object){
    object@tpcaResultTable
})

#### Get diffTpcaResultTable
setGeneric("diffTpcaResultTable", function(object)
    standardGeneric("diffTpcaResultTable"))

#' Extract diffTpcaResultTable
#' @aliases diffTpcaResultTable
#' @param object and object of class tpcaResult
#' @return a data frame containing the results
#' from a diffTpca analysis
#' @export
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
#' diffTpcaResultTable(tpcaObj)
setMethod("diffTpcaResultTable", "tpcaResult", function(object){
    object@diffTpcaResultTable
})

