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
                          DistMat = "matrix", 
                          ContrastDistMat = "matrix",
                          ComplexAnnotation = "data.frame",
                          ComplexBackgroundDistributionList = "list",
                          PPiAnnotation = "data.frame",
                          PPiRocTable = "data.frame", 
                          PPiRocTableAnno = "data.frame", 
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
              
              cat('Slot "ObjList":', names(ObjList(object)), 
                  "of class", class(ObjList(object)), 
                  "and length", length(ObjList(object)), "\n")
              
              cat('Slot "ContrastList":', names(ContrastList(object)), 
                  "of class", class(ContrastList(object)),
                  "and length", length(ContrastList(object)),"\n")
              
              cat('Slot "DistMat" with dimension:', dim(DistMat(object)), "\n") 

              cat('Slot "ContrastDistMat" with dimension:', 
                  dim(ContrastDistMat(object)), "\n") 
              
              cat('Slot "ComplexAnnotation" of class:', 
                  class(ComplexAnnotation(object)), 
                  'with dim:', dim(ComplexAnnotation(object)), "\n")
              
              cat('Slot "ComplexBackgroundDistributionList" of class:', 
                  class(ComplexBackgroundDistributionList(object)), 
                  'with length:', length(
                      ComplexBackgroundDistributionList(object)), "\n")
              
              cat('Slot "PPiAnnotation" of class:', 
                  class(PPiAnnotation(object)), 
                  'with dim:', dim(PPiAnnotation(object)), "\n")
              
              cat('Slot "PPiRocTable" of class:', 
                  class(PPiRocTable(object)), 
                  'with dim:', dim(PPiRocTable(object)), "\n")
              
              cat('Slot "PPiRocTableAnno" of class:', 
                  class(PPiRocTableAnno(object)), 
                  'with dim:', dim(PPiRocTableAnno(object)), "\n")
              
              cat('Slot "ComplexRocTable" of class:', 
                  class(ComplexRocTable(object)), 
                  'with dim:', dim(ComplexRocTable(object)), "\n")
              
              cat('Slot "summaryMethod":', SummaryMethod(object), "\n")
              
              cat('Slot "distMethod":', 
                  ifelse(!is.na(DistMethod(object)), 
                         DistMethod(object), "not specified") , "\n")
              
              cat('Slot "tpcaResultTable" of class:', 
                  class(tpcaResultTable(object)), 
                  'with dim:', dim(tpcaResultTable(object)), "\n")
              
              cat('Slot "diffTpcaResultTable" of class:', 
                  class(diffTpcaResultTable(object)), 
                  'with dim:', dim(diffTpcaResultTable(object)), "\n")
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generics
###

#### Get ObjList
setGeneric("ObjList", function(object)
    standardGeneric("ObjList"))

#' Extract ObjList
#' @aliases ObjList
#' @param object an object of class tpcaResult
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

#### Get ContrastList
setGeneric("ContrastList", function(object)
    standardGeneric("ContrastList"))

#' Extract ContrastList
#' @aliases ContrastList
#' @param object an object of class tpcaResult
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
#' ContrastList(tpcaObj)
setMethod("ContrastList", "tpcaResult", function(object){
    object@ContrastList
})

#### Get CtrlCondName
setGeneric("CtrlCondName", function(object)
    standardGeneric("CtrlCondName"))

#' Extract CtrlCondName
#' @aliases CtrlCondName
#' @param object and object of class tpcaResult
#' @return a character string describing the control 
#' condition
#' @export
#' @examples
#' tpcaObj <- new("tpcaResult")
#' CtrlCondName(tpcaObj)
setMethod("CtrlCondName", "tpcaResult", function(object){
    object@CtrlCondName
})

#### Set CtrlCondName
setGeneric("SetCtrlCondName", function(object, name)
    standardGeneric("SetCtrlCondName"))

#' Set CtrlCondName
#' @aliases SetCtrlCondName
#' @param object an object of class tpcaResult
#' @param name a character string describing the control 
#' condition
#' @return an object of class tpcaResult
#' @export
#' @examples
#' tpcaObj <- new("tpcaResult")
#' SetCtrlCondName(tpcaObj, "DMSO")
setMethod("SetCtrlCondName", "tpcaResult", function(object, name){
    object@CtrlCondName <- name
    return(object)
})


#### Get ContrastCondName
setGeneric("ContrastCondName", function(object)
    standardGeneric("ContrastCondName"))

#' Extract ContrastCondName
#' @aliases ContrastCondName
#' @param object and object of class tpcaResult
#' @return a character string describing the contrast 
#' condition
#' @export
#' @examples
#' tpcaObj <- new("tpcaResult")
#' ContrastCondName(tpcaObj)
setMethod("ContrastCondName", "tpcaResult", function(object){
    object@ContrastCondName
})

#### Set ContrastCondName
setGeneric("SetContrastCondName", function(object, name)
    standardGeneric("SetContrastCondName"))

#' Set ContrastCondName
#' @aliases SetContrastCondName
#' @param object an object of class tpcaResult
#' @param name a character string describing the contrast 
#' condition
#' @return an object of class tpcaResult
#' @export
#' @examples
#' tpcaObj <- new("tpcaResult")
#' SetContrastCondName(tpcaObj, "DMSO")
setMethod("SetContrastCondName", "tpcaResult", function(object, name){
    object@ContrastCondName <- name
    return(object)
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

#### Get PPiRocTable
setGeneric("PPiRocTable", function(object)
    standardGeneric("PPiRocTable"))

#' Extract PPiRocTable
#' @aliases PPiRocTable
#' @param object an object of class tpcaResult
#' @return a data frame containing the results
#' from a tpca analysis
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
#' PPiRocTable(tpcaObj)
setMethod("PPiRocTable", "tpcaResult", function(object){
    object@PPiRocTable
})

#### Set PPiRocTable
setGeneric("SetPPiRocTable", function(object, df)
    standardGeneric("SetPPiRocTable"))

#' Set PPiRocTable
#' @aliases SetPPiRocTable
#' @param object an object of class tpcaResult
#' @param df data.frame containg PPiRocTable to set
#' @return an object of class tpcaResult
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
#' SetPPiRocTable(tpcaObj, data.frame(FPR = 1, TPR = 1))
setMethod("SetPPiRocTable", "tpcaResult", function(object, df){
    object@PPiRocTable <- df
    return(object)
})

#### Get PPiRocTableAnno
setGeneric("PPiRocTableAnno", function(object)
    standardGeneric("PPiRocTableAnno"))

#' Extract PPiRocTableAnno
#' @aliases PPiRocTableAnno
#' @param object an object of class tpcaResult
#' @return a data frame containing annotation 
#' information for PPiRocTable
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
#' PPiRocTableAnno(tpcaObj)
setMethod("PPiRocTableAnno", "tpcaResult", function(object){
    object@PPiRocTableAnno
})

#### Set PPiRocTable
setGeneric("SetPPiRocTableAnno", function(object, df)
    standardGeneric("SetPPiRocTableAnno"))

#' Set PPiRocTableAnno
#' @aliases SetPPiRocTableAnno
#' @param object an object of class tpcaResult
#' @param df data.frame containg PPiRocTable 
#' annotation to set
#' @return an object of class tpcaResult
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
#' SetPPiRocTableAnno(tpcaObj, data.frame(pair = "A:B"))
setMethod("SetPPiRocTableAnno", "tpcaResult", function(object, df){
    object@PPiRocTableAnno <- df
    return(object)
})


#### Get ComplexBackgroundDistributionList
setGeneric("ComplexBackgroundDistributionList", function(object)
    standardGeneric("ComplexBackgroundDistributionList"))

#' Extract ComplexBackgroundDistributionList
#' @aliases ComplexBackgroundDistributionList
#' @param object and object of class tpcaResult
#' @return a list of data frames containing distances
#' of random complexes with different number of members
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
#' ComplexBackgroundDistributionList(tpcaObj)
setMethod("ComplexBackgroundDistributionList", "tpcaResult", function(object){
    object@ComplexBackgroundDistributionList
})

#### Set ComplexBackgroundDistributionList
setGeneric("SetComplexBackgroundDistributionList", function(object, lt)
    standardGeneric("SetComplexBackgroundDistributionList"))

#' Set ComplexBackgroundDistributionList
#' @aliases SetComplexBackgroundDistributionList
#' @param object an object of class tpcaResult
#' @param lt a list of data frames containing distances
#' of random complexes with different number of members
#' @return an object of class tpcaResult
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
#' SetComplexBackgroundDistributionList(tpcaObj, 
#'   list('3' = data.frame(pair = "A:B")))
setMethod("SetComplexBackgroundDistributionList", "tpcaResult", 
          function(object, lt){
    object@ComplexBackgroundDistributionList <- lt
    return(object)
})

#### Get ComplexRocTable
setGeneric("ComplexRocTable", function(object)
    standardGeneric("ComplexRocTable"))

#' Extract ComplexRocTable
#' @aliases ComplexRocTable
#' @param object and object of class tpcaResult
#' @return a data frame containing a complex
#' analysis roc table
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
#' ComplexRocTable(tpcaObj)
setMethod("ComplexRocTable", "tpcaResult", function(object){
    object@ComplexRocTable
})

#### Set ComplexRocTable
setGeneric("SetComplexRocTable", function(object, df)
    standardGeneric("SetComplexRocTable"))

#' Set ComplexRocTable
#' @aliases SetComplexRocTable
#' @param object and object of class tpcaResult
#' @param df data.frame containg ComplexRocTable to set
#' @return a data frame containing the results
#' from a tpca analysis
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
#' SetComplexRocTable(tpcaObj, data.frame(FPR = 1, TPR = 1))
setMethod("SetComplexRocTable", "tpcaResult", function(object, df){
    object@ComplexRocTable <- df
    return(object)
})

#### Get DistMat
setGeneric("DistMat", function(object)
    standardGeneric("DistMat"))

#' Extract DistMat
#' @aliases DistMat
#' @param object an object of class tpcaResult
#' @return a matrix containing the distance matrix
#' of all pariwise protein-protein melting curve 
#' distances computed from a TPP experiment
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
#' DistMat(tpcaObj)
setMethod("DistMat", "tpcaResult", function(object){
    object@DistMat
})

#### Set DistMat
setGeneric("SetDistMat", function(object, mat)
    standardGeneric("SetDistMat"))

#' Set DistMat
#' @aliases SetDistMat
#' @param object an object of class tpcaResult
#' @param mat matrix containg distance matrix to set
#' @return a data frame containing the results
#' from a tpca analysis
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
#' SetDistMat(tpcaObj, matrix(c(0, 1, 0, 1), ncol = 2))
setMethod("SetDistMat", "tpcaResult", function(object, mat){
    object@DistMat <- mat
    return(object)
})

#### Get ContrastDistMat
setGeneric("ContrastDistMat", function(object)
    standardGeneric("ContrastDistMat"))

#' Extract ContrastDistMat
#' @aliases ContrastDistMat
#' @param object an object of class tpcaResult
#' @return a matrix containing the constrast 
#' distance matrix of all pariwise protein-protein 
#' melting curve distances computed from a TPP 
#' experiment
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
#' ContrastDistMat(tpcaObj)
setMethod("ContrastDistMat", "tpcaResult", function(object){
    object@ContrastDistMat
})

#### Set DistMat
setGeneric("SetContrastDistMat", function(object, mat)
    standardGeneric("SetContrastDistMat"))

#' Set ContrastDistMat
#' @aliases SetContrastDistMat
#' @param object and object of class tpcaResult
#' @param mat matrix containg contrast distance 
#' matrix to set
#' @return a data frame containing the results
#' from a tpca analysis
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
#' SetContrastDistMat(tpcaObj, matrix(c(0, 1, 0, 1), ncol = 2))
setMethod("SetContrastDistMat", "tpcaResult", function(object, mat){
    object@ContrastDistMat <- mat
    return(object)
})

#### Get PPiAnnotation
setGeneric("PPiAnnotation", function(object)
    standardGeneric("PPiAnnotation"))

#' Extract PPiAnnotation
#' @aliases PPiAnnotation
#' @param object and object of class tpcaResult
#' @return a data frame containing the results
#' from a tpca analysis
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
#' PPiAnnotation(tpcaObj)
setMethod("PPiAnnotation", "tpcaResult", function(object){
    object@PPiAnnotation
})

#### Get ComplexAnnotation
setGeneric("ComplexAnnotation", function(object)
    standardGeneric("ComplexAnnotation"))

#' Extract ComplexAnnotation
#' @aliases ComplexAnnotation
#' @param object and object of class tpcaResult
#' @return a data frame containing the complex
#' annotation
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
#' ComplexAnnotation(tpcaObj)
setMethod("ComplexAnnotation", "tpcaResult", function(object){
    object@ComplexAnnotation
})

#### Set ComplexAnnotation
setGeneric("SetComplexAnnotation", function(object, df)
    standardGeneric("SetComplexAnnotation"))

#' Set ComplexAnnotation
#' @aliases SetComplexAnnotation
#' @param object an object of class tpcaResult
#' @param df data frame containing complex 
#' annotation
#' @return an object of class tpcaResult
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
#' SetComplexAnnotation(tpcaObj, data.frame(id = "complex1"))
setMethod("SetComplexAnnotation", "tpcaResult", function(object, df){
    object@ComplexAnnotation <- df
    return(object)
})

#### Get distMethod
setGeneric("DistMethod", function(object, method)
    standardGeneric("DistMethod"))

#' Extract DistMethod
#' @aliases DistMethod
#' @param object and object of class tpcaResult
#' @return a character string of the dist method
#' @export
#' @examples
#' tpcaObj <- new("tpcaResult")
#' DistMethod(tpcaObj)
setMethod("DistMethod", "tpcaResult", function(object){
    object@distMethod
})

#### Set distMethod
setGeneric("SetDistMethod", function(object, method)
    standardGeneric("SetDistMethod"))

#' Set distMethod
#' @aliases SetDistMethod
#' @param object an object of class tpcaResult
#' @param method character string of dist method
#' @return an object of class tpcaResult
#' @export
#' @examples
#' tpcaObj <- new("tpcaResult")
#' SetDistMethod(tpcaObj, "euclidean")
setMethod("SetDistMethod", "tpcaResult", function(object, method){
    object@distMethod <- method
    return(object)
})

#### Get SummaryMethod
setGeneric("SummaryMethod", function(object, method)
    standardGeneric("SummaryMethod"))

#' Extract SummaryMethod
#' @aliases SummaryMethod
#' @param object and object of class tpcaResult
#' @return a character string of the summary method
#' @export
#' @examples
#' tpcaObj <- new("tpcaResult")
#' SummaryMethod(tpcaObj)
setMethod("SummaryMethod", "tpcaResult", function(object){
    object@summaryMethod
})

#### Set summaryMethod
setGeneric("SetSummaryMethod", function(object, method)
    standardGeneric("SetSummaryMethod"))

#' Set summaryMethod
#' @aliases SetSummaryMethod
#' @param object an object of class tpcaResult
#' @param method character string of summary method
#' @return an object of class tpcaResult
#' @export
#' @examples
#' tpcaObj <- new("tpcaResult")
#' SetSummaryMethod(tpcaObj, "median")
setMethod("SetSummaryMethod", "tpcaResult", function(object, method){
    object@summaryMethod <- method
    return(object)
})

#### Get tpcaResultTable
setGeneric("tpcaResultTable", function(object)
    standardGeneric("tpcaResultTable"))

#' Extract tpcaResultTable
#' @aliases tpcaResultTable
#' @param object an object of class tpcaResult
#' @return a data frame containing the results
#' from a tpca analysis
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
#' tpcaResultTable(tpcaObj)
setMethod("tpcaResultTable", "tpcaResult", function(object){
    object@tpcaResultTable
})

#### Set tpcaResultTable
setGeneric("SetTpcaResultTable", function(object, df)
    standardGeneric("SetTpcaResultTable"))

#' Set TpcaResultTable
#' @aliases SetTpcaResultTable
#' @param object an object of class tpcaResult
#' @param df a data frame containing the results
#' from a tpca analysis
#' @return an object of class tpcaResult
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
#' SetTpcaResultTable(tpcaObj, data.frame(pair = "A:B"))
setMethod("SetTpcaResultTable", "tpcaResult", function(object, df){
    object@tpcaResultTable <- df
    return(object)
})


#### Get diffTpcaResultTable
setGeneric("diffTpcaResultTable", function(object)
    standardGeneric("diffTpcaResultTable"))

#' Extract diffTpcaResultTable
#' @aliases diffTpcaResultTable
#' @param object an object of class tpcaResult
#' @return a data frame containing the results
#' from a diffTpca analysis
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
#' diffTpcaResultTable(tpcaObj)
setMethod("diffTpcaResultTable", "tpcaResult", function(object){
    object@diffTpcaResultTable
})

#### Set diffTpcaResultTable
setGeneric("SetDiffTpcaResultTable", function(object, df)
    standardGeneric("SetDiffTpcaResultTable"))

#' Set diffTpcaResultTable
#' @aliases SetDiffTpcaResultTable
#' @param object an object of class tpcaResult
#' @param df a data frame containing the results
#' from a differential tpca analysis
#' @return an object of class tpcaResult
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
#' SetDiffTpcaResultTable(tpcaObj, data.frame(pair = "A:B"))
setMethod("SetDiffTpcaResultTable", "tpcaResult", function(object, df){
    object@diffTpcaResultTable <- df
    return(object)
})