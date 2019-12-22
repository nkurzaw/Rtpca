#' Run the TPCA analysis
#' 
#' @param objList inout list of objects, e.g.
#' ExpressionSets retrieved after TPP data import
#' or matrices or data frames
#' @param complexAnno data frame annotating
#' known protein complexes of interest to test 
#' @param ppiAnno data frame annotation known
#' protein-protein interactions (PPI) to test 
#' @param rownameCol in case the input objects are tibbles
#' this parameter takes in the name (character) of the column 
#' specifying protein names or ids
#' @param summaryFUN function to use to summarize measurements
#' across replicates, default is median
#' @param distMethod method to use within dist function,
#' default is 'euclidean'
#' @param stringScore score cutoff to filter
#' PPIs retrieved from StringDb
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
#' rownames(m1) <- 1:3
#' rownames(m2) <- 2:4
#' rownames(m3) <- 2:5
#' 
#' colnames(m1) <- paste0("X", 1:4)
#' colnames(m2) <- paste0("X", 1:4)
#' colnames(m3) <- paste0("X", 1:4)
#' 
#' mat_list <- list(
#'     m1, m2, m3
#' )
#' 
#' ppi_anno <- tibble(
#'     gene_name1 = "2",
#'     gene_name2 = "3",
#'     combined_score = 700,
#'     pair = "2:3")
#' 
#' runTPCA(
#'     objList = mat_list,
#'     complexAnno = NULL,
#'     ppiAnno = ppi_anno
#' )
runTPCA <- function(objList, 
                    complexAnno = NULL,
                    ppiAnno = NULL,
                    rownameCol = NULL,
                    summaryFUN = median,
                    distMethod = "euclidean",
                    stringScore = 700){
    tpcaObj <- .checkAnnoArguments(
        objList = objList,
        complexAnno = complexAnno,
        ppiAnno = ppiAnno
    )
    tpcaObj <- .createDistMatTpcaObj(
        tpcaObj = tpcaObj, 
        rownameCol = rownameCol, 
        summaryFUN = summaryFUN,
        distMethod = distMethod
    )
    if(!is.null(ppiAnno)){
        tpcaObj <- .createPPiRocTable(
            tpcaObj = tpcaObj,
            stringScore = stringScore
        )
    }
    return(tpcaObj)
}

#' Create distance matrix of all vs all protein 
#' melting profiles
#' 
#' @param objList list of objects suitable for the analysis,
#' currently allowed classes of objects are: matrices,
#' data.frames, tibbles and ExpressionSets
#' @param rownameCol in case the input objects are tibbles
#' this parameter takes in the name (character) of the column 
#' specifying protein names or ids
#' @param summaryFUN function to use to summarize measurements
#' across replicates, default is median
#' @param distMethod method to use within dist function,
#' default is 'euclidean'
#' 
#' @return a distance matrix of all pairwise protein
#' melting profiles
#' 
#' @examples 
#' m1 <- matrix(1:12, ncol = 4)
#' m2 <- matrix(2:13, ncol = 4)
#' m3 <- matrix(c(2:10, 1:7), ncol = 4)
#' 
#' rownames(m1) <- 1:3
#' rownames(m2) <- 2:4
#' rownames(m3) <- 2:5
#' 
#' colnames(m1) <- paste0("X", 1:4)
#' colnames(m2) <- paste0("X", 1:4)
#' colnames(m3) <- paste0("X", 1:4)
#' 
#' mat_list <- list(
#'     m1, m2, m3
#' )
#' 
#' createDistMat(mat_list)
#' 
#' expr1 <- ExpressionSet(m1)
#' expr2 <- ExpressionSet(m2)
#' expr3 <- ExpressionSet(m3)
#' 
#' exprSet_list <- list(
#'     expr1, expr2, expr3
#' )
#' 
#' createDistMat(exprSet_list)
#' 
#' @export
#' @importFrom stats dist
createDistMat <- function(objList, rownameCol = NULL, 
                           summaryFUN = median,
                           distMethod = "euclidean"){
    common_rownames <-
        .getCommonRownames(objList, rownameCol = rownameCol)
    mat_list <- .getMatList(objList, 
                            commonRownames = common_rownames,
                            rownameCol = rownameCol)
    summarized_mat <- .getSummarizedMat(matList = mat_list,
                                        FUN = summaryFUN)
    dist_mat <- as.matrix(dist(summarized_mat,
                               method = distMethod))
    return(dist_mat)
}

#' Plot PPI ROC curve
#' @param tpcaObj tpcaResult object
#' 
#' @import ggplot2
plotPPiRoc <- function(tpcaObj){
    ggplot(tpcaObj@PPiRocTable, aes(FPR, TPR)) + 
        geom_line() +
        geom_line(aes(x, y),
                  color = "gray",
                  data = tibble(x = 0:1, y = 0:1),
                  linetype = "dashed")
}

.createDistMatTpcaObj <- function(tpcaObj, rownameCol = NULL,
                                   summaryFUN = median,
                                   distMethod = "euclidean"){
    tpcaObj@summaryFUN <- summaryFUN
    tpcaObj@distMethod <- distMethod
    tpcaObj@DistMat <- createDistMat(
        objList = tpcaObj@ObjList, 
        rownameCol = rownameCol, 
        summaryFUN = summaryFUN,
        distMethod = distMethod
    )
    return(tpcaObj)
}

#' @importFrom methods new
.checkAnnoArguments <- function(objList, complexAnno = NULL,
                                ppiAnno = NULL){
    if(is.null(complexAnno) & is.null(ppiAnno)){
        stop(paste(c("Neither complex annotation nor a PPI", 
                     "annotation were supplied!\n",
                     "One of them is required."),
                   collapse = " "))
    }else if(!is.null(complexAnno) & !is.null(ppiAnno)){
        tpcaObj <- new("tpcaResult",
                       ObjList = objList,
                       ComplexAnnotation = complexAnno)
        return(tpcaObj)
    }else if(!is.null(complexAnno)){
        tpcaObj <- new("tpcaResult",
                       ObjList = objList,
                       ComplexAnnotation = complexAnno)
        return(tpcaObj)
    }else if(!is.null(ppiAnno)){
        tpcaObj <- new("tpcaResult",
                       ObjList = objList,
                       PPiAnnotation = ppiAnno)
        return(tpcaObj)
    }
}

#' @import dplyr
#' @import tidyr
.createPPiRocTable <- function(tpcaObj, stringScore = 700){
    colname <- value <- rowname <- pair <- combined_score <- 
        annotated <- NULL
    
    dist_df <- tpcaObj@DistMat %>% 
        tbl_df %>% 
        mutate(rowname = rownames(tpcaObj@DistMat)) %>% 
        gather(colname, value, -rowname) %>% 
        rowwise() %>% 
        mutate(pair = paste(sort(c(rowname, colname)), 
                            collapse = ":")) %>% 
        ungroup %>% 
        filter(rowname != colname,
               !duplicated(pair)) %>% 
        arrange(value) %>% 
        mutate(annotated = pair %in% 
                   filter(tpcaObj@PPiAnnotation, 
                          combined_score >= stringScore)$pair) %>% 
        mutate(TPR = cumsum(as.numeric(annotated)) / 
                    sum(as.numeric(annotated)), 
               FPR = cumsum(as.numeric(annotated == FALSE)) / 
                   (n() - cumsum(as.numeric(value != 0)) +  
                        cumsum(as.numeric(annotated == FALSE))))
    
    tpcaObj@PPiRocTable <- dist_df
    return(tpcaObj)
}


#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase exprs
.getCommonRownames <- function(objList, rownameCol = NULL){
    if(class(objList[[1]])[1] == "ExpressionSet"){
        return(Reduce(intersect, lapply(objList, function(x) 
            rownames(exprs(x)))))
    }else if(class(objList[[1]])[1] %in% c("matrix", "data.frame")){
        return(Reduce(intersect, lapply(objList, function(x) 
            rownames(x))))
    }else if((class(objList[[1]])[1] %in% c("tbl_df")) & 
             !is.null(rownameCol)){
        return(Reduce(intersect, lapply(objList, function(x) 
            x[[rownameCol]])))
    }else{
        stop(paste("getCommonRownames: Supplied object is neither", 
                   "ExpressionSet nor matrix or data.frame!\n"))
    }
}

#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr one_of
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase exprs
.getMatList <- function(objList, commonRownames, rownameCol = NULL){
    if(class(objList[[1]])[1] == "ExpressionSet"){
        return(lapply(objList, function(x){
            e_mat <- exprs(x)
            e_mat[rownames(e_mat) %in% commonRownames,]
        }))
    }else if(class(objList[[1]])[1] %in% c("matrix", "data.frame")){
        return(lapply(objList, function(x){
            as.matrix(x[rownames(x) %in% commonRownames,])
        }))
    }else if((class(objList[[1]])[1] %in% c("tbl_df")) & 
             !is.null(rownameCol)){
        return(lapply(objList, function(x){
            x_fil <- filter(x, eval(parse(text = rownameCol)) %in% 
                                commonRownames) 
            x_mat <- as.matrix(dplyr::select(x_fil, 
                                             -one_of(rownameCol)))
            rownames(x_mat) <- x_fil[[rownameCol]]
            return(x_mat)
        }))
    }else{
        stop(paste("getMatList: Supplied object is neither", 
                   "ExpressionSet nor matrix or data.frame!\n"))
    }
}

.checkMatDims <- function(matList){
    if(length(Reduce(intersect, lapply(matList, dim))) != 2){
        stop("checkMatDims: unequal matrix dimensions!")
    }
}

#' @importFrom stats median
.getSummarizedMat <- function(matList, FUN = median){
    .checkMatDims(matList)
    arr <- array(unlist(matList), 
                 c(dim(matList[[1]]), length(matList)))
    sumMat <- apply(arr, seq_len(2), FUN)
    rownames(sumMat) <- rownames(matList[[1]])
    colnames(sumMat) <- colnames(matList[[1]])
    return(sumMat)
}
