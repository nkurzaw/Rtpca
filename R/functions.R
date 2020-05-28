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
#' @param summaryMethod character string indicating a method 
#' to use to summarize measurements across replicates, 
#' default is "median", other options are c("mean", "rbind")
#' @param distMethod method to use within dist function,
#' default is 'euclidean'
#' @param doRocAnalysis logical indicating whether a ROC 
#' analysis should be performed which can be used to assess 
#' the predictive power of the dataset for protein-protein
#' interactions / protein complexes based on distanc between
#' melting curves of protein interactions partners
#' @param minCount integer indicating how many subunits 
#' of a complex should be qunatified to inlucde it into the
#' analysis, default is 3
#' @param nSamp integer indicating the number of random samples
#' which should be performed to estimate empirical null 
#' distributions, default is 10000
#' @param p_adj_method character string indicating a valid 
#' method to be used for multiple testing adjusment, default 
#' is "BH" which makes p.adjust use benjamini-hochberg, for 
#' additional options check ?p.adjust
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
#'     x = "2",
#'     y = "3",
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
                    summaryMethod = "median",
                    distMethod = "euclidean",
                    doRocAnalysis = TRUE,
                    minCount = 3,
                    nSamp = 10000,
                    p_adj_method = "BH"
                    ){
    message("Checking input arguments. \n")
    tpcaObj <- .checkAnnoArguments(
        objList = objList,
        complexAnno = complexAnno,
        ppiAnno = ppiAnno
    )
    message("\nCreating distance matrices. \n")
    tpcaObj <- .createDistMatTpcaObj(
        tpcaObj = tpcaObj, 
        rownameCol = rownameCol, 
        summaryMethod = summaryMethod,
        distMethod = distMethod
    )
    if(!is.null(ppiAnno)){
        message("\nTesting for complex co-aggregation. \n")
        tpcaObj <- .testForPPiCoAggregation(
            tpcaObj = tpcaObj, 
            nSamp = nSamp,
            p_adj_method = p_adj_method)
        if(doRocAnalysis){
            message("\nPerforming PPi ROC analysis. \n")
            tpcaObj <- .createPPiRocTable(
                tpcaObj = tpcaObj
            )
        }
    }
    if(!is.null(complexAnno)){
        message("\nTesting for complex co-aggregation. \n")
        tpcaObj <- .testForComplexCoAggregation(
            tpcaObj = tpcaObj, 
            minCount = minCount, 
            nSamp = nSamp,
            p_adj_method = p_adj_method
        )
        if(doRocAnalysis){
            message("\nPerformning Complex ROC analysis. \n")
            tpcaObj <- .createComplexRocTable(
                tpcaObj = tpcaObj
            )
        }
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
#' @param summaryMethod character string indicating a method 
#' to use to summarize measurements across replicates, 
#' default is "median", other options are c("mean", "rbind")
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
                           summaryMethod = "median",
                           distMethod = "euclidean"){
    common_rownames <-
        .getCommonRownames(objList, rownameCol = rownameCol)
    mat_list <- .getMatList(objList, 
                            commonRownames = common_rownames,
                            rownameCol = rownameCol)
    summarized_mat <- .getSummarizedMat(matList = mat_list,
                                        sumMethod = summaryMethod)
    dist_mat <- as.matrix(dist(summarized_mat,
                               method = distMethod))
    return(dist_mat)
}

#' Plot PPI ROC curve
#' @param tpcaObj tpcaResult object
#' @param computeAUC logical parameter indicating
#' whether area under the ROC should be computed
#' and indicated in the lower right corner of the 
#' plot
#' 
#' @return ggplot object of a receiver operating
#' curve (ROC)
#' 
#' @export
#' 
#' @import ggplot2
#' @importFrom pROC auc
#' 
#' @examples 
#' rocTab = data.frame(
#'     TPR = c(0, 0.1, 0.2, 0.4, 0.5, 0.7, 0.9, 1),
#'     FPR = c(0, 0.05, 0.1, 0.2, 0.5, 0.7, 0.9, 1)
#' )
#' 
#' tpcaTest <- new(
#'     "tpcaResult",
#'     PPiRocTable = rocTab)
#' 
#' plotPPiRoc(tpcaTest)
#' 
plotPPiRoc <- function(tpcaObj, computeAUC = FALSE){
    rocTab <- tpcaObj@PPiRocTable
    if(computeAUC){
        rocTabAnno <- tpcaObj@PPiRocTableAnno
        aucText <- paste0("AUC = ", as.character(
            round(suppressMessages(
                auc(rocTabAnno$annotated ~
                        rocTab$eucl_dist)[1]), 3)))
    }
    p <- ggplot(rocTab, 
           aes(FPR, TPR)) + 
        geom_path() +
        geom_line(aes(x, y),
                  color = "gray",
                  data = tibble(x = 0:1, y = 0:1),
                  linetype = "dashed") +
        theme_bw()
    if(computeAUC){
        p <- p + geom_text(aes(x, y, label = label), 
                           data = tibble(x = 0.75, y = 0.1,
                                         label = aucText))
    }
    return(p)
        
}

.testForPPiCoAggregation <- function(tpcaObj, nSamp = 10000,
                                     p_adj_method = "BH"){
    tpcaObj@ComplexAnnotation <- 
        .adaptComplexAnnoFromPPi(
            tpcaObj@PPiAnnotation)
    
    tpcaObj <-.testForComplexCoAggregation(
        tpcaObj = tpcaObj, 
        minCount = 2, 
        nSamp = nSamp,
        p_adj_method = p_adj_method
    )
}

.adaptComplexAnnoFromPPi <- function(ppiAnno){
    complex_anno <- ppiAnno %>% 
        dplyr::select(pair, x, y) %>% 
        gather(key, value, -pair) %>% 
        arrange(pair) %>% 
        dplyr::select(id = pair, protein = value)
    return(complex_anno)
}

.testForComplexCoAggregation <- function(tpcaObj, 
                                        minCount = 3,
                                        nSamp = 10000,
                                        p_adj_method = "BH"){
    tpcaObj <- .intersectComplexAnnotation(
        tpcaObj, minCount = minCount)
    tpcaObj <- .setBackgroundDistribution(
        tpcaObj, nSamp = nSamp
    )
    tpcaObj <- .computeTpcaResultTable(
        tpcaObj, p_adj_method = p_adj_method
    )
}

.intersectComplexAnnotation <- function(tpcaObj, minCount = 3){
    caDf <- tpcaObj@ComplexAnnotation
    caDfFil <- caDf %>% 
        distinct() %>% 
        filter(
            protein %in% 
            tpcaObj@CommonFeatures) %>% 
        group_by(id) %>% 
        mutate(count = n()) %>% 
        ungroup %>% 
        filter(count >= minCount)
    tpcaObj@ComplexAnnotation <- caDfFil
    return(tpcaObj)
}

.setBackgroundDistribution <- function(tpcaObj, nSamp = 10000){
    uniqueMembersN <- unique(tpcaObj@ComplexAnnotation$count)
    tpcaObj@ComplexBackgroundDistributionList <- 
        .createBackgroundDistList(distMat = tpcaObj@DistMat,
                                  nMemVec = uniqueMembersN)
    return(tpcaObj)
}

.createBackgroundDistList <- function(distMat, nMemVec, nSamp = 10000){
    backgList <- lapply(nMemVec, function(x) {
        .sampleBackgroundDistribution(distMat, nMem = x, nSamp = nSamp)
        })
    names(backgList) <- nMemVec
    return(backgList)
}

.sampleBackgroundDistribution <- function(distMat, nMem = 3, nSamp = 10000){
    nCol <- ncol(distMat)
    samplesOfN <- vapply(seq_len(nSamp), function(i){
        ids <- sample(seq_len(nCol), size = nMem, replace = FALSE)
        mean(distMat[ids, ids][upper.tri(matrix(0, ncol = nMem, nrow = nMem))])
    }, FUN.VALUE = 1.0)
    return(samplesOfN)
}

.computeTpcaResultTable <- function(tpcaObj, p_adj_method = "BH"){
    tpca_tab <- .getMeanDistVals4Complexes(
        complexAnno = tpcaObj@ComplexAnnotation,
        distMat = tpcaObj@DistMat
    )
    tpcaObj@tpcaResultTable <- .computeDistPValue(
        tpcaDf = tpca_tab,
        backgList = tpcaObj@ComplexBackgroundDistributionList,
        p_adj_method = p_adj_method
    )
    return(tpcaObj)
}

#' @import dplyr
.getMeanDistVals4Complexes <- function(complexAnno, distMat){
    unique_complexes <- unique(complexAnno$id)
    complex_table <- bind_rows(lapply(unique_complexes, function(cm_id){
        x <- filter(complexAnno, id == cm_id)$protein
        x_dim <- length(x)
        ids <- rownames(distMat) %in% x
        mean_dist <- mean(distMat[ids, ids][
            upper.tri(matrix(0, ncol = x_dim, nrow = x_dim))])
        out_tab <- tibble(
            complex_name = cm_id,
            count = filter(complexAnno, id == cm_id)$count[1],
            mean_dist = mean_dist
        )
        return(out_tab)
    }))
    return(complex_table)
}

#' @import dplyr
#' @importFrom stats p.adjust
#' @importFrom fdrtool fdrtool
#' @importFrom stats na.omit
.computeDistPValue <- function(tpcaDf, backgList, 
                               p_adj_method = "BH"){
    tpcaDf <- tpcaDf %>% 
        na.omit() %>% 
        rowwise() %>% 
        mutate(p_value = length(which(backgList[[as.character(count)]] <= 
                                          mean_dist))/
                   length(backgList[[as.character(count)]])) %>% 
        ungroup %>% 
        within(p_value[p_value == 0] <- .Machine$double.eps)
    if(p_adj_method %in% c("BH", "bonferroni")){
        tpcaDf <- tpcaDf %>% 
            mutate(p_adj = p.adjust(p_value, method = p_adj_method))
    }else if(p_adj_method %in% c("fdrtool")){
        tpcaDf <- tpcaDf %>% 
            mutate(l_fdr = fdrtool(tpcaDf$p_value, 
                                   statistic = "pvalue")$lfdr)
    }
    
    return(tpcaDf)
}

#' Plot TPCA analysis results
#' 
#' @param tpcaObj a tpcaObj after having performed
#' a differential analysis, see \code{runDiffTPCA}
#' @param alpha significance level / FDR at which 
#' null hypothesis should be rejected
#' 
#' @return ggplot displaying a volcano plot
#' 
#' @import ggplot2
#' @import dplyr
#' @export
#' @examples 
#' library(dplyr)
#' library(Biobase)
#' 
#' m1 <- matrix(1:28, ncol = 4)
#' m2 <- matrix(2:25, ncol = 4)
#' m3 <- matrix(c(2:10, 1:19), ncol = 4)
#' 
#' rownames(m1) <- 1:7
#' rownames(m2) <- 3:8
#' rownames(m3) <- 2:8
#' 
#' mat_list <- list(
#'     m1, m2, m3
#' )
#' 
#' complex_anno <- tibble(
#'     protein = c("3", "4", "5", 
#'        "4", "5", "6", "7"),
#'     id = c(rep("1", 3), rep("2", 4)),
#'     count = c(rep(3, 3), rep(4, 4)))
#' 
#' tpca_result <- runTPCA(
#'   mat_list, complexAnno = complex_anno)
#' 
#' plotTpcaVolcano(tpca_result)
#' 
plotTpcaVolcano <- function(tpcaObj, alpha = 0.1){
    plot_df <- tpcaObj@tpcaResultTable
    
    p <- ggplot(plot_df, aes(x = mean_dist, -log10(p_value))) + 
        geom_point(color = "gray", alpha = 0.75) + 
        geom_point(data = filter(plot_df, p_adj < alpha)) + 
        theme_bw() +
        labs(x = "mean distance",
             y = expression('-log'[10]*'('*italic(p)*' value)'))
    
    return(p)
}



.createDistMatTpcaObj <- function(tpcaObj, rownameCol = NULL,
                                  summaryMethod = "median",
                                  distMethod = "euclidean",
                                  cond = "c1"){

    if(length(tpcaObj@ContrastList) == 0){
        tpcaObj@CommonFeatures <- .getCommonRownames(
            objList = tpcaObj@ObjList,
            rownameCol = rownameCol
        )  
    }else{
        tpcaObj@CommonFeatures <- .getCommonRownames(
            objList = c(tpcaObj@ObjList, tpcaObj@ContrastList),
            rownameCol = rownameCol
        )  
    }
    
    distMat <- createDistMat(
        objList = tpcaObj@ObjList, 
        rownameCol = rownameCol, 
        summaryMethod = summaryMethod,
        distMethod = distMethod
    )
    if(length(tpcaObj@ContrastList) != 0){
        contrastDistMat <- createDistMat(
            objList = tpcaObj@ContrastList, 
            rownameCol = rownameCol, 
            summaryMethod = summaryMethod,
            distMethod = distMethod
        )
    }
    # chunkDim <- ifelse(ncol(distMat) < 250, 
    #                    ncol(distMat), 250)
    tpcaObj@summaryMethod <- summaryMethod
    tpcaObj@distMethod <- distMethod
    tpcaObj@DistMat <- distMat 
    dimnames(tpcaObj@DistMat) <- dimnames(distMat)
    if(length(tpcaObj@ContrastList) != 0){
        tpcaObj@ContrastDistMat <- contrastDistMat
        dimnames(tpcaObj@ContrastDistMat) <- dimnames(contrastDistMat)
    }
    
    # writeHDF5Array(
    #     distMat,
    #     chunkdim = c(chunkDim, chunkDim)
    # )
    
    return(tpcaObj)
}

#' @importFrom methods new
.checkAnnoArguments <- function(objList, 
                                contrastList = list(),
                                ctrlCondName = "control",
                                contrastCondName = "treatment",
                                complexAnno = NULL,
                                ppiAnno = NULL){
    if(is.null(complexAnno) & is.null(ppiAnno)){
        stop(paste(c("Neither complex annotation nor a PPI", 
                     "annotation were supplied!\n",
                     "One of them is required."),
                   collapse = " "))
    }else if(!is.null(complexAnno) & !is.null(ppiAnno)){
        stop(paste(c("Both complex annotationand PPI", 
                     "annotation were supplied!\n",
                     "Please only supply one of them!"),
                   collapse = " "))
        # tpcaObj <- new("tpcaResult",
        #                ObjList = objList,
        #                ComplexAnnotation = complexAnno)
        # return(tpcaObj)
    }else if(!is.null(complexAnno)){
        tpcaObj <- new("tpcaResult",
                       ObjList = objList,
                       ContrastList = contrastList,
                       CtrlCondName = ctrlCondName,
                       ContrastCondName = contrastCondName,
                       ComplexAnnotation = complexAnno)
        return(tpcaObj)
    }else if(!is.null(ppiAnno)){
        tpcaObj <- new("tpcaResult",
                       ObjList = objList,
                       ContrastList = contrastList,
                       CtrlCondName = ctrlCondName,
                       ContrastCondName = contrastCondName,
                       PPiAnnotation = ppiAnno)
        return(tpcaObj)
    }
}

#' @import dplyr
#' @import tidyr
.createPPiRocTable <- function(tpcaObj){
    colname <- value <- rowname <- pair <- 
        annotated <- NULL
    
    distDf <- tpcaObj@DistMat %>% 
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
                   tpcaObj@PPiAnnotation$pair) %>% 
        mutate(TPR = cumsum(as.numeric(annotated)) / 
                    sum(as.numeric(annotated)), 
               FPR = cumsum(as.numeric(annotated == FALSE)) / 
                   (n() - cumsum(as.numeric(value != 0)) +  
                        cumsum(as.numeric(annotated == FALSE))))
    # chunDimsMat <- c(ifelse(nrow(distDf) < 500,
    #                         nrow(distDf), 500), 3)
    # chunDimsAnno <- c(ifelse(nrow(distDf) < 500, 
    #                          nrow(distDf), 500), 2)
    # rocMat <- as.matrix(
    #     distDf %>% dplyr::select(TPR, FPR, value))
    # tpcaObj@PPiRocTable <- writeHDF5Array(
    #     x = rocMat, chunkdim = chunDimsMat)
    tpcaObj@PPiRocTable <- distDf %>% 
        dplyr::select(TPR, FPR, eucl_dist = value)
    # colnames(tpcaObj@PPiRocTable) <- 
    #     c("TPR", "FPR", "eucl_dist")
    # annoDf <- as.data.frame(
    #     dplyr::select(distDf, pair, annotated))
    # tpcaObj@PPiRocTableAnno <- writeHDF5Array(
    #     x = annoDf, chunkdim = chunDimsAnno)
    tpcaObj@PPiRocTableAnno <- distDf %>% 
        dplyr::select(pair, annotated)
    # colnames(tpcaObj@PPiRocTableAnno) <- 
    #     c("pair", "annotated")
    return(tpcaObj)
}


#' @import dplyr
#' @import tidyr
.createComplexRocTable <- function(tpcaObj, 
                                   nPermuts = 5){
    perm_tpca_tab_list <- .getMeanDistVals4RandComplexes(
        tpcaObj, nPermuts = nPermuts)
    
    tpcaObj@ComplexRocTable <- .computeComplexRocTable(
        tpcaObj = tpcaObj,
        tpca_tab = tpcaObj@tpcaResultTable, 
        perm_tpca_tab_list = perm_tpca_tab_list
    )
    
    return(tpcaObj)
}

#' @import dplyr
.permuteComplexAnno <- function(complexAnno){
    permComplexAnno <- complexAnno %>% 
        mutate(protein = sample(protein)) %>% 
        group_by(id) %>% 
        filter(!duplicated(protein)) %>% 
        ungroup
    return(permComplexAnno)
}

.getMeanDistVals4RandComplexes <- function(tpcaObj, nPermuts = 5){
    perm_tpca_tab_list <- lapply(seq_len(nPermuts), function(np){
        perm_complex_anno <- 
            .permuteComplexAnno(tpcaObj@ComplexAnnotation)
        perm_tpca_tab <- .getMeanDistVals4Complexes(
            complexAnno = perm_complex_anno,
            distMat = tpcaObj@DistMat
        )
        return(perm_tpca_tab)
        
    })
    return(perm_tpca_tab_list)
}


.computeComplexRocTable <- function(tpcaObj, tpca_tab, perm_tpca_tab_list){
    complex_roc_df <- bind_rows(lapply(seq_len(length(perm_tpca_tab_list)), function(i){
                                           
        perm_tpca_tab_i <-  .computeDistPValue(
            tpcaDf = perm_tpca_tab_list[[i]],
            backgList = tpcaObj@ComplexBackgroundDistributionList,
            p_adj_method = "BH"
        )         
                                           
        roc_df <- bind_rows(
            tpca_tab %>% 
                mutate(category = "TP"),
            perm_tpca_tab_i %>% 
                mutate(category = "FP")) %>% 
            arrange(p_value, mean_dist) %>% 
            mutate(rank = seq_len(n()),
                   tp_count = cumsum(as.numeric(category == "TP")),
                   fp_count = cumsum(as.numeric(category == "FP"))) %>% 
            mutate(TPR = tp_count / nrow(tpca_tab),
                   FPR = (fp_count/nrow(perm_tpca_tab_i))) 
    })) %>% 
        group_by(rank) %>% 
        summarize(TPR = mean(TPR, na.rm = TRUE),
                  FPR = mean(FPR, na.rm = TRUE))
    return(complex_roc_df)
}

#' Plot Complex ROC curve
#' @param tpcaObj tpcaResult object
#' @param computeAUC logical parameter indicating
#' whether area under the ROC should be computed
#' and indicated in the lower right corner of the 
#' plot
#' 
#' @return ggplot object of a receiver operating
#' curve (ROC)
#' 
#' @export
#' 
#' @import ggplot2
#' @importFrom pROC auc
#' 
#' @examples 
#' rocTab = data.frame(
#'     TPR = c(0, 0.1, 0.2, 0.4, 0.5, 0.7, 0.9, 1),
#'     FPR = c(0, 0.05, 0.1, 0.2, 0.5, 0.7, 0.9, 1)
#' )
#' 
#' tpcaTest <- new(
#'     "tpcaResult",
#'     ComplexRocTable = rocTab)
#' 
#' plotComplexRoc(tpcaTest)
#' 
plotComplexRoc <- function(tpcaObj, computeAUC = FALSE){
    rocTab <- tpcaObj@ComplexRocTable
    if(computeAUC){
        aucText <- paste0("AUC = ", as.character(
            round(.simpleAuc(rocTab$FPR, rocTab$TPR), 3)))
    }
    p <- ggplot(rocTab, 
                aes(FPR, TPR)) + 
        geom_path() +
        geom_line(aes(x, y),
                  color = "gray",
                  data = tibble(x = 0:1, y = 0:1),
                  linetype = "dashed") +
        theme_bw() 
    if(computeAUC){
        p <- p + geom_text(aes(x, y, label = label), 
                           data = tibble(x = 0.75, y = 0.1,
                                         label = aucText))
    }
    return(p)
    
}

#' @importFrom utils head
#' @importFrom utils tail
.simpleAuc <- function(x, y){
    sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
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
    dim_list <- lapply(matList, dim)
    if(!all(vapply(dim_list, identical, dim_list[[1]], 
            FUN.VALUE = TRUE))){
        stop("checkMatDims: unequal matrix dimensions!")
    }
}

#' @importFrom stats median
.getSummarizedMat <- function(matList, sumMethod = "median"){
    .checkMatDims(matList)
    if(sumMethod == "median"){
        arr <- array(unlist(matList), 
                     c(dim(matList[[1]]), length(matList)))
        sumMat <- apply(arr, seq_len(2), median)
        rownames(sumMat) <- rownames(matList[[1]])
        colnames(sumMat) <- colnames(matList[[1]])
    }
    return(sumMat)
}

.filterDistMat <- function(dist_mat, ppi_anno){
    mat_nrow <- length(which(rownames(dist_mat) %in% ppi_anno$x))
    dist_row_names_in_ppi_anno <- rownames(dist_mat) %in% ppi_anno$x
    dist_col_names_in_ppi_anno <- colnames(dist_mat) %in% ppi_anno$y
    filt_mat <- dist_mat[dist_row_names_in_ppi_anno,
                         dist_col_names_in_ppi_anno]
    filt_matrix <- matrix(filt_mat, nrow = mat_nrow)
    rownames(filt_matrix) <- rownames(dist_mat)[
        dist_row_names_in_ppi_anno]
    colnames(filt_matrix) <- colnames(dist_mat)[
        dist_col_names_in_ppi_anno]
    
    return(filt_matrix)
}

#' @import dplyr
#' @import tidyr
.distMat2AnnotatedPPiDf <- function(dist_mat, ppi_anno){
    dist_mat %>% 
        tbl_df %>% 
        mutate(rowname = rownames(dist_mat)) %>% 
        gather(key, value, -rowname) %>% 
        rowwise %>% 
        mutate(pair = paste(sort(c(rowname, key)), 
                            collapse = ":")) %>% 
        ungroup %>% 
        filter(rowname != key, !duplicated(pair), 
               pair %in% ppi_anno$pair)
}

#' @import dplyr
.combineCondDistMatsFstat <- function(dist_df_c1, 
                                      dist_df_c2){
    combo_df <- left_join(
        dist_df_c1 %>% 
            dplyr::select(pair, valueC1 = value), 
        dist_df_c2 %>% 
            dplyr::select(pair, valueC2 = value), 
        by = "pair") %>% 
        mutate(rssC1 = valueC1^2,
               rssC2 = valueC2^2) %>% 
        mutate(rssC1_rssC2 = rssC1 - rssC2) %>% 
        rowwise() %>% 
        mutate(min_rssC1_rssC2 = min(c(rssC1, rssC2))) %>% 
        ungroup %>% 
        mutate(f_stat = abs(rssC1_rssC2)/min_rssC1_rssC2)
    
    return(combo_df)
}

.getRandomProteinPairs <- function(common_rownames, n = 10000){
    random_prots <- 
        head(tibble(
            x = sample(common_rownames, n + 100, replace = TRUE),
            y = sample(common_rownames, n + 100, replace = TRUE)
        ) %>% 
            filter(x!=y), n = n) %>% 
        rowwise() %>% 
        mutate(pair = paste(sort(c(x, y)), collapse = ":")) %>% 
        ungroup 
    
    return(random_prots)
}

.compareConditions <- function(tpcaObj, prot_pairs){
    distMatC1 <- .filterDistMat(
        tpcaObj@DistMat,
        prot_pairs)
    
    distMatDfC1 <- .distMat2AnnotatedPPiDf(
        dist_mat = distMatC1,
        ppi_anno = prot_pairs
    )
    
    distMatC2 <- .filterDistMat(
        tpcaObj@ContrastDistMat,
        prot_pairs)
    
    distMatDfC2 <- .distMat2AnnotatedPPiDf(
        dist_mat = distMatC2,
        ppi_anno = prot_pairs
    )
    
    comboDf <- .combineCondDistMatsFstat(
        distMatDfC1,
        distMatDfC2
    )
    return(comboDf)
}

.compareConditionsRandom <- function(tpcaObj, n = 10000){
    commonRownames <- .getCommonRownames(
        c(tpcaObj@ObjList, tpcaObj@ContrastList)
    )
    
    random_prot_pairs <- .getRandomProteinPairs(
        common_rownames = commonRownames,
        n = n
    )
    
    comboRandDf <- .compareConditions(
        tpcaObj,
        prot_pairs = random_prot_pairs)
    
    return(comboRandDf)
}

#' @importFrom stats p.adjust
#' @importFrom fdrtool fdrtool
#' @importFrom stats na.omit
.computeEmpiricalPValue <- function(
    combo_df, combo_rand_df, p_adj_method = "BH"){
    
    empirical_p_df <- combo_df %>% 
        na.omit() %>% 
        rowwise() %>% 
        mutate(p_value = length(which(combo_rand_df$f_stat >= f_stat))/
                   nrow(combo_rand_df)) %>% 
        ungroup %>% 
        within(p_value[p_value == 0] <- .Machine$double.eps)
    if(p_adj_method %in% c("BH", "bonferroni")){
        empirical_p_df <- empirical_p_df %>% 
            mutate(p_adj = p.adjust(p_value, method = p_adj_method))
    }else if(p_adj_method %in% "fdrtool"){
        empirical_p_df <- empirical_p_df %>% 
            mutate(l_fdr = fdrtool(empirical_p_df$p_value, 
                                   statistic = "pvalue")$lfdr)
    }
    
    return(empirical_p_df)
}

#' Run differential TPCA analysis
#' 
#' @param objList input list of objects, e.g.
#' ExpressionSets retrieved after TPP data import
#' or matrices or data frames
#' @param contrastList input list of objects for
#' comparison at e.g. different treatment condtion,
#' same file formats work as for objList
#' @param ctrlCondName character string indicating the name
#' of the control condition, default is "control"
#' @param contrastCondName character string indicating the name
#' of the contrast condition, default is "treatment"
#' @param ppiAnno data frame annotation known
#' protein-protein interactions (PPI) to test 
#' @param complexAnno data frame annotating
#' known protein complexes of interest to test 
#' @param rownameCol in case the input objects are tibbles
#' this parameter takes in the name (character) of the column 
#' specifying protein names or ids
#' @param summaryMethod character string indicating a method 
#' to use to summarize measurements across replicates, 
#' default is "median", other options are c("mean", "rbind")
#' @param distMethod method to use within dist function,
#' default is 'euclidean'
#' @param n number of random protein pair draws to obtain
#' empirical p-value, default is 10000
#' @param p_adj_method method to be used for multiple
#' testing adjusment, default is "BH"
#' 
#' @return an object of class tpcaResult
#' with the following slots:
#' 1) ObjList: containing the supplied list of
#' objects
#' 
#' 
#' @export
#' @examples 
#' library(dplyr)
#' library(Biobase)
#' 
#' m1 <- matrix(1:28, ncol = 4)
#' m2 <- matrix(2:25, ncol = 4)
#' m3 <- matrix(c(2:10, 1:19), ncol = 4)
#' 
#' rownames(m1) <- 1:7
#' rownames(m2) <- 3:8
#' rownames(m3) <- 2:8
#' 
#' mat_list <- list(
#'     m1, m2, m3
#' )
#' 
#' c1 <- matrix(29:2, ncol = 4)
#' c2 <- matrix(26:3, ncol = 4)
#' c3 <- matrix(c(11:3, 20:2), ncol = 4)
#' 
#' rownames(c1) <- 1:7
#' rownames(c2) <- 3:8
#' rownames(c3) <- 2:8
#' 
#' contrast_list <- list(
#'     c1, c2, c3
#' )
#' 
#' ppi_anno <- tibble(
#'     x = c("3", "3"),
#'     y = c("5", "7"),
#'     pair = c("3:5", "3:7"))
#' 
#' ref_df <- tibble(
#'     pair = c("3:5", "3:7"),
#'     valueC2 = c(4, 8)
#' )
#' 
#' diff_tpca <- Rtpca:::runDiffTPCA(
#'   mat_list, contrast_list, ppiAnno = ppi_anno)
#' 
runDiffTPCA <- function(objList, 
                        contrastList,
                        ctrlCondName = "control",
                        contrastCondName = "treatment",
                        ppiAnno = NULL,
                        complexAnno = NULL,
                        rownameCol = NULL,
                        summaryMethod = "median",
                        distMethod = "euclidean",
                        n = 10000,
                        p_adj_method = "BH"){
    message("Checking input arguments. \n")
    tpcaObj <- .checkAnnoArguments(
        objList = objList,
        contrastList = contrastList,
        ctrlCondName = ctrlCondName,
        contrastCondName = contrastCondName,
        complexAnno = complexAnno,
        ppiAnno = ppiAnno
    )
    message("Creating distance matrices. \n")
    tpcaObj <- .createDistMatTpcaObj(
        tpcaObj = tpcaObj, 
        rownameCol = rownameCol, 
        summaryMethod = summaryMethod,
        distMethod = distMethod
    )
    if(!is.null(ppiAnno)){
        message("Comparing annotated protein-pairs across conditions. \n")
        combo_df <- .compareConditions(
            tpcaObj = tpcaObj, 
            prot_pairs = ppiAnno
        )
        message("Comparing random protein-pairs across conditions. \n")
        combo_rand_df <- .compareConditionsRandom(
            tpcaObj = tpcaObj, 
            n = n
        )
    }else if(!is.null(complexAnno)){
        NULL
    }
    message("Generating result table. \n")
    tpcaObj@diffTpcaResultTable <- .computeEmpiricalPValue(
        combo_df, 
        combo_rand_df, 
        p_adj_method = p_adj_method)
    return(tpcaObj)
}

#' Plot differential TPCA analysis results
#' 
#' @param tpcaObj a tpcaObj after having performed
#' a differential analysis, see \code{runDiffTPCA}
#' @param setXLim logical determining whether
#' x-axis limits should be set according to xlimit
#' @param xlimit numeric vector with two elements
#' determing the x-axis limits, only is implemented
#' if setXLim is set to TRUE
#' @param alpha significance level / FDR at which 
#' null hypothesis should be rejected
#' 
#' @return ggplot displaying a volcano plot
#' 
#' @import ggplot2
#' @import dplyr
#' @export
#' @examples 
#' library(dplyr)
#' library(Biobase)
#' 
#' m1 <- matrix(1:28, ncol = 4)
#' m2 <- matrix(2:25, ncol = 4)
#' m3 <- matrix(c(2:10, 1:19), ncol = 4)
#' 
#' rownames(m1) <- 1:7
#' rownames(m2) <- 3:8
#' rownames(m3) <- 2:8
#' 
#' mat_list <- list(
#'     m1, m2, m3
#' )
#' 
#' c1 <- matrix(29:2, ncol = 4)
#' c2 <- matrix(26:3, ncol = 4)
#' c3 <- matrix(c(11:3, 20:2), ncol = 4)
#' 
#' rownames(c1) <- 1:7
#' rownames(c2) <- 3:8
#' rownames(c3) <- 2:8
#' 
#' contrast_list <- list(
#'     c1, c2, c3
#' )
#' 
#' ppi_anno <- tibble(
#'     x = c("3", "3"),
#'     y = c("5", "7"),
#'     pair = c("3:5", "3:7"))
#' 
#' ref_df <- tibble(
#'     pair = c("3:5", "3:7"),
#'     valueC2 = c(4, 8)
#' )
#' 
#' diff_tpca <- runDiffTPCA(
#'   mat_list, contrast_list, ppiAnno = ppi_anno)
#' 
#' plotDiffTpcaVolcano(diff_tpca)
#' 
plotDiffTpcaVolcano <- function(tpcaObj, 
                                alpha = 0.1,
                                setXLim  = FALSE, 
                                xlimit = c(-0.75, 0.75)){
    plot_df <- tpcaObj@diffTpcaResultTable
    controlCond <- tpcaObj@CtrlCondName
    constrastCond <- tpcaObj@ContrastCondName
    p <- ggplot(plot_df, aes(x = sqrt(valueC1) - sqrt(valueC2), -log10(p_value))) + 
        geom_point(color = "gray", alpha = 0.75) + 
        geom_point(data = filter(plot_df, p_adj < alpha)) + 
        theme_bw() +
        labs(x = bquote(sqrt(italic(d)[.(controlCond)])~ ' - ' 
                            ~sqrt(italic(d)[.(constrastCond)])),
             y = expression('-log'[10]*'('*italic(p)*' value)'))
    
    if(setXLim){
        p <- p + coord_cartesian(xlim = xlimit)
    }
    
    return(p)
}

#' Plot thermal profile of protein pairs
#' 
#' @param tpcaObj a tpcaObj after having performed
#' a differential analysis, see \code{runDiffTPCA}
#' @param pair character vector of one or more
#' protein names
#' @param splinesDf numeric, degree of freedom of 
#' the spline fit to the melting curves
#' 
#' @return ggplot displaying the thermal profile 
#' of a protein pair across conditions
#' 
#' @import ggplot2
#' @import dplyr
#' @importFrom splines ns
#' @export
#' @examples 
#' set.seed(12)
#' m1 <- matrix(rnorm(50), ncol = 10)
#' m2 <- matrix(rnorm(50), ncol = 10)
#' 
#' rownames(m1) <- letters[1:5]
#' rownames(m2) <- letters[1:5]
#' 
#' colnames(m1) <- paste("fc", 1:10, sep = "_")
#' colnames(m2) <- paste("fc", 1:10, sep = "_")
#' 
#' pheno <- data.frame(
#'     temperature = seq(37, 67, length.out = 10))
#' rownames(pheno) <- paste("fc", 1:10, sep = "_")
#' 
#' eset1 <- ExpressionSet(
#'     assayData = m1,
#'     phenoData = AnnotatedDataFrame(pheno)
#' )
#' 
#' eset2 <- ExpressionSet(
#'     assayData = m2,
#'     phenoData = AnnotatedDataFrame(pheno)
#' )
#' 
#' tpcaObj <- new("tpcaResult",
#'                ObjList = list(eset1),
#' 
#'ContrastList = list(eset2),
#'                 CtrlCondName = "control",
#'                 ContrastCondName = "treatment")
#' 
#' plotPPiProfiles(tpcaObj, pair = c("b", "d"))
#' 
plotPPiProfiles <- function(tpcaObj, pair, splinesDf = 4){
    plot_df <- .getDf4PlotProfiles(tpcaObj, pair)
    ggplot(plot_df, aes(temperature, rel_value)) +
        geom_point(aes(color = gene_name)) +
        geom_smooth(method = "lm", 
                    aes(group = gene_name, color = gene_name),
                    formula = y ~ splines::ns(x, df = splinesDf),
                    se = FALSE, size = 0.5) +
        facet_wrap(~condition) +
        scale_color_brewer("protein", palette = "Set1") +
        labs(y = "fraction non-denatured") +
        theme_bw()
        
}

.getDf4PlotProfiles <- function(tpcaObj, pair){
    if(class(tpcaObj@ObjList[[1]])[1] == "ExpressionSet"){
        full_cond_df <- bind_rows(lapply(tpcaObj@ObjList, function(eset){
            sub_mat <- exprs(eset)[rownames(eset) %in% pair,]
            sub_cond1_df <- .gatherSubMat(
                sub_mat = sub_mat,
                temperature_anno = eset$temperature
            )
        })) %>% mutate(condition = tpcaObj@CtrlCondName)
        if(length(tpcaObj@ContrastList) != 0){
            if(class(tpcaObj@ContrastList[[1]])[1] == "ExpressionSet"){
                cond2_df <- bind_rows(lapply(tpcaObj@ContrastList, function(eset){
                    sub_mat <- exprs(eset)[rownames(eset) %in% pair,]
                    sub_cond2_df <- .gatherSubMat(
                        sub_mat = sub_mat,
                        temperature_anno = eset$temperature
                    )
                })) %>% mutate(condition = tpcaObj@ContrastCondName)
                full_cond_df <- bind_rows(
                    full_cond_df, cond2_df
                )
                return(full_cond_df) 
            }else{
                stop("ObjList is of class 'ExpressionSet', 
                     while ContrastList ist not!\n")
            }
        }else{
           return(full_cond_df) 
        }
    }
    else if(class(tpcaObj@ObjList[[1]])[1] %in% c("matrix", "data.frame")){
        full_cond_df <- bind_rows(lapply(tpcaObj@ObjList, function(mat){
            sub_mat <- mat[rownames(mat) %in% pair,]
            sub_cond1_df <- .gatherSubMat(
                sub_mat = sub_mat,
                temperature_anno = attributes(mat)$temperature
            )
        })) %>% mutate(condition = tpcaObj@CtrlCondName)
        if(length(tpcaObj@ContrastList) != 0){
            if(class(tpcaObj@ContrastList[[1]])[1] %in% c("matrix", "data.frame")){
                cond2_df <- bind_rows(lapply(tpcaObj@ContrastList, function(mat){
                    sub_mat <- mat[rownames(mat) %in% pair,]
                    sub_cond2_df <- .gatherSubMat(
                        sub_mat = sub_mat,
                        temperature_anno = attributes(mat)$temperature
                    )
                })) %>% mutate(condition = tpcaObj@ContrastCondName)
                full_cond_df <- bind_rows(
                    full_cond_df, cond2_df
                )
                return(full_cond_df) 
            }else{
                stop("ObjList is of class 'matrix' or 'data.frame', 
                     while ContrastList ist not!\n")
            }
        }
    }
    # }else if((class(tpcaObj@ObjList[[1]])[1] %in% c("tbl_df")) &
    #          !is.null(rownameCol)){
    #     
    # }
    # cond1_df <- tpcaObj@objList
} 

.gatherSubMat <- function(sub_mat, temperature_anno = NULL){
    if(!is.null(temperature_anno)){
        colnames(sub_mat) <- as.character(temperature_anno)
    }
    cond_df <- sub_mat %>%
        tbl_df %>%
        mutate(gene_name = rownames(sub_mat)) %>%
        gather(temperature, rel_value, -gene_name) %>% 
        mutate(temperature = as.numeric(temperature))
    
    return(cond_df)
}

