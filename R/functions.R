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
                    distMethod = "euclidean"){
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
            tpcaObj = tpcaObj
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
#' @param computeAUC logical parameter indicating
#' whether areau under the ROC should be computed
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
                  linetype = "dashed")
    if(computeAUC){
        p <- p + geom_text(aes(x, y, label = label), 
                           data = tibble(x = 0.75, y = 0.1,
                                         label = aucText))
    }
    return(p)
        
}

# testForComplexCoAggregation <- function(tpcaObj){
#     
# }
# 
.intersectComplexAnnotation <- function(tpcaObj, minCount = 3){
    caDf <- tpcaObj@ComplexAnnotation
    caDfFil <- filter(
        caDf, protein %in% 
            tpcaObj@CommonFeatures) %>% 
        group_by(id) %>% 
        mutate(count = n()) %>% 
        ungroup %>% 
        filter(count >= minCount)
    tpcaObj@ComplexAnnotation <- caDfFil
    return(tpcaObj)
}

.setBackgroundDistribution <- function(tpcaObj, nSamp = 10000){
    uniqueMembersN <- unique(tpcaObj@ComplexAnnotation$n)
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
    samplesOfN <- sapply(seq_len(nSamp), function(i){
        ids <- sample(seq_len(nCol), size = nMem, replace = FALSE)
        mean(distMat[ids, ids][upper.tri(matrix(0, ncol = nMem, nrow = nMem))])
    })
    return(samplesOfN)
}

# plotComplexRoc <- function(tpcaObj){
#     
# }

.createDistMatTpcaObj <- function(tpcaObj, rownameCol = NULL,
                                  summaryFUN = median,
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
        summaryFUN = summaryFUN,
        distMethod = distMethod
    )
    if(length(tpcaObj@ContrastList) != 0){
        contrastDistMat <- createDistMat(
            objList = tpcaObj@ContrastList, 
            rownameCol = rownameCol, 
            summaryFUN = summaryFUN,
            distMethod = distMethod
        )
    }
    # chunkDim <- ifelse(ncol(distMat) < 250, 
    #                    ncol(distMat), 250)
    tpcaObj@summaryFUN <- summaryFUN
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
                                complexAnno = NULL,
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
                       ContrastList = contrastList,
                       ComplexAnnotation = complexAnno)
        return(tpcaObj)
    }else if(!is.null(ppiAnno)){
        tpcaObj <- new("tpcaResult",
                       ObjList = objList,
                       ContrastList = contrastList,
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
    if(!all(sapply(dim_list, identical, dim_list[[1]]))){
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

.filterDistMat <- function(dist_mat, ppi_anno){
    dist_mat[rownames(dist_mat) %in% ppi_anno$x,
             colnames(dist_mat) %in% ppi_anno$y]
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
        mutate(rssC1 = sqrt(valueC1),
               rssC2 = sqrt(valueC2)) %>% 
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

.computeEmpiricalPValue <- function(
    combo_df, combo_rand_df, p_adj_method = "BH"){
    
    empirical_p_df <- combo_df %>% 
        na.omit() %>% 
        rowwise() %>% 
        mutate(p_value = length(which(combo_rand_df$f_stat >= f_stat))/
                   nrow(combo_rand_df)) %>% 
        ungroup %>% 
        mutate(p_adj = p.adjust(p_value, method = p_adj_method))
    
    return(empirical_p_df)
}

runDiffTPCA <- function(objList, 
                        contrastList,
                        ppiAnno = NULL,
                        rownameCol = NULL,
                        summaryFUN = median,
                        distMethod = "euclidean",
                        n = 10000,
                        p_adj_method = "BH"){
    tpcaObj <- .checkAnnoArguments(
        objList = objList,
        contrastList = contrastList,
        complexAnno = NULL,
        ppiAnno = ppiAnno
    )
    tpcaObj <- .createDistMatTpcaObj(
        tpcaObj = tpcaObj, 
        rownameCol = rownameCol, 
        summaryFUN = summaryFUN,
        distMethod = distMethod
    )
    combo_df <- .compareConditions(
        tpcaObj = tpcaObj, 
        prot_pairs = ppiAnno
    )
    
    combo_rand_df <- .compareConditionsRandom(
        tpcaObj = tpcaObj, 
        n = n
    )
    
    tpcaObj@diffTpcaResultTable <- .computeEmpiricalPValue(
        combo_df, 
        combo_rand_df, 
        p_adj_method = p_adj_method)
    return(tpcaObj)
}

plotDiffTpcaVolcano <- function(tpcaObj){
    plot_df <- tpcaObj@diffTpcaResultTable
    
    p <- ggplot(plot_df, aes(x = rssC1_rssC2, -log10(p_value))) + 
        geom_point(color = "gray", alpha = 0.25) + 
        geom_point(data = filter(plot_df, p_adj < 0.1)) + 
        theme_bw() +
        labs(x = expression('RSS'^c1* ' - RSS'^c2*''),
             y = expression('-log'[10]*'('*italic(p)*' value)'))
}
