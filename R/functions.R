# createEuclDistMat <- function(objList){
#     common_rownames <- 
#         getCommonRownames(objList)
#     exprs_list <- lapply(trData, function(x){
#             e_mat <- exprs(x)
#             e_mat[rownames(e_mat) %in% common_rownames,]
#         })
# }

#' @import Biobase 
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
            x_fil <- filter(x, eval(parse(text = rownameCol)) %in% commonRownames) 
            x_mat <- as.matrix(dplyr::select(x_fil, -one_of(rownameCol)))
            rownames(x_mat) <- x_fil[[rownameCol]]
            return(x_mat)
        }))
    }else{
        stop(paste("getMatList: Supplied object is neither", 
                   "ExpressionSet nor matrix or data.frame!\n"))
    }
}

# .getSummarizedMat <- function(matList, FUN = median){
#     
# }