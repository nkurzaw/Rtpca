# createEuclDistMat <- function(objList){
#     common_rownames <- 
#         getCommonRownames(objList)
#     exprs_list <- lapply(trData, function(x){
#             e_mat <- exprs(x)
#             e_mat[rownames(e_mat) %in% common_rownames,]
#         })
# }

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
