#' @title Data frame of annotated protein complexes by Ori 
#' et al.
#' @name ori_et_al_complexes_df
#' @docType data
#' @description data frame assigning proteins to 
#' annotated protein complexes 
#' @references Ori et al. (2016), Genome Biology, 17, 47
#' @format data frame with columns ensembl_id, protein
#' and id (complex identifier)
#' @examples data("ori_et_al_complexes_df")
#' @usage data("ori_et_al_complexes_df")
"ori_et_al_complexes_df"

#' @title Data frame of eukaryotic protein-protein
#' interactions inferred from annotated protein 
#' complexes by Ori et al. and StringDB interations
#' with a combined score of at least 900
#' @name ori_et_al_complex_ppis
#' @docType data
#' @description data frame assigning proteins to (in)directly 
#' interacting proteins within protein complexes
#' @references Ori et al. (2016), Genome Biology, 17, 47;
#' Jensen et al. (2009), Nucleic Acids Research, 
#' 37, D412–D416
#' @format data frame with columns complex_name, x, y, 
#' pair (unique pair id)
#' @examples data("ori_et_al_complex_ppis")
#' @usage data("ori_et_al_complex_ppis")
"ori_et_al_complex_ppis"

#' @title Data frame of annotated human protein-protein
#' interactions retrieved from stringDB with a combined
#' interaciton score euqal or higher than 700
#' @name string_ppi_df
#' @docType data
#' @description data frame assigning proteins to 
#' interacting proteins 
#' @references Jensen et al. (2009), Nucleic Acids Research, 
#' 37, D412–D416
#' @format data frame with columns x, y (gene symbol of 
#' interactors), combined_score, pair (unique pair id)
#' @examples data("string_ppi_df")
#' @usage data("string_ppi_df")
"string_ppi_df"