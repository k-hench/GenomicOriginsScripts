#' Get gxp data
#'
#' \code{get_gxp_both_models} imports gxp data
#'
#' @family Suppl Figure 9
#'
#' @export
get_gxp_both_models <- function(file, trait, model_type, path){
  vroom::vroom(str_c(path,file), delim = "\t") %>%
    left_join(hypo_chrom_start)  %>%
    mutate(trait = trait,
           model_type = model_type,
           gpos = GSTART + MID_POS)
}
