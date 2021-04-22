#' Get gxp data
#'
#' \code{get_gxp_both_models} imports gxp data
#'
#' @param file       input file
#' @param trait      string, focal trait
#' @param model_type string ("lm" /"lmm")
#' @param path       path to folder contain the GxP results
#'
#' @export
get_gxp_both_models <- function(file, trait, model_type, path){
  vroom::vroom(stringr::str_c(path, file), delim = "\t") %>%
    dplyr::left_join(hypo_chrom_start) %>%
    dplyr::mutate(trait = trait,
                  model_type = model_type,
                  gpos = GSTART + MID_POS)
}
