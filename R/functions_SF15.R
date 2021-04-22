#' Import Fst permutation tests results
#'
#' @param file string, path to fst permutation test results
#'
#' @export
get_random_fst <- function(file){
  nm <- file %>%
    stringr::str_remove(pattern = ".*\\/") %>%
    stringr::str_remove("_random_fst.tsv.gz")
  rn <- stringr::str_remove(nm, "_.*")
  sub_type <- stringr::str_remove(nm, "[a-z]{6}-[a-z]{6}_")
  vroom::vroom(file = file,
               delim = "\t",
               col_types = "dcdd") %>%
    dplyr::mutate(group = nm,
           run = rn,
           subset_type = sub_type)
}

#' Compute the percentile of original Fst copared to fst distribution
#'
#' @param data tibble
#'
#' @export
get_percentile <- function(data){
  ran <- data$weighted_fst[data$type == "random"]
  real_fst <- data$weighted_fst[data$type == "real_pop"]
  sum(ran < real_fst) / length(ran)
}

#' Count number of permutation Fst greather than the original Fst
#'
#' @param data tibble
#'
#' @export
get_n_above <- function(data){
  ran <- data$weighted_fst[data$type == "random"]
  real_fst <- data$weighted_fst[data$type == "real_pop"]
  sum(ran > real_fst)
}

#' Count number of fst permutations
#'
#' @param data tibble
#'
#' @export
get_n_total <- function(data){
  ran <- data$weighted_fst[data$type == "random"]
  length(ran)
}
