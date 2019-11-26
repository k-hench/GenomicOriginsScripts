#' Compute summary stats an absolute for Fst outlier
#'
#' \code{get_fst_percentile} computes summary stats for Fst outlier for a given Fst threshold.
#'
#' @family General functions
#'
#' @export
get_fst_fixed <- function(file, run, fst_threshold,...){

  data <- hypo_import_windows(file,...) %>%
    mutate(rank = rank(WEIGHTED_FST,ties.method = "random"))%>%
    mutate(thresh = fst_threshold) %>%
    mutate(outl = (WEIGHTED_FST>thresh) %>% as.numeric()) %>%
    filter(outl == 1 )

  if(nrow(data) == 0){
    return(tibble(run = run, n = 0, avg_length = NA, med_length = NA, min_length = NA, max_length = NA,
                  sd_length = NA, overal_length = NA, threshold_value = fst_threshold))
    } else {
  data %>%
    # next, we want to collapse overlapping windows
    group_by(CHROM) %>%
    # we check for overlap and create 'region' IDs
    mutate(check = 1-(lag(BIN_END,default = 0)>BIN_START),
           ID = str_c(CHROM,'_',cumsum(check))) %>%
    ungroup() %>%
    # then we collapse the regions by ID
    group_by(ID) %>%
    summarise(run = run[1],
              run = run[1],
              treshold_value = thresh[1],
              CHROM = CHROM[1],
              BIN_START = min(BIN_START),
              BIN_END = max(BIN_END)) %>%
    mutate(PEAK_SIZE = BIN_END-BIN_START) %>%
    summarize(run = run[1],
              run = run[1],
              n = length(ID),
              avg_length = mean(PEAK_SIZE),
              med_length = median(PEAK_SIZE),
              min_length = min(PEAK_SIZE),
              max_length = max(PEAK_SIZE),
              sd_length = sd(PEAK_SIZE),
              overal_length = sum(PEAK_SIZE),
              threshold_value = treshold_value[1])
    }
}

#' Collapse all Fst peaks of a given run
#'
#' \code{collapse_peaks} collapses all Fst peaks of a given run
#'
#' @family Figure SX
#'
#' @export
collapse_peaks <- function(x){
  table(x) %>%
    str_c(names(.),.,sep = ':',collapse = ', ')
  }

#' Reformat run name
#'
#' \code{reformat_run_name} formats run to "specloc-specloc" format.
#'
#' @family Figure SX
#'
#' @export
reformat_run_name <- function(run){
  run %>%
    str_remove(pattern = '^.*/') %>%
    str_replace(pattern = '([a-y]{3})-([a-y]{3})-([a-y]{3})','\\2\\1-\\3\\1')
  }
