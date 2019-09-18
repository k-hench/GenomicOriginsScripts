#' Import a single topology weighting data set
#'
#' \code{get_twisst_data} a single topology weighting data set
#'
#' The topology weighting file containg the weights is read,
#' normalized and filtered for missing data.
#' Then the file containg the window positions is read,
#' filtered according to the weights file and its columns
#' are renamed.
#'
#' Optionally, the data is smoothed, then the weights and
#' the position information is merged, reformated and the
#' genomic positions are added.
#'
#' @family Figure 2
#'
#' @export
get_twisst_data <- function(loc,w_in,d_in,smooth = FALSE, span = 0.01){
  ########## read data ##################
  weights <- read.table(str_c(w_path, w_in), header = T)
  #normalise rows so weights sum to 1
  weights <- weights / apply(weights, 1, sum)
  #exclude any rows where data is missing
  good_rows <- which(is.na(apply(weights,1,sum)) == F)
  weights <- weights[good_rows,]

  #retrieve the names of the topologies
  window_data <- read.table(str_c(d_path, loc, '/', d_in), header = T)

  window_data <- window_data[good_rows,] %>%
    set_names(., nm = c('CHROM','BIN_START','BIN_END','BIN_MID','N_SITES','lnL'))

  if(smooth == TRUE){
    weights <- smooth.weights(window_positions = window_data$BIN_MID, weights_dataframe = weights,
                              span = span, window_sites = window_data$N_SITES)
  }

  weights  %>%
    mutate(rn = row_number()) %>%
    bind_cols(window_data) %>%
    gather(key = 'topo', value = 'weight', -rn:-lnL) %>%
    left_join(hypo_karyotype) %>%
    mutate(GPOS = BIN_MID + GSTART,
           topo_nr = str_remove(topo, 'topo') %>%
             as.numeric()) %>%
    as_tibble()
}

#' Import all topology weighting data sets for a location
#'
#' \code{match_twisst_files} imports all topology weighting data sets for a location.
#'
#' First the weight and position files of all topology weighting data
#' sets (the different linkage groups) for a given location are collected.
#'
#' The topology weighting data import funtion is iterated over all
#' data sets and combined to create a full data set.
#'
#' Then, the full data set is summarized to then be able to rank
#' the topologies based on their overall weighting.
#'
#' The topology ranks are added to the full data set and a window
#' column caontaing the location is added for faceting.
#'
#' @family Figure 2
#'
#' @export
match_twisst_files <- function(loc, window_size = 50,panel = 'a'){
  weight_files <- dir(w_path, pattern = loc) %>% .[grepl(str_c('w',window_size),.)]
  data_files <- dir(str_c(d_path, loc, '/'), pattern = 'LG.*data.tsv') %>% .[grepl(str_c('w',window_size),.)]

  data <- tibble(w_in = weight_files,
                 d_in = data_files) %>%
    purrr::pmap(get_twisst_data, smooth = FALSE,loc=loc) %>%
    bind_rows()

  topo_summary <- data %>%
    group_by(topo_nr) %>%
    summarise(mean = mean(weight),
              median = median(weight),
              sum = sum(weight),
              sd = sd(weight)) %>%
    mutate(rank = 106-rank(mean, ties.method = 'random')) %>%
    gather(key = 'stat', value = 'val', mean:sd)

  topo_rank <- topo_summary %>%
    filter(stat == 'mean') %>%
    select(topo_nr, rank)

  data <- data %>%
    left_join(topo_rank) %>%
    mutate(topo3 = str_pad(topo_nr, width = 3, pad = '0'),
           topo_rel = topo_nr/max(topo_nr))

  data <- data %>%
    mutate(window = str_c(project_case(panel),':~weighting[',loc,']'))
  return(data)
}
