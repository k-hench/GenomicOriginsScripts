#' summarize model statistics
#'
#' \code{summarise_model} summarizes model statistics in a tibble row.
#'
#' @param data tibble with column 'mod', containing models to be summarized
#'
#' @family Suppl Figure 3
#'
#' @export
summarise_model <- function(data){
  data$mod %>%
    purrr::map_dfr(broom::glance) %>%
    dplyr::bind_cols(.,
              data$mod %>%
                purrr::map_dfr(broom::tidy) %>%
                dplyr::mutate(grp = (dplyr::row_number()+1)%/%2) %>%
                dplyr::select(grp, term, estimate) %>%
                tidyr::spread(key = 'term', value = 'estimate') %>%
                dplyr::select(-grp) %>%
                purrr::set_names(., nm = c('intercept', 'slope')))
}

#' Plot two fishes with location
#'
#' \code{plot_fishes_location} creates a ggplot annotation layer with two hamlets and a flag.
#'
#' This is used to annotate plots.
#'
#' @param left  string, three letter species identifier
#' @param right string, three letter species identifier
#' @param loc   string, sample location (bel [Belize]/ hon [Honduras]/ pan [Panama])
#'
#' @family Suppl Figure 3
#'
#' @export
plot_fishes_location <- function (left, right, loc) {
  nr_left <- which(str_sub(hypoimg::hypo_img$spec, start = 1, end = 3) == left)
  nr_right <- which(str_sub(hypoimg::hypo_img$spec, start = 1, end = 3) == right)
  nr_loc <- which(str_sub(hypoimg::hypo_flag$geo, start = 1, end = 3) == loc)

  p <- ggplot2::ggplot() +
    ggplot2::annotation_custom(hypoimg::hypo_flag$flag[[nr_loc]], xmin = -0.28, xmax = 0.28, ymin = -Inf, ymax = Inf)+
    ggplot2::coord_fixed(xlim = c(-1, 1)) +
    ggplot2::theme_void() +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::scale_y_continuous(limits = c(-0.4, 0.38)) +
    ggplot2::annotation_custom(hypoimg::hypo_img$r[[nr_right]], xmin = -1, xmax = -0.05, ymin = -Inf, ymax = Inf) +
    ggplot2::annotation_custom(hypoimg::hypo_img$l[[nr_left]], xmin = 0.05, xmax = 1, ymin = -Inf, ymax = Inf)

  tibble::tibble(run = str_c(left, loc, '-', right, loc),
         grob = list(p %>% ggplot2::ggplotGrob()))
}

#' The custom grob geom
#'
#' This geom builds heavily on the answer by baptiste on the
#' tidiverse github forum:
#' https://github.com/tidyverse/ggplot2/issues/1399
#'
#' @param mapping      Aesthetic mapping; ggplot2::aes()
#' @param data         data frame
#' @param stat         ggplot2 stat
#' @param position     string, "identity"
#' @param na.rm        logical
#' @param show.legend  logical
#' @param inherit.aes  logical
#' @param ...          arguments passed to ggplot2::layer()
#'
#' @family Suppl Figure 3
#'
#' @export
geom_hypo_grob2 <- function(mapping = NULL,
                            data = NULL,
                            stat = "identity",
                            position = "identity",
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = FALSE,
                            ...) {
  ggplot2::layer(
    geom = hypo_geom_grob_custom2,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

hypo_geom_grob_custom2 <- ggplot2::ggproto(
  "hypo_geom_grob_custom2",
  Geom,
  setup_data = function(self, data, params) {
    data <- ggplot2::ggproto_parent(Geom, self)$setup_data(data, params)
    data
  },

  draw_group = function(data, panel_scales, coord) {
    vp <- grid::viewport(x=data$rel_x,
                         y=data$rel_y,
                         h = data$height,
                         width = data$width,
                         angle = data$angle)
    g <- grid::editGrob(data$grob[[1]], vp=vp)
    ggplot2:::ggname("geom_hypo_grob2", g)
  },

  required_aes = c("grob","rel_x","rel_y"),
  default_aes = list(height = 1, width = 1, angle = 0)
)

#' Compute summary stats an absolute for Fst outlier
#'
#' \code{get_fst_percentile} computes summary stats for Fst outlier for a given Fst threshold.
#'
#' @param file          path to fst file
#' @param run           string, species pair
#' @param fst_threshold numeric, threshold for fst outliers
#' @param ...           arguments passed to hypogen::hypo_import_windows()
#'
#' @family  Suppl Figure 3
#'
#' @export
get_fst_fixed <- function(file, run, fst_threshold,...){

  data <- hypogen::hypo_import_windows(file, ...) %>%
    dplyr::mutate(rank = rank(WEIGHTED_FST, ties.method = "random"))%>%
    dplyr::mutate(thresh = fst_threshold) %>%
    dplyr::mutate(outl = (WEIGHTED_FST > thresh) %>% as.numeric()) %>%
    dplyr::filter(outl == 1 )

  if(nrow(data) == 0){
    return(tibble(run = run, n = 0, avg_length = NA, med_length = NA, min_length = NA, max_length = NA,
                  sd_length = NA, overal_length = NA, threshold_value = fst_threshold))
    } else {
      data %>%
        # next, we want to collapse overlapping windows
        dplyr::group_by(CHROM) %>%
        # we check for overlap and create 'region' IDs
        dplyr::mutate(check = 1-(lag(BIN_END,default = 0)>BIN_START),
                      ID = str_c(CHROM,'_',cumsum(check))) %>%
        dplyr::ungroup() %>%
        # then we collapse the regions by ID
        dplyr::group_by(ID) %>%
        dplyr::summarise(run = run[1],
                         run = run[1],
                         treshold_value = thresh[1],
                         CHROM = CHROM[1],
                         BIN_START = min(BIN_START),
                         BIN_END = max(BIN_END)) %>%
        dplyr::mutate(PEAK_SIZE = BIN_END-BIN_START) %>%
        dplyr::summarize(run = run[1],
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
#' @param x string vector with outlier peak ids
#'
#' @family  Suppl Figure 3
#'
#' @export
collapse_peaks <- function(x){
  table(x) %>%
    stringr::str_c(names(.),.,sep = ':', collapse = ', ')
  }

#' Reformat run name
#'
#' \code{reformat_run_name} formats run to "specloc-specloc" format.
#'
#' @param run string, species pair
#'
#' @family  Suppl Figure 3
#'
#' @export
reformat_run_name <- function(run){
  run %>%
    stringr::str_remove(pattern = '^.*/') %>%
    stringr::str_replace(pattern = '([a-y]{3})-([a-y]{3})-([a-y]{3})','\\2\\1-\\3\\1')
  }
