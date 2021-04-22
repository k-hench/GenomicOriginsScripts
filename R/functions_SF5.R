#' crate global fst bar table
#'
#' \code{fst_bar_row_run} crates a re-scaled global fst bar table.
#'
#' @param fst global fst value
#' @param run string, species pair
#'
#'
#' @export
fst_bar_row_run <- function(fst, run){
  tibble::tibble(xmin = rescale_fst(0),
         xmax = rescale_fst(fst),
         xmin_org = 0,
         xmax_org = fst,
         ymin = 0,
         ymax= .15,
         run = run)
}

#' normalize global fst
#'
#' \code{rescale_fst} normalize global fst to 0 - 1.
#'
#' @param fst numeric, global fst value
#'
#' @export
rescale_fst <- function(fst){
  start <- 0
  end <- 1
  fst_max <- max(globals$weighted)

  scales::rescale(fst, from = c(0, fst_max), to = c(start, end))
}

#' reorder species pairs
#'
#' \code{refactor_run} reorder species pairs according to fst.
#'
#' @param self    string vector, species pairs
#' @param globals reference table, order of species pairs according to reference table
#'
#' @export
refactor_run <- function(self, globals){
  factor(as.character(self$run),
         levels = c(levels(globals$run)))
}

#' create hamlet pair annotations
#'
#' \code{plot_pair_run} prepares the hamlet pair annotations.
#'
#' @param run   string, species pair
#' @param loc   string, sample location (bel [Belize]/ hon [Honduras]/ pan [Panama])
#' @param left  string, left species
#' @param right string, right species
#'
#' @export
plot_pair_run <- function(run,loc,left,right){
  tibble::tibble( run = run, loc = loc,
                  grob = list(anno_pair_flag(loc, left, right) %>%
                                ggplot2::ggplotGrob()))

}

#' create hamlet pair annotations
#'
#' \code{plot_pair_run} creates the hamlet pair annotations.
#'
#' @param loc                string, sample location (bel [Belize]/ hon [Honduras]/ pan [Panama])
#' @param left               string, left species
#' @param right              string, right species
#' @param circle_color       string, color
#' @param circle_fill_left   string, color
#' @param circle_fill_right  string, color
#' @param circle_lwd         numeric
#' @param plot_names         logical, should species labels be plotted
#' @param plot_name_size     numeric, size of labels
#' @param font_family        string, font family
#' @param ...                catch all for excess parameters for purrr::pmap
#'
#' @export
anno_pair_flag <- function (loc,left, right, circle_color = NA, circle_fill_left = "white",
                            circle_fill_right = "lightgray", circle_lwd = 0.5, plot_names = FALSE,
                            plot_name_size = 3, font_family = "sans", ...) {
  nr_left <- which(hypoimg::hypo_img$spec %>% str_sub(.,1,3) == left)
  nr_right <- which(hypoimg::hypo_img$spec %>% str_sub(.,1,3) == right)
  nr_flag <- which((hypoimg::hypo_flag$geo %>% str_sub(.,1,3)) == loc)
  p <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-0.4, 0.38)) +
    ggplot2::annotation_custom(hypoimg::hypo_flag$flag[[nr_flag]], xmin = -0.3, xmax = .3, ymin = -Inf, ymax = Inf) +
    ggplot2::annotation_custom(hypoimg::hypo_img$l[[nr_right]], xmin = 0.05, xmax = 1, ymin = -Inf, ymax = Inf) +
    ggplot2::annotation_custom(hypoimg::hypo_img$r[[nr_left]], xmin = -1, xmax = -0.05, ymin = -Inf, ymax = Inf)+
    ggplot2::coord_cartesian(xlim = c(-1.1,1.1))
  return(p)
}
