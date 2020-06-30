#' crate global fst bar table
#'
#' \code{fst_bar_row_run} crates a rescaled global fst bar table.
#'
#' @family Suppl Figure 2
#'
#' @export
fst_bar_row_run <- function(fst,run){
  tibble(xmin = rescale_fst(0),
         xmax = rescale_fst(fst),
         xmin_org = 0,
         xmax_org = fst,
         ymin = 0,
         ymax= .15, run = run)
}

#' normalize global fst
#'
#' \code{rescale_fst} normalize global fst to 0 - 1.
#'
#' @family Suppl Figure 2
#'
#' @export
rescale_fst <- function(fst){
  start <- 0
  end <- 1
  fst_max <- max(globals$weighted)

  scales::rescale(fst,from = c(0,fst_max), to = c(start,end))
}

#' reorder species pairs
#'
#' \code{refactor_run} reorder species pairs according to fst.
#'
#' @family Suppl Figure 2
#'
#' @export
refactor_run <- function(self,globals){
  factor(as.character(self$run),
         levels = c(levels(globals$run)))
}

#' create hamlet pair annotations
#'
#' \code{plot_pair_run} prepares the hamlet pair annotations.
#'
#' @family Suppl Figure 2
#'
#' @export
plot_pair_run <- function(run,loc,left,right){
  tibble( run = run,loc = loc,
           grob = list(anno_pair_flag(loc,left,right) %>%
                         ggplotGrob()))

}

#' create hamlet pair annotations
#'
#' \code{plot_pair_run} creates the hamlet pair annotations.
#'
#' @family Suppl Figure 2
#'
#' @export
anno_pair_flag <- function (loc,left, right, circle_color = NA, circle_fill_left = "white",
                            circle_fill_right = "lightgray", circle_lwd = 0.5, plot_names = FALSE,
                            plot_name_size = 3, font_family = "sans", ...) {
  nr_left <- which(hypo_img$spec %>% str_sub(.,1,3) == left)
  nr_right <- which(hypo_img$spec %>% str_sub(.,1,3) == right)
  nr_flag <- which((hypo_flag$geo %>% str_sub(.,1,3)) == loc)
  p <- ggplot() +
    theme_void() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.4, 0.38)) +
    annotation_custom(hypo_flag$flag[[nr_flag]], xmin = -0.3, xmax = .3, ymin = -Inf, ymax = Inf) +
    annotation_custom(hypo_img$l[[nr_right]], xmin = 0.05, xmax = 1, ymin = -Inf, ymax = Inf) +
    annotation_custom(hypo_img$r[[nr_left]], xmin = -1, xmax = -0.05, ymin = -Inf, ymax = Inf)+
    coord_cartesian(xlim = c(-1.1,1.1))
  return(p)
}
