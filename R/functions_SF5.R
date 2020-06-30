#' hamlet pair annotation
#'
#' \code{fish_plot} provides hamlet pair annotation.
#'
#' @family Suppl. Figure SX
#'
#' @export
fish_plot <- function(spec){

  p <- ggplot()+
    hypoimg::hypo_anno_flag(geo = loc_names[str_sub(spec,4,6)] %>% str_to_lower(),ymax = 0)+
    hypoimg::hypo_anno_l(species = sp_names[str_sub(spec,1,3)],ymin = 0)+
    scale_y_continuous(limits = c(-.8,1))+
    theme_void()

  tibble(spec = spec, grob = list(p %>% ggplotGrob()))
}

#' hamlet pair annotation 2
#'
#' \code{fish_plot2} provides hamlet pair annotation with alternative axis specifications.
#'
#' @family Suppl. Figure SX
#'
#' @export
fish_plot2 <- function(spec){

  p <- ggplot()+
    hypoimg::hypo_anno_flag(geo = loc_names[str_sub(spec,4,6)] %>% str_to_lower(),xmax = 0)+
    hypoimg::hypo_anno_r(species = sp_names[str_sub(spec,1,3)],xmin = 0)+
    scale_y_continuous(limits = c(-1,1))+
    scale_x_continuous(limits = c(-.4,1))+
    theme_void()

  tibble(spec = spec, grob = list(p %>% ggplotGrob()))
}
