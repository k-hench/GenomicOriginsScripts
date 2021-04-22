#' hamlet pair annotation
#'
#' \code{fish_plot} provides hamlet pair annotation.
#'
#' @param spec string six letter population name
#'
#' @export
fish_plot <- function(spec){

  p <- ggplot()+
    hypoimg::hypo_anno_flag(geo = loc_names[stringr::str_sub(spec, 4, 6)] %>% stringr::str_to_lower(), ymax = 0)+
    hypoimg::hypo_anno_l(species = sp_names[stringr::str_sub(spec, 1, 3)], ymin = 0)+
    ggplot2::scale_y_continuous(limits = c(-.8,1))+
    ggplot2::theme_void()

  tibble::tibble(spec = spec, grob = list(p %>% ggplot2::ggplotGrob()))
}

#' hamlet pair annotation 2
#'
#' \code{fish_plot2} provides hamlet pair annotation with alternative axis specifications.
#'
#' @param spec string six letter population name
#'
#' @export
fish_plot2 <- function(spec){

  p <- ggplot()+
    hypoimg::hypo_anno_flag(geo = loc_names[stringr::str_sub(spec,4,6)] %>% stringr::str_to_lower(), xmax = 0)+
    hypoimg::hypo_anno_r(species = sp_names[stringr::str_sub(spec,1,3)], xmin = 0)+
    ggplot2::scale_y_continuous(limits = c(-1, 1))+
    ggplot2::scale_x_continuous(limits = c(-.4, 1))+
    ggplot2::theme_void()

  tibble::tibble(spec = spec, grob = list(p %>% ggplot2::ggplotGrob()))
}
