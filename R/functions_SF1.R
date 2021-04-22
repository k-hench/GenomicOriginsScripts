#' PCA plot without fish annotations
#'
#' @param loc string, three letter location abbreviation
#' @param mode string, type of pca data
#' @param pc1 integer, number of PC for x axis
#' @param pc2 integer, number of PC for y axis
#'
#' @export
pca_plot_no_fish <- function(loc, mode, pc1 = 1, pc2 = 2){
  if(loc == ""){ loc2 <- "all" } else { loc2 <- loc }
  titles <- c(loc_names,"all") %>% purrr::set_names(nm = c("bel.", "hon.", "pan.", "flo", "all"))
  clr_pca <- c(clr_loc, "black") %>% purrr::set_names(nm = c("bel.", "hon.", "pan.", "flo", "all"))

  set.seed(42)
  evs <- stringr::str_c(pca_dir, loc , mode, ".exp_var.txt.gz") %>%
    readr::read_tsv()
  stringr::str_c(pca_dir, loc , mode,".scores.txt.gz") %>%
    readr::read_tsv() %>%
    dplyr::mutate(spec = str_sub(id, -6,-4)) %>%
    ggplot2::ggplot(ggplot2::aes_string(x = stringr::str_c("EV0", pc1),
                                        y = stringr::str_c("EV0", pc2),
                                        fill = "spec"))+
    ggforce::geom_mark_ellipse(aes(color = spec),
                               fill = "transparent",
                               linetype = 3,
                               size = .3,
                               expand = unit(5, "pt"))+
    ggplot2::geom_point(shape = 21,
                        aes(color = ggplot2::after_scale(prismatic::clr_darken(fill))), size = .7) +
    ggplot2::labs(x = str_c("PC",pc1," (", sprintf("%.1f",evs$exp_var[[ pc1 ]]), " %)"),
         y = str_c("PC",pc2," (", sprintf("%.1f",evs$exp_var[[ pc2 ]]), " %)"))+
    ggplot2::scale_fill_manual(values = clr)+
    ggplot2::scale_color_manual(values = clr_alt %>%
                         prismatic::clr_alpha(alpha = .7) %>%
                         purrr::set_names(nm = names(clr_alt)))+
    ggplot2::labs(title = stringr::str_c(titles[[ loc2 ]])) +
    ggplot2::theme_minimal()+
    ggplot2::theme(legend.position = "none",
          panel.grid = ggplot2::element_blank(),
          plot.background = ggplot2::element_blank(),
          panel.background = ggplot2::element_rect(color = clr_pca[[ loc2 ]],
                                          fill = "transparent"),
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          text = ggplot2::element_text(size = plot_text_size),
          plot.title = ggplot2::element_text(color = clr_pca[[ loc2 ]]))
}
