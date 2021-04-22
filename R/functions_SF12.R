#' plot region tree (with Serranus outgroup)
#'
#' @param tree tree object
#' @param show_legend logical
#'
#' @export
plot_outl_tree_s <- function(tree, show_legend = TRUE){
  p <- ggtree::open_tree(ggtree::ggtree(tree,
                        layout = lyout,
                        aes(color = spec), size = .3), 180) +
    ggtree::geom_tiplab(aes(color = lab2spec(label),
                    label = stringr::str_sub(label, -6, -1)),
                size = 1.7,
                hjust = -.1)+
    ggtree::geom_treescale(width = .01,
                           x = .001,
                           y = 158,
                           offset = -1,
                           fontsize = 1.5,
                           color = clr_neutral) +
    ggtree::geom_nodepoint(aes(fill = support_class, size = support_class),
                           shape = 21) +
    scale_color_manual(values = c(ungrouped = clr_neutral,
                                  clr2),
                       guide = FALSE) +
    ggplot2::scale_fill_manual(values = c(`(0,50]` = "transparent",
                                 `(50,70]` = "white",
                                 `(70,90]` = "gray",
                                 `(90,100]` = "black"),
                      drop = FALSE) +
    ggplot2::scale_size_manual(values = c(`(0,50]` = 0,
                                 `(50,70]` = 1,
                                 `(70,90]` = 1,
                                 `(90,100]` = 1),
                      na.value = 0,
                      drop = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Node Support Class", title.position = "top", ncol = 2),
           size = ggplot2::guide_legend(title = "Node Support Class", title.position = "top", ncol = 2)) +
    ggplot2::theme_void()

  if(show_legend){
    p <- p +
      ggplot2::theme(legend.position = c(.5, .75),
            legend.justification = c(.5, 1),
            legend.text = ggplot2::element_text(size = plot_text_size))
  } else {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  ggplot2::ggplot() +
    ggplot2::coord_equal(xlim = c(0, 1),
                ylim = c(-.01, .54),
                expand = 0) +
    ggplot2::annotation_custom(grob = ggplot2::ggplotGrob(p),
                      ymin = -.6, ymax = .6,
                      xmin = -.1, xmax = 1.1) +
    ggplot2::theme_void()
}
