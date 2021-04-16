#' Import revPoMo tree
#'
#' @param file tree file
#'
#' @export
import_tree <- function(file){
  tag = file %>%
    stringr::str_remove(".*/") %>%
    stringr::str_remove("_pop.cf.treefile")
  tibble::tibble(outlier = stringr::str_sub(tag, 1, 6) %>%
                   stringr::str_to_upper() %>%
                   stringr::str_replace("\\.","_"),
                 type = stringr::str_sub(tag, 8, -1) ,
                 tree = list(read.tree(file)))
}

#' Root hamlet tree
#'
#' Just for the flexible outgroup.
#'
#' @param tree string, input tree file
#' @param type string, tree type
#'
#' @export
root_hamlets <- function(tree, type){
  outgr <- switch (type,
                   "155" = "floflo",
                   "hyS" = c("tabhon", "torpan")
  )
  ape::root(phy = tree, outgroup = outgr)
}

#' Conditional Tree Rotation
#'
#' @param tree       tree object
#' @param rotate     logical, should the tree be rotated
#' @param root_node  node, where to rotate the tree
#'
#' @export
conditional_rotate <- function(tree, rotate = TRUE, root_node){
  if(rotate){ tree %>% ggtree::rotate(root_node) } else{ tree }
}

#' Conditional Tree Highlighting
#'
#' @param tree          tree object
#' @param highl         logical, should there be highlighting
#' @param higl_node     node, where to start highlighting
#' @param support_guide logical, show legend for node support
#'
#' @export
conditional_highlight <- function(tree,
                                  highl = TRUE,
                                  higl_node,
                                  support_guide = FALSE ){
  if(highl){
    p <- tree %>%
      ggtree::groupClade(.node = higl_node, group_name =  "clade") %>%
      ggtree::ggtree(layout = "circular", aes(color = clade), size = plot_lwd)
  } else {
    p <- tree %>%
      ggtree::ggtree(layout = "circular", size = plot_lwd)
  }
  p <- p +
    ggtree::geom_nodepoint(aes(fill = support_class,
                               size = support_class),
                           shape = 21) +
    ggplot2::scale_fill_manual(values = c(`(0,50]` = "transparent",
                                          `(50,70]` = "white",
                                          `(70,90]` = "gray",
                                          `(90,100]` = "black"),
                               drop = FALSE,
                               guide = FALSE) +
    ggplot2::scale_size_manual(values = c(`(0,50]` = 0,
                                          `(50,70]` = .9,
                                          `(70,90]` = .9,
                                          `(90,100]` = .9),
                               na.value = 0,
                               drop = FALSE,
                               guide = FALSE)

  if(support_guide){
    p <- p +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Support", row = 1, keywidth = unit(5, "pt")),
                      size = ggplot2::guide_legend(title = "Support", row = 1, keywidth = unit(5, "pt")))
  }
  p
}

#' Plot Tree
#'
#' @param tree          tree object
#' @param angle_in      numeric, tree opening angle
#' @param color         string, highlight color
#' @param higl_clade    logical, should there be highlighting
#' @param higl_node     node, where to start highlighting
#' @param root_node     node, where to rotate the tree
#' @param rotate_root   logical, should the tree be rotated
#' @param xlim          numeric vector, x limits
#' @param support_guide logical, show legend for node support
#' @param open_angle    numeric, tree opening angle
#'
#' @export
plot_tree <- function(tree, angle_in = 0,
                      color = "red",
                      higl_clade = TRUE,
                      higl_node = NA,
                      root_node = NA,
                      rotate_root = FALSE,
                      xlim = c(-2,4.5),
                      support_guide = FALSE,
                      open_angle = 180){
  p <- ( tree %>%
          conditional_highlight(highl = higl_clade, higl_node = higl_node, support_guide = support_guide) %>%
          conditional_rotate(rotate = rotate_root, root_node = root_node) %>%
          ggtree::open_tree(angle = open_angle) %>%
          ggtree::rotate_tree(angle = angle_in)) +
    ggtree::geom_tiplab(offset = diff(xlim) * .07, size = plot_text_size / ggplot2::.pt *.7) +
    ggplot2::scale_color_manual(values  = c(`0` = "black", `1` = color),
                                guide = FALSE) +
    ggplot2::scale_x_continuous(limits = xlim) +
    ggplot2::theme_void()

  ggplot2::ggplot() +
    ggplot2::coord_equal(xlim = c(-.1, .95),
                         ylim = c(-.01, .52),
                         expand = 0) +
    ggplot2::annotation_custom(grob = ggplot2::ggplotGrob(p),
                               ymin = 0 , ymax = 1,
                               xmin = 0, xmax = 1) +
    ggplot2::theme_void()
}
