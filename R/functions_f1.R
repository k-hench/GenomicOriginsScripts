#' Import and summarize Fst data
#'
#' \code{summarize_fst} imports and summarizes Fst data.
#'
#' Fst data is imported, filtered for non-overlapping windows
#' and boxplot statistics are created.
#'
#' @family Figure 1
#'
#' @export
summarize_fst <- function(file){
  get_fst(file) %>%
    filter(BIN_START %% 50000 == 1) %>%
    select(GPOS, WEIGHTED_FST, run) %>%
    summarise(`mean_weighted-fst` = mean(WEIGHTED_FST),
              `median_weighted-fst` = median(WEIGHTED_FST),
              `lower_weighted-fst` = quantile(WEIGHTED_FST,.25),
              `upper_weighted-fst` = quantile(WEIGHTED_FST,.75),
              `sd_weighted-fst` = sd(WEIGHTED_FST),
              `lowpoint_weighted-fst` = whisker_points(WEIGHTED_FST)["low"],
              `highpoint_weighted-fst` = whisker_points(WEIGHTED_FST)["high"],
              run = run[1])
}

#' Import and summarize dxy data
#'
#' \code{summarize_dxy} imports and summarizes dxy data.
#'
#' dxy data is imported, filtered for non-overlapping windows
#' and boxplot statistics are created.
#'
#' @family Figure 1
#'
#' @export
summarize_dxy <- function(file){
  get_dxy(file) %>%
    filter(BIN_START %% 50000 == 1) %>%
    summarise(mean_dxy = mean(dxy),
              median_dxy = median(dxy),
              lower_dxy = quantile(dxy,.25),
              upper_dxy = quantile(dxy,.75),
              sd_dxy = sd(dxy),
              lowpoint_dxy = whisker_points(dxy)["low"],
              highpoint_dxy = whisker_points(dxy)["high"],
              run = run[1])
}

#' Whisker points
#'
#' \code{whisker_points} calculates the extent of boxplot whiskers.
#'
#' This is used to plot boxplots without outliers.
#'
#' @family Figure 1
#'
#' @export
whisker_points <- function(x){
  q1 <- quantile(x,.25)
  q2 <- quantile(x,.75)
  iqr <- q2-q1
  lp <- min(x[x > q1-(1.5*iqr)])
  hp <- max(x[x < q2+(1.5*iqr)])
  return(c(low = lp, high = hp))
}

# anno_pair_flag <- function (loc,left, right, circle_color = NA, circle_fill_left = "white",
#                             circle_fill_right = "lightgray", circle_lwd = 0.5, plot_names = FALSE,
#                             plot_name_size = 3, font_family = "sans", ...) {
#   nr_left <- which(hypo_img$spec %>% str_sub(.,1,3) == left)
#   nr_right <- which(hypo_img$spec %>% str_sub(.,1,3) == right)
#   nr_flag <- which((hypo_flag$geo %>% str_sub(.,1,3)) == loc)
#   p <- ggplot() +
#     theme_void() +
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(limits = c(-0.4, 0.38)) +
#     annotation_custom(hypo_flag$flag[[nr_flag]], xmin = -0.3, xmax = .3, ymin = -Inf, ymax = Inf) +
#     annotation_custom(hypo_img$l[[nr_right]], xmin = 0.05, xmax = 1, ymin = -Inf, ymax = Inf) +
#     annotation_custom(hypo_img$r[[nr_left]], xmin = -1, xmax = -0.05, ymin = -Inf, ymax = Inf)+
#     coord_cartesian(xlim = c(-1.1,1.1))
#   return(p)
# }

#' Plot single fish
#'
#' \code{plot_fish} creates a ggplot annotations layer with a single hamlet.
#'
#' This is used to annotate plots.
#'
#' @family Figure 1
#'
#' @export
plot_fish <- function(short, x = 0, y = 0, height = .75, width = .75){

  hypo_anno_l(sp_names[short],
              xmin = x - .5 * width, xmax = x + .5 * width,
              ymin = y - .5 * height, ymax = y + .5 * height)
}

#' make a faint version of a color
#'
#' \code{make_faint_clr} creates a faint version of a given color.
#'
#' @family Figure 1
#'
#' @export
make_faint_clr <- function(loc){
  colorRampPalette(scales::colour_ramp(c('white',clr_loc[loc]))(c(.5,1)))(2)
  }


#' crate a network layout for pairwise species comparisons
#'
#' \code{network_layout} creates a network layout for pairwise species comparisons.
#'
#' @family Figure 1
#'
#' @export
network_layout <- function(n, rotate = 0, label = NULL, weight = 2, loc = NA){
  tau <- seq(0,2*pi,length.out = n+1)[1:n] + rotate
  nodes <- tibble(idx = 1:n,
                  x = sin(tau),
                  y = cos(tau))

  edges <- cross_df(tibble(pop1 = nodes$idx,
                           pop2 = nodes$idx), .filter = `>=`) %>%
    arrange(pop1) %>%
    left_join(., nodes %>% select(idx, x,y) %>% set_names(., nm = c('pop1','x','y'))) %>%
    left_join(., nodes %>% select(idx, x,y) %>% set_names(., nm = c('pop2','xend','yend'))) %>%
    mutate(idx= row_number(),
           xmid_shift = (x*weight + xend)/(1*weight+1),
           ymid_shift = (y*weight + yend)/(1*weight+1),
           xmid = (x+xend)/2,
           ymid = (y+yend)/2) %>%
    select(idx, pop1:ymid)

  if(!is.null(label)){
    nodes <- nodes %>%
      mutate(label = label)

    edges <- edges %>% mutate(lab1 = label[pop1],
                              lab2 = label[pop2],
                              run = str_c(lab1,'-',lab2)) %>%
      left_join(.,run_ord)
  }

  tibble(loc = loc, nodes = list(nodes), edges = list(edges))
}


#' Plot pariwise species comparisons per location
#'
#' \code{plot_network} creates a ggplot with a network of species comparisons.
#'
#' This creates a network for all pairwise species comparisons within a location.
#'
#' @family Figure 1
#'
#' @export
plot_network <- function(loc, nodes, edges,asp = .8, sep = 0, node_lab_shift = 0){

  clr_prep <- scales::colour_ramp(c('black',clr_loc[loc]))(c(.4,1))
  clrs <- colorRampPalette(clr_prep)(max(edges$idx))

  p <- nodes %>%
    ggplot(aes(x, y))+
    coord_fixed(ratio = asp)+
    geom_segment(data = edges, aes(xend = xend, yend = yend),
                 color = clr_loc[loc], size = plot_lwd)

  for (k in nodes$idx) {
    p <- p + plot_fish(short = str_sub(nodes$label[k],1,3),
                       x = nodes$x[k],
                       y = nodes$y[k])
  }

  p + geom_label(data = edges,
                 aes(x = xmid_shift + sign(xmid_shift)*sep,
                                   y = ymid_shift + sign(ymid_shift)*sep*asp,
                                   label = run_ord),
                 color = clr_loc[loc],
                 label.padding = unit(1,'pt'),
                 label.size = 0,
                 size = (5/14) * plot_text_size)+
    scale_fill_manual(values = clr_loc, guide = FALSE) +
    scale_x_continuous(limits = c(-1.3,1.3),expand = c(0,0))+
    scale_y_continuous(expand = c(0,.1))+
    theme_void()
}
