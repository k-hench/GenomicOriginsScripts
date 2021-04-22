#' Import and summarize Fst data
#'
#' \code{summarize_fst} imports and summarizes Fst data.
#'
#' Fst data is imported, filtered for non-overlapping windows
#' and boxplot statistics are created.
#'
#' @param file input file
#'
#' @family Figure 1
#'
#' @export
summarize_fst <- function(file){
  get_fst(file) %>%
    dplyr::filter(BIN_START %% 50000 == 1) %>%
    dplyr::select(GPOS, WEIGHTED_FST, run) %>%
    dplyr::summarise(`mean_weighted-fst` = mean(WEIGHTED_FST),
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
#' @param file input file
#'
#' @family Figure 1
#'
#' @export
summarize_dxy <- function(file){
  get_dxy(file) %>%
    dplyr::filter(BIN_START %% 50000 == 1) %>%
    dplyr::summarise(mean_dxy = mean(dxy),
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
#' @param x input vector
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

#' Plot single fish
#'
#' \code{plot_fish} creates a ggplot annotations layer with a single hamlet.
#'
#' This is used to annotate plots.
#'
#' @param short three letter abbreviation of hamlet species
#' @param x horizontal center of grob
#' @param y vertical center of grob
#' @param height grob height
#' @param width grob width
#'
#' @family Figure 1
#'
#' @export
plot_fish <- function(short, x = 0, y = 0, height = .75, width = .75){
  hypoimg::hypo_anno_l(sp_names[short],
              xmin = x - .5 * width, xmax = x + .5 * width,
              ymin = y - .5 * height, ymax = y + .5 * height)
}

#' make a faint version of a color
#'
#' \code{make_faint_clr} creates a faint version of a given color.
#'
#' @param loc three letter location abbreviation
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
#' @param n nomber of nodes
#' @param rotate rotation angle
#' @param label node labels
#' @param weight link weighting
#' @param loc location (Belize/ Honduras/ Panama)
#'
#' @family Figure 1
#'
#' @export
network_layout <- function(n, rotate = 0, label = NULL, weight = 2, loc = NA){
  tau <- seq(0, 2*pi, length.out = n+1)[1:n] + rotate
  nodes <- tibble::tibble(idx = 1:n,
                          x = sin(tau),
                          y = cos(tau))

  edges <- purrr::cross_df(tibble::tibble(pop1 = nodes$idx,
                                          pop2 = nodes$idx), .filter = `>=`) %>%
    dplyr::arrange(pop1) %>%
    dplyr::left_join(., nodes %>% select(idx, x,y) %>% purrr::set_names(., nm = c('pop1','x','y'))) %>%
    dplyr::left_join(., nodes %>% select(idx, x,y) %>% purrr::set_names(., nm = c('pop2','xend','yend'))) %>%
    dplyr::mutate(idx = row_number(),
                  xmid_shift = (x*weight + xend)/(1*weight+1),
                  ymid_shift = (y*weight + yend)/(1*weight+1),
                  xmid = (x+xend)/2,
                  ymid = (y+yend)/2) %>%
    dplyr::select(idx, pop1:ymid)

  if(!is.null(label)){
    nodes <- nodes %>%
      dplyr::mutate(label = label)

    edges <- edges %>%
      dplyr::mutate(lab1 = label[pop1],
                    lab2 = label[pop2],
                    run = stringr::str_c(lab1,'-',lab2)) %>%
      dplyr::left_join(., run_ord)
  }

  tibble::tibble(loc = loc, nodes = list(nodes), edges = list(edges))
}


#' Plot pariwise species comparisons per location
#'
#' \code{plot_network} creates a ggplot with a network of species comparisons.
#'
#' This creates a network for all pairwise species comparisons within a location.
#'
#' @param loc location (Belize/ Honduras/ Panama)
#' @param nodes tibble containing network nodes
#' @param edges tibble containing network edges
#' @param asp plot aspect ration
#' @param sep buffer network segments from nodes
#' @param node_lab_shift relative label position on edge
#'
#' @family Figure 1
#'
#' @export
plot_network <- function(loc, nodes, edges, asp = .8, sep = 0, node_lab_shift = 0){

  loc_edge <- c(bel = .68, hon = .66, pan = .82) -.03
  clr_prep <- (scales::colour_ramp(c("black",
                                     clr_loc[loc])))(c(0.4, 1))
  clrs <- colorRampPalette(clr_prep)(max(edges$idx))
  p <- nodes %>% ggplot(aes(x, y)) +
    coord_fixed(ratio = asp) +
    geom_segment(data = edges,
                 aes(xend = xend, yend = yend),#, size = median),
                 size = .1,
                 color = clr_loc[loc]) +#, size = plot_lwd)
    # scale_size(limits = c(0, 1), range = c(.1, 4))+
    scale_size(limits = c(0, 1), range = c(.1, 2))+
    scale_color_manual(values = clr)
  for (k in nodes$idx) {
    p <- p + plot_fish_lwd(short = str_sub(nodes$label[k], 1,  3),
                           x = nodes$x[k], y = nodes$y[k], height = .7, width = .7)
  }
  p + geom_label(data = edges, aes(x = xmid_shift + sign(xmid_shift) * sep,
                                   y = ymid_shift + sign(ymid_shift) * sep * asp,
                                   label = run_ord),
                 color = clr_loc[loc],
                 label.padding = unit(1, "pt"),
                 label.size = 0,
                 size = plot_text_size * .5 /ggplot2::.pt) +
    scale_fill_manual(values = clr_loc,
                      guide = FALSE) +
    scale_x_continuous(limits = c(-1.3, 1.3),
                       expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0.1)) +
    theme_void()
}


#' Plot Hamlet drawing (alt)
#'
#' Plot hamlet drawing with adjustable line width for line elements of drawing.
#'
#' @param name       string, species name
#' @param x          numeric, horizontal center of hamlet image
#' @param y          numeric, vertical center of hamlet image
#' @param height     numeric, height of hamlet image
#' @param width      numeric, width of hamlet image
#' @param lwd        numeric, line width of hamlet drawing line elements
#' @param line_color string, color code, color of hamlet drawing line elements
#'
#' @export
plot_fish2 <- function (name, x = 0, y = 3, height = 4, width = 4, lwd = .15, line_color = "black") {
  spec <- stringr::str_remove(string = name, "Hypoplectrus_")
  hypo_anno_l_lwd(spec, xmin = x - 0.5 * width, xmax = x +
                    0.5 * width, ymin = y - 0.5 * height, ymax = y + 0.5 *
                    height, lwd = lwd, line_color = line_color)
}

#' Hamlet Annotation (alt)
#'
#' Add ggplot hamlet drawing  annotation with adjustable line width for line elements of drawing.
#'
#' @param species       string, species name
#' @param xmin          numeric, left border of hamlet image
#' @param xmax          numeric, right border of hamlet image
#' @param ymin          numeric, lower border of hamlet image
#' @param ymax          numeric, upper border of hamlet image
#' @param lwd        numeric, line width of hamlet drawing line elements
#' @param line_color string, color code, color of hamlet drawing line elements
#'
#' @export
hypo_anno_l_lwd <- function (species, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, lwd = .15, line_color = "black") {
  stopifnot(length(species) == 1)
  stopifnot(is.character(species))
  stopifnot(species %in% hypo_img$spec)
  nr_species <- which(hypo_img$spec == species)
  annotation_custom(grid::editGrob(grob = hypo_img$l[[nr_species]],
                                   gPath = "GRID.picComplexPath.*", grep = TRUE,
                                   gp = grid::gpar( lwd = lwd, col = line_color
                                   ),
                                   global = TRUE, strict = FALSE) ,
                    xmin = xmin,
                    xmax = xmax, ymin = ymin, ymax = ymax)
}


#' Hamlet Annotation (alt)
#'
#' Add ggplot hamlet drawing  annotation with adjustable line width for line elements of drawing.
#'
#' @param short       string, three letter species name abbreviation
#' @param x          numeric, horizontal center of hamlet image
#' @param y          numeric, vertical center of hamlet image
#' @param height     numeric, height of hamlet image
#' @param width      numeric, width of hamlet image
#' @param lwd        numeric, line width of hamlet drawing line elements
#' @param line_color string, color code, color of hamlet drawing line elements
#'
#' @export
plot_fish_lwd <- function (short, x = 0, y = 3, height = 5, width = 5, lwd = .15, line_color = "transparent") {
  hypo_anno_l_lwd(sp_names[short], xmin = x - 0.5 * width, xmax = x +
                    0.5 * width, ymin = y - 0.5 * height, ymax = y + 0.5 *
                    height, lwd = lwd, line_color = line_color)
}

#' Automated PCA plotting
#'
#' @param loc string, three letter location abbreviation
#'
#' @export
pca_plot <- function(loc){
  set.seed(42)
  evs <- stringr::str_c("2_analysis/pca/", loc ,".exp_var.txt.gz") %>%
    readr::read_tsv()
  stringr::str_c("2_analysis/pca/", loc ,".scores.txt.gz") %>%
    readr::read_tsv() %>%
    dplyr::mutate(spec = stringr::str_sub(id, -6,-4)) %>%
    ggplot2::ggplot(aes(x = EV01, y = EV02, fill = spec))+
    ggforce::geom_mark_ellipse(aes(color = spec),
                               fill = "transparent",
                               linetype = 3,
                               size = .3,
                               expand = unit(5, "pt"))+
    ggplot2::geom_point(shape = 21, aes(color = ggplot2::after_scale(prismatic::clr_darken(fill))), size = .7) +
    (pca_fish_pos$data[[which(pca_fish_pos$loc == loc)]] %>%
       purrr::pmap(plot_fish_lwd))+
    labs(x = str_c("PC1 (", sprintf("%.1f",evs$exp_var[[1]]), " %)"),
         y = str_c("PC2 (", sprintf("%.1f",evs$exp_var[[2]]), " %)"))+
    ggplot2::scale_fill_manual(values = clr)+
    ggplot2::scale_color_manual(values = clr_alt %>%
                                  prismatic::clr_alpha(alpha = .7) %>%
                                  purrr::set_names(nm = names(clr_alt)))+
    ggplot2::labs(title = loc_names[[loc]])+
    ggplot2::theme_minimal()+
    ggplot2::theme(legend.position = "none",
                   panel.grid = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(color = clr_loc[[loc]],
                                                            fill = "transparent"),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = plot_text_size),
                   plot.title = ggplot2::element_text(color = clr_loc[[loc]]))
}


#' A modded version of BAMMtools::plot.bammdata
#'
#' @param x              An object of class bammdata.
#' @param tau            A numeric indicating the grain size for the calculations. See documentation for dtRates.
#' @param method         A character string indicating the method for plotting the phylogenetic tree. method = "phylogram" (default) plots the phylogenetic tree using rectangular coordinates. method = "polar" plots the phylogenetic tree using polar coordinates.
#' @param xlim           A numeric vector of coordinates for the x-axis endpoints. Defaults to NULL, in which case they are calculated automatically. The x-axis corresponds to time when the phylogeny is plotted facing to the left or to the right. The time at the root equals zero.
#' @param ylim           A numeric vector of coordinates for the y-axis endpoints. Defaults to NULL, in which case they are calculated automatically. Tips are plotted at integer values beginning from zero and stepping by one when the phylogeny is plotted facing to the left or to the right.
#' @param vtheta         A numeric indicating the angular separation (in degrees) of the first and last terminal nodes. Ignored if method = "phylogram".
#' @param rbf            A numeric indicating the length of the root branch as a fraction of total tree height. Ignored if method = "phylogram".
#' @param show           A logical indicating whether or not to plot the tree. Defaults to TRUE.
#' @param labels         A logical indicating whether or not to plot the tip labels. Defaults to FALSE.
#' @param legend         A logical indicating whether or not to plot a legend for interpreting the mapping of evolutionary rates to colors. Defaults to FALSE.
#' @param spex           A character string indicating what type of macroevolutionary rates should be plotted. "s" (default) indicates speciation rates, "e" indicates extinction rates, and "netdiv" indicates net diversification rates. Ignored if ephy$type = "trait".
#' @param lwd            A numeric specifying the line width for branches.
#' @param cex            A numeric specifying the size of tip labels.
#' @param pal            A character string or vector of mode character that describes the color palette. See Details for explanation of options.
#' @param mask           An optional integer vector of node numbers specifying branches that will be masked with mask.color when plotted.
#' @param mask.color     The color for the mask.
#' @param colorbreaks    A numeric vector of percentiles delimiting the bins for mapping rates to colors. If NULL (default) bins are calculated from the rates that are passed with the bammdata object.
#' @param logcolor       Logical. Should colors be plotted on a log scale.
#' @param breaksmethod   Method used for determining color breaks. See help file for assignColorBreaks.
#' @param color.interval Min and max value for the mapping of rates. If NULL, then min and max are inferred from the data. NA can also be supplied for one of the two values. See details.
#' @param JenksSubset    If breaksmethod = "jenks", the number of regularly spaced samples to subset from the full rates vector. Only relevant for large datasets. See help file for assignColorBreaks.
#' @param par.reset      A logical indicating whether or not to reset the graphical parameters when the function exits. Defaults to FALSE.
#' @param direction      A character string. Options are "rightwards", "leftwards", "upwards", and "downwards", which determine the orientation of the tips when the phylogeny plotted.
#' @param labelcolor     string, color code, a vector of colors specifying the label colors
#' @param ...            Further arguments passed to par.
#'
#' @export
bammplot_k <- function (x, tau = 0.01, method = "phylogram", xlim = NULL,
                        ylim = NULL, vtheta = 5, rbf = 0.001, show = TRUE, labels = FALSE,
                        legend = FALSE, spex = "s", lwd = 1, cex = 1, pal = "RdYlBu",
                        mask = integer(0), mask.color = gray(0.5), colorbreaks = NULL,
                        logcolor = FALSE, breaksmethod = "linear", color.interval = NULL,
                        JenksSubset = 20000, par.reset = FALSE, direction = "rightwards",
                        labelcolor = "black",
                        ...)
{
  if (inherits(x, "bammdata")) {
    if (attributes(x)$order != "cladewise") {
      stop("Function requires tree in 'cladewise' order")
    }
    phy <- BAMMtools:::as.phylo.bammdata(x)
  }
  else stop("Object ephy must be of class bammdata")
  if (!spex %in% c("s", "e", "netdiv")) {
    stop("spex must be 's', 'e' or 'netdiv'.")
  }
  if (length(pal) == 1 && !pal %in% names(get("palettes",
                                              envir = .colorEnv)) && pal != "temperature" && pal !=
      "terrain")
    pal <- rep(pal, 3)
  else if (length(pal) == 2)
    pal <- c(pal, pal[2])
  if (breaksmethod == "linear" & !is.null(color.interval)) {
    if (length(color.interval) != 2) {
      stop("color.interval must be a vector of 2 numeric values.")
    }
  }
  if (!is.binary.phylo(phy)) {
    stop("Function requires fully bifurcating tree")
  }
  if (any(phy$edge.length == 0)) {
    warning("Tree contains zero length branches. Rates for these will be NA and coerced to zero")
  }
  if (!("dtrates" %in% names(x))) {
    x <- dtRates(x, tau)
  }
  NCOLORS <- 64
  if (!is.null(color.interval)) {
    if (x$type == "trait") {
      ratesRange <- range(x$dtrates$rates)
    }
    else if (x$type == "diversification") {
      if (tolower(spex) == "s") {
        ratesRange <- range(x$dtrates$rates[[1]])
      }
      else if (tolower(spex) == "e") {
        ratesRange <- range(x$dtrates$rates[[2]])
      }
      else if (tolower(spex) == "netdiv") {
        ratesRange <- range(x$dtrates$rates[[1]] - x$dtrates$rates[[2]])
      }
    }
    if (all(!is.na(color.interval))) {
      brks <- seq(min(color.interval[1], ratesRange[1]),
                  max(color.interval[2], ratesRange[2]), length.out = (NCOLORS +
                                                                         1))
      intervalLength <- length(which.min(abs(color.interval[1] -
                                               brks)):which.min(abs(color.interval[2] - brks)))
    }
    else if (is.na(color.interval[1])) {
      brks <- seq(ratesRange[1], max(color.interval[2],
                                     ratesRange[2]), length.out = (NCOLORS + 1))
      intervalLength <- length(1:which.min(abs(color.interval[2] -
                                                 brks)))
    }
    else if (is.na(color.interval[2])) {
      brks <- seq(min(color.interval[1], ratesRange[1]),
                  ratesRange[2], length.out = (NCOLORS + 1))
      intervalLength <- length(which.min(abs(color.interval[1] -
                                               brks)):length(brks))
    }
    NCOLORS <- round((NCOLORS^2)/intervalLength)
  }
  if (is.null(colorbreaks)) {
    colorbreaks <- assignColorBreaks(x$dtrates$rates, NCOLORS,
                                     spex, logcolor, breaksmethod, JenksSubset)
  }
  if (x$type == "trait") {
    colorobj <- BAMMtools:::colorMap(x$dtrates$rates, pal, colorbreaks,
                                     logcolor, color.interval)
  }
  else if (x$type == "diversification") {
    if (tolower(spex) == "s") {
      colorobj <- BAMMtools:::colorMap(x$dtrates$rates[[1]], pal,
                                       colorbreaks, logcolor, color.interval)
    }
    else if (tolower(spex) == "e") {
      colorobj <- BAMMtools:::colorMap(x$dtrates$rates[[2]], pal,
                                       colorbreaks, logcolor, color.interval)
    }
    else if (tolower(spex) == "netdiv") {
      colorobj <- BAMMtools:::colorMap(x$dtrates$rates[[1]] - x$dtrates$rates[[2]],
                                       pal, colorbreaks, logcolor, color.interval)
    }
  }
  else {
    stop("Unrecognized/corrupt bammdata class. Type does not equal 'trait' or 'diversification'")
  }
  edge.color <- colorobj$cols
  tH <- max(x$end)
  phy$begin <- x$begin
  phy$end <- x$end
  tau <- x$dtrates$tau
  if (method == "polar") {
    ret <- setPolarTreeCoords(phy, vtheta, rbf)
    rb <- tH * rbf
    p <- BAMMtools:::mkdtsegsPolar(ret$segs[-1, ], tau, x$edge)
  }
  else if (method == "phylogram") {
    ret <- BAMMtools:::setPhyloTreeCoords(phy)
    p <- BAMMtools:::mkdtsegsPhylo(ret$segs[-1, ], tau, x$edge)
  }
  else {
    stop("Unimplemented method")
  }
  x0 <- c(ret$segs[1, 1], p[, 1])
  x1 <- c(ret$segs[1, 3], p[, 2])
  y0 <- c(ret$segs[1, 2], p[, 3])
  y1 <- c(ret$segs[1, 4], p[, 4])
  offset <- table(p[, 5])[as.character(unique(p[, 5]))]
  if (length(mask)) {
    edge.color[p[, 5] %in% mask] <- mask.color
  }
  arc.color <- c(edge.color[1], edge.color[match(unique(p[,
                                                          5]), p[, 5]) + offset])
  edge.color <- c(edge.color[1], edge.color)
  if (show) {
    op <- par(no.readonly = TRUE)
    if (length(list(...))) {
      par(...)
    }
    if (legend) {
      par(mar = c(5, 4, 4, 5))
    }
    plot.new()
    ofs <- 0
    if (labels) {
      if (method == "phylogram")
        ofs <- max(nchar(phy$tip.label) * 0.03 * cex *
                     tH)
      else ofs <- max(nchar(phy$tip.label) * 0.03 * cex)
    }
    if (method == "polar") {
      if (is.null(xlim) || is.null(ylim)) {
        if (is.null(xlim))
          xlim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs)
        if (is.null(ylim))
          ylim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs)
      }
      plot.window(xlim = xlim, ylim = ylim, asp = 1)
      segments(x0, y0, x1, y1, col = edge.color, lwd = lwd,
               lend = 2)
      arc(0, 0, ret$arcs[, 1], ret$arcs[, 2], c(rb, rb +
                                                  phy$end/tH), border = arc.color, lwd = lwd)
      if (labels) {
        for (k in 1:length(phy$tip.label)) {
          text(ret$segs[-1, ][phy$edge[, 2] == k, 3],
               ret$segs[-1, ][phy$edge[, 2] == k, 4], phy$tip.label[k],
               cex = cex, srt = (180/pi) * ret$arcs[-1,
               ][phy$edge[, 2] == k, 1], adj = c(0, NA),
               col = labelcolor[[k]])
        }
      }
    }
    if (method == "phylogram") {
      direction <- match.arg(direction, c("rightwards",
                                          "leftwards", "downwards", "upwards"))
      if (direction == "rightwards") {
        bars <- BAMMtools:::redirect(cbind(x0, y0, x1, y1), 0)
        arcs <- BAMMtools:::redirect(ret$arcs, 0)
        bars[, c(1, 3)] <- tH * bars[, c(1, 3)]
        arcs[, c(1, 3)] <- tH * arcs[, c(1, 3)]
        ret$segs[-1, c(1, 3)] <- tH * ret$segs[-1, c(1,
                                                     3)]
      }
      else if (direction == "leftwards") {
        bars <- BAMMtools:::redirect(cbind(x0, y0, x1, y1), pi)
        bars[, c(2, 4)] <- abs(bars[, c(2, 4)])
        arcs <- BAMMtools:::redirect(ret$arcs, pi)
        arcs[, c(2, 4)] <- abs(arcs[, c(2, 4)])
        bars[, c(1, 3)] <- tH * bars[, c(1, 3)]
        arcs[, c(1, 3)] <- tH * arcs[, c(1, 3)]
        ret$segs[-1, c(1, 3)] <- -tH * ret$segs[-1,
                                                c(1, 3)]
      }
      else if (direction == "downwards") {
        bars <- BAMMtools:::redirect(cbind(x0, y0, x1, y1), -pi/2)
        arcs <- BAMMtools:::redirect(ret$arcs, -pi/2)
        bars[, c(2, 4)] <- tH * bars[, c(2, 4)]
        arcs[, c(2, 4)] <- tH * arcs[, c(2, 4)]
        ret$segs <- BAMMtools:::redirect(ret$segs, -pi/2)
        ret$segs[, c(2, 4)] <- tH * ret$segs[, c(2,
                                                 4)]
      }
      else if (direction == "upwards") {
        bars <- BAMMtools:::redirect(cbind(x0, y0, x1, y1), pi/2)
        bars[, c(1, 3)] <- abs(bars[, c(1, 3)])
        arcs <- BAMMtools:::redirect(ret$arcs, pi/2)
        arcs[, c(1, 3)] <- abs(arcs[, c(1, 3)])
        bars[, c(2, 4)] <- tH * bars[, c(2, 4)]
        arcs[, c(2, 4)] <- tH * arcs[, c(2, 4)]
        ret$segs <- BAMMtools:::redirect(ret$segs, pi/2)
        ret$segs[, c(1, 3)] <- abs(ret$segs[, c(1, 3)])
        ret$segs[, c(2, 4)] <- tH * ret$segs[, c(2,
                                                 4)]
      }
      if (is.null(xlim) && direction == "rightwards")
        xlim <- c(0, tH + ofs)
      if (is.null(xlim) && direction == "leftwards")
        xlim <- c(-(tH + ofs), 0)
      if (is.null(ylim) && (direction == "rightwards" ||
                            direction == "leftwards"))
        ylim <- c(0, phy$Nnode)
      if (is.null(xlim) && (direction == "upwards" ||
                            direction == "downwards"))
        xlim <- c(0, phy$Nnode)
      if (is.null(ylim) && direction == "upwards")
        ylim <- c(0, tH + ofs)
      if (is.null(ylim) && direction == "downwards")
        ylim <- c(-(tH + ofs), 0)
      plot.window(xlim = xlim, ylim = ylim)
      segments(bars[-1, 1], bars[-1, 2], bars[-1, 3],
               bars[-1, 4], col = edge.color[-1], lwd = lwd,
               lend = 2)
      isTip <- phy$edge[, 2] <= phy$Nnode + 1
      isTip <- c(FALSE, isTip)
      segments(arcs[!isTip, 1], arcs[!isTip, 2], arcs[!isTip,
                                                      3], arcs[!isTip, 4], col = arc.color[!isTip],
               lwd = lwd, lend = 2)
      if (labels) {
        if (direction == "rightwards")
          text(ret$segs[isTip, 3], ret$segs[isTip, 4],
               phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex,
               pos = 4, offset = 0.25,
               col = labelcolor)
        else if (direction == "leftwards")
          text(ret$segs[isTip, 3], ret$segs[isTip, 4],
               phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex,
               pos = 2, offset = 0.25,
               col = labelcolor)
        else if (direction == "upwards")
          text(ret$segs[isTip, 3], ret$segs[isTip, 4],
               phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex,
               pos = 4, srt = 90, offset = 0,
               col = labelcolor)
        else if (direction == "downwards")
          text(ret$segs[isTip, 3], ret$segs[isTip, 4],
               phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex,
               pos = 2, srt = 90, offset = 0,
               col = labelcolor)
      }
    }
  }
  index <- order(as.numeric(rownames(ret$segs)))
  if (show) {
    if (method == "phylogram") {
      assign("last_plot.phylo", list(type = "phylogram",
                                     direction = direction, Ntip = phy$Nnode + 1,
                                     Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index,
                                                                                       3], yy = ret$segs[index, 4], pp = par(no.readonly = TRUE)),
             envir = .PlotPhyloEnv)
    }
    else if (method == "polar") {
      assign("last_plot.phylo", list(type = "fan", Ntip = phy$Nnode +
                                       1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index,
                                                                                            3], yy = ret$segs[index, 4], theta = ret$segs[index,
                                                                                                                                          5], rb = rb, pp = par(no.readonly = TRUE)),
             envir = .PlotPhyloEnv)
    }
    if (legend) {
      addBAMMlegend(x = list(coords = ret$segs[-1, ],
                             colorbreaks = colorobj$breaks, palette = colorobj$colpalette,
                             colordens = colorobj$colsdensity), location = "right")
    }
  }
  if (par.reset) {
    par(op)
  }
  invisible(list(coords = ret$segs[-1, ], colorbreaks = colorobj$breaks,
                 palette = colorobj$colpalette, colordens = colorobj$colsdensity))
}
