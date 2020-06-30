#' Manage plot panels for a single outlier window
#'
#' \code{plot_curtain} calls all individual plotting functions needed for a single outlier zoom
#'
#' The panels containg the gene models, fst, dxy delta dxy,
#' genotype x phenotype association and topology weighting are called.
#'
#' Then they are combined using \code{cowplot::plot_grid} and
#' depending on the column position (first column/ other column),
#' the y axis titles are removed.
#'
#' @family Figure 5
#'
#' @export
plot_curtain <- function(loc = 'bel', outlier_id, outlier_nr, lg, start, end,
                         cool_genes,text = TRUE, label, trait, ...){
  p_g <- plot_panel_anno(lg = lg, outlier_id = outlier_id, label = label,
                         start = start,end = end, genes = cool_genes)
  p_fst <- plot_panel_fst(lg = lg, start = start,end = end)
  p_dxy <- plot_panel_dxy(lg = lg, start = start,end = end)
  p_delta_dxy <- plot_panel_delta_dxy(lg = lg, start = start,end = end)
  p_gxp <- plot_panel_gxp(lg = lg, start = start,end = end, trait = trait)
  p_t1 <- plot_panel_twisst(loc = loc, lg = lg, start = start,end = end, window_size = 200,
                            neighbours = tibble(left_neighbour = 'ind', right_neighbour = 'may', palette = twisst_clr['Blue']),
                            xlab = FALSE, highlight_mode = 'pair')
  p_t2 <- plot_panel_twisst(loc = loc, lg = lg, start = start,end = end, window_size = 200,
                            neighbours = tibble(left_neighbour = 'ind', right_neighbour = 'pue', palette = twisst_clr['Bars']),
                            xlab = FALSE, highlight_mode = 'pair')
  p_t3 <- plot_panel_twisst(loc = loc, lg = lg, start = start,end = end, window_size = 200,
                            neighbours = tibble(left_neighbour = 'uni', pops = list(pops_bel), palette = twisst_clr['Butter']),
                            xlab = FALSE, highlight_mode = 'isolation' )

  if(text){
    p_curtain <- cowplot::plot_grid(p_g,
                                    p_gxp,
                                    p_fst, p_dxy,
                                    p_delta_dxy,
                                    p_t1,p_t2,p_t3,
                                    ncol = 1, align = 'v',
                                    rel_heights = c(1, rep(.85, 7)))
  } else {
    p_curtain <- cowplot::plot_grid(p_g + no_title(),
                                    p_gxp + no_title(),
                                    p_fst + no_title(axis.text.y = element_blank(),
                                                     axis.ticks.y = element_blank()),
                                    p_dxy + no_title(axis.text.y = element_blank(),
                                                     axis.ticks.y = element_blank()),
                                    p_delta_dxy + no_title(axis.text.y = element_blank(),
                                                           axis.ticks.y = element_blank()),
                                    p_t1 + no_title(axis.text.y = element_blank(),
                                                    axis.ticks.y = element_blank()),
                                    p_t2 + no_title(axis.text.y = element_blank(),
                                                    axis.ticks.y = element_blank()),
                                    p_t3 + no_title(axis.text.y = element_blank(),
                                                    axis.ticks.y = element_blank()),
                                    ncol = 1, align = 'v',
                                    rel_heights = c(1, rep(.85, 7)))
  }
  p_curtain
}

#' Create the topology weighting panel
#'
#' \code{plot_panel_twisst} plots the topolgy weighting panel of a selected outlier window
#'
#' The topology weighting data for the selected location and
#' LG is fetched.
#' Then the number of possible topolgies is determined based
#' on the location and all topolgies are imported.
#' Dependend on th highlighting mode (pair/ isolation),
#' the topologies of interest are determined.
#' The topologies of interest are labelled in the dataset and
#' the color scheme for the highlighting is created.
#'
#' Finally, the topologie weighting panel is plotted.
#'
#' @family Figure 5
#'
#' @export
plot_panel_twisst <- function(loc, lg, start, end, window_size = twisst_size ,neighbours,xlab = TRUE, highlight_mode, ...){
  # fetch and subset data set for location
  data <- data_tables[[loc]]%>%
    filter(CHROM == lg ) %>%
    mutate(topo3 = str_pad(topo_nr, width = 3, pad = '0'))

  # determine number of topologies
  ntopo <- data$topo %>% unique() %>% length()

  # import topologies
  topo_plots <- read_lines(file = str_c(w_path,loc,'.LG01.w', twisst_size, '.phyml_bionj.weights.tsv.gz'),n_max = ntopo) %>%
    str_remove_all('#') %>%
    tibble(pre=.) %>%
    separate(pre,
             into = c('topo','tree'),sep = ' ') %>%
    rowwise() %>%
    mutate(phylo = list(read.tree(text = tree) %>% unroot()),
           topo_nr  = str_remove(topo,'topo') %>% as.numeric()) %>%
    ungroup()

  # highlight topologies of interest based on highlighting mode (pair / isolation)
  if(highlight_mode == 'pair'){
    topo_highlight <- neighbours %>%
      dplyr::select(-palette) %>%
      purrr::pmap(get_neighbour_topos,topo_plots = topo_plots) %>%
      set_names(.,nm = LETTERS[(length(neighbours$left_neighbour)+1):2]) %>%
      bind_cols() %>%
      gather(key = 'prefix',value = 'topo_nr')
  }  else if (highlight_mode == 'isolation') {
    topo_highlight <- neighbours %>%
      dplyr::select(-palette) %>%
      purrr::pmap(get_isolated_topos,topo_plots = topo_plots) %>%
      set_names(.,nm = LETTERS[(length(neighbours$left_neighbour)+1):2]) %>%
      bind_cols() %>%
      gather(key = 'prefix',value = 'topo_nr')
  }

  # label highlighted topologies
  data <- data %>%
    left_join(.,topo_highlight) %>%
    mutate(prefix = ifelse(is.na(prefix),"A",prefix),
           topo4 = str_c(prefix,topo3))

  # small highlighting function (needs to be defined here for the '<<-' to work)
  highlighter <- function(prefix,palette,...){
    high_n <- length(cols_topo[grepl(prefix,names(cols_topo))])
    high_base <- c(darken(palette), lighten(palette,factor = .1))
    high_clr <- colorRampPalette(high_base)(high_n)
    cols_topo[grepl(prefix,names(cols_topo))] <<- high_clr
  }

  # neutral topology color scheme
  cols_topo <- colorRampPalette(c("#00000082","#DEDEDE82"))(ntopo)
  names(cols_topo) <- data$topo4 %>% factor() %>% levels()

  # override colors of highlighted topologies
  neighbours %>%
    mutate(prefix = LETTERS[2:(length(neighbours$left_neighbour)+1)]) %>%
    purrr::pmap(highlighter)

  # plot topology weighting panel ----------------------------------
  p1 <- data %>%
    ggplot(aes(x = BIN_MID,
               y = weight,
               color = topo4,
               fill = topo4,
               group = str_c(topo4,CHROM, sep = '_')))+
    # add weighting results
    geom_area(position = 'stack',size = .2)+
    # add outlier area
    geom_rect(inherit.aes = FALSE, data = tibble(start = start, end=end),
              aes(xmin = start, xmax = end),ymin = -Inf, ymax = Inf,
              fill=rgb(1,1,1,.3), color = rgb(1,1,1,.9))+
    # use same boundaries for all panels
    coord_cartesian(xlim = c(start-window_buffer,end+window_buffer))+
    # use highlighter color scheme
    scale_color_manual(values = cols_topo, guide = FALSE)+
    scale_fill_manual(values = alpha(cols_topo,.5), guide = FALSE)+
    # layout x ayis
    scale_x_continuous(limits = c(start-window_buffer*1.25,end+window_buffer*1.25),expand = c(0,0),
                       labels = ax_scl) +
    # layout y ayis
    scale_y_continuous(expression(bolditalic(w)), expand = c(0.01, 0.01))+
    # use same plot appreance for all panels
    theme_panels()

  # if asked for it, drop x axis labels and title
  if(!xlab){p1 <- p1+theme(axis.text.x = element_blank())}

  p1
}

#' Create the fst panel
#'
#' \code{plot_panel_fst} plots the fst panel of a selected outlier window
#'
#' @family Figure 5
#'
#' @export
plot_panel_fst <- function(lg, start, end, xlab = TRUE, ...){
  ggplot2::ggplot() +
    # add outlier area
    geom_rect(inherit.aes = FALSE, data = tibble(start = start, end=end),
              aes(xmin = start, xmax = end),ymin = -Inf, ymax = Inf,
              fill=rgb(.9,.9,.9,.3), color = rgb(.9,.9,.9,.9)) +
    # add subsetted fst data
    geom_line(data = fst_data %>%
                filter(CHROM == lg, BIN_MID>start-window_buffer*1.25,BIN_MID<end+window_buffer*1.25) ,
              aes(x = BIN_MID, y = WEIGHTED_FST,
                  color = weighted_fst,#as.numeric(run),
                  group = run),size = .2) +
    # use same boundaries for all panels
    coord_cartesian(xlim = c(start-window_buffer,end+window_buffer))+
    # define color scheme
    scale_color_gradientn(name = expression(global~weighted~italic(F[ST])),
                          colours = hypogen::hypo_clr_LGs[1:24])+
    # layout x ayis
    scale_x_continuous(name = lg, expand = c(0,0), position = 'top') +
    # layout y ayis
    scale_y_continuous(name = expression(bolditalic(F[ST])),
                       expand = c(0,0),
                       limits = c(-0.07, 0.92))+
    # legend styling
    guides(color=guide_colorbar(barheight = unit(5,'pt'),barwidth = unit(200,'pt')))+
    # use same plot appreance for all panels
    theme_panels()
}

#' Create the dxy panel
#'
#' \code{plot_panel_dxy} plots the dxy panel of a selected outlier window
#'
#' The dxy data is subsetted to mathch the outlier area,
#' then the panel is plotted.
#'
#' @family Figure 5
#'
#' @export
plot_panel_dxy <- function(lg, start, end, ...){
  plot_data  <- dxy_data %>%
    select(CHROM, BIN_MID, dxy, weighted_fst, run) %>%
    filter(CHROM == lg,
           BIN_MID > (start-window_buffer),
           BIN_MID < (end+window_buffer))

  ggplot2::ggplot() +
    # add outlier area
    geom_rect(inherit.aes = FALSE, data = tibble(start = start, end=end),
              aes(xmin = start, xmax = end),ymin = -Inf, ymax = Inf,
              fill=rgb(.9,.9,.9,.3),color = rgb(.9,.9,.9,.9)) +
    # add dxy data
    geom_line(data = plot_data ,
              aes(x = BIN_MID, y = dxy, color = weighted_fst,
                  group = run),size = .2) +
    # use same boundaries for all panels
    coord_cartesian(xlim = c(start-window_buffer,end+window_buffer))+
    # define color scheme
    scale_color_gradientn(name = expression(global~weighted~italic(F[ST])),
                          colours = hypogen::hypo_clr_LGs[1:24])+
    # layout x ayis
    scale_x_continuous(name = lg, expand = c(0,0),position = 'top') +
    # layout y ayis
    scale_y_continuous(name = expression(bolditalic(d[XY])),
                       expand = c(0,0),
                       limits = c(0.0009, 0.0075))+
    # legend styling
    guides(color=guide_colorbar(barheight = unit(5,'pt'),barwidth = unit(300,'pt')))+
    # use same plot appreance for all panels
    theme_panels()
}

#' Create the genotype x phenotype panel
#'
#' \code{plot_panel_gxp} plots the genotype x phenotype panel of a selected outlier window
#'
#'
#' @family Figure 5
#'
#' @export
plot_panel_gxp <- function(lg, start, end, trait, ...){
  ggplot2::ggplot() +
    # add outlier area
    geom_rect(inherit.aes = FALSE,
              data = tibble(start = start, end=end),
              aes(xmin = start, xmax = end),
              ymin = -Inf, ymax = Inf,
              fill=rgb(.9,.9,.9,.3),color = rgb(.9,.9,.9,.9)) +
    # add subsetted gxp data
    geom_line(data = gxp_data %>%
                filter(CHROM == lg,
                       MID_POS>start-window_buffer*1.25,
                       MID_POS<end+window_buffer*1.25) ,
              aes(x = MID_POS, y = AVG_p_wald,
                  color = trt),size = .6) +
    hypoimg::geom_hypo_grob(data = tibble(grob = list(trait_grob[[trait]])),
                              #tibble(grob = hypoimg::hypo_trait_img$grob_circle[hypoimg::hypo_trait_img$trait == trait]),
                            aes(grob = grob, angle = 0, height = .65),
                            inherit.aes = FALSE, x = .9, y = 0.65)+
    # use same boundaries for all panels
    coord_cartesian(xlim = c(start-window_buffer,end+window_buffer))+
    # define color scheme
    scale_color_manual(name = 'GxP Trait',
                       values = gxp_clr)+
    # layout x ayis
    scale_x_continuous(name = lg, expand = c(0,0),position = 'top') +
    # layout y ayis
    scale_y_continuous(name = expression(bolditalic(-log[10](p))), expand = c(0,0))+
    # legend styling
    guides(color=guide_legend(keyheight =  unit(5,'pt'),
                              keywidth = unit(50,'pt'),
                              override.aes = list(size = 2))) +
    # use same plot appreance for all panels
    theme_panels()
}

#' Create the delta dxy panel
#'
#' \code{plot_panel_delta_dxy} plots the dxy panel of a selected outlier window
#'
#' The delta dxy data is subsetted to mathch the outlier area,
#' then the panel is plotted.
#'
#' @family Figure 5
#'
#' @export
plot_panel_delta_dxy <- function(lg, start, end, ...){
  dxy_summary <- data_dxy_summary %>%
    filter(scaffold == lg, mid>start-window_buffer*1.25,mid<end+window_buffer*1.25)

  ggplot2::ggplot() +
    # add delta dxy data
    geom_area(data = dxy_summary, aes(x = mid, y = delta_dxy),
              fill = 'lightgrey', color = 'darkgray')+
    # add outlier area
    geom_rect(inherit.aes = FALSE, data = tibble(start = start, end=end),
              aes(xmin = start, xmax = end),ymin = -Inf, ymax = Inf,
              fill=rgb(1, 1, 1, .3), color = rgb(1, 1, 1, .9)) +
    # use same boundaries for all panels
    coord_cartesian(xlim = c(start-window_buffer,end+window_buffer))+
    # layout x ayis
    scale_x_continuous(name = lg,
                       limits = c(start-window_buffer*1.25,
                                  end+window_buffer*1.25),
                       expand = c(0, 0),
                       position = 'top') +
    # layout y ayis
    scale_y_continuous(name = expression(bolditalic(Δ~d[XY])),
                       expand = c(0, 0),
                       breaks = c(0, 0.0025, 0.005),
                      # labels = c("0", "", "0.005"),
                       limits = c(0, 0.0051))+
    # use same plot appreance for all panels
    theme_panels()
}

#' Custom annotation base plot
#'
#' \code{custom_annoplot} extracts and plots hypoplectrus annotation data
#'
#' This is a modification of hypogen::hypo_annotation_baseplot().
#' The highlighting of the outlier area is added and the coloration
#' of genes is simplified to ony a single color.
#'
#' @family Figure 5
#'
#' @export
custom_annoplot <- function (..., searchLG, xrange, genes_of_interest = c(), genes_of_sec_interest = c(),
                             anno_rown = 3, width = 0.1, gene_color = 'darkgray', start, end) {
  # import annotation data of this LG and range
  df_list <- hypo_annotation_get(searchLG = searchLG, xrange = xrange,
                                 genes_of_interest = genes_of_interest,
                                 genes_of_sec_interest = genes_of_sec_interest,
                                 anno_rown = anno_rown)

  ggplot2::ggplot() +
    # add exons
    ggplot2::geom_rect(data = df_list[[2]],
                       aes(xmin = ps, xmax = pe, ymin = yl - (width/2), ymax = yl + (width/2), group = Parent),
                       fill = alpha(gene_color,.6), col = gene_color, lwd = 0.1) +
    # add outlier area
    geom_rect(inherit.aes = FALSE, data = tibble(start = start, end = end),
              aes(xmin = start, xmax = end),
              ymin = -Inf, ymax = Inf,
              fill=rgb(1,1,1,.3),color = rgb(1,1,1,.9)) +
    # add gene direction if known
    ggplot2::geom_segment(data = (df_list[[1]] %>% filter(strand %in% c("+", "-"))),
                          aes(x = ps, xend = pe,
                              y = yl, yend = yl, group = Parent),
                          lwd = 0.2, arrow = arrow(length = unit(2,"pt"), type = "closed"),
                          color = "black") +
    # add gene extent if direction is unknown
    ggplot2::geom_segment(data = (df_list[[1]] %>% filter(!strand %in% c("+", "-"))),
                          aes(x = ps, xend = pe,
                              y = yl, yend = yl, group = Parent), lwd = 0.2, color = "black") +
    # add gene label
    ggplot2::geom_text(data = df_list[[1]] %>% filter(label %in% genes_of_interest), size = 3,
                       aes(x = labelx, label = gsub("hpv1g000000", ".", label), y = yl - 0.5))
}

#' Create annotation panel
#'
#' \code{plot_panel_anno} plots the gene models of a selected outlier window
#'
#' The overal column title is created and the annotation data is plotted.
#'
#' @family Figure 5
#'
#' @export
plot_panel_anno <- function(outlier_id, label, lg, start, end, genes = c(),...){
  ttle <- str_sub(outlier_id,1,4) %>%
    str_c(.,' (',project_inv_case(label),')')

  p <- custom_annoplot(searchLG = lg,
                       xrange = c(start-window_buffer*1.25,end+window_buffer*1.25),
                       genes_of_interest = genes,
                       anno_rown = 6, start = start, end = end)+
    # layout x ayis
    scale_x_continuous(name = ttle,
                       position = 'top',
                       expand = c(0,0),
                       limits = c(start-window_buffer*1.25,end+window_buffer*1.25),
                       labels = ax_scl)+
    # layout y ayis
    scale_y_continuous(name = expression(bolditalic(Genes)), expand = c(0,.4))+
    # use same boundaries for all panels
    coord_cartesian(xlim = c(start-window_buffer,end+window_buffer))+
    # special panel layout for annotation panel
    theme_hypo()+
    theme(
      panel.background = element_rect(fill = rgb(.9,.9,.9),color = rgb(1,1,1,0)),
      legend.position = 'none',
      axis.title.x = element_text(),
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank())

  # use correct greec symbols in labels if needed
  if(outlier_id == 'LG17_1'){
    p$layers[[5]]$data$label <- p$layers[[5]]$data$label %>%
      str_replace(.,'alpha',"\u03B1")%>%
      str_replace(.,'beta',"\u03B2")
  }
  p
}

#' Plot fish for topology legend
#'
#' \code{plot_fish_zoom} creates a hamlet annotation from  a predefined table.
#'
#' This funtion is used in the creation of the topologie legend.
#'
#' @family Figure 5
#'
#' @export
plot_fish_zoom <- function(anno_tab,idx, height = .1, width = 5){
  annotation_custom(grob = anno_tab$grob[[idx]] %>% ggdraw() %>% cowplot::as_grob(),
                    xmin = anno_tab$x[[idx]] - width/2,
                    xmax = anno_tab$x[[idx]] + width/2,
                    ymin = anno_tab$y[[idx]] - height/2,
                    ymax = anno_tab$y[[idx]] + height/2 )
}

#' Create topology legend
#'
#' \code{plot_leg} creates the elements of the topology highlighting legend.
#'
#' This funtion is used to create the topologie legend of a single a single
#' highlighing scheme.
#' Depending on the highlighting mode (pair/ isolation), the background geometry
#' is created and the hamlet annotations are collected.
#'
#' Then the geometry is plotted an the hamlets are added.
#'
#' @family Figure 5
#'
#' @export
plot_leg <- function(spec1 = 'puella', spec2 = 'maya', color = "#084082ff", size = .5, mode = 'pair'){
  if(mode == 'pair'){
    # geometry of collapsed part of topology
    tri <- tibble(x = c(0,3,3,0),
                  y = c(0,-.5,.5,0))

    # branches of the pairs
    lns <- tibble(x = rep(c(-1,0),2),
                  y = c(.75,0,-.75,0),
                  grp = rep(c('a','b'),
                            each = 2))

    # collect hamlet pair annotations and generic hamlet
    ann <- tibble(grob = list(hypoimg::hypo_img$r[[which(hypoimg::hypo_img$spec == spec1)]],
                              hypoimg::hypo_img$r[[which(hypoimg::hypo_img$spec == spec2)]],
                              generic_hamlet_img),
                  x = c(-1,-1,2.45),
                  y = c(.75,-.75,0))

    # plot legened element
    ggplot(tri, aes(x,y))+
      geom_polygon(data = tri,color = color,fill = alpha(color,.4), size = size )+
      geom_line(data = lns,aes(group = grp),color = color, size = size ) +
      plot_fish_zoom(anno_tab = ann,idx = 1,width = 1,height = .5)+
      plot_fish_zoom(anno_tab = ann,idx = 2,width = 1,height = .5)+
      plot_fish_zoom(anno_tab = ann,idx = 3,width = 1,height = .75)+
      coord_equal(xlim = c(-1.4,3),
                  ylim = c(-1,1))+
      theme_void()
  } else  if(mode == 'isolation'){
    # geometry of collapsed part of topology
    tri <- tibble(x = c(0,3,3,0),
                  y = c(0,-.5,.5,0))

    # branches of the species
    lns <- tibble(x = c(-1,0),
                  y = c(0,0),
                  grp = rep('a',each = 2))

    # collect hamlet species annotation and generic hamlet
    ann <- tibble(grob = list(hypoimg::hypo_img$r[[which(hypoimg::hypo_img$spec == spec1)]],
                              generic_hamlet_img),
                  x = c(-1,2.45),
                  y = c(0,0))

    # plot legened element
    ggplot(tri, aes(x,y))+
      geom_polygon(data = tri,color = color,fill = alpha(color,.4), size = size )+
      geom_line(data = lns,aes(group = grp),color = color, size = size ) +
      plot_fish_zoom(anno_tab = ann,idx = 1,width = 1,height = .5)+
      plot_fish_zoom(anno_tab = ann,idx = 2,width = 1,height = .75)+
      coord_equal(xlim = c(-1.4,3),
                  ylim = c(-1,1))+
      theme_void()
  }
}

#' Create the fst population tree panel
#'
#' \code{plot_fst_poptree} plots the fst population tree panel of a selected outlier window
#'
#' @family Figure 5
#'
#' @export
plot_fst_poptree <- function(gid, data_nj, ...){

  data_nj$edge.length[data_nj$edge.length < 0 ] <- 0

  p <- data_nj %>%
    tidygraph::as_tbl_graph(layout = "dendrogram",
                            length = length,
                            directed = TRUE) %>%
    mutate(leaf = node_is_leaf()) %>%
    ggraph(layout = "dendrogram",
           length = length,
           circular = TRUE) +
    geom_edge_link()+
    geom_node_point(aes(filter = leaf,
                        fill = str_sub(name, 1, 3),
                        shape = str_sub(name, 4, 6)),
                    size = 2.5,)+
    scale_fill_manual("Species",
                      values = clr,
                      labels = sp_names %>%
                        str_c("*H. ",.,"*") %>%
                        set_names(nm = names(sp_names)),
                      guide = guide_legend(override.aes = list(shape = 21),
                                           nrow = 1))+
    scale_shape_manual("Location",
                       labels = loc_names,
                       values = 21:23) +
    theme_void()+
    theme(legend.position = "none",
          legend.text = element_markdown()) +
    coord_equal()+
    scale_y_reverse()

  p
}
