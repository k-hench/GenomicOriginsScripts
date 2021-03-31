#' Prepate topology weighting data for zoom plot
#'
#' \code{prep_data} imports the twisst data of a selected LG
#'
#' Based on the focal outlier names, the linkage group is determined.
#' The topology weighting results and positions files are located and
#' filtered for the focal LG and location.
#'
#' Then, the selected files are imported.
#'
#'
#' @param loc sample location (bel [Belize]/ hon [Honduras]/ pan [Panama])
#' @param window_size size of twisst windows (bp)
#' @param ... catch-all parameter for purrr::pmap to ignore excessive parameters
#'
#' @family Figure 5
#'
#' @export
prep_data <- function(loc, window_size = twisst_size, ...){
  LG_select <- outlier_pick %>% str_sub(.,1,4) %>% unique()

  weight_files_prep <- dir(w_path, pattern = loc) %>% .[grepl(str_c('w',window_size),.)]
  weight_files_select <- weight_files_prep %>% str_sub(.,5,8)
  weight_files <- weight_files_prep[weight_files_select %in% LG_select]

  data_files_prep <- dir(str_c(d_path, loc, '/'), pattern = 'LG.*data.tsv') %>% .[grepl(str_c('w',window_size),.)]
  data_files_select <- data_files_prep %>% str_sub(.,5,8)
  data_files <- data_files_prep[data_files_select %in% LG_select]

  data <- tibble(w_in = weight_files,
                 d_in = data_files) %>%
    purrr::pmap(get_twisst_data, smooth = FALSE, loc = loc) %>%
    bind_rows()
}

#' Defin topology weighing base color
#'
#' \code{get_clr} gets the most intense color of a given colorbrewer palette
#'
#'
#' @param palette string - RColorBrewer palette identifier
#'
#' @family Figure 5
#'
#' @export
get_clr <- function(palette){
  palette %>%
    purrr::map(.f = function(palette){RColorBrewer::brewer.pal(5,palette)[5]}) %>%
    unlist()
  }

#' Scale bp to Mb
#'
#' \code{ax_scl} scales the x axis of ggplots
#'
#' This is a scaling function to transform bp positions to Mb.
#'
#' @param x continuous numerical vector
#'
#' @family Figure 5
#'
#' @export
ax_scl <- function(x){ x/(10^6) }


#' ggplot layout template
#'
#' \code{theme_panels} provides a template for the panels of Figure 3
#'
#' @param ... parameters to funnel through to ggplot2::theme()
#'
#' @family Figure 5
#'
#' @export
theme_panels <- function(...){
  hypogen::theme_hypo() +
    hypogen::theme_hypo_anno_extra()+
    ggplot2::theme(text = ggplot2::element_text(size = plot_text_size),
          legend.position = 'none',
          axis.title.y = ggplot2::element_text(angle = 90),
          axis.line.y = ggplot2::element_line(size = plot_lwd),
          axis.ticks.y = ggplot2::element_line(size = plot_lwd),
          axis.title.x = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(t = 1, r = 1, b = 3, l = 1),
          ...)
}

#' drop ggplot y axis titles
#'
#' \code{no_title} provides further styling for panels of Figure 3
#'
#' This function removes the  y axis title and grid elements
#' of second and third column panels in Figure 3.
#'
#' @param ... parameters to funnel through to ggplot2::theme()
#'
#' @family Figure 5
#'
#' @export
no_title <- function(...){
  theme(axis.title.y = element_blank(),
        panel.grid.major.x = element_line(colour = hypogen::hypo_clr_lg),
        panel.grid.minor.x = element_line(colour = hypogen::hypo_clr_lg),
        ...)
}

#' Import genotype x phenotype association in long format
#'
#' \code{get_gxp_long} imports genotype x phenotype association data in long format.
#'
#' The phenotype trait name is extracted from the file name.
#' Then, the data is imported, the genomic poistion and the trait name is added as columns.
#'
#' @param file gpx results file
#' @param kb gxp window size (matching file for string manipulation)
#'
#' @family Figure 5
#'
#' @export
get_gxp_long <- function(file, kb = 10){
  trt <- file %>%
    str_remove('^.*/') %>%
    str_remove(str_c('.lm.',kb, 'k.',kb/10, 'k.txt.gz'))

  data <- file %>%
    vroom::vroom(delim = "\t") %>%
    left_join(.,hypogen::hypo_chrom_start) %>%
    add_gpos()%>%
    dplyr::mutate(trt = trt)
  data
}
