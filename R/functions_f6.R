#' Import Geva data (Allele Age)
#'
#' @param chrom string, linkage group identifier
#' @param start numeric, start position (bp)
#' @param end   numeric, end position (bp)
#' @param ...   catch all parameter, for usage with purrr::pmap() and unexpected parameters
#'
#' @export
import_geva_data <- function(chrom, start, end, ...){
  vroom::vroom(file = stringr::str_c(geva_path, chrom,".sites.txt.gz"), delim = " ") %>%
    dplyr::left_join(vroom::vroom(stringr::str_c(geva_path, chrom,".marker.txt.gz"), delim = " ")) %>%
    dplyr::arrange(Position) %>%
    dplyr::filter(between(Position, start, end)) %>%
    dplyr::select(Chromosome, Position, MarkerID, Clock, Filtered:PostMedian) %>%
    dplyr::mutate(CHROM = str_c("LG", stringr::str_pad(Chromosome, width = 2, pad = "0")))
}

#' Import windowed gxp data for focal trait
#'
#' @param trait string, gxp focal trait
#' @param chrom string, linkage group identifier
#' @param start numeric, start position (bp)
#' @param end   numeric, end position (bp)
#'
#' @export
gxp_importer <- function(trait, chrom, start, end ){
  vroom::vroom(stringr::str_c(gxp_path, trait,".lm.GxP.txt.gz"), delim = "\t") %>%
    dplyr::filter(CHROM == chrom, dplyr::between(POS, start, end) ) %>%
    dplyr::select(CHROM, POS, p_wald) %>%
    purrr::set_names(nm = c("CHROM", "Position", trait))
}

#' Import windowed gxp data for all traits
#'
#' @param chrom string, linkage group identifier
#' @param start numeric, start position (bp)
#' @param end   numeric, end position (bp)
#' @param ...   catch all parameter, for usage with purrr::pmap() and unexpected parameters
#'
#' @export
import_gxp_data <- function(chrom, start, end, ...){
  bar_data <- gxp_importer("Bars", chrom, start, end)
  peduncle_data <- gxp_importer("Peduncle", chrom, start, end)
  snout_data <- gxp_importer("Snout", chrom, start, end)

  bar_data %>%
    dplyr::left_join(peduncle_data, by = c(CHROM = "CHROM", Position = "Position")) %>%
    dplyr::left_join(snout_data, by = c(CHROM = "CHROM", Position = "Position"))
}

#' Combine allele age and gxp data import
#'
#' @param gid   string, focal outlier id
#' @param chrom string, linkage group identifier
#' @param start numeric, start position (bp)
#' @param end   numeric, end position (bp)
#' @param ...   catch all parameter, for usage with purrr::pmap() and unexpected parameters
#'
#' @export
get_gxp_and_geva <- function(gid, chrom, start, end, ...){
  import_geva_data(chrom, start, end) %>%
    dplyr::left_join(import_gxp_data(chrom, start, end)) %>%
    dplyr::mutate(gid = gid)
}

#' Custom reverse log transformation
#'
#' @param base numeric, delog base
#'
#' @export
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                    log_breaks(base = base),
                    domain = c(1e-100, Inf))
}

#' Single grob annotation for specific x range
#'
#' @param g   grob oject
#' @param x1  numeric, left border
#' @param x2  numeric, right border
#' @param ... catch all parameter, for usage with purrr::pmap() and unexpected parameters
#'
#' @export
grid_piece_x <- function(g, x1, x2, ...){
  ggplot2::annotation_custom(grob = g, xmin = x1, xmax = x2, ...)
}

#' Tiled Grob background for ggplot
#'
#' @param xlim_in numeric vector, x range of ggplot
#' @param colors  string vector, colors to use for background gradients
#' @param fun     function, x transformation function
#' @param step    numeric, step width for x breaks (after transformation)
#' @param ...     catch all parameter, for usage with purrr::pmap() and unexpected parameters
#'
#' @export
theme_gradient_bg_x <- function(xlim_in, colors = c("#7B0664", "#E32219"), fun = identity, step = 1,...){
  g <- grid::rasterGrob(cbind(colors[[1]], colors[[2]]), width = unit(1, "npc"),
                        height = unit(1, "npc"),
                        interpolate = TRUE)

  xstart <- floor(fun(min(xlim_in)))
  xend <- ceiling(fun(max(xlim_in)))

  breaks <- seq(xstart, xend, by = step)

  purrr::map2(.x = breaks[1:(length(breaks)-1)],
              .y = breaks[2:length(breaks)],
              .f = grid_piece_x,
              g = g, ...)
}
