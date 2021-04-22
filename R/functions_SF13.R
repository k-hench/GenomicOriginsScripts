#' import sc data sheet
#'
#' \code{read_sc} imports a sc data sheet.
#'
#' @param file input file
#'
#'
#' @export
read_sc <-function(file){
  length_comment <- readr::read_lines(file) %>% stringr::str_which(.,'^#') %>% max()
  length_format <- readr::read_lines(file) %>% stringr::str_which(.,'^format*') %>% max()
  pos_frame <- readr::read_lines(file) %>% stringr::str_which(.,'^frame*') %>% max()

  head_length <- max(length_comment, length_format, pos_frame)

  readr::read_delim(file, delim = ' ', skip = head_length,
                    col_names = c('type', 'cell', 'delim', 'value')) %>%
    dplyr::filter(!(type == "goto"), !(is.na(value))) %>%
    dplyr::select(-type,-delim) %>%
    dplyr::mutate(cell = str_replace(cell,"^([A-Z])([0-9])","\\1_\\2")) %>%
    tidyr::separate(cell, into = c('column','row'),sep = '_') %>%
    dplyr::mutate(row = as.numeric(row)) %>%
    tidyr::spread(key = column, value = value) %>%
    setNames(.,nm = .[1,] %>% as.character()) %>%
    dplyr::filter(`0` != 0) %>%
    dplyr::select(-`0`) %>%
    dplyr::mutate_all(readr::parse_guess)
}

#' Get adimxture data
#'
#' \code{data_amdx} imports the admixture data
#'
#' @param gid       string, identifier of fst outlier ID (eg "LG04_1")
#' @param k         numeric, number of clusters
#' @param admx_path path to folder with admixture results
#'
#' @export
data_amdx <- function(gid, k, admx_path){
  q_file <- stringr::str_c(admx_path, "hapmap.",gid,".",k,".Q")
  ind_file <- stringr::str_c(admx_path, "pop.",gid,".",k,".txt")

  vroom::vroom(q_file, delim = " ",
               col_names = stringr::str_c("bin", stringr::str_pad(1:k, width = 2,pad = 0))) %>%
    dplyr::bind_cols(vroom::vroom(ind_file,
                                  col_names = c("id","spec","loc"))) %>%
    dplyr::mutate(id_nr = stringr::str_sub(id, 1, -7),
                  id_order = stringr::str_c(loc,spec,  id_nr),
                  gid = gid) %>%
    tidyr::pivot_longer(names_to = "bin", values_to = "prop", cols = contains("bin"))
}

#' Plot adimxture data
#'
#' \code{adm_plot} plots the admixture data
#'
#' @param gid_in string, identifier of fst outlier ID (eg "LG04_1")
#' @param data   tibble, input data
#'
#' @export
adm_plot <- function(data, gid_in){
  data %>%
    dplyr::left_join(sample_order) %>%
    dplyr::left_join(pheno_facet) %>%
    dplyr::filter(gid == gid_in) %>%
    ggplot2::ggplot(aes(x = factor(ord_nr))) +
    ggplot2::geom_bar(stat = "identity",
                      position = "stack",
                      aes(y = prop, fill = bin)) +
    ggplot2::geom_point(data = pheno_plot_data %>%
                          dplyr::filter(trait == gid_traits[[gid_in]]),
                        aes(y = -0.1,  fill = factor(phenotype)),
                        shape = 21) +
    ggplot2::scale_y_continuous(gid_labels[[gid_in]],
                                expand = c(0, 0),
                                breaks = c(-0.1, 0, 0.5, 1),
                                labels = c(trait_icons[[gid_in]], 0, 0.5, 1),
                                limits = c(-0.125, 1)) +
    ggplot2::scale_x_discrete(breaks = sample_order$ord_nr,
                              labels = sample_order$id) +
    ggplot2::scale_fill_manual(values = c(RColorBrewer::brewer.pal(4, "Greys")[2:3] %>%
                                            purrr::set_names(nm = c("bin01", "bin02")),
                                          `0` = "white", `1` = "black"), na.value = "gray") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 9),
                   legend.position = "none",
                   axis.title.y = ggplot2::element_text(vjust = -6),
                   axis.title.x = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggtext::element_markdown())
}

#' Add fish annotation
#'
#' \code{add_spec_drawing} adds a fish annotation to a plot
#'
#'
#' @param spec hamlet species
#' @param pos position on ggplot
#'
#' @export
add_spec_drawing <- function(spec, pos){
  plot_fish(short = spec, x = pos, y = .5,height = .9,width = 10)
}
