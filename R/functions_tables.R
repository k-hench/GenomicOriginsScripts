#' Collect hamlet genes
#'
#' \code{get_genes} collects the genes within a outlier window.
#'
#' @param gid   string, identifier of fst outlier ID (eg "LG04_1")
#' @param chrom string, linkage group identifier (eg "LG04")
#' @param start numeric, start position of window (bp)
#' @param end   numeric, end position of window (bp)
#' @param ...   catch-all parameter to allow excessive parameters through purrr::pmap
#'
#' @family Table functions
#'
#' @export
get_genes <- function(gid, chrom, start, end, ...){
  xrange <- c(start, end)
  gfffile <- system.file("extdata",
                         stringr::str_c("HP.annotation.named.",
                                        chrom, ".gff.gz"),
                         package = "hypogen")
  gff_filter <- list(seqid = chrom)
  data <- as.data.frame(rtracklayer::readGFF(gzfile(gfffile), filter = gff_filter)) %>%
    dplyr::mutate(Parent = as.character(Parent))

  mRNAs <- data %>%
    dplyr::filter(type == "mRNA", end > xrange[1],
                  start < xrange[2]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(yl = dplyr::row_number()%%4 +  2) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(gid = gid,
                  label = unlist(strsplit(tolower(Parentgenename), "_"))[1] %>%
                    stringr::str_replace(pattern = 'hpv1g000000', replacement = 'hp...') %>%
                    stringr::str_c('\\textit{', ., '}')) %>%
    dplyr::select(gid, label)
  mRNAs
}

#' Export latex tables
#'
#' \code{export_2_latex} exports a formatted latex table.
#'
#' @param table tibble, input table
#' @param name  string, name for exported file
#'
#' @family General functions
#'
#' @export
export_2_latex <- function(table, name){
  table %>%
    dplyr::mutate(`\\\\\\hline` = '\\\\') %>%
    readr::write_delim(., path = '.tmp.tex',delim = '&')
  # clean last column
  n <- names(table) %>% length()

  # open latex table
  readr::write_lines(stringr::str_c('\\begin{tabular}{',
                                    stringr::str_c(rep(' c',n), collapse = ''),
                                    ' }'), name)
  # add table body
  readr::read_lines('.tmp.tex') %>%
    stringr::str_replace(.,'&\\\\$','\\\\') %>%
    readr::write_lines(name, append = TRUE)

  # close latex table
 readr::write_lines('\\end{tabular}', name, append = TRUE)
  # remove temporary file
  file.remove('.tmp.tex')
  message(stringr::str_c('Exportet latex table to "',name,'".'))
}

#' Wrap content in multirow
#'
#' \code{as_multirow} formats multirow content for latex tables.
#'
#' @param x string, cell content
#' @param n integer, number of rows to combine
#'
#' @family Table functions
#'
#' @export
as_multirow <- function(x, n){
  stringr::str_c('\\multirow{', n ,'}{*}{', x ,'}')
}
