#' Collect hamlet genes
#'
#' \code{get_genes} collects the genes within a outlier window.
#'
#' @family Table functions
#'
#' @export
get_genes <- function(gid, chrom, start, end, ...){
  xrange <- c(start, end)
  gfffile <- system.file("extdata", stringr::str_c("HP.annotation.named.",
                                                   chrom, ".gff.gz"), package = "hypogen")
  gff_filter <- list(seqid = chrom)
  data <- as.data.frame(rtracklayer::readGFF(gzfile(gfffile), filter = gff_filter)) %>%
    dplyr::mutate(Parent = as.character(Parent))

  mRNAs <- data %>% dplyr::filter(type == "mRNA", end > xrange[1],
                                  start < xrange[2]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(yl = row_number()%%4 +  2) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(gid = gid,
                  label = unlist(strsplit(tolower(Parentgenename), "_"))[1] %>%
                    str_replace(pattern = 'hpv1g000000', replacement = 'hp...') %>%
                    str_c('\\textit{', ., '}')) %>%
    dplyr::select(gid, label)
  mRNAs
}

#' Export latex tables
#'
#' \code{export_2_latex} exports a formated latex table.
#'
#' @family General functions
#'
#' @export
export_2_latex <- function(table, name){
  table %>%
    mutate(`\\\\\\hline` = '\\\\') %>%
    write_delim(.,path = '.tmp.tex',delim = '&')
  # clean last column
  n <- names(table) %>% length()

  # open latex table
  write_lines(str_c('\\begin{tabular}{',
                    str_c(rep(' c',n),collapse = ''),
                    ' }'), name)
  # add table body
  read_lines('.tmp.tex') %>%
    str_replace(.,'&\\\\$','\\\\') %>%
    write_lines(name, append = TRUE)

  # close latex table
  write_lines('\\end{tabular}', name, append = TRUE)
  # remove temporary file
  file.remove('.tmp.tex')
  message(str_c('Exportet latex table to "',name,'".'))
}



#' Wrap content in multirow
#'
#' \code{as_multirow} formats multirow content for latex tables.
#'
#' @family Table functions
#'
#' @export
as_multirow <- function(x,n){
  str_c('\\multirow{', n ,'}{*}{', x ,'}')
}
