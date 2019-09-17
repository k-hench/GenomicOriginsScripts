#' Process script input
#'
#' \code{process_input} processes the script input and echos the received parameters.
#'
#' @family General functions
#'
#' @export
process_input <- function(script_name, args){
  cat(str_c(cli::rule( left = str_c(crayon::bold('Script: '),crayon::red(script_name))),'\n'))
  args_out <- args[7:length(args)]
  cat(str_c(crayon::green('Parameters read:'),'\n'))
  cat(' ')
  cat(str_c(crayon::green(cli::symbol$star),' ', 1:length(args_out),': ',crayon::green(args_out),'\n'))
  cat(str_c(cli::rule(right = getwd()),'\n'))
  return(args_out)
}

#' Add genomic position
#'
#' \code{add_gpos} adds the genomic position to a tibble of genomic stats.
#'
#' The imported data typically comes with one column specifieing the linkage group (LG)
#' and one specifieng the position of that LG.
#' This function translates the postion into the genomic postion where
#' all LGs have been concatinated (eg. LG02, POS:1 => length(LG01)+1)
#'
#' @export
add_gpos <- function(tib, ...){
  tib %>%
    left_join(hypogen::hypo_chrom_start) %>%
    mutate(GPOS = GSTART+(BIN_START+BIN_END)/2)
}

#' Import Fst data
#'
#' \code{get_fst} imports Fst data, adds the genomic position and the run.
#'
#' The run name (the pair wise species comparison) is extracted from the file name,
#' the data is imported, column names are standardized and genomic position and
#' run name are added.
#'
#' @family General functions
#'
#' @export
get_fst <- function(file){
  run <- file %>%
    str_extract(.,'[a-z]{3}-[a-z]{3}-[a-z]{3}.50k') %>%
    str_remove(.,'.50k') %>%
    str_replace(.,'([a-z]{3})-([a-z]{3})-([a-z]{3})','\\2\\1-\\3\\1')

  read_tsv(file) %>%
    add_gpos() %>%
    mutate(run = run)
}

#' Import dxy data
#'
#' \code{get_dxy} imports dxy data, adds the genomic position and the run.
#'
#' The run name (the pair wise species comparison) is extracted from the file name,
#' the data is imported, column names are standardized and genomic position and
#' run name are added.
#'
#' @family General functions
#'
#' @export
get_dxy <- function(file){
  run <- file %>%
    str_extract(.,'[a-z]{6}-[a-z]{6}.50kb') %>%
    str_remove(.,'.50kb')

  read_csv(file) %>%
    setNames(., nm = c("CHROM","BIN_START","BIN_END", "BIN_MID", "N_SITES",
                       "PI_POP1", "PI_POP2", "dxy", "Fst") ) %>%
    add_gpos() %>%
    mutate(run = run)
}

#' Import genotype x phenotype association
#'
#' \code{get_gxp} imports genotype x phenotype association data and adds the genomic position.
#'
#' The phenotype trait name is extracted from the file name and
#' converted into a label for plotting.
#' Then, the data is imported and reduced to the columns containg
#' the genomic poistion and the p values of the wald test.
#' The column containing the p values is renamed according to the trait.
#'
#' @family General functions
#'
#' @export
get_gxp <- function(file){

  trait <- file %>%
    str_remove('^.*/') %>%
    str_remove('\\..*')

  trait_label <- str_c(trait_panels[trait],':italic(p)~(GxP[',trait,'])')

  file %>%
    read_tsv() %>%
    add_gpos() %>%
    select(GPOS, AVG_p_wald) %>%
    setNames(., nm = c('GPOS',trait_label))
}

#' Iterate left_join over list
#'
#' \code{join_list} binds a list of tables according to a common column.
#'
#' @family General functions
#'
#' @export
join_list <- function(lst){
  lst %>% purrr::reduce(left_join)
}
