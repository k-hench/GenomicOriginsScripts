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
#' @param tib input tibble
#' @param ... catchall input for purrr::pmap
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
get_fst <- function(file, kb = "50k"){
  run <- file %>%
    str_extract(., str_c('[a-z]{3}-[a-z]{3}-[a-z]{3}.',kb)) %>%
    str_remove(., str_c('.',kb)) %>%
    str_replace(.,'([a-z]{3})-([a-z]{3})-([a-z]{3})','\\2\\1-\\3\\1')

  #read_tsv(file) %>%
  vroom::vroom(file, delim = "\t") %>%
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
get_dxy <- function(file, kb = "50k"){
  run <- file %>%
    str_extract(., str_c('[a-z]{6}-[a-z]{6}.',kb)) %>%
    str_remove(., str_c('.',kb))

  #read_csv(file) %>%
  vroom::vroom(file, delim = ",") %>%
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

  trait_label <- str_c(trait_panels[trait],':italic(p)[',trait,']')

  vroom::vroom(file, delim = "\t") %>%
    add_gpos() %>%
    select(GPOS, AVG_p_wald) %>%
    setNames(., nm = c('GPOS',trait_label))
}


#' Import diversity data
#'
#' \code{get_pi} imports diversity data and adds the genomic position.
#'
#' The data is imported and genomic position and species name are added.
#'
#' @family General functions
#'
#' @export
get_pi <- function(file){
  vroom::vroom(file, delim = "\t") %>%
    left_join(., hypogen::hypo_chrom_start) %>%
    mutate(gpos = (BIN_START + BIN_END)/2 + GSTART,
           spec = file %>% str_remove('^.*/') %>% str_sub(.,1,6))
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

#' Darken a given color
#'
#' \code{darken} makes a darker version of a given color.
#'
#' @family General functions
#'
#' @export
darken <- function(color, factor=.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

#' Lighten a given color
#'
#' \code{lighten} makes a lighter version of a given color.
#'
#' @family General functions
#'
#' @export
lighten <- function(color, factor=.4){
  col <- col2rgb(color)
  col <- col + (255-col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

#' Reorder factor levels
#'
#' \code{refactor} redorders levels of a factor levels based on an external table.
#'
#' @family General functions
#'
#' @export
refactor <- function(self,globals){
  factor(as.character(self$run),
         levels = c(levels(globals$run)))
}

#' Import fst summary data for outlier windows
#'
#' \code{get_fst_summary_data} imports pair wise fst summary data for outlier windows.
#'
#'
#' @family General functions
#'
#' @export
get_fst_summary_data <- function(gid){
  file <- str_c(fst_summary_path, gid, ".fst.all.tsv")
  data <- read_tsv(str_c(file)) %>%
    mutate(weighted_fst = ifelse(weighted_fst < 0, 0, weighted_fst))

  data_w <- data %>%
    select(run, weighted_fst) %>%
    separate(run, into = c("pop1", "pop2"), sep = "-") %>%
    bind_rows(tibble( pop1 = c("abehon", "unipan"), pop2 = pop1)) %>%
    arrange(pop1, pop2) %>%
    pivot_wider(names_from = pop2,
                values_from = weighted_fst)

  data_m <- data_w %>% select(abehon:unipan) %>% as.matrix()
  rownames(data_m) <- data_w$pop1

  data_m_t <- data_m %>% t()

  data_m2 <- data_m
  data_m2[lower.tri(data_m2)] <- data_m_t[lower.tri(data_m_t)]
  data_m2[is.na(data_m2)] <- 0
  data_nj <- nj(data_m2)

  tibble(gid = gid,
         data = list(data),
         data_m2 = list(data_m2),
         data_nj = list(data_nj))
}
