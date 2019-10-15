#' import sc data sheet
#'
#' \code{read_sc} imports a sc data sheet.
#'
#' @family Suppl. Figure SX2
#'
#' @export
read_sc <-function(file){
  length_comment <- read_lines(file) %>% str_which(.,'^#') %>% max()
  length_format <- read_lines(file) %>% str_which(.,'^format*') %>% max()
  pos_frame <- read_lines(file) %>% str_which(.,'^frame*') %>% max()

  head_length <- max(length_comment, length_format, pos_frame)

  read_delim(file,delim = ' ',skip = head_length,
             col_names = c('type','cell','delim','value')) %>%
    filter(!(type == "goto"), !(is.na(value))) %>%
    select(-type,-delim) %>%
    mutate(cell = str_replace(cell,"^([A-Z])([0-9])","\\1_\\2")) %>%
    separate(cell, into = c('column','row'),sep = '_') %>%
    mutate(row = as.numeric(row)) %>%
    spread(key = column, value = value) %>%
    setNames(.,nm = .[1,] %>% as.character()) %>%
    filter(`0` != 0) %>%
    select(-`0`) %>%
    mutate_all(parse_guess)
}
