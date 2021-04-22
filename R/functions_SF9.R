#' Convert tree labels
#'
#' @param label string, sample id
#'
#' @export
lab2spec <- function(label){
  x <- stringr::str_sub(label, start = -6, end = -4) %>% stringr::str_remove(.,"[0-9.]{1,3}$") %>% stringr::str_remove(.," ")
  ifelse(x == "",'ungrouped', x)
}
