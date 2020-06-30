#' Import msmc data
#'
#' \code{get_msmc} imports msc data.
#'
#' The run name (the pair wise species comparison) is extracted from the file name,
#' the data is imported and the run name is added.
#'
#' @family Figure 2
#'
#' @export
get_msmc <- function(file, msmc_path, mu = 3.7e-8, gen = 1){
  vroom::vroom(str_c(msmc_path, file), delim = '\t') %>%
    gather("Side", "time_value", 2:3) %>%
    arrange(time_index) %>%
    mutate(YBP=time_value/mu*gen,
           Ne = (1/lambda)/mu,
           run_nr = str_replace(file, pattern = 'run([0-9]*).*', replacement = '\\1') %>% as.numeric(),
           spec = str_replace(file, pattern = 'run[0-9]*\\.([a-z]*).*', replacement = '\\1'),
           loc = str_replace(file, pattern = 'run[0-9]*\\.[a-z]{3}\\.([a-z]*).*', replacement = '\\1'),
           run = str_c(spec, loc))
}

#' Import cross coalescence data
#'
#' \code{get_msmc} imports cross coalescence data and extracts according run information.
#'
#' The data is imported and run components are extracted from the a reference table,
#' and added to the data.
#'
#' @family Figure 2
#'
#' @export
get_cc <- function(file, cc_groups, cc_path, mu = 3.7e-8, gen = 1){
  cc_run <- str_replace(file, pattern = 'cc_run\\.([0-9]*).*', replacement = '\\1') %>%
    as.numeric()
  specs <- c(cc_groups$spec_1[cc_groups$run_nr == cc_run],
             cc_groups$spec_2[cc_groups$run_nr == cc_run]) %>%
    sort()
  loc <- cc_groups$geo[cc_groups$run_nr == cc_run]

  vroom::vroom(str_c(cc_path, file), delim = '\t') %>%
    gather("Side", "time_value", 2:3) %>%
    arrange(time_index) %>%
    mutate(YBP=time_value/mu*gen,
           Cross_coal = 2 * lambda_01 / (lambda_00 + lambda_11),
           run_nr = cc_run,
           spec_1 = specs[1],
           spec_2 = specs[2],
           loc = loc,
           run = str_c(spec_1,loc,'-',spec_2,loc))
}
