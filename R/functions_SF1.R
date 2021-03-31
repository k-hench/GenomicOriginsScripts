#' Get hybridization data
#'
#' \code{getPofZ} imports hybridization data
#'
#' @param base_dir path containing the newhybrids results
#' @param folder   sub-folder of interest within the results folder
#'
#' @family Suppl Figure 1
#'
#' @export
getPofZ <- function(base_dir, folder){
  runname <- folder %>% str_remove("newHyb.") %>% str_remove(".80SNPs.txt_Results")

  pops <- c(str_sub(runname,1,6),str_sub(runname,-6,-1))
  result_dir <- str_c(base_dir, folder,"/")

  pofz <- dir(result_dir, pattern = "PofZ.txt")
  inds <- dir(result_dir, pattern = "_individuals.txt")

  colN <- c("P1", "P1_bc", "P2", "P2_bc")

  NHres <- vroom::vroom(str_c(result_dir, pofz),
                        delim = '\t',
                        skip = 1,
                        col_names  = c('indNR', 'IndivName', colN[1], colN[3],
                                       'F1', 'F2', colN[2], colN[4])) %>%
    mutate(IndivName = vroom::vroom(str_c(result_dir, inds),
                                    delim = ',',
                                    col_names =  c('IndivName'))[,1] %>%
             unname() %>%
             unlist() )

  bin_tib <- tibble(bin_generic = c(colN, "F1", "F2"),
                    bin = c(paste0(c(pops[1],pops[1],pops[2],pops[2]),
                                   c('_pure','_BC','_pure','_BC')), "F1", "F2"))

  data <- NHres %>%
    pivot_longer(names_to = 'bin_generic',
                 values_to= "post_prob",
                 cols = c(-indNR, -IndivName)) %>%
    left_join(bin_tib) %>%
    mutate(run = runname,
           loc = str_sub(run,-3,-1),
           ind_order = str_c(str_sub(IndivName,-6,-1),"_",str_sub(IndivName,1,-7)))

  return(data)
}

#' Plot hybridization data
#'
#' \code{getPofZ} plot hybridization data of a location
#'
#' @param loc sample location (bel [Belize]/ hon [Honduras]/ pan [Panama])
#'
#' @family Suppl Figure 1
#'
#' @export
plot_loc <- function(loc){
  data <- map_dfr(.x = folders[str_detect(folders, loc)],
                  .f = getPofZ,
                  base_dir = base_dir)

  is_hybr <- data %>%
    filter(!(grepl(pattern = "_pure", bin))) %>%
    filter(post_prob > .99) %>%
    .$ind_order %>% unique()

  data <- data %>%
    mutate(ind_label = ifelse(ind_order %in% is_hybr, str_c("**", ind_order, "**"), ind_order),
           run2 = str_c("*H. ", sp_names[str_sub(string = run,1,3)],"* - *H. ", sp_names[str_sub(string = run,8,10)],"*"))

  data_labs <- data %>% filter(!duplicated(ind_order)) %>% select(ind_order, ind_label)

  lvls <- c("P1", "P1_bc","F1", "F2", "P2_bc", "P2")

  clr <- paletteer_c("ggthemes::Red-Green-Gold Diverging",3) %>%
    c(.,clr_lighten(.)) %>% color() %>% .[c(1,4,2,5,6,3)] %>%
    set_names(nm = lvls)

  data  %>%
    ggplot(aes(x = ind_order, y = post_prob, fill = bin_generic))+
    geom_bar(position = 'stack', stat = "identity")+
    scale_fill_manual(values = clr) +
    scale_x_discrete(breaks = data_labs$ind_order, labels = data_labs$ind_label)+
    labs(y = "Posterior probability", title = loc_names[loc])+
    facet_grid(run2~.)+
    theme_minimal()+
    theme(legend.position = "bottom",
          strip.text.x  = element_markdown(),
          strip.text.y  = element_markdown(),
          axis.text.x = element_markdown(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 4))
}

#' Custom ggplot theme
#'
#'\code{theme_hyb} ggplot theme for hybridization plots
#'
#' @param legend.position string, according to ggplot2::theme()
#' @param ...             arguments funneled through to ggplot2::theme()
#'
#' @family Suppl Figure 1
#'
#' @export
theme_hyb <-  function(legend.position = "none",...){
  list(scale_y_continuous(breaks = c(0,.5,1)),
       theme(legend.position = legend.position,
             legend.background = element_rect(fill = "white",colour = rgb(1,1,1,0)),
             legend.direction = "horizontal",
             legend.justification = c(1,1),
             strip.text.y = element_markdown(angle = 0,hjust = 0),
             ...))
}

#' Adjust plot sizes
#'
#' \code{label_spacer} adjusts plot sizes with respect to presence of labels
#'
#' @param x    label position
#' @param plus label offset
#'
#' @family Suppl Figure 1
#'
#' @export
label_spacer <- function(x, plus = 1.1){x + plus}
