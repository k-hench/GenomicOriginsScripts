#' Get adimxture data
#'
#' \code{data_amdx} imports the admixture data
#'
#'
#' @family Figure SX
#'
#' @export
data_amdx <- function(gid, k, admx_path){
  q_file <- str_c(admx_path, "hapmap.",gid,".",k,".Q")
  ind_file <- str_c(admx_path, "pop.",gid,".",k,".txt")

  vroom::vroom(q_file, delim = " ",
               col_names = str_c("bin", str_pad(1:k,width = 2,pad = 0))) %>%
    bind_cols(vroom::vroom(ind_file,
                           col_names = c("id","spec","loc"))) %>%
    mutate(id_nr = str_sub(id,1,-7),
           id_order = str_c(loc,spec,  id_nr),
           gid = gid) %>%
    pivot_longer(names_to = "bin", values_to = "prop", cols = contains("bin"))
}

#' Plot adimxture data
#'
#' \code{adm_plot} plots the admixture data
#'
#'
#' @family Figure SX
#'
#' @export
adm_plot <- function(gid_in){
  data %>%
    left_join(sample_order) %>%
    left_join(pheno_facet) %>%
    filter(gid == gid_in) %>%
    ggplot(aes(x = factor(ord_nr)))+
    geom_bar(stat = 'identity',
             position = "stack",
             aes (y = prop, fill = bin))+
    geom_point(data = pheno_plot_data %>%
                 filter(trait == gid_traits[[gid_in]]),
               aes(y = -.075,
                   fill = factor(phenotype)),
               shape = 21)+
    scale_y_continuous(gid_labels[[gid_in]],
                       expand = c(0,0),
                       breaks = c(-.075,0,.5,1),
                       labels = c(gid_traits[[gid_in]],0,.5,1),
                       limits = c(-.125,1))+
    scale_x_discrete(breaks = sample_order$ord_nr,
                     labels = sample_order$id) +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(4,"Greys")[2:3] %>%
                                   set_names(nm = c("bin01", "bin02")),
                                 `0` = "white", `1` = "black"),
                      na.value = "gray")+
    theme_minimal()+
    theme(plot.title = element_text(size = 9),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank())
}

#' Add fish annotation
#'
#' \code{add_spec_drawing} adds a fish annotation to a plot
#'
#'
#' @family Figure SX
#'
#' @export
add_spec_drawing <- function(spec, pos){
  plot_fish(short = spec,x = pos, y = .5,height = .9,width = 10)
}
