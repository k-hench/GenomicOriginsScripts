#' summarise model statistics
#'
#' \code{summarise_model} sumamrizes model statistics in a tibble row.
#'
#' @family Suppl. Figure SX
#'
#' @export
summarise_model <- function(data){
  data$mod %>%
    purrr::map(broom::glance) %>%
    bind_rows() %>%
    bind_cols(.,
              data$mod %>%
                purrr::map(broom::tidy) %>%
                bind_rows() %>%
                mutate(grp = (row_number()+1)%/%2) %>%
                select(grp,term,estimate) %>%
                spread(key = 'term',value = 'estimate') %>%
                select(-grp) %>%
                set_names(., nm = c('intercept', 'slope')))
}

#' Plot two fishes with location
#'
#' \code{plot_fishes_location} creates a ggplot annotation layer with two hamlets and a flag.
#'
#' This is used to annotate plots.
#'
#' @family Suppl. Figure SX
#'
#' @export
plot_fishes_location <- function (left, right,loc) {
  nr_left <- which(str_sub(hypo_img$spec,start = 1,end = 3) == left)
  nr_right <- which(str_sub(hypo_img$spec,start = 1,end = 3) == right)
  nr_loc <- which(str_sub(hypo_flag$geo,start = 1,end = 3) == loc)

  p <- ggplot() +
    annotation_custom(hypo_flag$flag[[nr_loc]], xmin = -0.28, xmax = 0.28, ymin = -Inf, ymax = Inf)+
    coord_fixed(xlim = c(-1, 1)) +
    theme_void() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(limits = c(-0.4, 0.38)) +
    annotation_custom(hypo_img$r[[nr_right]], xmin = -1, xmax = -0.05, ymin = -Inf, ymax = Inf) +
    annotation_custom(hypo_img$l[[nr_left]], xmin = 0.05, xmax = 1, ymin = -Inf, ymax = Inf)

  tibble(run = str_c(left, loc, '-', right, loc),
         grob = list(p %>% ggplotGrob()))
}

#' The custom grob geom
#'
#' This geom builds heavily on the answer by baptiste on the
#' tidiverse github forum:
#' https://github.com/tidyverse/ggplot2/issues/1399
#'
#' @family Suppl. Figure SX
#'
#' @export
geom_hypo_grob2 <- function(mapping = NULL,
                            data = NULL,
                            stat = "identity",
                            position = "identity",
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = FALSE,
                            ...) {
  layer(
    geom = hypo_geom_grob_custom2,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

hypo_geom_grob_custom2 <- ggproto(
  "hypo_geom_grob_custom2",
  Geom,
  setup_data = function(self, data, params) {
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    data
  },

  draw_group = function(data, panel_scales, coord) {
    vp <- grid::viewport(x=data$rel_x,
                         y=data$rel_y,
                         h = data$height,
                         width = data$width,
                         angle = data$angle)
    g <- grid::editGrob(data$grob[[1]], vp=vp)
    ggplot2:::ggname("geom_hypo_grob2", g)
  },

  required_aes = c("grob","rel_x","rel_y"),
  default_aes = list(height = 1, width = 1, angle = 0)
)
