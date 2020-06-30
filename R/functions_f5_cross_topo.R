#' Cross species pairs
#'
#' \code{cross_spec} creates a cross of all possible species pairs
#'
#' @family Figure 5
#'
#' @export
cross_spec <- function(grph){
  test <- names(V(grph)) %>% .[!str_detect(.,'Node.*')]
  flter <- function(x, y) x >= y
  purrr::cross_df(tibble(x=test,y=test), .filter = flter) %>% arrange(x)
}

#' Create distance of species pair
#'
#' \code{get_dist} creates a tibble containg distance of a single species pair within a topology
#'
#' @family Figure 5
#'
#' @export
get_dist <- function(grph,x,y){
  tibble( out = distances(grph,v = x, to = y)[1] ) %>%
    set_names(.,nm=str_c(x,'-',y))
}


#' Create tibble of pairwise species distances
#'
#' \code{dist_tibble} creates a tibble of all pairwise species distances for a topology
#'
#' @family Figure 5
#'
#' @export
dist_tibble <- function(phyl){
  grph <- igraph::as.igraph(phyl %>% ape::unroot())
  cross_spec(grph) %>%
    purrr::pmap(get_dist,grph=grph) %>%
    bind_cols()
}

#' Get pairwise highlighing topologies
#'
#' \code{get_neighbour_topos} finds all topologies where two species are sister species
#'
#' First the distances of all species within all topologies are created.
#' This is subset for the topologies where the distance between the two species equals 2 (2 edges).
#'
#' @family Figure 5
#'
#' @export
get_neighbour_topos <- function(topo_plots,left_neighbour,right_neighbour){
  topo_plots$phylo %>%
    purrr::map(dist_tibble) %>%
    bind_rows() %>%
    bind_cols(topo_plots,.) %>%
    filter( .data[[str_c(left_neighbour,'-',right_neighbour)]] == 2) %>%
    select(topo_nr) %>%
    unlist() %>%
    unname()
}

#' Get distance matrix of a topology
#'
#' \code{distances_tree} unroots a topology and creates a distance matrix from it
#'
#' @family Figure 5
#'
#' @export
distances_tree <- function(tree){
  tree_dist <-  tree %>%
    ape::unroot() %>%
    as.igraph() %>%
    distances() %>%
    as.matrix()
}

#' Find minimum distance to other species in topology
#'
#' \code{min_dist} finds the minimum distance of a focal species to all other species in a topology
#'
#' @family Figure 5
#'
#' @export
min_dist <- function(dist, pop, pops){
  other_pops <- pops[!(pops == pop)]
  tibble(min_dist = min(dist[pop,other_pops]))
}

#' Determine if species is isolated in topology
#'
#' \code{is_isolated} determines if a species has a sister species within a topology
#'
#' @family Figure 5
#'
#' @export
is_isolated <- function(x, pop, pops){
  purrr::map(.x = x,
             .f = min_dist,
             pop = pop,
             pops = pops) %>%
    bind_rows() %>%
    .$min_dist == 3
}

#' get isolated topologies
#'
#' \code{is_isolated} finds all topologies where a species is isolated
#'
#' First the distance matrixesof all topologies are created and attached
#' to the topologies.
#' Then, the isolation status of a species is determined based on the
#' distance matrixes and the topologies are filtered for the isolation cases.
#'
#' The index number of those topologies is exported.
#'
#' @family Figure 5
#'
#' @export
get_isolated_topos <- function(topo_plots, left_neighbour, pops){
  topo_plots$phylo %>%
    purrr::map(distances_tree) %>%
    tibble(dist = .) %>%
    bind_cols(topo_plots,.) %>%
    mutate(isolated = is_isolated(dist, pop = left_neighbour, pops = pops)) %>%
    filter(isolated) %>%
    select(topo_nr) %>%
    unlist() %>%
    unname()
}
