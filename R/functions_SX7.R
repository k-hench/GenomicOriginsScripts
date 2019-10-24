#' reroot tree manually
#'
#' \code{root_manual} reroot a tree manually.
#'
#' @family Suppl. Figure SX
#'
#' @export
root_manual <- function(object, outgroup, node = NULL, resolve.root = TRUE, ...) {
  tree <- ape::root.phylo(object, outgroup = outgroup,
                          node = node, resolve.root = resolve.root, ...)

  if (Nnode(tree) != Nnode(object)) {
    return(tree)
  }

  attr(tree, "reroot") <- TRUE
  node_map <- reroot_node_mapping(object, tree)
  attr(tree, "node_map") <- node_map
  return(tree)
}

#' setup phylo tree data
#'
#' \code{setup_tree_data} sets up phylo tree data.
#'
#' @family Suppl. Figure SX
#'
#' @export
setup_tree_data <- function(data) {
  if (nrow(data) == length(unique(data$node)))
    return(data)
  data[match(unique(data$node), data$node),]
}

#' import phylogenetic tree
#'
#' \code{get_tree} imports a phylogenetic tree.
#'
#' @family Suppl. Figure SX
#'
#' @export
get_tree <- function(loc, mode, pttrn, tree_dir, ...){
  tree <- read.tree(stringr::str_c(tree_dir, 'all_samples.',loc,'.',mode,'.whg.SNP.tree'))

  if(loc == 'all'){
    pttrn <- '.*floflo'
    tree <- root_manual(tree,outgroup = which(tree$tip.label %in% c(tree$tip.label[stringr::str_detect(pattern = pttrn,string = tree$tip.label)])))
  }
  tibble::tibble(loc = loc,tree = list(tree))
}

#' ggplot layer for node support
#'
#' \code{StatTreeData} provides a custom ggplot layer for node support.
#'
#' @family Suppl. Figure SX
#'
#' @export
StatTreeData <-  ggproto("StatTreeLabel", Stat,
                         required_aes = "node",
                         compute_group = function(data, scales) {
                           setup_tree_data(data)
                         }
)

#' backbone for ggplot layer for node support
#'
#' \code{geom_nodepoint_support} provides the backbone for a custom ggplot layer for node support.
#'
#' @family Suppl. Figure SX
#'
#' @export
geom_nodepoint_support <- function (mapping = NULL, data = NULL, position = "identity",
                                    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...)
{
  library(magrittr)
  self_mapping <- aes_(node = ~node, subset = ~(!isTip))
  if (is.null(mapping)) {
    mapping <- self_mapping
  }
  else {
    if (!is.null(mapping$subset)) {
      self_mapping <- aes_string(node = "node", subset = paste0(as.expression(get_aes_var(mapping,
                                                                                          "subset")), "&!isTip"))
    }
    mapping %<>% modifyList(self_mapping)
  }
  geom_point_support(mapping, data, position, na.rm, show.legend,
                     inherit.aes, stat = StatTreeData, ...)
}

