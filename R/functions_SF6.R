#' reroot tree manually
#'
#' \code{root_manual} reroot a tree manually.
#'
#' @param object       tree object
#' @param outgroup     vector of nodes to define outgroup
#' @param node         node to root tree at
#' @param resolve.root logical
#' @param ...          catch-all parameter to allow excessive parameters through purrr::pmap
#'
#' @family Suppl Figure 6
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
#' @param data tree object
#'
#' @family Suppl Figure 6
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
#' @param loc      string, sample location (bel [Belize]/ hon [Honduras]/ pan [Panama])
#' @param file     input file
#' @param tree_dir directory with phylogenetic results
#' @param ...      catch-all parameter to allow excessive parameters through purrr::pmap
#'
#' @family Suppl Figure 6
#'
#' @export
get_tree <- function(loc, file, tree_dir, ...){
  tree <- read.tree(stringr::str_c(tree_dir, file))

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
#' @family Suppl Figure 6
#'
#' @export
StatTreeData <-  ggplot2::ggproto("StatTreeLabel", Stat,
                         required_aes = "node",
                         compute_group = function(data, scales) {
                           setup_tree_data(data)
                         }
)

#' backbone for ggplot layer for node support
#'
#' \code{geom_nodepoint_support} provides the backbone for a custom ggplot layer for node support.
#'
#' @param mapping     Aesthetic mapping; ggplot2::aes()
#' @param data        data frame
#' @param position    string, "identity"
#' @param na.rm       logical
#' @param show.legend logical
#' @param inherit.aes logical
#' @param ...         arguments passed to geom_point_support()
#'
#' @family Suppl Figure 6
#'
#' @export
geom_nodepoint_support <- function (mapping = NULL, data = NULL, position = "identity",
                                    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...)
{
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
