% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_f5_cross_topo.R
\name{get_isolated_topos}
\alias{get_isolated_topos}
\title{get isolated topologies}
\usage{
get_isolated_topos(topo_plots, left_neighbour, pops)
}
\arguments{
\item{topo_plots}{tibble containing the topologies (in column 'phylo')}

\item{left_neighbour}{string (species identifier for left neighbor)}

\item{pops}{all populations within the topology}
}
\description{
\code{is_isolated} finds all topologies where a species is isolated
}
\details{
First the distance matrixes of all topologies are created and attached
to the topologies.
Then, the isolation status of a species is determined based on the
distance matrixes and the topologies are filtered for the isolation cases.

The index number of those topologies is exported.
}
\seealso{
Other Figure 5: 
\code{\link{ax_scl}()},
\code{\link{cross_spec}()},
\code{\link{custom_annoplot}()},
\code{\link{dist_tibble}()},
\code{\link{distances_tree}()},
\code{\link{get_clr}()},
\code{\link{get_dist}()},
\code{\link{get_gxp_long}()},
\code{\link{get_neighbour_topos}()},
\code{\link{is_isolated}()},
\code{\link{min_dist}()},
\code{\link{no_title}()},
\code{\link{plot_curtain}()},
\code{\link{plot_fish_zoom}()},
\code{\link{plot_fst_poptree}()},
\code{\link{plot_leg}()},
\code{\link{plot_panel_anno}()},
\code{\link{plot_panel_delta_dxy}()},
\code{\link{plot_panel_dxy}()},
\code{\link{plot_panel_fst}()},
\code{\link{plot_panel_gxp}()},
\code{\link{plot_panel_twisst}()},
\code{\link{prep_data}()},
\code{\link{theme_panels}()}
}
\concept{Figure 5}
