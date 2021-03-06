% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_f1.R
\name{hypo_anno_l_lwd}
\alias{hypo_anno_l_lwd}
\title{Hamlet Annotation (alt)}
\usage{
hypo_anno_l_lwd(
  species,
  xmin = -Inf,
  xmax = Inf,
  ymin = -Inf,
  ymax = Inf,
  lwd = 0.15,
  line_color = "black"
)
}
\arguments{
\item{species}{string, species name}

\item{xmin}{numeric, left border of hamlet image}

\item{xmax}{numeric, right border of hamlet image}

\item{ymin}{numeric, lower border of hamlet image}

\item{ymax}{numeric, upper border of hamlet image}

\item{lwd}{numeric, line width of hamlet drawing line elements}

\item{line_color}{string, color code, color of hamlet drawing line elements}
}
\description{
Add ggplot hamlet drawing  annotation with adjustable line width for line elements of drawing.
}
