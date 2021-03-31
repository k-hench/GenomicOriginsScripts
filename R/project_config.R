#' Hamlet colors (standard)
#'
#' \code{clr} defines the standard colors for hamlets species within plots.
#'
#' This is used in most color scales thoughout the paper
#' where individual species are indicated by color.
#'
#' @export
clr <- c(
  abe = "#E5E5A1",
  flo = "#ABA7C4",
  gum = "#E3A258",
  ind = "#1C6FCE",
  may = "#7EA7C2",
  nig = "#000000",
  pue = "#E17366",
  ran = "#7EBDB3",
  tab = "#A1C75E",
  tor = "#E2B8CE",
  uni = "#ffffff"
)

#' Hamlet colors (alternative)
#'
#' \code{clr2} defines the alternative colors for hamlets species within plots.
#'
#' This is used in color scales where cloring H. unicolor
#' white does not work.
#'
#' @export
clr2 <- c(
  abe = "#E5E5A1",
  flo = "#ABA7C4",
  gum = "#E3A258",
  ind = "#1C6FCE",
  may = "#7EA7C2",
  nig = "#000000",
  pue = "#E17366",
  ran = "#7EBDB3",
  tab = "#A1C75E",
  tor = "#E2B8CE",
  uni = "#CCCCCC"
)

#' Default plot color
#'
#' \code{plot_clr} is the default color of data in plots.
#'
#' This color is used to plot neutral data thoughout the figures.
#'
#' @export
plot_clr <- rgb(.2,.2,.2)

#' Default background bar color
#'
#' \code{clr_below} is the default background bar color.
#'
#' @export
clr_below <- 'lightgray'#rgb(.8,.8,.8)

#' Default gene color
#'
#' \code{clr_genes} is the default color of genes zoom plot.w
#'
#' @export
clr_genes <-rgb(.4,.4,.4)

#' Default plot size
#'
#' \code{plot_size} sets the default size of data points in plots.
#'
#' This sets the default size of data points thoughout the figures.
#'
#' @export
plot_size <- .2


#' Default plot linewidth
#'
#' \code{plot_lwd} sets the linewidth in plots.
#'
#' This sets the default linewidth thoughout the figures.
#'
#' @export
plot_lwd <- .25

#' Default plot font size
#'
#' \code{plot_text_size} sets the default font size of text in plots.
#'
#' This sets the default font size of text thoughout the figures.
#'
#' @export
plot_text_size <- 7

#' Small plot font size
#'
#' \code{plot_text_size_small} sets the small font size of text in plots.
#'
#' This sets the small font size of text for the figures (legends).
#'
#' @export
plot_text_size_small <- 6

#' Outlier color
#'
#' \code{outlr_clr} defines the color used to highlight Fst outlier windows in plots.
#'
#' This defines the color used to highlight Fst outlier windows in plots.
#'
#' @export
outlr_clr <- rgb(1,0,0,.2)

#' Location names
#'
#' \code{loc_names} stores the geographic names of the sample abrreviations.
#'
#' This is used to translate between the 3-letter sample abrreviations and
#' the actual geographic names of the sample locations.
#'
#' @export
loc_names <- c(
	bel = "Belize",
	hon = "Honduras",
	pan = "Panama",
	flo = "Florida"
)

#loc_labs <- c( bel = "Belize",
#               flo = "Florida",
#               hon = "Honduras",
#               pan = "Panama")

#' Location colors
#'
#' \code{clr_loc} defines the colors for sample locations within plots.
#'
#' This is used in color scales thoughout the paper
#' where the locations are indicated by color.
#'
#' @export
clr_loc = c(
	bel = "#E41A1C",
	hon = "#377EB8",
	pan = "#4DAF4A",
	flo = "#984EA3"
)


#' Location shapes
#'
#' \code{shps} defines the shapes for sample locations within plots.
#'
#' This is used in shape scales thoughout the paper
#' where the locations are indicated by shape.
#'
#' @export
shps <- c(bel = 21,
          hon = 23,
          pan = 22,
          flo = 24)

#' Fill function for adxiture classes
#'
#' \code{fll_fun} creates a fill scale for admixture classes.
#'
#' This is used to keep a consistent color scheme thoughout
#' various admixture plots of different k.
#'
#' @param n integer, number of clusters
#'
#' @export
fll_fun <- function(n){viridis::inferno(n)}

#' Named fill function for adxiture classes
#'
#' \code{fll_n} creates a named fill scale for admixture classes.
#'
#' This is used to keep a consistent color scheme thoughout
#' various admixture plots of different k.
#' The named variant is used to assign the colors consistently.
#'
#' @param n integer, number of clusters
#'
#' @export
fll_n <- function(n){ fll_fun(n) %>% setNames(., nm = stringr::str_c("pop_", 1:n)) }

#' Species labels
#'
#' \code{sp_labs} stores the scientific names of the hamlet abrreviations.
#'
#' This is used to translate between the 3-letter species abrreviations and
#' the actual formated scientific names of the hamlet species.
#'
#' @export
sp_labs <- c(
  abe = expression(italic(H.~aberrans)),
  flo = expression(italic(H.~floridae)),
  gum = expression(italic(H.~gummigutta)),
  ind = expression(italic(H.~indigo)),
  may = expression(italic(H.~maya)),
  nig = expression(italic(H.~nigricans)),
  pue = expression(italic(H.~puella)),
  ran = expression(italic(H.~randallorum)),
  tab = expression(italic(S.~tabacarius)),
  tor = expression(italic(S.~tortugarum)),
  uni = expression(italic(H.~unicolor))
)

#' Species names
#'
#' \code{sp_names} stores the species names of the hamlet abrreviations.
#'
#' This is used to translate between the 3-letter species abrreviations and
#' the actual species names of the hamlet species.
#'
#' @export
sp_names <- c(
  abe = "aberrans",
  flo = "floridae",
  gum = "gummigutta",
  ind = "indigo",
  may = "maya",
  nig = "nigricans",
  pue = "puella",
  ran = "randallorum",
  tab = "tabacarius",
  tor = "tortugarum",
  uni = "unicolor"
)

#' Population order
#'
#' \code{pop_levels} provides the hierachic order of the populations.
#'
#' The sampled populatinos are alphabetically ordered, first
#' by location, then by species.
#'
#' @export
pop_levels <- c("indbel", "maybel", "nigbel", "puebel", "unibel",
                "abehon", "gumhon", "nighon", "puehon", "ranhon", "unihon",
                "nigpan", "puepan", "unipan")

#' Project case
#'
#' \code{project_case} manages the format of labels for figure panels.
#'
#' Currently all figure sub panels are labeled using lower case.
#' This function potentially needs to be changed in case of
#' resubmission to a different journal.
#'
#' @param x string
#'
#' @export
project_case <- function(x){
	str_to_lower(x)
}

#' Inverted project case
#'
#' \code{project_inv_case} manages the format of labels within figures.
#'
#' Labels within figures (outlier highlights) should be in inverse case
#' compared to figure panel label to avoid confusion.
#' This function potentially needs to be changed in case of
#' resubmission to a different journal.
#'
#' @param x string
#'
#' @export
project_inv_case <- function(x){
  str_to_upper(x)
}


#' Default figure width (double column)
#'
#' \code{f_width} defines the full width for figures (mm).
#'
#' 7 inch = 178mm
#'
#' @export
f_width <- 7

#' Default figure width (single column)
#'
#' \code{f_width_half} defines the half width for figures (mm).
#'
#' 3.43 inch = 87mm
#'
#' @export
f_width_half <- 3.43
