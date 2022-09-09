#' Multiple STR score ECDF plots
#' 
#' Allows plotting the STR scores of many loci with legend.
#' Plots may be written directly to PDF files. 
#' 
#' @param strscore exstra_score object.
#' @param prefix File name prefix.
#' @param dir Output directory for PDF plots. 
#' @param plot_cols If a named vector, colours to use for each sample (names) for all loci. 
#'                  Unspecified samples are treated as controls. 
#'                  May be a named (loci) list with named vectors as above for specific
#'                  locus colours.  
#' @param color_only List (named with loci) of a character vectors of samples to color at each locus.
#' @param plot_types vector; Specify the plot output device. May have multiple of the 
#'          following options if these are desired. 
#'          Set to: 1 for current device, 
#'          2 for output to a single PDF file with multiple pages (when >1 locus), and
#'          3 for a PDF file per locus. 
#' @param alpha_control Transparency alpha value for control samples. 
#' @param alpha_case Transparency alpha value for case samples. 
#' @param legend logical; if TRUE, include a legend.
#' @param legend_control logical; if TRUE, include generic control specification in legend.
#' @param controls_label character; a generic label in the legend for controls. 
#' @param cases_label character; a generic label for case samples in the legend. 
#' @param custom_legend A named (sample) vector of colors, for the legend.
#' @param ... Further arguments to \code{\link{plot.exstra_score}}
#' 
#' @examples 
#' data(exstra_wgs_pcr_2)
#' # Plot 4 loci
#' par(mfrow = c(2, 2))
#' plot_multi(exstra_wgs_pcr_2[c("HD", "SCA6", "FRDA", "SCA1")], alpha_case = 0.2)
#' 
#' @importFrom graphics par
#' @export
plot_multi <- function(strscore,
  prefix = "exSTRa_plot",
  dir = "images/", 
  plot_cols = NULL, # list of the color of samples of each locus. Or just a named vector of colours to use for all loci. 
  color_only = NULL, # list (name = locus) of samples to color at each locus
  plot_types = 1, # 1 to current device, 2 to a single PDF file with multiple pages, 3 to individual PDF files per locus.
  alpha_control = 0.5, 
  alpha_case = NULL,
  legend = TRUE, 
  legend_control = TRUE, 
  controls_label = "controls", 
  cases_label = "cases",
  custom_legend = NULL, # a named (sample) vector of colors, for the legend
  ... # further plot arguments
) {
  
  # plot_many_str_score(mixed_cohorts_wgs, "STR-score-WGS-01-to-03", plot_cols, 
  # plot_types = 3, alpha_case = 0.5, alpha_control = 0.3, legend = FALSE)
  
  plot_many_str_score(strscore, 
    typename = prefix, 
    plot_cols, 
    loci = NULL, 
    color_only = color_only, 
    plot_types = plot_types, 
    dirbase = paste0(dir, "/"), 
    alpha_control = alpha_control, 
    alpha_case = alpha_case,
    legend = legend, 
    legend_control = legend_control, 
    controls_label = controls_label, 
    cases_label = cases_label,
    custom_legend = custom_legend, 
    ...)
  
}
