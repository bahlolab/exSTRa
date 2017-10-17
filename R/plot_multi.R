#' Multiple STR score ECDF plots
#' 
#' Allows plotting the STR scores of many loci, including direct to files.
#' 
#' 
#' @export
plot_multi <- function(strscore,
  prefix = "exSTRa_plot",
  dir = "images/", 
  plot_cols = NULL, # list of the color of samples of each locus. Or just a named vector of colours to use for all loci. 
  color_only = NULL, # list (name = locus) of samples to color at each locus
  plottypes = 1, # 1 to current device, 2 to a single PDF file with multiple pages, 3 to individual PDF files per locus.
  alpha_control = 0.5, 
  alpha_case = NULL,
  legend = TRUE, 
  legend_control = TRUE, 
  controls_label = "controls", 
  custom_legend = NULL, # a named (sample) vector of colors, for the legend
  ... # further plot arguments
) {
  
  
  
  # plot_many_str_score(mixed_cohorts_wgs, "STR-score-WGS-01-to-03", plot_cols, 
  # plottypes = 3, alpha_case = 0.5, alpha_control = 0.3, legend = FALSE)
  
  plot_many_str_score(strscore, 
    typename = prefix, 
    plot_cols, 
    loci = NULL, 
    color_only = color_only, 
    plottypes = plottypes, 
    dirbase = paste0(dir, "/"), 
    alpha_control = alpha_control, 
    alpha_case = alpha_case,
    legend = legend, 
    legend_control = legend_control, 
    controls_label = controls_label, 
    custom_legend = custom_legend, 
    ...)
  
}
