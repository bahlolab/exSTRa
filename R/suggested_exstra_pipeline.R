# create graphs of the suggested pipeline

#' Create a graph of the suggested pipeline
#' 
#' @param file Output SVG file name, if required.
#' @param show.with.file If output to a file, determines whether to display in R also. Has no effect when file is NULL.
#' @return The loci of \code{x}.
#' @examples 
#' \dontrun{
#' suggested_exstra_pipeline()
#' }
#' @seealso \pkg{DiagrammeR}
#' @seealso \pkg{DiagrammeRsvg}
#' @importFrom methods show
#' @export
suggested_exstra_pipeline <- function(file = NULL, show.with.file = FALSE) {
  gv <- DiagrammeR::grViz(system.file("extdata", "exSTRa_pipeline_graph.gv", package = "exSTRa"))
  if(!is.null(file)) {
    write(DiagrammeRsvg::export_svg(gv), file)
    if(show.with.file) {
      show(gv)
    }
  } else {
    show(gv)
  }
}



