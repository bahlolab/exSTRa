#' Read repeat score data produced from Bio::STR::exSTRa
#' 
#' \code{read_score} reads output from the Perl module Bio::STR::exSTRa.
#' 
#' @param file Path of file to read.
#' @param database Path of database file (character) or exstra database (exstra_db object).
#' @param groups.regex Sets the sample groups with a regular expression. 
#'                     This is a named character vector of R regular expressions,
#'                     where names correspond to the group and regular expressions to 
#'                     match the sample names. 
#'                     The first match in the vector takes priority over other matches, 
#'                     such that it may be sensible to have the final regex as the empty
#'                     string: "". 
#'                     Group names should only be 'case', 'control' or 'null'. 
#' @param groups.samples
#'            Sets the sample groups with a named list of character vectors, with names the group
#'            and vectors of sample names. 
#'            Group names should only be 'case', 'control' or 'null'. 
#' @param filter.low.counts If TRUE, apply automatic filtering of counts lower than
#'            expected by chance (assuming independent and uniform DNA base distribution).
#' 
#' @references 
#' Bio::STR:exSTRa Perl module
#' \url{https://github.com/bahlolab/Bio-STR-exSTRa}
#' 
#' @examples 
#' str_score <- read_score (
#'     file = system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt.gz", package = "exSTRa"), 
#'     database = system.file("extdata", "repeat_expansion_disorders_hg19.txt", package = "exSTRa"), 
#'     groups.regex = c(control = "^WGSrpt_0[24]$", case = ""), 
#'     filter.low.counts = TRUE
#' )
#' str_score
#' str_score$samples
#' 
#' # Plot the HD locus only:
#' plot(str_score["HD"])
#' 
#' 
#' # Defining cases by sample name directly:
#' str_score_HD_cases <- read_score (
#'     file = system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"), 
#'     database = system.file("extdata", "repeat_expansion_disorders_hg19.txt", package = "exSTRa"),
#'     groups.samples = list(case = c("WGSrpt_10", "WGSrpt_12")),
#'     filter.low.counts = TRUE
#' )
#' str_score_HD_cases
#' str_score_HD_cases$samples
#' 
#' plot(str_score_HD_cases["HD"])
#' 
#' # for greater control, use object from read_exstra_db() instead
#' str_db <- read_exstra_db(
#'             system.file("extdata", "repeat_expansion_disorders_hg19.txt", package = "exSTRa")
#'           )
#' str_score <- read_score (
#'     file = system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"), 
#'     database = str_db, 
#'     groups.regex = c(control = "^WGSrpt_0[24]$", case = ""), 
#'     filter.low.counts = TRUE
#' )
#' str_score
#' str_score$samples
#' 
#' @export
read_score <- function(file, database, groups.regex = NULL, groups.samples = NULL, filter.low.counts = TRUE) {
  # Load the STR counts
  # Groups should be named null, control and case
  if(is.character(database)) {
    # as database is presumbly a file, try to read from it
    database <- read_exstra_db(database)
  }
  out <- strs_read_(file, database, groups.regex, groups.samples, this.class = "exstra_score")
  # TODO: checks for this input
  out$data$prop <- with(out$data, rep / mlength)
  #strscore <- exstra_score_new_(out$data, out$db)
  strscore <- exstra_score_new_(out$data, database)
  if(filter.low.counts) {
    # Filter low counts, assumed wanted by default
    strscore <- filter_low_scores(strscore)
  }
  # Remove any loci without data:
  strscore$db <- strscore$db[locus %in% strscore$data$locus]
  return(strscore)
}
