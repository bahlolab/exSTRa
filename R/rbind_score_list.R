# Combine multiple exstra_score objects, checking for sample name clashes

# TODO: can roxygen make this the one man page?

# old name: rbind.exstra_score.list
#' Combine multiple exSTRa data objects.
#' 
#' Allows data from multiple \code{exstra_score} objects to be combined, similarly to
#' \code{\link[data.table]{rbindlist}} (data.table) and \code{\link[base]{rbind}}.
#' The \code{exstra_score} objects may either be given directly to \code{rbind_score}, or
#' as a list to \code{rbind_score_list} that is easier to use when the number of objects
#' to combine is not known in advance. 
#' 
#' @param strscore_list A list containing \code{exstra_score} objects.
#' @param ... \code{exstra_score} objects to combine.
#' @param idcol The name of the column giving the names of \code{strscore_list} in the output.
#' @param allow_sample_clash If TRUE, allows a sample name to occur multiple times in the 
#'         \code{exstra_score} objects, otherwise duplicates cause an error.
#' @param fill If TRUE, missing columns are filled with NAs. Default FALSE. 
#' @return An exstra_score object.
#' 
#' @seealso \code{\link[data.table]{rbindlist}}, \code{\link[base]{rbind}}
#' 
#' @examples
#' # create a list of distinct samples
#' score_list <- list(group1 = exstra_wgs_pcr_2[, "WGSrpt_20"], 
#'                    controls = exstra_wgs_pcr_2[, c("WGSrpt_02", "WGSrpt_04")],
#'                    group2 = exstra_wgs_pcr_2[, c("WGSrpt_15", "WGSrpt_16")]
#'                   )
#' combined_scores <- rbind_score_list(score_list)
#' combined_scores$samples
#' 
#' # Without using a list:
#' rbind_score(group1 = exstra_wgs_pcr_2[, "WGSrpt_20"], 
#'             controls = exstra_wgs_pcr_2[, c("WGSrpt_02", "WGSrpt_04")],
#'             group2 = exstra_wgs_pcr_2[, c("WGSrpt_15", "WGSrpt_16")]
#'             )$samples
#' 
#' # Combining with a repeated sample name, possibly due to multiple experiments.
#' # Splitting data for this example:
#' score_1 <- exstra_wgs_pcr_2["DM1", "WGSrpt_15"]
#' score_2 <- exstra_wgs_pcr_2["HD", "WGSrpt_15"]
#' # Combine with the sample name repeated across exstra_score objects:
#' rbind_score_list(list(score_1, score_2), allow_sample_clash = TRUE)
#' # or
#' rbind_score(score_1, score_2, allow_sample_clash = TRUE)
#' @export
rbind_score_list <- function(strscore_list, idcol = "data_group", 
  allow_sample_clash = FALSE, fill = FALSE) {
  assert("strscore_list must be a list", inherits(strscore_list, "list"))
  if(length(strscore_list) == 0) {
    stop("List is empty")
  }
  assert("Not all elements are rep score data", is.exstra_score(strscore_list[[1]]))
  if(length(strscore_list) == 1) {
    return(strscore_list[[1]])
  }
  for(i in seq_along(strscore_list)) {
    assert(paste("Element", i, "is not exstra_score"), is.exstra_score(strscore_list[[i]]))
    assert("STR database is of mixed types", strscore_list[[1]]$input_type == strscore_list[[i]]$input_type)
  }
  
  # Could be written much better, all in one go here instead, rather than recursion
  
  data.new <- rbindlist(lapply(strscore_list, function(x) { x$data }), idcol = idcol, fill = fill)
  db.new.db <- rbindlist(lapply(strscore_list, function(x) { x$db }), fill = fill)
  setkey(db.new.db, locus)
  db.new.db <- unique(db.new.db)
  db.new <- exstra_db_new_(db.new.db, input_type = strscore_list[[1]]$input_type)
  new_strscore <- exstra_score_new_(data.new, db.new)
  new_strscore$samples <- rbindlist(lapply(strscore_list, function(x) { x$samples }), idcol = idcol, fill = TRUE)
  setkey(new_strscore$samples, sample)
  if(!allow_sample_clash) {
    test <- table (new_strscore$samples$sample)
    if(max(test) > 1) {
      stop("A sample name is duplicated in inputs, for sample names: ", 
        paste(names(which(test > 1)), collapse = ", "), 
        "\nSet allow_sample_clash = TRUE if this is ok. "
      )
    }
  }
  return(new_strscore)
}

# convinient version of rbind_exstra_score_list() without the use of lists
#' @rdname rbind_score_list
#' @export
rbind_score <- function(..., idcol = "data_group", allow_sample_clash = FALSE, fill = FALSE) {
  rbind_score_list(list(...), idcol = idcol, allow_sample_clash = allow_sample_clash, fill = fill)
}
