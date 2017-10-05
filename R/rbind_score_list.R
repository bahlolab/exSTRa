# Combine multiple exstra_score objects, checking for sample name clashes

# TODO: can roxygen make this the one man page?

# old name: rbind.exstra_score.list
#' @export
#' Combine multiple exSTRa data objects
#' 
#' @param strscore_list A list containing exstra_score objects.
#' @param idcol The name of the column giving the names of strscore_list in the output.
#' @param allow_sample_clash If TRUE, allows a sample name to occur multiple times in the exstra_score objects, otherwise duplicates cause an error. 
#' @param ... Additional arguments to rbindlist().
#' @return An exstra_score object.
#' @examples
rbind_score_list <- function(strscore_list, idcol = "data_group", allow_sample_clash = FALSE, ...) {
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
  
  data.new <- rbindlist(lapply(strscore_list, function(x) { x$data }), idcol = idcol, ...)
  db.new.db <- rbindlist(lapply(strscore_list, function(x) { x$db }), ...)
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
#' @export
rbind_score <- function(..., idcol = "data_group", allow_sample_clash = FALSE, fill = FALSE) {
  rbind_score_list(list(...), idcol = idcol, allow_sample_clash = allow_sample_clash, fill = fill)
}
