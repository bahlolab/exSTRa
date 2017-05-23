# TODO:
# this function is from http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
# so may need a rewrite
# old name: add.alpha
#' @export
add_alpha_ <- function(col, alpha = 1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  if(length(col) == 0) {
    stop("col is of length 0")
  }
  apply(sapply(col, col2rgb)/255, 2, 
    function(x) 
      rgb(x[1], x[2], x[3], alpha=alpha))  
}