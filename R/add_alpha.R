# this function is taken from
# http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
# old name: add.alpha
#' Add transparency to a color value
#'
#' Adapted from 
#' Markus Gesmann (Apr 30, 2013) How to change the alpha value of colours in R. 
#' Retrieved from 
#' https://magesblog.com/post/2013-04-30-how-to-change-alpha-value-of-colours-in/
#' 
#' Used in exSTRa with permission.
#' 
#' @param col Input color character.
#' @param alpha Transparency level, a numeric from 0 to 1 inclusive.
#' @return A character vector.
#' @importFrom grDevices col2rgb rgb
#' @export
add_alpha <- function(col, alpha = 1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  if(length(col) == 0) {
    return(col)
  }
  apply(sapply(col, col2rgb)/255, 2, 
    function(x) 
      rgb(x[1], x[2], x[3], alpha=alpha))  
}