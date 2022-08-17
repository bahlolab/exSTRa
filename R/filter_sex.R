# old name: str_filter_sex
#' filter rep_score_data by sex
#' 
#' @param strscore exstra_score object
#' @param sex the sex to filter, can be any of:
#'            - "all":     no filtering
#'            - "male":    only male samples
#'            - "female":  only female samples
#'            - "missing": only missing samples
#'            - "known":   only samples with sex assigned
#' @param safe When TRUE, missing sex assignments with cause an error when \code{sex} is any of "all", "male" or "female".
#' @export
filter_sex <- function(strscore, sex = "known", safe = TRUE) {
  # filter rep_score_data by sex
  # sex can be:
  #   "all":     no filtering
  #   "male":    only male samples
  #   "female":  only female samples
  #   "missing": only missing samples
  #   "known":   only samples with sex assigned
  if(sex %in% c("all", "male", "female")) {
    if(safe) {
      # Check that no data is missing
      if(sum(is.na(strscore$samples$sex)) != 0) {
        stop("In str_filter_sex(), some samples have not been assigned a sex.")
      }
    }
    if(sex == "all") {
      return(strscore)
    } else if(sex == "male") {
      return(strscore[, sex == "male"])
    } else if(sex == "female") {
      return(strscore[, sex == "female"])
    }
  } else if (sex == "missing") {
    strscore[, is.na(sex)]
  } else if (sex == "known") {
    strscore[, !is.na(sex)]
  } else {
    stop("Bad sex assignment")
  }
}