# Algorithm for imputing missing data using quantiles from:
# Munoz JF, Rueda M (2009) New imputation methods for missing data using quantiles. Journal of Computational and Applied Mathematics 232:305-317. doi: 10.1016/j.cam.2009.06.011
#
# Only algorithm AL1 implemented here

library(testit)

munoz_rueda_al1 <- function(x, extra.missing = 0, names = TRUE) {
    # implement the imputation from J.F. Munoz, M. Rueda 2009
    # Only give imputed values
    assert("x should be numeric vector", is.vector(x), is.numeric(x))
    assert("extra.missing should be >= 0", is.vector(extra.missing), is.numeric(extra.missing), extra.missing >= 0)
    ## Step 1
    m <- sum(is.na(x)) + extra.missing
    x <- x[!is.na(x)]
    r <- length(x)
    if(m == 0) {
        return(numeric(0))
    }
    if(m == 1) {
        return("50%" = quantile(x, 0.5))
    }
    ## Step 2
    if(r / (m + 1) > 0.5) {
        S <- 1 / round(r / (m + 1)) # S = 1 / B
    } else {
        S <- 1
    }
    
    # Looped steps
    w <- seq(0, 1, S) # vector may be set at the start
    k <- 0 # first use k=1 due to loop
    i <- seq_len(m) # always goes to m
    d2 <- rep(NA, m)
    d2_func <- function(y_m, y_r) {
        y_bar_m <- mean(y_m)
        y_bar_2_m <- sum(y_m ^ 2) / length(y_m) # here this is an imputed value
        y_bar_r <- mean(y_r)
        y_bar_2_r <- sum(y_r ^ 2) / length(y_r) # here this is the actual value
        abs((y_bar_m - y_bar_r) / y_bar_r) + abs((y_bar_2_m - y_bar_2_r) / y_bar_2_r)
    }
    while(k == 0 || w[k] < 1){
        k <- k + 1 # part of step 6
        ## Step 3
        alpha <- w[k] * i / (m + 1) + (1 - w[k]) * (i - 1) / (m - 1)
        ## Step 4
        y_star <- quantile(x, alpha, names = FALSE, type = 1) # may not be correct
        ## Step 5
        # with criterion (7)
        d2[k] <- d2_func(y_star, x)
        ## Step 6 is
        # already coded into loop and w
    }
    ## Step 7
    k_max <- which.max(d2)
    ## Step 8
    quantile(x, w[k_max] * i / (m + 1) + (1 - w[k_max]) * (i - 1) / (m - 1), names = names, type = 1)
}

munoz_rueda_al1_include <- function(x, extra.missing = 0, ...) {
    # Perform munoz_rueda_al1 but also include original values
    # If extra.missing = 0, then replaces NAs with imputed values randomly placed
    # If extra.missing != 0, does not maintain order
    y <- munoz_rueda_al1(x = x, extra.missing = extra.missing, names = FALSE, ...)
    if(extra.missing == 0){
        if(length(y) == 1) {
            x[is.na(x)] <- y
        } else {
            x[is.na(x)] <- sample(y, length(y))
        }
        return(x)
    } else {
        return(c(x[!is.na(x)], y))
    }
}

# munoz_rueda_al1(c(1, 5, 7, 5, 4, 8, 3))
# 
# munoz_rueda_al1(c(1, 5, 7, 5, 4, NA, 8, 3)) 
# 
# munoz_rueda_al1(c(1, 5, 7, NA, 5, NA, 4, 8, 3))
# 
# munoz_rueda_al1(c(1, 5, 7, NA, 5, NA, 4, 8, 3), 10)
# munoz_rueda_al1(c(1, 5, 7, NA, 5, NA, 4, 8, 3), 10, names = FALSE)
# 
# munoz_rueda_al1_include(c(1, 5, 7, NA, 5, NA, 4, 8, 3), 10)
# 
# munoz_rueda_al1_include(c(1, 5, 7, NA, 5, NA, 4, 8, 3))
# 
# munoz_rueda_al1_include(c(NA, NA, 1, 5, 7, NA, 5, NA, 4, 8, NA, 3))
# 
# munoz_rueda_al1_include(numeric(0), 15)

# Example where different k is chosen
# munoz_rueda_al1(c(1, 5, 7, NA, 5, NA, 4, 8, 3))
# munoz_rueda_al1(c(1, 5, 7, NA, 5, NA, 8, 3))


# y <- c(31L, 31L, 34L, 26L, 16L, 27L, 27L, 34L, 27L, 33L, 24L, 11L, 
#     27L, 21L, 19L, 27L, 31L, 24L, 26L, 28L, 27L, 30L, 34L, 35L, 27L, 
#     25L, 31L, 34L, 26L, 17L, 26L, 34L, 27L, 21L, 17L, 34L, 26L, 17L, 
#     26L, 27L, 31L, 34L, 31L, 22L, 31L, 38L)
# munoz_rueda_al1_include(c(NA, y), 0)
# munoz_rueda_al1_include(c(y), 0)
# munoz_rueda_al1_include(c(y), 1)
