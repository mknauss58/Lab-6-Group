#' Implements simple linear regression by gradient descent
#'
#' @param dat A data frame
#' @param response The name of a response variable in the data frame (unquoted)
#' @param explanatory The name of the explanatory variable in the data frame (unquoted)
#'
#' @return A data frame of coefficients
#'
#' @import dplyr
#'
#' @export
slr_gd <- function(dat, response, explanatory){

  ### Compute coefficients by gradient descent
  ### Return a data frame of the same form as in the `simple_linear_regression`
  beta_1 = 0
  beta_0 = 0

  learn = 0.001
  numgd = 10000
  n = nrow(dat)

  x <- dat %>% pull({{explanatory}})
  x <- scale(x)
  y <- dat %>% pull({{response}})
  y <- scale(y)

  explan_name <- dat %>%
    select({{explanatory}}) %>%
    names()

  for (i in 1:numgd){
    pred = x*beta_1 + beta_0

    deriv_beta1 = (-2/n) * sum(x * (y-pred))
    deriv_beta0 = (-2/n) * sum(y-pred)

    beta_1 <- beta_1 - (learn * deriv_beta1)
    beta_0 <- beta_0 - (learn * deriv_beta0)
  }

  results <- tibble::tibble(
    Intercept = beta_0,
    Slope = beta_1
  )

  names(results)[2] <- explan_name
  return(results)
}


#' Implements linear regression with many predictors by gradient descent
#'
#' This function computes coefficients for multiple regression by gradient descent
#' All columns of the provided data frame are used as predictors, except the
#' one specified as a response.
#'
#' No interaction terms are included.
#'
#'
#' @param dat A data frame
#' @param response The name of a response variable in the data frame (unquoted)
#'
#' @return A data frame of coefficients
#'
#' @import dplyr
#'
#'@export
mlr_gd <- function(dat, response) {

  ### Compute coefficients by gradient descent
  ### Return a data frame of the same form as in the `multiple_linear_regression`
  betas = matrix(0,ncol(dat))

  learn = 0.001
  numgd = 10000
  n = nrow(dat)

  x <- dat %>% select(-{{response}})
  x <- scale(x)
  x <- as.matrix(x)
  int <- as.matrix(rep(1, nrow(x)), ncol = 1)
  x <- cbind(int, x)
  y <- dat %>% select({{response}})
  y <- scale(y)
  y <- as.matrix(y)

  for (i in 1:numgd){
    gradient <- (t(x) %*% (y - (x %*% betas)))
    betas <- betas - (learn * -2*gradient)
  }

  results <- as.data.frame(t(betas))
  names(results)[1] <- "Intercept"

  return(results)

}

#' Implements linear regression with many predictors by matrix decomposition
#'
#' This function computes coefficients for multiple regression by QR matrix decomposition
#' All columns of the provided data frame are used as predictors, except the
#' one specified as a response.
#'
#' No interaction terms are included.
#'
#'
#' @param dat A data frame
#' @param response The name of a response variable in the data frame (unquoted)
#'
#' @return A data frame of coefficients
#'
#' @import dplyr
#'
#'@export
mlr_qr <- function(dat, response) {

  x <- dat %>% select(-{{response}})
  xm <- as.matrix(x)
  xm <- cbind(1,xm)
  y <- dat %>% pull({{response}})
  ym <- as.matrix(y)

  decomp <- qr(xm)

  Q <- qr.Q(decomp)
  R <- qr.R(decomp)

  results <- solve(R) %*% t(Q) %*% ym

  results <- as.data.frame(t(results))
  names(results)[1] <- "Intercept"

  return(results)
}
