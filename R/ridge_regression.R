#' Implements ridge regression with many predictors
#'
#' This function computes coefficients for ridge regression
#' All columns of the provided data frame are used as predictors, except the
#' one specified as a response.
#'
#' No interaction terms are included.
#'
#'
#' @param dat A data frame
#' @param response The name of a response variable in the data frame (unquoted)
#' @param lambda A vector of penalty terms to try
#'
#' @return A data frame of coefficients
#'
#' @import dplyr
#'
#' @export
ridge_regression <- function(dat, response, lambda) {


  x <- dat %>% select(-{{response}})
  x <- as.matrix(x)
  x <- scale(x)
  int <- as.matrix(rep(1, nrow(x)), ncol = 1)
  x <- cbind(int, x)
  y <- dat %>% select({{response}})
  y <- as.matrix(y)

  results <- solve(crossprod(x) + diag(lambda,nrow = ncol(x))) %*% (t(x) %*% y)

  results <- as.data.frame(t(results))
  results = cbind(results, lambda)
  names(results)[1] <- "Intercept"

  return(results)

}

#' Determines the best penalty term from a set of options
#'
#' This function uses a randomly chosen test and training set
#'
#' No interaction terms are included.
#'
#'
#' @param train_dat A data frame to construct the model from
#' @param test_dat A data frame to test the model on
#' @param response The name of a response variable in the data frame (unquoted)
#' @param lambda A vector of penalty terms to try
#'
#' @return A data frame of penalty terms and resulting errors
#'
#' @import dplyr
#'
#' @export
find_best_lambda <- function(train_dat, test_dat, response, lambdas) {


  lambda_errors = matrix(, ncol = 2, nrow = length(lambdas))
  for (i in 1:length(lambdas)){

    coefs <- ridge_regression({{train_dat}}, {{response}}, {{lambdas}}[i])
    coefs <- as.matrix(coefs[1:ncol(coefs)-1])
    predictions <- predict_from_coefs({{test_dat}}, {{response}}, coefs)
    reponse <- test_dat %>% select({{response}})
    SSE <- sum(reponse - predictions)

    lambda_errors[i,1] = lambdas[i]
    lambda_errors[i,2] = SSE

  }

  lambda_errors <- as.data.frame(lambda_errors)
  names(lambda_errors) <- c("lambda", "error")

  return(lambda_errors)

}
