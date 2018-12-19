
# Estimate parameters for a single sample --------------------------------------

# EstimateParameters -----------------------------------------------------------
#' Computes the fitted parameters for a data sample and a model.
#'
#' @param data: compatible with the specified model.
#' @param model: A string specifying the model to compute the log-likehood for.
#'
#' @return A vector/matrix of the fitted parameters.
#'
#' @examples
#' data = GenerateData(25, "norm", c(mu = 0.3, sd = 1.2))
#' models = c("norm0", "norm1", "norm")
#' EstimateParameters(data)
#' @export
EstimateParameters <- function(data, model){
  do.call(paste("EstimateParameters.", model, sep = ""),
          list(data = data, model = model))
}

# EstimateParameters.norm ------------------------------------------------------
#' @export
EstimateParameters.norm <- function(data, model){
  # Computes the fitted parameters for a data sample and a norm(mu, sd) model.

  n <- length(data)
  mu <- sum(data)/n
  sd <- (crossprod(data-mu)/n)**.5
  return(c(mu, sd))
}

# EstimateParameters.norm0 -----------------------------------------------------------

#' @export
EstimateParameters.norm0 <- function(data, model){
  # Returns empty for norm(0,1) model.

  return()
}

# EstimateParameters.norm1 -----------------------------------------------------------
#' @export
EstimateParameters.norm1 <- function(data, model){
  # Computes the fitted parameters for a data sample and a norm(mu, 1) model.

  n <- length(data)
  mu <- sum(data)/n
  return(mu)
}

# EstimateParameters.rwalk -----------------------------------------------------------
#' @export
EstimateParameters.rwalk <- function(data, model){
  # Computes fitted parameters for a data sample and random walk (mu, sd) model.

  n <- length(data)
  steps <- c(data[1], diff(data))
  EstimateParameters.norm(steps, "norm")
  mu <- sum(steps)/n
  sd <- (crossprod(steps-mu)/n)**.5
  return(c(mu, sd))
}

# EstimateParameters.rwalk0 -----------------------------------------------------------
#' @export
EstimateParameters.rwalk0 <- function(data, model){
  # Computes fitted parameters for a data sample and random walk (0, sd) model.

  n <- length(data)
  steps <-c(data[1], diff(data))
  (crossprod(steps)/n)^.5
}

# EstimateParametersMulti -------------------------------------
#' Computes the fitted parameters for multiple data samples and a model.
#'
#' @param data: A matrix of data samples compatible with the specified model.
#' @param model: A string specifying the model to compute the log-likehood for.
#'
#' @return A vector/matrix of the fitted parameters.
#'
#' @examples
#' data = GenerateDataMulti(25, "norm", c(mu = 0.3, sd = 1.2), 10)
#' models = c("norm0", "norm1", "norm")
#' EstimateParameters(data)
#'
#' @export
EstimateParametersMulti <- function(data, model, flat = F){
  # Computes the fitted parameters for multiple data samples and a model.
  #
  # Args:
  #   data: compatible with the specified model.
  #   model: A string specifying the model to compute the log-likehood for
  #
  # Returns:
  #   A vector/matrix of the fitted parameters.
  fnc <- paste("EstimateParametersMulti.", model, sep = "")
  args <- list(data = data, model = model, flat = flat)
  do.call(fnc, args)
}

#' @export
EstimateParametersMulti.norm <- function(data, model, flat = F, ...){
  # Computes fitted parameters for multiple data samples and norm(mu, sd) model.

  n <- nrow(data)
  N <- ncol(data)
  mus <- crossprod(rep(1/n, n), data)
  r <-crossprod(data - matrix(rep(mus, each = n), ncol = N))
  sds <- {diag(r)/n}**.5
  out = if (flat) {
    c(mus,sds)
  } else {
    matrix(c(mus, sds), nrow = 2, byrow = TRUE)
  }
  return(out)
}

#' @export
EstimateParametersMulti.norm1 <- function(data, model, ...){
  # Computes fitted parameters for multiple data samples and norm(mu, sd) model.

  n <- nrow(data)
  N <- ncol(data)
  mus <- crossprod(rep(1/n, n), data)
  mus
}

#' @export
EstimateParametersMulti.norm0 <- function(data, model, ...){
  # Computes fitted parameters for multiple data samples and norm(mu, sd) model.

  return()
}

