# loglik -----------------------------------------------------------------------
#' Computes log-likelihood for data and a specified model.
#'
#' @param data: A data sample compatible witih the model.
#' @param model: A string containing the model name.
#' @param ic: String giving  information criterion to use (e.g. 'AICc', 'BIC')
#'
#' @return: A list containing the following:.
#'
#' @examples
#'
#' @export
loglik <-  function(data, model, compress = F){
  # Computes the log-likelihood for a dataset and a model.
  #
  # Args:
  #   data: A vector/matrix data sample compatible with the specified model.
  #   model: A string specifying the model to compute the log-likehood for.
  #   compress: A boolean specifying if output should be of list or vector type
  #
  # Returns:
  #   loglik: The log-likelihood for the given data sample
  #   param: The fitted model parameters, which are an intermediate step.
  do.call(paste("loglik.", model, sep = ""),
          list(data = data, model = model))
}

#' @export
loglik.norm <- function(data, model, compress = F ){
  # Computes the log-likelihood for a dataset and a norm(mu, sd) model.
  #
  # Args:
  #   data: A vector/matrix data sample compatible with the specified model.
  #   model: A string specifying the model to compute the log-likehood for.
  #   compress: A boolean specifying if output should be of list or vector type
  #
  # Returns:
  #   loglik: The log-likelihood for the given data sample
  #   param: The fitted model parameters, which are an intermediate step.
  param <- EstimateParameters.norm(data, model)
  logl <- sum(log(dnorm(data, param[1], param[2])))
  out <- c(logl, param) %>% setNames(c("log.likelihood", rep("est.parameters", 2)))
  return(out)
}

#' @export
loglik.norm1 <- function(data, model, compress = F){
  # Computes the log-likelihood for a dataset and a norm(mu, 1) model.
  #
  # Args:
  #   data: A vector/matrix data sample compatible with the specified model.
  #   model: A string specifying the model to compute the log-likehood for.
  #   compress: A boolean specifying if output should be of list or vector type
  #
  # Returns:
  #   loglik: The log-likelihood for the given data sample.
  #   param: The fitted model parameters, which are an intermediate step.
  param <- EstimateParameters.norm1(data, model)
  logl <- sum(log(dnorm(data, param[1],1)))
  out <- c(logl, param) %>% setNames(c("log.likelihood", "est.parameters"))
  return(out)
}

#' @export
loglik.norm0 <- function(data, model, compress = F){
  # Computes the log-likelihood for a dataset and a norm(0, 1) model.
  #
  # Args:
  #  data: A vector/matrix data sample compatible with the specified model.
  #  model: A string specifying the model to compute the log-likehood for
  #  compress: boolean specifying if output should be of list or vector type
  #
  # Returns:
  #   The log-likelihood for the given data sample
  logl <- sum(log(dnorm(data, 0,1)))
  return(logl)
}


# loliks ----------------------------------------------------------------------
#' Computes log-likelihood for many data samples and specified model.
#'
#' @param data: A matrix of data samples compatible with the model.
#' @param model: A string containing the model name.
#' @param ic: String giving information criterion to use (e.g. 'AICc', 'BIC')
#'
#' @return: A list containing the following:.
#'
#' @examples
#' data = GenerateDataMulti(25, "norm", c(0.3, 1.2), 100)
#' logliks(data, "norm")
#'
#' @export
logliks <- function(data, model, compress = F){
  # Computes the log-likelihood for many datasets and a model.
  #
  # Args:
  #   data: A matrix/array of data samples compatible with the specified model.
  #   model: A string specifying the model to compute the log-likelihoods for.
  #   compress: A boolean specifying if output should be of list or vector type.
  #
  # Returns:
  #   loglik: The log-likelihood for the given data samples.
  #   param: The fitted model parameters, which are an intermediate step.
  do.call(paste("logliks.", model, sep = ""),
          list(data = data, model = model))
}

#' @export
logliks.norm <- function(data, model, compress = F){
  # Computes the log-likelihood for many datasets and a norm(mu, sd) model.
  #
  # Args:
  #   data: A matrix/array of data samples compatible with the specified model.
  #   model: A string specifying the model to compute the log-likelihoods for.
  #   compress: A boolean specifying if output should be of list or vector type
  #
  # Returns:
  #   loglik: The log-likelihood for the given data samples.
  #   param: The fitted model parameters, which are an intermediate step.
  n <- nrow(data)
  N <- ncol(data)
  k = 2
  params <- EstimateParametersMulti(data, "norm", T)
  liks <- crossprod(rep(1, n),
                    log(dnorm(data,
                              rep(params[1:N], each = n),
                              rep(params[(N+1):(2*N)], each = n))))
  out <- matrix(c(liks, params) , ncol = 1+k )
  'colnames<-'(out, c("log.likelihood", rep("parameter", k)))
}

#' @export
logliks.norm1 <- function(data, model, compress = F){
  # Computes the log-likelihoods for many datasets and a norm(mu, 1) model.
  #
  # Args:
  #   data: A matrix/array of data samples compatible with the specified model.
  #   model: A string specifying the model to compute the log-likelihoods for.
  #   compress: A boolean specifying if output should be of list or vector type.
  #
  # Returns:
  #   loglik: The log-likelihood for the given data samples.
  #   param: The fitted model parameters, which are an intermediate step.
  n <- nrow(data)
  N <- ncol(data)
  params <- EstimateParametersMulti(data, "norm1")
  liks <- crossprod(rep(1, n), log(dnorm(data, rep(params[1,], each = n), 1)))
  out <- matrix(c(liks, c(t(params))),ncol = 2)
  'colnames<-'(out, c("log.likelihood", "parameter"))
}


logliks.norm0 <- function(data, model, compress = F){
  # Computes the log-likelihood for many datasets and a norm(0, 1) model.
  #
  # Args:
  #   data: A matrix/array of data samples compatible with the specified model.
  #   model: A string specifying the model to compute the log-likelihoods for.
  #   compress: A boolean specifying if output should be of list or vector type
  #
  # Returns:
  #   The log-likelihood for the given data samples.
  n <- nrow(data)
  N <- ncol(data)
  params <- EstimateParametersMulti(data, "norm")
  liks <- crossprod(rep(1, n), log(dnorm(data, 0, 1))) %>% t
  'colnames<-'(liks, "log.likelihood")
}


