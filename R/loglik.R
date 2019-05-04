# logLik -----------------------------------------------------------------------
#' Computes log-likelihood for data and a specified model.
#'
#' @param model: An ecicModel object.
#' @param data: A data sample compatible with the model.
#'
#'
#' @return: A list containing the following:.
#'
#' @examples
#'
#' @export
logLik.ecicModel = function(model, data, compress = F){
  NextMethod()
}
#' @export
logLik.character = function(model, data, compress = F){
  model1 = tryCatch(ecicModel(model),
                    error = function(e){
                      stop(paste(model, "is not a valid ecicModel type."),
                           call. = FALSE)
                    })
  logLik(model1, data)
}
#' @export
logLik.norm <- function(model, data, compress = F ){
  parameters <- suppressMessages(EstimateParameters.norm(model, data)$parameters)
  logl <- sum(log(dnorm(data, parameters[1], parameters[2])))


  out <- list(log.likelihood = logl, parameters = parameters)

  return(out)
}
#' @export
logLik.rwalk <- function(model, data, compress = F ){
  ep <- suppressMessages(EstimateParameters.rwalk(model, data))
  parameters <- ep$parameters
  steps <- ep$steps
  logl <- sum(log(dnorm(steps, parameters[1], parameters[2])))

  out <- list(log.likelihood = logl, parameters = parameters)

  return(out)
}
#' @export
logLik.lmECIC <- function(model, data, compress = F){
  ep <- EstimateParameters(model, data)
  parameters <- ep$parameters
  resid <- ep$resid
  logl <- sum(log(dnorm(resid, 0, parameters['sd'])))

  out <- list(log.likelihood = logl, parameters = parameters)

  return(out)
}
# logLikMulti ----------------------------------------------------------------------
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
#' logLikMulti(data, "norm")
#'
#' @export
logLikMulti <- function(model, data, compress = F){
  # Computes the log-likelihood for many datasets and a model.
  #
  # Args:
  #   data: A matrix/array of data samples compatible with the specified model.
  #   model: A string specifying the model to compute the log-likelihoods for.
  #   compress: A boolean specifying if output should be of list or vector type.
  #
  # Returns:
  #   logLik: The log-likelihood for the given data samples.
  #   param: The fitted model parameters, which are an intermediate step.
  UseMethod("logLikMulti", model)
}
#' @export
logLikMulti.ecicModel <- function(model, data, compress = F){
  NextMethod()
}

#' @export
logLikMulti.norm <- function(model, data, compress = F){
  # Computes the log-likelihood for many datasets and a norm(mu, sd) model.
  #
  # Args:
  #   data: A matrix/array of data samples compatible with the specified model.
  #   model: A string specifying the model to compute the log-likelihoods for.
  #   compress: A boolean specifying if output should be of list or vector type
  #
  # Returns:
  #   logLik: The log-likelihood for the given data samples.
  #   param: The fitted model parameters, which are an intermediate step.

  parameters <- suppressMessages( EstimateParametersMulti(model, data)$parameters )
  n <- nrow(data)
  N <- ncol(data)
  logl <- crossprod(rep(1, n),
                    log(dnorm(data,
                              rep(parameters[1,], each = n),
                              rep(parameters[2,], each = n))))
  out <- list(log.likelihood = logl, parameters =  parameters)
  return(out)
}
#' @export
logLikMulti.rwalk <- function(model, data, compress = F){
  # Computes the log-likelihood for many datasets and a norm(mu, sd) model.
  #
  # Args:
  #   data: A matrix/array of data samples compatible with the specified model.
  #   model: A string specifying the model to compute the log-likelihoods for.
  #   compress: A boolean specifying if output should be of list or vector type
  #
  # Returns:
  #   logLik: The log-likelihood for the given data samples.
  #   param: The fitted model parameters, which are an intermediate step.

  ep <- suppressMessages(EstimateParametersMulti(model, data))
  parameters <- ep$parameters
  steps <- ep$steps
  n <- nrow(data)
  N <- ncol(data)
  k = 2
  likelihood = dnorm(steps,
                     rep(parameters[1, ], each=n),
                     rep(parameters[2,], each=n))
  logl <- crossprod(rep(1, n),
                    log(dnorm(steps,
                              rep(parameters[1, ], each = n),
                              rep(parameters[2, ], each = n))))

  out <- list(log.likelihood = logl, parameters =  parameters)
  return(out)
}
#' @export
logLikMulti.lmECIC = function(model, data, compress = F){
  ep <- EstimateParametersMulti(model, data)
  parameters <- ep$parameters
  k = nrow(parameters)
  resid <- ep$resid
  N = ncol(resid)
  logl <- sapply(1:N, function(x) sum(log(dnorm(resid[,x], 0, parameters['sd',x]))))
  out <- list(log.likelihood = logl, parameters =  parameters)

}

