# IC ---------------------------------------------------------------------------
#' Computes a given information criterion for data and a specified model.
#'
#' @param data: A data sample compatible witih the model.
#' @param model: A string containing the model name.
#' @param ic: String giving  information criterion to use (e.g. 'AICc', 'BIC')
#'
#' @return: A list containing the following:.
#'
#' @examples
#'
#'
#' @export
IC <- function(model, data, ic = 'AIC', compress = F){
  # Computes the given information criterion for a single dataset and model.
  #
  # Args:
  #   data: A vector/matrix data sample compatible with the specified model.
  #   model: A string specifying the model to fit and compute the score for
  #   verbose: If TRUE, prints sample covariance; if not, not. Default is TRUE.
  #   compress: boolean specifying if output should be of list or vector type
  #
  # Returns:
  #   ic: the criterion value for the given data
  #   param: the fitted model parameters, which are an intermediate step
  #
  # Error handling
  class(ic) = ic
  if (is.character(model)){
    model1 = tryCatch(ecicModel(model),
                      error = function(e){
                        stop(paste(model, "is not a valid ecicModel type."),
                             call. = FALSE)
                      })
    IC(model1, data, ic, compress)
  } else if (inherits(model, "ecicModel"){
  UseMethod("IC", ic)
  } else {
    stop("Unexpected input for model in IC function.", call. = FALSE)
  }
}

#' @export
IC.AIC <- function(model, data, ic = 'AIC', compress = F){
  logl <- logLik(model, data)
  k = model$k
  IC = -2*logl$log.likelihood + 2*k
  out = list(ic=IC, parameters = logl$parameters)
  #names(out) = c('AIC', model$parameter.names)
  return(out)
}

# ICMulti ----------------------------------------------------------------------
#' Computes an information criterion for many data samples and specified model.
#'
#' @param n: Sample size.
#' @param model: A string containing the model name.
#' @param ic: String giving information criterion to use (e.g. 'AICc', 'BIC')
#'
#' @return: A list containing the following:.
#'
#' @examples
#'
#'
#' @export
ICMulti <- function(model, data, ic = "AIC"){
  class(ic) = ic
  if (is.character(model)){
      model1 = tryCatch(ecicModel(model),
                        error = function(e){
                          stop(paste(model, "is not a valid ecicModel type."),
                               call. = FALSE)
                        })
      ICMulti(model1, data, ic)
    } else if (inherits(model, "ecicModel"){
      UseMethod("ICMulti", ic)
    } else {
      stop("Unexpected input for model in ICMulti function.", call. = FALSE)
    }
}

#' @export
ICMulti.AIC <- function(model, data, ic){
  # Computes the Akaike information criterion for many datasets and a model.
  #
  # Args:
  #   data: A matrix/array of data samples compatible with the specified model.
  #   model: A string specifying the model to fit and compute the score for
  #   ic: which information criterion (e.g. AIC, AICc, BIC)
  #   compress: boolean specifying if output should be of list or vector type
  #
  # Returns:
  #   ic: the criterion value for the given data samples
  #   param: the fitted model parameters, which are an intermediate step
  if (inherits(data,"matrix"){
  x <- logLikMulti(model, data)
  liks = x$log.likelihood
  parameters = x$parameters
  ic <- (-2*liks)+2*model$k
  out = list(ic = ic, parameters = parameters)
  return(out)
}
  if (inherits(data,"list"){
    if (inherits(data[[1]],"paleoTS"){
      fits <- logLikMulti(model, data)
      out = c()
      for (fit in fits){
        fit$ic = -2*fit$log.likelihood + 2*model$k
        out = c(out, list(fit))
      }
      out
    }
  }
}

# ICMultiMulti------------------------------------------------------------------
#' Computes an information criterion for many data samples and list of models.
#'
#' @param n: Sample size.
#' @param model: A vector of strings containing the model names.
#' @param ic: String giving information criterion to use (e.g. 'AICc', 'BIC')
#'
#' @return: A list containing the following:.
#'
#' @examples
#'
#'
#' @export
ICMultiMulti = function(models, data, ic = "AIC", compress = F){
  if (inherits(data,"matrix"){
    models = ecicModelList(models)
    p = length(models)
    data.type = models[[1]]$data.type

    if(data.type==1){
    n = nrow(data)
    N = ncol(data)
    x = lapply(models, function(x) ICMulti(x, data, ic))

    ics = sapply(x, function(y) y$ic)
    parameters = lapply(x, function(y) y$parameters)

    out = list(ic = ics, parameters = parameters)
    return(out)
    }
  }
    if (inherits(data,"list"){
      if (inherits(data,"paleoTS"){
        out = lapply(models, function(model) ICMulti(model, data))
        return(out)
      }
    }



}

#' @export
AICweights = function(scores){
  weights = function(row){
    ex = sapply(row, function(x) exp(-0.5*x))
    return(ex/sum(ex))
  }
  if (inherits(scores,"matrix"){
    return(t(apply(scores, 1, weights)))
  } else {
    return(weights(scores))
  }

}
