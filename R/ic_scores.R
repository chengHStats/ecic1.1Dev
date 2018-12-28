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
  } else if (class(model)[1] == "ecicModel"){
  UseMethod("IC", ic)
  } else {
    stop("Unexpected input for model in IC function.", call. = FALSE)
  }
}

#' @export
IC.AIC <- function(model, data, ic = 'AIC', compress = F){
  logl <- logLik(model, data)
  k = model$k
  out = c(-2*logl[1] + 2*k, logl[-1])
  names(out) = c('AIC', model$parameter.names)
  out
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
    } else if (class(model)[1] == "ecicModel"){
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
  liks <- logLikMulti(model, data)
  ic <- (-2*liks[,1])+2*model$k
  out = cbind(ic, liks[,-1])
  colnames(out) = c("AIC", colnames(liks)[-1])
  out
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
  models = ecicModelList(models)
  p = length(models)
  n = nrow(data)
  N = ncol(data)
  ics = lapply(models, function(x) ICMulti(x, data, ic))
  ncols = sapply(ics, ncol)
  ic.ix = c(1, cumsum(ncols[-p])+1)
  ics1 = lapply(1:p, function(ix){
    x = ics[[ix]]
    'colnames<-'(x, paste(names(models)[ix], colnames(x), sep = '.'))
  } )
  out = do.call(cbind, ics1)%>%data.frame
  params = lapply(ics, function(x) x[,-1])
  if(!compress){
    out = list(ic = out[,ic.ix], parameters = params)
    colnames(out$ic) = names(models)
  }
  return(out)
}


