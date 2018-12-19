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
IC <- function(data, model, ic = 'AIC', compress = F){
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

   do.call(paste("IC.",ic,sep = ""),list(data = data, model = model,ic = ic))
}

#' @export
IC.AIC <- function(data, model, ic = 'AIC', compress = F){
  logl <- loglik(data, model)
  out = c(-2*logl[1] + 2*length(logl[-1]), logl[-1])
  names(out) = NULL
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
ICMulti <- function(data, model, ic = "AIC"){
  fnc <- paste("ICMulti.", ic, sep = "")
  args <- list(data = data, model = model, ic = ic)
  do.call(fnc, args)
}

#' @export
ICMulti.AIC <- function(data, model, ic){
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
  liks <- logliks(data, model)
  k <- ncol(liks)-1
  ic <- (-2*liks[,1])+2*k
  out = cbind(ic, liks[,-1])
  colnames(out) = c("ic", colnames(liks)[-1])
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
ICMultiMulti = function(data, models, ic, compress = F){
  names(models) <- models
  p = length(models)
  n = nrow(data)
  N = ncol(data)
  liks = sapply(models, function(x) ICMulti(data, x, ic))
  ncols = sapply(liks, ncol)
  ic.ix = c(1, cumsum(ncols[-p])+1)
  liks1 = sapply(1:p, function(ix){
    x = liks[[ix]]
    'colnames<-'(x, paste(models[ix], colnames(x), sep = '.'))
  } )
  out = do.call(cbind, liks1)%>%data.frame
  params = sapply(liks, function(x) x[,-1])
  if(!compress){
    out = list(ic = out[,ic.ix], parameters = params)
    colnames(out$ic) = models
  }
  return(out)
}


