# ModelFrequencies--------------------------------------------------------------
#'  Generate Bootstrap Information Criterion Frequencies
#'
#' @param n: Sample size.
#' @param true: A string giving the name of the model to generate data from.
#' @param param: A vector of parameters specifying a fitted true model.
#' @param models: A string vector of model names for selecting between.
#' @param N: Number of bootstrap samples.
#' @param ic: String giving which information criterion to use (e.g. 'AICc', 'BIC')
#'
#' @return: A list containing the following:.
#'
#' @examples
#' models = c("norm0", "norm1", "norm")
#' ModelFrequencies(25, "norm", c(0.3, 1.2), models, 1000, "AIC")
#'
#' @export
ModelFrequencies <- function(n, true, param, models, N = 1000, ic = 'AIC'){
  # Bootstrap IC frequencies for a model with specified param and sample size.
  #
  # Args:
  #   n: Sample size.
  #   true: A string giving the name of the model to generate data from.
  #   param: A vector of model parameters.
  #   models: A string vector of model names for selecting between.
  #   N: Number of bootstrap samples.
  #   ic: String giving which information criterion to use (e.g. 'AICc', 'BIC')
  #
  # Returns:
  #   A vector of the frequencies each model attained the lowest IC value.
  names(models) <- models
  data <- GenerateDataMulti(n, true, param, N)
  icmm <- ICMultiMulti(data, models, ic)
  scores = icmm$ic

  freqs = table(
    factor(apply(scores,1,function(x) models[which.min(x)]),
           levels = models))/N
  return(list(frequencies = freqs, data = data, ic.scores = scores, parameters = icmm$parameters))
}
