# ModelFrequencies--------------------------------------------------------------
#'  Generate Bootstrap Information Criterion Frequencies
#'
#' @parameters n: Sample size.
#' @parameters true: A string giving the name of the model to generate data from.
#' @parameters parameters: A vector of parameterseters specifying a fitted true model.
#' @parameters models: A string vector of model names for selecting between.
#' @parameters N: Number of bootstrap samples.
#' @parameters ic: String giving which information criterion to use (e.g. 'AICc', 'BIC')
#'
#' @return: A list containing the following:.
#'
#' @examples
#' models = c("norm0", "norm1", "norm")
#' ModelFrequencies(25, "norm", c(0.3, 1.2), models, 1000, "AIC")
#'
#' @export
ModelFrequencies <- function(n, true, parameters, models, N = 1000, ic = 'AIC'){
  models = ecicModelList(models)
  data <- GenerateDataMulti(n, true, parameters, N)
  icmm <- ICMultiMulti(models, data, ic)
  scores = icmm$ic

  freqs = table(
    factor(apply(scores,1,function(x) names(models)[which.min(x)]),
           levels = names(models)))/N
  return(structure(list(frequencies = freqs,
                        data = data,
                        ic.scores = scores,
                        parameters = icmm$parameters,
                        true = true),
                   class = "ecicFrequencies"))
}
