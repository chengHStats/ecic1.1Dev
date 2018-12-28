#' Corrects for bias in parameter estimates due to model uncertainty.
#'
#' @param n: Sample size.
#' @param true: A string giving the name of the model to generate data from.
#' @param models: A string vector of model names for selecting between.
#' @param param: A vector of parameters specifying a fitted true model.
#' @param N: Number of bootstrap samples.
#' @param ic: String giving which information criterion to use (e.g. 'AICc', 'BIC')
#
#' @return: The corrected model parameters.
#'
#' @examples
#' data = GenerateData(25, "norm", c(mu = 0.3, sd = 1.2))
#' models = c("norm0", "norm1", "norm")
#' ecic = ECIC(data, models, alpha = 0.05, N = 100, ic = "AIC", correct = T)
#'
#' @export
BiasCorrect <- function(n, true, parameters, models, N = 1000, ic = 'AIC'){
  # Corrects biased parameter estimates stemming from partitioning by best model
  #
  parameters.full = parameterCheck(true, parameters)
  models = ecicModelList(models)
  p = length(true$parameter.names)
  true.ix = which(names(models) == true$ID)
  if(p > 0){
    newdata <- suppressMessages(GenerateDataBest(n, true, parameters, true, models, N))
    params.boot = suppressMessages(EstimateParametersMulti(true, newdata)$parameters)

    if(p == 1){
      param.boot <- params.boot %>% mean
    } else{
      param.boot <- params.boot%>% rowMeans
    }
    param.corrected = 2*unlist(parameters.full)[true$parameter.names]-param.boot
    return(structure(list(parameters = param.corrected, data = newdata),
                     class = 'ecicBiasCorrect'))
  }
}
