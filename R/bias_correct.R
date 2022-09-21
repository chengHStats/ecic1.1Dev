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
BiasCorrect <- function(n, true.model, parameters, models, N = 1000, ic = 'AIC', genBest = TRUE)
{
  # Corrects biased parameter estimates stemming from partitioning by best model
  # checks for consistency for parameters and model input
  parameters.full = parameterCheck(true.model, parameters)
  models = ecicModelList(models)
  if(true.model$data.type==1)
  {
    p = length(true.model$parameter.names)
    true.model.ix = which(names(models) == true.model$ID)
    if(p > 0)
    {
      if (genBest)
        newdata <- GenerateDataBest(n, true.model, parameters, true.model, models, N)
      else
        newdata <- suppressMessages(GenerateDataMulti(n, true.model, parameters, true.model, models, N))
      params.boot = suppressMessages(EstimateParametersMulti(true.model, newdata)$parameters)
      if(p == 1)
        param.boot <- params.boot %>% mean
      else
        param.boot <- params.boot%>% rowMeans
      param.corrected = 2*unlist(parameters.full)[true.model$parameter.names]-param.boot
      return(structure(list(parameters = param.corrected, data = newdata),class = 'ecicBiasCorrect'))
    }
  }
  if(true.model$data.type=="paleoTS")
  {
    p = true.model$k
    parameters.full = parameterCheck(true.model, parameters)
    if(p > 0)
    {
      if(genBest)
        newdata <- suppressMessages(GenerateDataBest(n, true.model, parameters, true.model, models, N))
      else 
        newdata <- suppressMessages(GenerateDataMulti(n, true.model, parameters, true.model, models, N))
      fits.boot <- suppressMessages(EstimateParametersMulti(true.model, newdata))
      if(inherits(true.model,"paleoGRW"))
      {
        if("ms" %in% names(true.model$fixed.parameters))
        {
          params.model <- sapply(c("anc", "vs"), function(x) parameters.full[[x]])
          params.boot <- matrix(unlist(sapply(fits.boot,function(x) x[c('anc', 'vs')])),
                              nrow = 2,dimnames=list(c('anc', 'vs'), c()))
        }
        else
        {
          params.model <- sapply(c("anc", "ms", "vs"), function(x) parameters.full[[x]])
          params.boot <- matrix(unlist(sapply(fits.boot,function(x) x[c('anc', 'ms', 'vs')])),
                                nrow = 3,dimnames=list(c('anc', 'ms', 'vs'), c()))
        }
      }
      if (inherits(true.model,"paleoStasis"))
      {
        params.model <- sapply(c("theta", "omega"), function(x) parameters.full[[x]])
        params.boot <- matrix(unlist(sapply(fits.boot,function(x) x[c('theta', 'omega')])),
                                    nrow = 2,dimnames=list(c('theta','omega'),c()))
      }
      if(p == 1)
        param.boot <- params.boot %>% mean
      else
      {
        param.boot <- params.boot%>% rowMeans
        if(inherits(true.model,"paleoStasis"))
          param.boot['omega'] = min(param.boot['omega'], 0.01)
        if(inherits(true.model,"paleoGRW"))
          param.boot['vs'] = min(param.boot['vs'], 0.01)
      }
    param.corrected =  as.list(params.model-param.boot)
    param.corrected$vp = parameters.full$vp
    param.corrected$nn = parameters.full$nn
    param.corrected$ns = parameters.full$ns
    param.corrected$tt = parameters.full$tt
    if("ms" %in% names(true.model$fixed.parameters)) 
      param.corrected$ms=0
    return(structure(list(parameters = param.corrected, data = newdata),
                     class = 'ecicBiasCorrect'))
    }
  }
}
