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

BiasCorrect <- function(n, true, param, models, N = 1000, ic = 'AIC'){
  # Corrects biased parameter estimates stemming from partitioning by best model
  #
  # Args:
  #   data: A matrix/array of data samples compatible with the specified model.
  #   true: A string giving the name of the model to generate data from.
  #   models: A string vector of model names for selecting between.
  #   param: A vector of model parameters.
  #   N: Number of bootstrap samples.
  #   ic: String giving which information criterion to use (e.g. 'AICc', 'BIC')
  #   compress: A boolean specifying if output should be of list or vector type
  #
  # Returns:
  #   param: The corrected model parameters.
  names(models) <- models
  p = length(param)
  true.ix = which(models == true)
  if(p > 0){
  newdata <- GenerateDataBest(n, true, param, true, models, N)
  params.boot = EstimateParametersMulti(newdata, true)

  if(p == 1){
  param.boot <- params.boot %>% mean
  } else{
  param.boot <- params.boot%>% colMeans
  }
  param.corrected = 2*param-param.boot
  return(list(parameters = param.corrected, data = newdata))
  }
}

ICDists = function(n, true, param,best, models, N = 1000, ic = 'AIC'){
  # Estimates distributions of IC differences between models from a fitted model.
  #
  # Args:
  #   n: Sample size.
  #   true: A string giving the name of the model to generate data from.
  #   models: A string vector of model names for selecting between.
  #   param: A vector of model parameters.
  #   best: A string given the name of the model that should be chosen as best.
  #   N: Number of bootstrap samples.
  #   ic: String giving which information criterion to use (e.g. 'AICc', 'BIC')
  #
  # Returns:
  #   differences: A vector giving the empirical distribution of IC Differences.
  #   est.parameters: Parameter estimates for true model for each sample.
  #   scores: Matrix of information criterion values for each model and sample.
  p <-  length(models)

  names(models) <- models
  params.true <- NULL
  best.ix <- which(models==best)
  mf1 = ModelFrequencies(n, true, param, models, N, ic)

  mins = apply(mf1$ic.scores, 1, function(x) which.min(x))
  params.true = mf1$parameters[[true]] %>% as.matrix


  n.best <- sum(mins==best.ix)
  which.best <- mins==best.ix

  scores.keep = mf1$ic.scores[which.best,]
  data.keep = mf1$data[, which.best]

  if(!is.null(dim(params.true))){

    params.keep = params.true[which.best,]
  }

  data2 <- GenerateDataBest(n, true, param, best, models, N-(n.best), ic)
  icmm <- ICMultiMulti(data2, models, ic)
  scores2 = icmm$ic
  scores.full <- rbind(scores.keep,scores2)

  params.true2 = icmm$parameters[[true]] %>% as.matrix
  rownames(scores.full) <- NULL
  dif <- sort(scores.full[,best.ix] - apply(scores.full[,-best.ix], 1, min),
             method = 'radix')

  data.full = cbind(data.keep, data2)

  #if(!is.null(params.full)){
    params.full = NULL
    if(!is.null(params.true)){
      params.full = rbind(params.true, params.true2)
    }
    # if(k == 1){
    #   params.true = c(estscores1.keep[[true]][-1,], estscores2[[true]][-1,])
    #
    # }else{
    #   params.true = cbind(estscores1.keep[[true]][-1,], estscores2[[true]][-1,])
    #
    # }
  #}
  return(list(differences = dif, frequencies = mf1$frequencies,
              est.parameters = params.full, scores = scores.full,
              data = data.full, data.bc = mf1$data))
}




