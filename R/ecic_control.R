#' Bootstrap distributions of Information Criterion Differences
#'
#' @param n: Sample size.
#' @param true: A string giving the name of the model to generate data from.
#' @param param: A vector of parameters specifying a fitted true model.
#' @param best: A string giving the name of the model to be given as best.
#' @param models: A string vector of model names for selecting between.
#' @param N: Number of bootstrap samples.
#' @param ic: String giving which information criterion to use (e.g. 'AICc', 'BIC')
#
#' @return An ecicControl object containing the following:
#'
#'   differences:
#'
#'   frequencies: a modelFrequencies object containing bootstrap frequencies,
#'   fitted parameters for the bootstrap samples, scores for the bootstrap samples
#'   and the samples themselves.
#'
#' @examples
#' my.models = c("norm0" "norm1", "norm")
#' ecicControl(n = 25, true = "norm", param = c(0.3, 1.2),
#'             my.models, N = 1000, ic = 'AIC')
#'
#' @export
ecicControl = function(n, true, parameters,best, models, N = 1000, ic = 'AIC'){
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

  models = ecicModelList(models)
  params.true <- NULL
  true.ix <- which(names(models)==true$ID)
  best.ix <- which(names(models)==best$ID)
  mf1 = suppressMessages(ModelFrequencies(n, true, parameters, models, N, ic))

  if(mf1$frequencies[best.ix] == 0){
    return(structure(list(differences = NA,
                          est.parameters = NA, scores = NA,
                          data.control = NA, frequencies = mf1),
                     class = 'ecicControl'))
  }

  mins = apply(mf1$ic.scores, 1, function(x) which.min(x))
  if(true$data.type!="paleoTS"){
  params.true = mf1$parameters[[true.ix]] %>% as.matrix



  n.best <- sum(mins==best.ix)
  which.best <- mins==best.ix

  scores.keep = mf1$ic.scores[which.best,]
  data.keep = mf1$data[, which.best]

  if(!is.null(dim(params.true))){
    params.keep = params.true[, which.best]
  }

  data2 <- suppressMessages(GenerateDataBest(n, true, parameters, best, models, N-(n.best), ic))
  icmm <- ICMultiMulti(models, data2, ic)
  scores2 = icmm$ic
  scores.full <- rbind(scores.keep,scores2)



  data.full = cbind(data.keep, data2)
  if(min(is.na(data.full))==1){
    return(structure(list(differences = c(0),
                          est.parameters = NA, scores = NA,
                          data.control = NA, frequencies = mf1),
                     class = 'ecicControl'))
  }
  params.true2 = icmm$parameters[[true.ix]] %>% as.matrix
  rownames(scores.full) <- NULL
  dif <- sort(scores.full[,best.ix] - apply(as.matrix(scores.full[,-best.ix]), 1, min),
              method = 'radix')

  weights.full <- AICweights(scores.full)
  ratio <- sort(weights.full[,best.ix] / apply(as.matrix(weights.full[,-best.ix]), 1, max),
         method = 'radix')

  #if(!is.null(params.full)){
  params.full = NULL
  if(!is.null(params.true)){
    params.full = cbind(params.true, params.true2)
  }
  # if(k == 1){
  #   params.true = c(estscores1.keep[[true]][-1,], estscores2[[true]][-1,])
  #
  # }else{
  #   params.true = cbind(estscores1.keep[[true]][-1,], estscores2[[true]][-1,])
  #
  # }
  #}
  #}
  }else{
    mins = apply(mf1$ic.scores, 1, function(x) which.min(x))

    params.true = mf1$parameters[true.ix][[1]]



    n.best <- sum(mins==best.ix)
    which.best <- mins==best.ix

    scores.keep = mf1$ic.scores[which.best,]
    data.keep = mf1$data[which.best]

    if(!is.null(params.true)){

      params.keep = params.true[which.best]
    }

    data2 <- suppressMessages(GenerateDataBest(n, true, parameters, best, models, N-(n.best), ic))
    icmm <- ICMultiMulti(models, data2, ic)
    scores2 = sapply(icmm, function(x) sapply(x, function(y) y$ic))
    scores.full <- rbind(scores.keep,scores2)

    params.true2 = lapply(icmm[[true.ix]], function(x) x$parameters)
    rownames(scores.full) <- NULL
    dif <- sort(scores.full[,best.ix] - apply(as.matrix(scores.full[,-best.ix]), 1, min),
                method = 'radix')

    data.full = c(data.keep, data2)

    #if(!is.null(params.full)){
    params.full = NULL
    if(!is.null(params.true)){
      params.full = c(params.keep, params.true2)

    }
  }
  return(structure(list(differences = dif, ratios = ratio,
              est.parameters = params.full, scores = scores.full,
              data.control = data.full, frequencies = mf1),
              class = 'ecicControl'))

}
