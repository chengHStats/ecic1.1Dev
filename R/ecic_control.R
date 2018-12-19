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
ecicControl = function(n, true, param,best, models, N = 1000, ic = 'AIC'){
  # Estimates distributions of IC differences between models from a fitted model.
  #
  # Args:
  #   n: Sample size.adsfasdfdfs
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
              data.control = data.full, data.freq = mf1$data))
}
