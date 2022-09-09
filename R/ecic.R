length.paleoTS = function(data){
  length(data$mm)
}
#' Perform the Error Control for Information Criteria (ECIC) Procedure.
#'
#' @param data: A vector/matrix data sample compatible with the specified models.
#' @param models: A string vector specifying the models to use.
#' @param alpha: A vector of desired type 1 error rates 
#'
#' @return An ecic object with the following attributes:

#'    \code{decisions:  }A dataframe containing the results of the ECIC procedure along
#'    with the decisions given by the best scoring model and the Burnham &
#'    Anderson (BA) rule (accept best model if delta < -10)
#'
#'    \code{observed:  } An ecicObserved object containing the observed IC delta,
#'    the information criterion values for each model given the data, fitted
#'    parameters for each model, and the observed data itself
#'
#'    \code{bias.correction:  }An ecicBias object containing the corrected model
#'    parameters and the data used for the bias correction procedure
#'
#'    \code{error.control:  }A list of ecicControl objects corresponding to each
#'    alternative model where each ecicControl object contains information related
#'    to the error control bootstrap procedure: a vector of model bootstrap
#'    frequencies, a vector of bootstrap IC deltas, estimated parameters
#'    corresponding to each bootstrap sample, a matrix of model scores, and the
#'    bootstrap data itself
#'
#'
#' @examples
#' my.data = GenerateData(25, "norm", c(mu = 0.3, sd = 1.2))
#' my.models = c("norm0", "norm1" "norm")
#' ECIC(my.data, my.models, alpha = 0.05, N = 1000, ic = "AIC")
#'
#' @export
ECIC = function(models, data, alpha = c(0.01, 0.05, 0.1), N = 1000, ic = 'AIC', genBest = TRUE){
  # create an ecicModelList object out of the model set provided as input
  models = ecicModelList(models)
  # store the number of models in the model set
  p <- length(models)
  # store the number of data sets inputted
  n <- length(data)
  # compute the information criterion value for each model... outputs list of IC value and MLE parameters 
  obs = lapply(models, function(x) IC(x, data, ic)) #observed
  # what does inherits do exactly?? (ask Beckett)
  if (!inherits(data, "paleoTS")) 
  {
    params.obs = lapply(obs, function(x) x$parameters) # extract the MLE parameters for each model in the model set
    scores.obs = sapply(obs, function(x) x$ic) # extract the IC value for each model in the model set
  } 
  else 
  {
    n = length(data$nn) ## potential bug here since n is already set to length(data) before the if-else condition
    params.obs = lapply(obs, function(x) x$parameters)
    scores.obs = sapply(obs, function(x) x$ic)
  }
  # compute AIC weights (for Burnham and Anderson's rule of thumb approach) for each model in the model set
  weights.obs = AICweights(scores.obs)
  # determine which index holds the model with the lowest IC value
  best.ix = which.min(scores.obs)
  # use index to store the model with the lowest IC value 
  best = models[[best.ix]]
  # create a new model set without the model with the lowest IC value
  alt.models = models[-best.ix]
  # take the score of the best model and subtract the best remaining score (where smaller is better)
  dif.obs = scores.obs[best.ix]-min(scores.obs[-best.ix]) 
  # take the ratio of the observed best model's Akaike weight to the next best Akaike weight
  ratio.obs = weights.obs[best.ix]/max(weights.obs[-best.ix]) 
  # store the name of the best model (as selected by using the lowest IC value) as the name of dif.obs
  names(dif.obs) = best$ID

  bc = lapply(alt.models, function(x)
    BiasCorrect(n, x, params.obs[[x$ID]],
                models, N, ic, genBest))

  icd = lapply(alt.models, function(x) ecicControl(n, x, bc[[x$ID]]$parameters, best, models, N, ic))


  best.freq = sapply(icd, function(x) x$frequencies$frequencies[best.ix])
  alpha.primes = lapply(alpha,
                        function(a)  sapply(best.freq, function(x) ifelse(x==0, 1, a/x)))
  alpha.primes = lapply(alpha.primes, function(x) sapply(x, function(a) min(a, 1)))
  alpha.primes.N = lapply(alpha.primes, function(x) round(x * N))
  differences = lapply(icd, function(x) x$differences)
  ratios = lapply(icd, function(x) x$ratios)

  ecic.thresholds = lapply(alpha.primes, function(alpha.prime) {
    thresh = sapply(1:(p-1), function(x){
    try({
    # dif = differences[[x]]
    # dif.N = length(dif)
    # out = dif[ceiling(alpha.prime[x] * dif.N)]
    # if(alpha.prime[x] == 1) out = 0
    # out
    ratio = ratios[[x]]
    ratio.N = length(ratio)
    out = ratio[floor(ratio.N - (alpha.prime[x] * ratio.N))]
    if(alpha.prime[x] == 1) out = 1
    out


    })
    })
    names(thresh) = names(alt.models)
    thresh
    }
  )
  ecic.thresholds = lapply(ecic.thresholds, function(x)  sapply(x, function(y) ifelse(is.na(y), 1, y)))
  ecic.thresholds = lapply(ecic.thresholds, function(x)  sapply(x, function(y) ifelse(length(y) == 0, 1, y)))

  names(ecic.thresholds) = alpha

  decision.ecic = sapply(ecic.thresholds, function(thresh) ifelse(dif.obs > max(thresh),
                          best$ID, "No Decision")%>%unname)

  decision.ba = ifelse(ratio.obs > 2.7,
                       best$ID, "No Decision")%>%unname


  return(
    list(
      decisions = list(
        ecic = decision.ecic,
        ba = decision.ba,
        ic = best$ID,
        thresholds = ecic.thresholds),
      observed = list(delta = dif.obs,
        scores = scores.obs,
        ratio = ratio.obs,
        parameters = lapply(bc, function(x) x$parameters),
        data = data
      )%>% structure(class = c("ecicDecisions")),
      bias.correction = bc ,
      error.control = icd)
  ) %>% structure(class = "ecic")
  }

#' @export
summary.ecic = function(ecic){
  print("Summary Method")
}
#' @export
plot.ecic = function(ecic){
  print("Plot Method")
}

plot.ecicBiasCorrect = function(bc){

}
plot.ecicFrequencies = function(ec){


}



