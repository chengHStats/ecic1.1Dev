length.paleoTS = function(data){
  length(data$mm)
}
#' Perform the Error Control for Information Criteria (ECIC) Procedure.
#'
#' @param data: A vector/matrix data samplecompatible with the specified models.
#' @param model: A string vector specifying the models to use.
#'
#' @return An ecic object with the following attributes:
#'    \code{decisions:  }A dataframe containing the results of the ECIC procedure, along
#'    with the decisions given by the best scoring model and the Burnham &
#'    Anderson (BA) rule (accept best model if delta < -10).
#'
#'    \code{observed:  } an ecicObserved object containing the observed IC delta,
#'    the information criterion values for each model given the data, fitted
#'    parameters for each model, and the observed data itself.
#'
#'    \code{bias.correction:  }An ecicBias object containing the corrected model
#'    parameters and the data used for the bias correction procedure.
#'
#'    \code{error.control:  } A list of ecicControl objects corresponding to each
#'    alternative model. Each ecicControl object contains information related
#'    to the error control bootstrap procedure: a vector of model bootstrap
#'    frequencies, a vector of bootstrap IC deltas, estimated parameters
#'    corresponding to each bootstrap sample, a matrix of model scores, and the
#'    bootstrap data itself.
#'
#'
#' @examples
#' my.data = GenerateData(25, "norm", c(mu = 0.3, sd = 1.2))
#' my.models = c("norm0", "norm1" "norm")
#' ECIC(my.data, my.models, alpha = 0.05, N = 1000, ic = "AIC")
#'
#' @export
ECIC = function(models, data, alpha = 0.05, N = 1000, ic = 'AIC'){
  models = ecicModelList(models)
  p <- length(models)
  n <- length(data)

  obs = lapply(models, function(x) IC(x, data, ic)) #observed

  if (class(data)!= "paleoTS"){
  params.obs = lapply(obs, function(x) x$parameters)
  scores.obs = sapply(obs, function(x) x$ic)
  } else {
    params.obs = lapply(obs, function(x) x$parameters)
    scores.obs = sapply(obs, function(x) x$ic)

  }


  best.ix = which.min(scores.obs)
  best = models[[best.ix]]

  alt.models = models[-best.ix]

  dif.obs = scores.obs[best.ix]-min(scores.obs[-best.ix])
  names(dif.obs) = best$ID
  bc = suppressMessages(lapply(alt.models, function(x)
    BiasCorrect(n, x, params.obs[[x$ID]],
                models, N, ic)))

  icd = lapply(alt.models, function(x)
    ecicControl(n, x, bc[[x$ID]]$parameters, best, models, N, ic))

  best.freq = sapply(icd, function(x) x$frequencies$frequencies[best.ix])
  alpha.primes = sapply(best.freq,
                        function(x) ifelse(x==0, 1, alpha/x))
  alpha.primes = sapply(alpha.primes, function(x) min(x, 1))
  alpha.primes.N = round(alpha.primes * N)
  differences = sapply(icd, function(x) x$differences)
  ecic.thresholds = sapply(1:(p-1), function(x){
    try({
    dif = differences[[x]]
    dif.N = length(dif)
    out = dif[round(alpha.primes[x] * dif.N)]
    if(alpha.primes[x] == 1) out = 0
    out
    })
  }
    )
  ecic.thresholds = sapply(ecic.thresholds, function(x) ifelse(is.na(x), 0, x))
  names(ecic.thresholds) = names(alt.models)
  decision.ecic = ifelse(dif.obs < min(ecic.thresholds),
                          best$ID, "No Decision")%>%unname

  decision.ba = ifelse(dif.obs < -2,
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



