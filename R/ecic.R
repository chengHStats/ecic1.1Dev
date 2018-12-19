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
ECIC = function(data, models, alpha = 0.05, N = 1000, ic = 'AIC'){
  names(models) <- models
  p <- length(models)
  n <- length(data)

  obs = sapply(models, function(x) IC(data, x, ic)) #observed
  params.obs = sapply(obs, function(x) x[-1])
  scores.obs = sapply(obs, function(x) x[1])

  best.ix = which.min(scores.obs)
  best = models[best.ix]
  alt.models = models[-best.ix]

  dif.obs = scores.obs[best.ix]-min(scores.obs[-best.ix])
  names(dif.obs) = best
  bc = lapply(alt.models, function(x)
    BiasCorrect(n, x, params.obs[[x]],
                models, N, ic))

  icd = lapply(alt.models, function(x)
    ecicControl(n, x, bc[[x]]$parameters, best, models, N, ic))

  best.freq = sapply(icd, function(x) x$frequencies[best.ix])
  alpha.primes = sapply(best.freq,
                        function(x) ifelse(x==0, 1, alpha/x))
  alpha.primes.N = round(alpha.primes * N)
  differences = sapply(icd, function(x) x$differences)
  ecic.thresholds = sapply(1:(p-1), function(x)
    differences[alpha.primes.N[x],x]
    )
  names(ecic.thresholds) = alt.models
  decision.ecic = ifelse(dif.obs < min(ecic.thresholds),
                          best, "No Decision")%>%unname
  decision.ba = ifelse(dif.obs < -10,
                       best, "No Decision")%>%unname


  return(
    list(
      decisions = list(
        ecic = decision.ecic,
        ba = decision.ba,
        ic = best,
        ecic.thresholds = ecic.thresholds
      ),
      observed = list(delta = dif.obs,
        scores = scores.obs,
        parameters = params.obs,
        data = data
      )%>% structure(class = c("ecic", "decisions")),
      bias.correction = bc %>% structure(class = c("ecic", "biascorrection")),
      error.control = icd
    )%>% structure(class = c("ecic", "errorcontrol"))
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



