
# Estimate parameters for a single sample --------------------------------------


# EstimateParameters -----------------------------------------------------------
#' Computes the fitted parameters for a data sample and a model.
#'
#' @param data: compatible with the specified model.
#' @param model: A string specifying the model to compute the log-likehood for.
#'
#' @return A vector/matrix of the fitted parameters.
#'
#' @examples
#' data = GenerateData(25, "norm", c(mu = 0.3, sd = 1.2))
#' models = c("norm0", "norm1", "norm")
#' EstimateParameters(data)
#' @export
EstimateParameters <- function(model, data){
  UseMethod("EstimateParameters", model)
  # do.call(paste("EstimateParameters.", model, sep = ""),
  #         list(data = model, data = model))
}
#' @export
EstimateParameters.character <- function(model, data){
 model2 = tryCatch(ecicModel(model),
                   error = function(e){
                     stop(paste(model, "is not a valid ecicModel type."), call. = FALSE)
                   })
 EstimateParameters(model2, data)
}
#' @export
EstimateParameters.ecicModel <- function(model, data){

  # Check that the data is appropriate to the model.
  if (model$data.type == 1){
    assert_that(is.numeric(data),
                (is.vector(data) | (is.matrix(data) && ncol(data) == 1 ) ))
  }
  NextMethod("EstimateParameters", model)
}
# EstimateParameters.norm ------------------------------------------------------
#' @export
EstimateParameters.norm <- function(model, data){
  # Computes the fitted parameters for a data sample and a norm(mu, sd) model.
  n <- length(data)
    if(is.null(model$fixed.parameters$mean)){
    mean <- sum(data)/n
  } else{
    mean <- model$fixed.parameters$mean
  }
  if(is.null(model$fixed.parameters$sd)){
    sd <- (crossprod(data-mean)/n)**.5
  } else{
    sd <- model$fixed.parameters$sd
  }
  return(list(parameters = setNames(c(mean, sd), c("mean", "sd"))))
}

# EstimateParameters.rwalk -----------------------------------------------------------
#' @export
EstimateParameters.rwalk <- function(model, data){
  # Computes fitted parameters for a data sample and random walk (mu, sd) model.
  if(length(model$fixed.parameters) < 2){
    n <- length(data)
    steps <- c(data[1], diff(data))
    if(is.null(model$fixed.parameters$step.mean)){
      step.mean <- sum(steps)/n
    }
    if(is.null(model$fixed.parameters$step.sd)){
      step.sd <- (crossprod(steps-step.mean)/n)**.5
    }
  } else {
    step.mean = model$fixed.parameters$step.mean
    step.sd = model$fixed.parameters$step.sd
  }
  return(list(parameters = setNames(c(step.mean, step.sd),
                                    c("step.mean", "step.sd")),
              steps = steps))
}
#' @export
EstimateParameters.lmECIC = function(model, data){
  if(length(data) != model$n){
    stop(paste("Data is not correct length for ", model$ID,
               ", expecting length ", n, ", but input is length ",
               length(data), ".", sep = ""))
  }
  coef = crossprod(model$model.matrix) %>% solve %>%
    crossprod(crossprod(model$model.matrix,data))
  fit = model$model.matrix %*% coef
  resid = data - fit
  sd = sd(resid)
  out = c(coef, sd) %>% setNames(model$parameter.names)
  return(list(parameters = out, resid = resid))
}

# EstimateParametersMulti -------------------------------------
#' Computes the fitted parameters for multiple data samples and a model.
#'
#' @param data: A matrix of data samples compatible with the specified model.
#' @param model: A string specifying the model to compute the log-likehood for.
#'
#' @return A vector/matrix of the fitted parameters.
#'
#' @examples
#' data = GenerateDataMulti(25, "norm", c(mu = 0.3, sd = 1.2), 10)
#' models = c("norm0", "norm1", "norm")
#' EstimateParameters(data)
#'
#' @export
EstimateParametersMulti <- function(model, data, flat = F){
  # Computes the fitted parameters for multiple data samples and a model.
  #
  # Args:
  #   data: compatible with the specified model.
  #   model: A string specifying the model to compute the log-likehood for
  #
  # Returns:
  #   A vector/matrix of the fitted parameters.
  UseMethod("EstimateParametersMulti", model)
}
#' @export
EstimateParametersMulti.character <- function(model, data){
  model2 = tryCatch(ecicModel(model),
                    error = function(e){
                      stop(paste(model, "is not a valid ecicModel type."), call. = FALSE)
                    })
  EstimateParametersMulti(model2, data)
}
EstimateParametersMulti.ecicModel <- function(model, data, flat = F){
  if (model$data.type == 1){
    assert_that(is.numeric(data), is.matrix(data))
  }
  NextMethod()
}
#' @export
EstimateParametersMulti.norm <- function(model, data, flat = F, ...){
  # Computes fitted parameters for multiple data samples and norm(mu, sd) model.

  n <- nrow(data)
  N <- ncol(data)
  if(is.null(model$fixed.parameters$mean)){
    means <- crossprod(rep(1/n, n), data)
  } else {
    means = rep(model$fixed.parameters$mean, N)
  }
  if(is.null(model$fixed.parameters$sd)){
  r <-crossprod(data - matrix(rep(means, each = n), ncol = N))
  sds <- {diag(r)/n}**.5
  } else {
    sds = rep(model$fixed.parameters$sd, N)
  }
  out = if (flat) {
    c(means,sds)
  } else {
    matrix(c(means, sds), nrow = 2, byrow = TRUE, dimnames = list(c("mean", "sd"), NULL))
  }
  return(list(parameters = out))
}
#' @export
EstimateParametersMulti.rwalk <- function(model, data, flat = F, ...){
  # Computes fitted parameters for multiple data samples and norm(mu, sd) model.


  if(is.null(model$fixed.parameters$step.mean) |
     is.null(model$fixed.parameters$step.sd)){
    steps = rbind(data[1,], diff(data))
  }
  if(is.null(model$fixed.parameters$step.mean)){
    step.means <- crossprod(rep(1/n, n), steps)
  } else {
    step.means = rep(model$fixed.parameters$step.mean, N)
  }
  if(is.null(model$fixed.parameters$step.sd)){
    r <-crossprod(steps - matrix(rep(step.means, each = n), ncol = N))
    step.sds <- {diag(r)/n}**.5
  } else {
    step.sds = rep(model$fixed.parameters$step.sd, N)
  }
  out = if (flat) {
    c(step.means,step.sds)
  } else {
    matrix(c(step.means, step.sds), nrow = 2, byrow = TRUE,
           dimnames = list(c("step.mean", "step.sd"), NULL))
  }
  return(parameters = out, steps = steps)
}
#' @export
EstimateParametersMulti.lmECIC = function(model, data){
  n <- nrow(data)
  N <- ncol(data)
  if(n != model$n){
    stop(paste("Data is not correct length for ", model$ID,
               ", expecting length ", model$n, ", but input is length ",
               n, ".", sep = ""))
  }
  x = as.matrix(model$model.matrix[,-1])
  m = crossprod(model$model.matrix) %>% solve
  coefs = apply(data, 2, function(x) crossprod(m, crossprod(model$model.matrix,x)))
  fit = model$model.matrix %*% coefs
  resid = data-fit
  sds = apply(resid,2,sd)
  out = rbind(coefs, sds)
  rownames(out) = model$parameter.names
  return(list(parameters = out, residuals = resid))
}


