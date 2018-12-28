#' @export
ecicModel <- function(model.name, ID = model.name, fix = list()){
  UseMethod("ecicModel", model.name)
}
#' @export
ecicModel.character <- function(model.name, ID = model.name, fix = list()){
  assert_that(is.character(model.name))
  assert_that(is.character(ID))
  assert_that(is.list(fix))
  class(model.name) = model.name
  UseMethod("ecicModel", model.name)
}
#' @export
ecicModel.norm <- function(model.name, ID = model.name, fix = list()) {
  parameter.names <- c("mean", "sd") # Define parameters for model.

  # Error handling for fixed parameters.
  if(length(fix) > 0){
  # Check if fixed parameters have correct names.
    for (parameter in names(fix) ){
        if ( ! ( parameter %in% parameter.names ) ) {
          message(paste("Model", model.name, "has parameters",
                      paste(parameter.names, collapse = ", "),
                      "so the argument", parameter, "was ignored."))
          fix$parameter <- NULL
        }
    }
    # Assert that sd is a single numeric value greater than zero.
    if ("sd" %in% names(fix)){
     assert_that(is.numeric(fix$sd), length(fix$sd) == 1, fix$sd > 0)
    }
    # Assert that mean is a single numeric value.
    if("mean" %in% names(fix)){
      assert_that(is.numeric(fix$mean), length(fix$mean) == 1)
    }
  } # End error handling for fixed parameters.

  k <- 2 - length(fix)

  return(structure(list(ID = ID, name = "norm", parameter.names = parameter.names,
                        fixed.parameters = fix,
                        k = k, data.type = 1), class = c("ecicModel", "norm")))
}
#' @export
ecicModel.rwalk = function(model.name, ID = model.name, fix = list()){
  parameter.names <- c("step.mean", "step.sd") # Define parameters for model.
  # Error handling for fixed parameters.
  if(length(fix) > 0){
    # Check if fixed parameters have correct names.
    for (parameter in names(fix) ){
      if ( ! ( parameter %in% parameter.names ) ) {
        print(paste("Model", model.name, "has parameters",
                    paste(parameter.names, collapse = ", "),
                    "so the argument", parameter, "was ignored."))
        fix$parameter <- NULL
      }
    }
    # Assert that sd is a single numeric value greater than zero.
    if ("step.sd" %in% names(fix)){
      assert_that(is.numeric(fix$step.sd), length(fix$step.sd) == 1, fix$sd > 0)
    }
    # Assert that mean is a single numeric value.
    if("step.mean" %in% names(fix)){
      assert_that(is.numeric(fix$step.mean), length(fix$step.mean) == 1)
    }
  } # End error handling for fixed parameters.

  k <- 2 - length(fix)

  return(structure(list(ID = ID, name = "rwalk", parameter.names = parameter.names,
                        fixed.parameters = fix,
                        k = k, data.type = 1), class = c("ecicModel", "rwalk")))

}
#' @export
ecicModel.lm <- function(model.name, ID, fix = list()){
  parameter.names = c(names(coef(model.name)), "sd")
  k = length(parameter.names)-1
  n = length(model.name$fitted.values)
  return(structure(list(ID = ID,
                        name = "lm",
                        parameter.names = parameter.names,
                        fixed.parameters = list(),
                        k = k,
                        n = n,
                        model.matrix = model.matrix(model.name),
                        data.type = 1
                        ), class = c("ecicModel", "lmECIC")))
}
#' @export
ecicModelList = function(models = list()){
  models = lapply(models, function(x){
    if (is.character(x)){
      tryCatch(ecicModel(x, ID = x),
               error = function(e){
                 stop(paste(x, "is not a valid ecicModel type."),
                      call. = FALSE)
               })
    } else if (class(x)[1] == "ecicModel"){
      x
    }
    else{
      stop("Unexpected input in model list.")
    }
  })
  # Check if all models have different IDs.
  IDs = sapply(models, function(x) x$ID)
  if (length(unique(IDs)) < length(models)){
    stop("Duplicate IDs among models. Change the IDs such that each is unique.")
  }
  names(models) = IDs
  # Check if all models have the same data type.
  data.types = sapply(models, function(x) x$data.type)
  if(diff(range(data.types)) != 0){
    stop("Specified models have incompatible data types.")
  }
  ns = lapply(models, function(x) x[["n"]])

  # Check if all models have same data space.

  return(models)
}
