#' @export
parameterCheck = function(model, parameters){
parameters = as.list(parameters)

  # Check for extraneous parameters.
  for (parameter in names(parameters) ){
    if ( ! ( parameter %in% model$parameter.names ) ) {
      message(paste("Model", model$name, "has parameters",
                    paste(model$parameter.names, collapse = ", "),
                    "so the argument", parameter, "was ignored."))
      parameters[[parameter]] <- NULL
    }
  }
  # Check for duplicate parameters.
  for (parameter in names(parameters)){
    if (parameter %in% names(model$fixed.parameters)){
      message(paste("The given model has a fixed value of", parameter, "=",
                    model$fixed.parameters[[parameter]], "so the argument",
                    parameter,"=", parameters[parameter], "was ignored."))
      parameters[[parameter]] <- NULL
    }
  }
  # Check if all parameters are supplied
  parameters <- c(parameters, model$fixed.parameters)

parameters
}
