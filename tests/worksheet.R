# Specify model set (norm(0,1), norm(mu, 1), norm(mu, sd)).
models = c("norm0", "norm1", "norm")

# Generate a sample of size 25 from norm(0.3, 1.2).
data = GenerateData(25, "norm", c(0.3, 1.2))

# Estimate model parameters for norm(mu, 1).
EstimateParameters(data, "norm1")

# Perform the ECIC procedure on this data with alpha = 0.05 and N = 100.
test.ecic = ECIC(data, models, 0.05, 100)

# Get the decisions from the ECIC procedure.
test.ecic$decisions

# Get the observed score difference and thresholds.
test.ecic$observed$delta
test.ecic$thresholds

# Plot the Difference Distributions from the error control step
test.ecic$error.control$norm0$differences %>% density %>% plot
test.ecic$error.control$norm1$differences %>% density %>% lines

# Methods for the ecic class
summary(test.ecic)
plot(test.ecic)
