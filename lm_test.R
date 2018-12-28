x.lm <- 1:10
noise.lm <- rnorm(10, 0, 5.4)
y.lm <- 1.4*x.lm + 0.02*(x.lm^2)+ 0.06*(x.lm^3) + noise.lm
lm1 <- lm(y.lm ~ x.lm)
lm2 <- lm(y.lm ~ poly(x.lm, 2, raw = TRUE))
lm3 <- lm(y.lm ~ poly(x.lm, 3, raw = TRUE))

model1 <- ecicModel(lm1, "mylm1")
model2 <- ecicModel(lm2, "mylm2")
model3 <- ecicModel(lm3, "mylm3")
models.lm <- list(model1, model2, model3)

parameters <- EstimateParameters(model2, y.lm)$parameters

GenerateData(25, model2, parameters) # Returns error
data = GenerateData(10, model2, parameters)

ecic.lmtest = ECIC(models.lm, data, N = 1000)
