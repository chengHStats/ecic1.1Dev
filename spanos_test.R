spanos = ECIC:::datasets$spanos

lm.spanos.1 = lm(y ~ x, data = spanos)
lm.spanos.2 = lm(y ~ poly(x, 2, raw = TRUE), data = spanos)
lm.spanos.3 = lm(y ~ poly(x, 3, raw = TRUE), data = spanos)

spanos.1 = ecicModel(lm.spanos.1, "spanos.1")
spanos.2 = ecicModel(lm.spanos.2, "spanos.2")
spanos.3 = ecicModel(lm.spanos.3, "spanos.3")

models.spanos = ecicModelList(list(spanos.1, spanos.3))

parameters = EstimateParameters(spanos.1, spanos.y)$parameters
data.spanos = GenerateData(35, spanos.1, parameters)


ECIC(models.spanos, data.spanos, N = 1000)
system.time(ECIC(models.spanos, data.spanos, N = 1000))

spanos.sim = lapply(1:1000, function(x){
  data.spanos = GenerateData(35, spanos.1, parameters)
  ECIC(models.spanos, data.spanos, N = 1000)
  })
spanos.sim2 = lapply(1:4000, function(x){
  data.spanos = GenerateData(35, spanos.1, parameters)
  ECIC(models.spanos, data.spanos, N = 1000)
 })
sapply(spanos.sim, function(x) x$decisions)
sapply(spanos.sim, function(x) x$observed)

# ECIC Decision
sapply(spanos.sim, function(x) x$decisions[1]) %>% unlist %>% table

# BA Decision
sapply(spanos.sim, function(x) x$decisions[2]) %>% unlist %>% table

# Naive AIC Implementation
sapply(spanos.sim, function(x) x$decisions[3]) %>% unlist %>% table
