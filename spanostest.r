library(ECIC)
grw = ecicModel("rwalk", ID = "grw")
urw = ecicModel("rwalk",
                ID="urw",
                fix=list(step.mean=0))
stasis = ecicModel("norm",
                   ID="stasis")
models = ecicModelList(list(
  grw,
  urw,
  stasis
))

true.model = grw
true.param = list(step.mean=0.1, step.sd = 1.1)

data = GenerateData(n, true.model, true.param)
ECIC(models,
     data,
     alpha = 0.05,
     N = 1000,
     ic = "AIC")

