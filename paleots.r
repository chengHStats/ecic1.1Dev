library(paleoTS)
library(ECIC)
grw = ecicModel("paleoGRW", ID="grw")
urw = ecicModel("paleoGRW", ID="urw", fix = list(ms = 0))
stasis = ecicModel("paleoStasis", ID="stasis")

models = ecicModelList(list(grw, urw, stasis))

paleoData = sim.GRW(ms = .1,vs=0.8, ns = 10)


parameters = EstimateParameters(grw, paleoData)
ModelFrequencies(10, grw, parameters, models, 100)$frequencies
data = GenerateData(6, grw, parameters)


  library(tictoc)
  tic()
  paleoECICtest = ECIC(models, data, N=50)
  toc()

