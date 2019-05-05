library(ECIC)
library(paleoTS)


grw = ecicModel("paleoGRW", ID="grw")
urw = ecicModel("paleoGRW", ID="urw", fix = list(ms = 0))
stasis = ecicModel("paleoStasis", ID="stasis")

models = ecicModelList(list(grw, urw, stasis))

paleoData = sim.Stasis(ns = 6)
parameters = EstimateParameters(stasis, paleoData)

EstimateParameters(grw, paleoData)
EstimateParameters(urw, paleoData)
EstimateParameters(stasis, paleoData)


ECIC(models, paleoData, N = 100)


paleoData = sim.Stasis(ns = 6)
parameters = EstimateParameters(stasis, paleoData)


d = GenerateDataBest(6, stasis, parameters, urw, models, 100, 'AIC')


testmm = c( -0.338391311, -0.179746973, -0.039437045,  0.290798380,  0.001312343, -0.513058254)
testvv = rep(1, 6)
testnn = rep(20, 6)
testtt = 0:5
testMM = rep(-0.09465211, 6)
testdata = as.paleoTS(testmm, testvv, testnn, testtt, testMM)
seed= 1944
set.seed(seed)
print(testdata)
fitSimple(testdata, "GRW")
sapply(1:100, function(x) EstimateParameters(urw, testdata))
