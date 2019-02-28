fitSimple(y, model = )
data = sim.GRW(ms = 1)
paleoGRW
paleoURW
myfit.grw = fitSimple(data, "GRW")
myfit.urw = fitSimple(data, "URW")
myfit.stasis = fitSimple(data, "Stasis")
parameters = list(ns = 20, )
GenerateData.paleoTS = function(n, model, parameters){
  model.type = strsplit(model$label, ".", fixed = TRUE)[[1]][2]
  class(model) = model.type
  NextMethod("GenerateData", model)
}

GenerateData.GRW = function(n, model, parameters){
 ns <- n
 ms <-parameters$ns
 vs <- parameters$vs
 vp <- parameters$vp
 nn <- parameters$nn
 tt <- parameters$tt
 sim.GRW(ns, ms, vs, vp, nn,  tt)
 }
GenerateData.URW = function(n, model, parameters){
  ns <- n
  ms <- 0
  vs <- parameters$vs
  vp <- parameters$vp
  nn <- parameters$nn
  tt <- parameters$tt
  sim.GRW(ns, ms, vs, vp, nn,  tt)
}
GenerateData.Stasis = function(n, model, parameters){
  ns <- n
  theta <- parameters$theta
  omega <- parameters$omega
  vp <- parameters$vp
  nn <- parameters$nn
  tt <- parameters$tt
  sim.Stasis(ns, theta, omega, vp, nn,  tt)
}
x<- sim.Stasis(ns=30, theta=10, omega=1)
s1<- fitSimple(x, model="URW")
s2<- fitSimple(x, model="Stasis")
s3<- fitSimple(x, model="StrictStasis")
