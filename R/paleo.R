
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
