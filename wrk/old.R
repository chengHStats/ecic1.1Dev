####VALUES USED FOR TESTING PROCESS
library(magrittr)
library(fitdistrplus)
library(parallel)


smoothdata = gen(70, "smooth", param.smooth)
est(smoothdata, "smooth")

est  = function(data, dist){
  if (dist == "rayleigh")
    return(fitdist(data, "rayleigh", start = 2)$estimate)
  if (dist == "weibull")
    return(fitdist(data, "weibull")$estimate)
  if (dist == "exp")
    return(fitdist(data, "exp")$estimate)
  if (dist == "norm2")
    return(c("mu" = mean(data), "sd" = sd(data)))
  if (dist == "norm1")
    return("mu" = mean(data))
  if (dist == "uwalk"){
    steps = c(data[1], diff(data))
    return(c("sd" = sd(steps)))
  }
  if (dist == "gwalk"){
    steps = c(data[1], diff(data))
    return(c("mu" = mean(steps), "sd" = sd(steps)))
  }
  if(dist == "seg"){
    # segs = apply(breaklist, 1, function(x) fit.seg(data, x))
    # ssrs = sapply(segs, function(x) x$ssr)
    # best = which.min(ssrs)
    # ssr = min(ssrs)
    # resid = unlist(segs[[best]]$resid)
    # sd = sqrt(min(ssrs)/(length(data)))
    fit = fit.seg(data, c(10,20))
    ssr = fit$ssr
    sd = sqrt(ssr/length(data))
    #breaks = setNames(breaklist[best,], c("T1", "T2"))
    coefs = fit$coefs #segs[[best]]$coefs
    return(unlist(c(coefs, c(10, 20), sd)))
  }
  if(dist == "seg2"){
    fit = segmented(lm(data~ ix), seg.Z = ~ix, psi = c(10, 20), control = seg.control(it.max = 0, maxit.glm = 0, h = 0, n.boot = 0))
    return(c(coefficients(fit), sqrt(sum(fit$residuals^2)/25)))
  }
  if (dist == "smooth"){
    sm = summary(nlsLM(as.vector(data) ~ a + b * exp(c*(1:length(data))), control = list(warnOnly = TRUE, maxiter = 20) , start = list(a = 0.2,b = 1, c =.01)))
    return(setNames(c(sm$coefficients[-2,1], sm$coefficients[2,1], "sd"  = sqrt((sm$sigma/length(data)))), c("a", "b", "c", "sd")))
  }
  if (dist == "smooth2"){
    x = exp(0.2 * 1:25)
    xmatrix = matrix(c(rep(1,25),x), ncol = 2)
    bhat = solve((t(xmatrix) %*% xmatrix)) %*% (t(xmatrix) %*% data)
    yhat = bhat[1] + bhat[2]*x
    resid = (data-yhat)
    sd = sqrt(sum(resid^2)/25)
    co = c(bhat, sd)
    sc = sum(log(dnorm(resid, mean(resid), sd(resid))))*-2 + 6
    return(co)
  }
  if(dist == "lm1"){
    bhat = solve(t(spanos.X1)%*%spanos.X1)%*%t(spanos.X1)%*%matrix(data, ncol = 1) %>%unname
    resid = data-spanos.X1%*%bhat
    sd = sqrt(t(resid)%*%(resid)/33)
    return(setNames(c(bhat, sd), c("a0", "a1", "sd")))}
  if(dist %in% c("GRW", "URW", "Stasis") )return(list("Parameters" = fitSimple(data, dist)$parameters, "vv" = data$vv, "nn" = data$nn))
  if(dist == "lm2"){
    bhat = solve(t(spanos.X2)%*%spanos.X2)%*%t(spanos.X2)%*%matrix(data, ncol = 1) %>%unname
    resid = data-spanos.X2%*%bhat
    sd = sqrt(t(resid)%*%(resid)/32)
    return(setNames(c(bhat, sd), c("b0", "b1", "b2", "b3", "sd")))}
  if(dist == "lm3"){
    bhat = solve(t(spanos.X3)%*%spanos.X3)%*%t(spanos.X3)%*%matrix(data, ncol = 1) %>%unname
    resid = data-spanos.X3%*%bhat
    sd = sqrt(t(resid)%*%(resid)/32)
    return(setNames(c(bhat, sd), c("b0", "b1", "b2", "sd")))}
}

gen= function(n, dist, param){
  if (dist == "rayleigh")
    return(rweibull(n, 2, param*sqrt(2)))
  if (dist == "weibull")
    return(rweibull(n, param[1], param[2]))
  if (dist == "exp")
    return(rexp(n, param))
  if (dist == "norm2")
    return(rnorm(n, param[1], param[2]))
  if (dist == "norm0")
    return(rnorm(n, 0, 1))
  if (dist == "norm1")
    return(rnorm(n, param[1], 1))
  if(dist == "uwalk")
    return(cumsum(rnorm(n, 0, param[1])))
  if(dist == "gwalk")
    return(cumsum(rnorm(n, param[1], param[2])))
  if (dist == "seg"){
    int = param[1]
    s1  = param[2]
    s2  = param[3]
    s3  = param[4]
    T1  = round(param[5])
    T2  = round(param[6])
    sd = param[7]

    y1 = int + (1:(T1) * s1)
    y2 = y1[T1] + (1:(T2-T1))*s2
    y3= y2[T2-T1] + (1:(n-T2)) * s3


    return(c(y1,y2,y3) + rnorm(n, 0, sd))
  }
  if (dist == "smooth"){
    a = param[1]
    b = param[2]
    c = param[3]
    sd = param[4]
    return(a + b*exp((1:n)*c) + rnorm(n, 0, sd))
  }
  if (dist == "GRW"){
    anc = param$Parameters[[1]]
    ms  = param$Parameters[[2]]
    vs  = param$Parameters[[3]]
    vv  = param$vv[1]
    nn  = param$nn
    return(sim.GRW(n, ms, vs, vv, nn))
  }
  if (dist == "URW") {
    anc = param$Parameters[[1]]
    vs  = param$Parameters[[2]]
    vv  = param$vv[1]
    nn  = param$nn
    return(sim.GRW(n, ms, vs, vv, nn))
  }
  if (dist == "Stasis"){
    theta  = param$Parameters[[1]]
    omega  = param$Parameters[[2]]
    vv  = param$vv[1]
    nn  = param$nn
    return(sim.Stasis(n, theta, omega, vv, nn))
  }
  if(dist == "lm1"){
    a0 = param[1]
    a1 = param[2]
    sd = param[3]
    x = spanos.x
    return(a0+(a1*x) + rnorm(n, 0, sd))}
  if(dist == "lm2"){
    b0 = param[1]
    b1 = param[2]
    b2 = param[3]
    b3 = param[4]
    sd = param[5]
    x = spanos.x
    return(b0 + (b1*x) + (b2 * x^2) + (b3 * x^3) + rnorm(n, 0, sd))}
  if(dist == "lm3"){
    b0 = param[1]
    b1 = param[2]
    b2 = param[3]
    sd = param[4]
    x = spanos.x
    return(b0 + (b1*x) + (b2 * x^2) + rnorm(n, 0, sd))}
  if (dist == "smooth2"){
    return(param[1] + param[2] * exp(.2 * 1:n) + rnorm(n, 0, param[3]))
  }
  if (dist == "seg2"){
    it = param[1]
    s1 = param[2]
    s2 = param[3]
    s3 = param[4]
    sd = param[5]

    y1 = it + s1*(1:10)
    y2 = y1[10] + s2*(1:10)
    y3 = y2[10] + s3*(1:5)

    return(c(y1,y2,y3) + rnorm(25, 0, sd))}
}
gens = function(n, dist, param, N){
  if (dist == "norm2")
    return(rnorm(n*N, param[1], param[2]) %>% matrix(ncol = N))
  if (dist == "norm0")
    return(rnorm(n*N, 0, 1) %>% matrix(ncol = N))
  if (dist == "norm1")
    return(rnorm(n*N, param[1], 1) %>% matrix(ncol = N))
  if(dist == "uwalk"){
    steps = rnorm(n*N, 0, param[1]) %>% matrix(ncol = N)
    return(apply(steps, 2, cumsum))
  }
  if(dist == "gwalk"){
    steps = rnorm(n*N, param[1], param[2]) %>% matrix(ncol = N)
    return(apply(steps, 2, cumsum))}
  if (dist == "GRW"){
    anc = param$Parameters[[1]]
    ms  = param$Parameters[[2]]
    vs  = param$Parameters[[3]]
    vv  = param$vv[1]
    nn  = param$nn
    return(lapply(1:N, function(x) sim.GRW(n, ms, vs, vv, nn)))
  }
  if (dist == "URW") {
    anc = param$Parameters[[1]]
    vs  = param$Parameters[[2]]
    vv  = param$vv[1]
    nn  = param$nn
    return(lapply(1:N, function(x) sim.GRW(n, ms, vs, vv, nn)))
  }
  if (dist == "Stasis"){
    theta  = param$Parameters[[1]]
    omega  = param$Parameters[[2]]
    vv  = param$vv[1]
    nn  = param$nn
    return(lapply(1:N, function(x) sim.Stasis(n, theta, omega, vv, nn)))
  }
  if (dist == "seg") return(sapply(1:N, function(x) gen(n, "seg", param)))
  if (dist == "seg2") return(sapply(1:N, function(x) gen(n, "seg2", param)))
  if (dist == "smooth") return(sapply(1:N, function(x) gen(n, "smooth", param)))
  if (dist == "smooth2") return(sapply(1:N, function(x) gen(n, "smooth2", param)))
  if (dist == "lm1") return(sapply(1:N, function(x) gen(n, "lm1", param)))
  if (dist == "lm2") return(sapply(1:N, function(x) gen(n, "lm2", param)))
  if (dist == "lm3") return(sapply(1:N, function(x) gen(n, "lm3", param)))

}
ests = function(data, dist){
  if (dist == "rayleigh")
    return(fitdist(data, "rayleigh", start = 2)$estimate)
  if (dist == "weibull")
    return(fitdist(data, "weibull")$estimate)
  if (dist == "exp")
    return(fitdist(data, "exp")$estimate)
  if (dist == "norm2")
    return( t(c(colMeans(data), apply(data, 2, sd))) %>% matrix(nrow = 2, byrow = TRUE, dimnames = list(c("mu", "sd"))))
  if (dist == "norm1")
    return("mu" = colMeans(data) %>% matrix(nrow = 1))
  if (dist == "uwalk")
    return(apply(data, 2, function(x) est(x, "uwalk")) %>% t)
  if (dist == "gwalk")
    return(apply(data, 2, function(x) est(x, "gwalk")))
  if(dist %in% c("GRW", "URW", "Stasis") )return(list("Parameters" = lapply(data,  function(x) fitSimple(x, dist)$parameters), "vv" = data[[1]]$vv, "nn" = data[[1]]$nn ))
  if (dist == "seg")
    return(apply(data, 2, function(x) est(x, "seg")))
  if (dist == "smooth")
    return(apply(data, 2, function(x) est(x, "smooth")))
  if (dist == "smooth2")
    return(apply(data, 2, function(x) est(x, "smooth2")))
  if (dist == "seg2")
    return(apply(data, 2, function(x) est(x, "seg2")))
  if (dist == "lm1")
    return(apply(data, 2, function(x) est(x, "lm1")))
  if (dist == "lm2")
    return(apply(data, 2, function(x) est(x, "lm2")))
  if (dist == "lm3")
    return(apply(data, 2, function(x) est(x, "lm3")))
}

distscore = function(data, dist, fun = "AIC") {
  score = 0
  n = length(data)
  if(fun == "AIC"){
    if (dist == "rayleigh") score = fitdist(data, "rayleigh", start = 2)$aic
    if (dist == "norm1"){ score = sum(log(dnorm(data, mean(data), 1)))*-2 + 2}
    if (dist == "norm0") {score = sum(log(dnorm(data, 0, 1)))*-2}
    if (dist == "norm2"){ score = sum(log(dnorm(data, mean(data), sd(data))))*-2 + 4}
    if(dist == "lm1"){
      bhat = solve( t(spanos.X1)%*%spanos.X1) %*%t(spanos.X1)%*%matrix(data, ncol = 1) %>% unname
      resid = data-spanos.X1%*%bhat
      sd = sqrt(t(resid)%*%(resid)/32)
      score = sum(log(dnorm(resid, 0, sd)))*-2 + 6
    }
    if(dist == "lm2"){
      bhat = solve(t(spanos.X2)%*%spanos.X2)%*%t(spanos.X2)%*%matrix(data, ncol = 1) %>%unname
      resid = data-spanos.X2%*%bhat
      sd = sqrt(t(resid)%*%(resid)/32)
      score = sum(log(dnorm(resid, 0, sd)))*-2 + 10
    }
    if(dist == "lm3"){
      bhat = solve(t(spanos.X3)%*%spanos.X3)%*%t(spanos.X3)%*%matrix(data, ncol = 1) %>%unname
      resid = data-spanos.X3%*%bhat
      sd = sqrt(t(resid)%*%(resid)/32)
      score = sum(log(dnorm(resid, 0, sd)))*-2 + 8
    }
    if (dist == "uwalk"){
      steps = c(data[1], diff(data))
      score = sum(log(dnorm(steps, 0, sd(steps))))*-2 + 2
    }
    if (dist == "gwalk"){
      steps = c(data[1], diff(data))
      param = c(mean(steps), sd(steps))
      score = sum(log(dnorm(steps, param[1], param[2]))) * -2 + 4
    }
    if (dist == "seg"){
      # segs = apply(breaklist, 1, function(x) fit.seg(data, x))
      # ssrs = sapply(segs, function(x) x$ssr)
      # best = which.min(ssrs)
      # ssr = min(ssrs)
      fit = fit.seg(data, c(10, 20))
      resid = unlist(fit$resid)

      logl = sum(log(dnorm(resid, 0, sqrt(sum(resid^2)/length(data)))))
      return(-2*logl + 10)
    }
    if (dist == "seg2"){
      return(AIC(segmented(lm(data~ ix), seg.Z = ~ix, psi = c(10, 20), control = seg.control(it.max = 0, maxit.glm = 0, h = 0, n.boot = 0))))

    }
    if (dist == "smooth"){
      return(as.numeric((-2 * logLik(nlsLM(as.vector(data) ~ a + exp(c*(1:length(data))), algorithm  = "plinear", control = list(warnOnly = TRUE, maxiter = 10) , start = list(a = 0, c =.2)))) + 8))
    }
    if (dist == "smooth2"){
      x = exp(0.2 * 1:25)
      xmatrix = matrix(c(rep(1,25),x), ncol = 2)
      bhat = solve((t(xmatrix) %*% xmatrix)) %*% (t(xmatrix) %*% data)
      yhat = bhat[1] + bhat[2]*x
      resid = (data-yhat)
      sd = sqrt(sum(resid^2)/25)
      score = sum(log(dnorm(resid, mean(resid), sd(resid))))*-2 + 6
      return(score)}
    return(score)
  }
  if (fun == "AICc"){
    if (dist == "rayleigh") score = fitdist(data, "rayleigh", start = 2)$aic
    if (dist == "norm1") score = sum(log(dnorm(data, mean(data), 1)))*-2 + (2/n-1)
    if (dist == "norm0") score = sum(log(dnorm(data, 0, 1)))*-2
    if (dist == "norm2"){
      score = sum(log(dnorm(data, mean(data), sd(data))))*-2 + 4 + (12/(n-3)) }
    if (dist == "uwalk"){
      steps = c(data[1], diff(data))
      score = sum(log(dnorm(steps, 0, sd(steps))))*-2 + 2 + (4/(n-2))
    }
    if (dist == "gwalk"){
      steps = c(data[1], diff(data))
      param = c(mean(steps), sd(steps))
      score = sum(log(dnorm(steps, param[1], param[2]))) * -2 + 4 + (12/(n-3))
    }}
  if(dist %in% c("GRW", "URW", "Stasis") ) score = fitSimple(data, dist)$AICc

  return(score)
}

scores = function(data, dists, fun = 'AIC'){
  sapply(dists, function(x) distscore(data, x, fun))
}
###distscores takes a matrix of data and for certain distributions processes the ho
distscores2 = function(data, dist, fun = 'AIC') {
  score = 0
  if(is.matrix(data)) n = length(data[,1])
  if (fun == 'AIC'){
    if (dist == "rayleigh") score = fitdist(data, "rayleigh", start = 2)$aic
    if (dist == "norm1") scores = colMeans(data) %>%(function(x) log(dnorm(data, rep(x,each =  n), 1))) %>% matrix(nrow = n) %>% colSums() *-2 +2
    if (dist == "norm0") scores = log(dnorm(data,0,1)) %>% matrix(nrow = n) %>% apply(2, sum)*-2
    if (dist == "norm2") scores = matrix(c(colMeans(data), apply(data, 2, sd)), ncol = 2) %>% (function(x) log(dnorm(data, rep(x[,1], each = n), rep(x[,2], each = n)))) %>% matrix(nrow = n) %>% colSums()*-2 + 4
    if (dist == "uwalk") scores = apply(data, 2, function(x) distscore(x, "uwalk"))
    if (dist == "gwalk") scores = apply(data, 2, function(x) distscore(x, "gwalk"))
    if (dist == "seg2") scores = apply(data, 2, function(x) distscore(x, "seg2"))
    if (dist == "smooth") scores = apply(data, 2, function(x) distscore(x, "smooth"))
    if (dist == "smooth2") scores = apply(data, 2, function(x) distscore(x, "smooth2"))

    if (dist == "lm1") scores = apply(data, 2, function(x) distscore(x, "lm1"))
    if (dist == "lm2") scores = apply(data, 2, function(x) distscore(x, "lm2"))
    if (dist == "lm3") scores = apply(data, 2, function(x) distscore(x, "lm3"))

    if(dist %in% c("GRW", "URW", "Stasis") ) scores = sapply(data, function(x) fitSimple(x, dist)$AICc)
  }
  if (fun == 'AICc'){
    if (dist == "norm2") scores = matrix(c(colMeans(data), apply(data, 2, sd)), ncol = 2) %>% (function(x) log(dnorm(data, rep(x[,1], each = n), rep(x[,2], each = n)))) %>% matrix(nrow = n) %>% colSums()*-2 + 4 + (12/(n-3))
    if (dist == "uwalk") scores = apply(data, 2, function(x) distscore(x, "uwalk", 'AIC')) +  (4/(n-2))
    if (dist == "gwalk") scores = apply(data, 2, function(x) distscore(x, "gwalk", 'AIC'))+ (12/(n-3))

  }
  return(scores)
}

scores2 = function(data, dists, fun = 'AIC'){
  as.matrix(t(sapply(dists, function(x) distscores2(data, x, fun))))
}
###FUNCTIONS FOR PERMUTATION TESTING OVER DISTRIBUTIONS
gen.best = function(n, true, param, best,dists, N, ic = 'AIC', ...){
  check = gen(n, true, param)
  if(is.numeric(check)){
    mins = 0
    best.ix = which(dists==best)
    it1 = 0
    nonzero = FALSE
    while(!(best.ix %in% mins) &  it1 < 1000){
      data = gens(n, true, param, N*3)
      mins = apply(scores2(data, dists, ic), 2, which.min)
      data.best2 = as.matrix(data[,mins==best.ix])
      if(best.ix %in% mins) nonzero = TRUE
      it1 = it1 + N
    }
    if(!nonzero)return(0)
    it = 0
    while(dim(data.best2)[2] < N & it < 4){
      data = gens(n, true, param, N*3)
      mins = apply(scores2(data, dists, ic), 2, which.min)
      data.best = data[,mins==best.ix]
      data.best2 = cbind(data.best2, data.best)
      it = it + 1
    }
    if(dim(data.best2)[2] > N) out = as.matrix(data.best2[,1:N])
    if(dim(data.best2)[2] <= N) out = as.matrix(data.best2)
    return(out)}


  if(is.list(check)){
    mins = 0
    best.ix = which(dists==best)
    while(!(best.ix %in% mins)){
      data = gens(n, true, param, N*3)
      mins = apply(scores2(data, dists),2, which.min)
      data.best2 = lapply((1:(N*3))[mins == best.ix], function(x) data[[x]])
    }
    while(length(data.best2) < N){
      data = gens(n, true, param, N*3)
      mins = apply(scores2(data, dists), 2, which.min)
      data.best = lapply((1:(N*3))[mins == best.ix], function(x) data[[x]])
      data.best2 = c(data.best2, data.best)
    }
    data.out = lapply(1:N, function(x) data.best2[[x]])
    if(N == 1)data.out = data.out[[1]]
    return(data.out)
  }
}

gen.best.perm = function(n, dists, param, N, ic = 'AIC'){f
  pv = perm.vec(length(dists))
  lists = apply(pv, 2, function(x) list(n = n, true = dists[x[1]], param = param,
                                        best = dists[x[2]], N = N))
  names(lists) = apply(pv, 2, function(x) paste(dists[x[1]], dists[x[2]] ))
  lapply(lists, function(x) do.call("gen.best", x))
}

###Gets the vector of indices for every permutation of the distributions
perm.vec = function(length){
  p = length
  toprow = rep((1:p), rep(p-1, p))
  botrow = sapply((1:p), function(x) (1:p)[-x])
  matrix(c(toprow, botrow), byrow = TRUE, nrow = 2)
}
###Tests all permutations of the dist list over a function
perm.test = function(fn, dists, ...){
  args = list(...)
  full.args = c(list(fn = fn, dists = dists), args)
  x = length(dists) %>% matrix(0, nrow = ., ncol = .)
  vec = perm.vec(length(dists))
  perm.dist = apply(vec, 2, function(x) dists[x]) %>% split(., rep(1:ncol(.), each = nrow(.)))
  out = lapply(perm.dist, function(x) do.call(fn, c(list("best" = x[1], "true" = x[2]),args)) )
  names(out) = lapply(perm.dist, function(x) paste("best = ", x[1], ", true = ",x[2], sep = ""))
  return(list(out = out, info = full.args))
}
###sets  up data to fill a matrix off the diagonal
matfill = function(data, p){
  mat = matrix(0, ncol = p, nrow = p)
  mat[row(mat) != col(mat)] = data
  mat
}
###takes the results form perm.test and puts in an array
perm.to.mat = function(results, mean = TRUE){
  dists = results$info$dists
  p = length(dists)
  outnames = results$out[[1]] %>% call(ifelse(is.vector(.), "names", "colnames"), .) %>% eval
  outs = length(outnames)
  means = results$out %>% lapply(function(x) call(ifelse(is.vector(x), "mean", "colMeans"), x) %>% eval) %>%
    do.call(rbind, .) %>% data.frame

  array(sapply(means, function(x) matfill(x, p)), c(p, p, outs),
        dimnames = list(paste(dists, "best"), paste(dists, "true"), outnames))
}

scores = function(data, dists) { sapply(dists, function(x) distscore(data,x)) }

####FIND PARAMETERS WITH EVEN AIC DISTRIBUTION

##### MLE'S / CONDITIONING+CENSORING / BIAS CORRECTION

ic.MLE.true = function(true, param, dists,  n=50, N = 1000, ic = 'AIC'){
  data = gens(n, true, param, N)
  if(is.matrix(data)){
    mins = scores2(data, dists, ic) %>% apply(2, which.min)
    mles = t(ests(data, true))
    newMLES = setNames(split(1:N, mins) %>% lapply(function(x) mles[x,]), names(dists)) %>%  (function(x) sapply(x, function(y) flexMean(y, 2)))
  }
  if(is.list(data)){
    mins = scores2(data, dists, ic) %>% apply(2, which.min)
    mles = matrix(unlist(ests(data, true)$Parameters), nrow= N, byrow = TRUE)
    newMLES =  setNames(split(1:N, mins) %>% lapply(function(x) mles[x,]), names(dists)) %>%  (function(x) sapply(x, function(y) flexMean(y, 2)))
  }
  return(newMLES)
}


###### NOT USING ANYMORE?


#performs our bias correction on a dataset - "full" gives mean and sd instead of just the parameters for the "true" dist
biascorrect = function(data, true, best, dists, N = 100, full = FALSE, ic = 'AIC'){
  # if(is.numeric(data)){
  n = length(data)
  if(!is.null(est(data, true))){
    fit = est(data, true)

    best.ind = which(dists==best)
    true.ind = which(dists==true)
    fit2 = c(mean(data), sd(data))

    newdata.best = gen.best(n, true, fit, best, dists, N, ic)
    if(!is.matrix(newdata.best))if(newdata.best == 0) return(0)
    N = dim(newdata.best)[2]

    if(full == FALSE){

      meanfit.best = matrix(apply(newdata.best, 2, function(x) est(x, true)), ncol = N) %>% flexMean(1)
      fit = fit*2 - meanfit.best
    }
    if(full == TRUE){
      meanfit.best = apply(newdata.best, 2, function(x) c(mean(x), sd(x))) %>% rowMeans
      fit = fit2*2 - meanfit.best
    }

    return(fit)
  }
}


#
# if(class(data) == "paleoTS"){
#   n = length(data$mm)
#   fit = est(data, true)
#   fit.param = fit$Parameters
#   fit.vv = fit$vv
#   fit.nn = fit$nn
#   newdata.best = gen.best(n, true, fit, best, dists, N)
#   meanfit.best = lapply(newdata.best, function(x) est(x, true)$Parameters) %>% unlist %>% matrix(nrow = N, byrow = TRUE) %>% colMeans
#   fit.param2 = fit.param*2 - meanfit.best
#   out = list("Parameters" = fit.param2, "vv" = fit.vv, "nn" = fit.nn)
#   return(out)
# }




###Gets MLEs for N datasets, N2 is length of parametric bootstrap for bias correction


# ####Run this to get MLES under all 6 permutations of distributions
# uncorrected.mles = perm.test("ic.MLE", dists, param = c(0.3, 1.2), n = 50, N= 1000, N2= 1000, FALSE) %>% perm.to.mat
# corrected.mles = perm.test("ic.MLE", dists, param = c(0.3, 1.2), n = 50, N= 1000,N2 = 1000, TRUE) %>% perm.to.mat


##### CONDITIONAL PROBABILITIES

###TRUE PROBABILITIES USING DIRECT SIMULATION
ic.probs.true = function(true, param, dists, n = 50, N = 10000, parallel = FALSE, ic = 'AIC'){
  require(parallel)
  if(parallel) clust = makeCluster(detectCores()-1)
  data = gens(n, true, param, N)
  if(is.matrix(data)){
    if(parallel){
      clusterExport(clust, c("data", "distscore", "%>%", "dists", "fit.seg"))
      mins = parApply(clust, data,2, function(x) sapply(dists, function(y) distscore(x,y, ic)) %>% which.min)
    }
    if(!parallel) mins = apply(data,2, function(x) sapply(dists, function(y) distscore(x,y, ic)) %>% which.min)
    out = rep(0, length(dists))
    names(out) = dists
    probs = table(dists[mins])/N
    out[names(probs)] = probs
    if(parallel) stopCluster(clust)}
  if(is.list(data)){
    if(class(data[[1]]) == "paleoTS"){
      mins = scores2(data, dists) %>% apply(2, which.min)
      out = rep(0, length(dists))
      names(out) = dists
      probs = table(dists[mins])/N
      out[names(probs)] = probs
    }
  }
  out
}

probs.test = function(dists, param, n = 50, N = 1000, ic = 'AIC'){
  out = as.matrix(sapply(dists, function(x) ic.probs.true(x, param, n, N, ic)))
  rownames(out) %<>% paste("best = ", .)
  out
}

probs.to.param = function(true, param, dists, prob = 1/length(dists), stepsize = 0.01, eps = 0.05, max.it = 100, n = 50, N = 10000){
  if (length(prob) == 1) prob %<>% rep(length(dists))
  k = length(param)
  param = matrix(param, nrow = k)
  prob.use = ic.probs.true(true, param,dists, n, N) %T>% print
  delta = sum(abs(prob.use - prob))
  i = 0
  mu = 0
  while(delta > eps & i < max.it){
    newparams = (rnorm(5*k, mu, stepsize)+c(param)) %>% matrix(nrow = k)
    prob.2 =  apply(newparams, 2 ,function(x) ic.probs.true(true, x, dists, n, N))
    delta.2 = apply(prob.2, 2, function(x) sum(abs(x-prob)))
    mu = c(mu, apply(newparams, 2, function(x) mean(param - x)))[which.min(c(delta, delta.2))] %T>% print
    param  = cbind(param, newparams)[,which.min(c(delta, delta.2))]
    prob.use = cbind(prob.use, prob.2)[,which.min(c(delta, delta.2))] %T>% print
    delta = sum(abs(prob.use-prob))
    i = (i + 1)  %T>% print
  }
  param
}

probs.to.param.2 = function(true, param, dists, prob = 1/length(dists), eps = c(0.05, 0.01), stepsize = c(0.01, 0.005), max.it = c(50,100), n = 50, N = c(1000,10000) ){
  param = probs.to.param(true, param, dists, prob, stepsize[1], eps[1], max.it[1], n, N[1] )
  param = probs.to.param(true, param, dists, prob, stepsize[2], eps[2], max.it[2], n, N[2])
}


###PARAMETRIC BOOTSTRAP


est.probs = function(data, dists, true, N = 1000, N2 = 1000, correct = TRUE, parallel = FALSE, ic = 'AIC'){
  if (parallel) {require(parallel)}
  if(!class(data) == "paleoTS") n = length(data)
  if(class(data) == "paleoTS") n = length(data$mm)
  best = dists[data %>% scores(dists, ic) %>% which.min]
  if(correct == FALSE){
    newdata = est(data, true) %>% gens(n, true, ., N)
  }
  if(correct == TRUE){
    fit = biascorrect(data, true, best, dists, N2, full = FALSE, ic)
    newdata = gens(n, true, fit, N)
  }

  mins = scores2(newdata, dists, ic) %>% apply(2, which.min)

  (table(factor(dists)[mins])/N)[dists]

}

test.probs = function(best, true, dists, param, n = 50, N= 100, N2 = 100, correct = TRUE){
  data = gen.best(n, true, param, best, dists, N)
  if(!is.matrix(data))if(data==0)return(0)
  t(apply(data, 2, function(x) est.probs(x, dists, true, N, N2, correct)))
}

test.probs.all = function(dists, param, n, N =100, N2 = 100, correct = TRUE, ic = 'AIC'){
  args = list(dists = dists,
              param = param,
              n = n,
              N = N,
              N2 = N2,
              correct = correct,
              ic = ic)
  vec = perm.vec(length(dists))
  perm.dist = apply(vec, 2, function(x) dists[x]) %>% split(., rep(1:ncol(.), each = nrow(.)))
  results = lapply(perm.dist, function(x) do.call("test.probs", c(list(best = x[1], true = x[2]) ,args)))
  names(results) = lapply(perm.dist, function(x) paste(x[1], "true,", x[2], "best"))
  results
}

#  test.probs.all(dists, param, n, N, N2, FALSE)
#  test.probs.all(dists, param, n, N, N2, TRUE)


###COMPARE NON-BIAS CORRECTED TO BIAS CORRECTED




##### ###GETTING PERCENTILES

###QUASI TRUE
test.thresh = function(thresh, best, true, dists, param, n, N){
  best.ix = which(names(dists)==best)
  true.ix = which(names(dists)==true)
  data = gen.best(n, true, param, best, dists, N)

  scores.data = apply(data, 2, function(x) scores(x, dists))
  difs = apply(scores.data, 2, function(x) x[best.ix] - min(x[-best.ix])) %>%  sort

  return( length(difs[difs < thresh])/N)
}

threshdirect = function(best, true, dists, param, n, N, alpha, ic = 'AIC'){
  best.ix = which(names(dists)==best)
  true.ix = which(names(dists)==true)
  data = gen.best(n, true, param, best,dists, N, ic)
  if(!is.matrix(data))if(data==0) return(0)
  scores.data = scores2(data, dists, ic)
  difs = apply(scores.data, 2, function(x) x[best.ix] - min(x[-best.ix])) %>%  sort
  thresh = difs[ceiling(alpha*N)]
  ifelse(is.na(thresh), 0, thresh)
  thresh
}

threshpar = function(data, best, true, dists, N, alpha, correct = TRUE, ic = 'AIC'){
  best.ix = which(names(dists)==best)
  true.ix = which(names(dists)==true)
  if(!class(data) == "paleoTS") n = length(data)
  if(class(data) == "paleoTS") n = length(data$mm)

  if(correct == TRUE){
    fit =  biascorrect(data, true, best,dists, N, full = FALSE, ic = ic)
    data = gen.best(n, true, fit, best, dists, N, ic)
    if(!is.matrix(data))if(data==0) return(0)
  }
  if(correct == FALSE){
    fit = est(data, true)
    data = gen.best(n, true, fit, best, dists, N,ic =  ic)
    if(!is.matrix(data))if(data==0) return(0)
  }
  scores.data = scores2(data, dists, ic)
  difs = apply(scores.data, 2, function(x) x[best.ix] - min(x[-best.ix])) %>%  sort
  thresh = difs[ceiling(alpha*N)]
  ifelse(is.na(thresh), 0, thresh)
  thresh
}


###PLOT DIFFERENCE DISTRIBUTIONS, SHOW THRESHOLDS

###SEE HOW THRESHOLDS CHANGE FOR VARYING SAMPLE SIZES (AND PARAMETERS?)

#ESTIMATION

###COMPARE PARAMETRIC TO TWO-STAGE BOOTSTRAP

###LOOK AT THRESHOLD VALUES OBTAINED AS WELL AS THE PERCENTILE THEY ARE

###ACTUAL PERCENTILE VALUES AS A FUNCTION OF BOOTSTRAP SIZE?

#####ESTIMATING THRESHOLDS USING ESTIMATIONS OF MODEL PROBABILITIES

###HERE WE DON'T KNOW THE A1 VALUE AND WE ARE ESTIMATING IT

threshpar.est = function(data, best, true, dists, N,N2, alpha, correct = TRUE, ic = 'AIC'){
  best.ix = which(names(dists)==best)
  true.ix = which(names(dists)==true)
  if(!class(data) == "paleoTS") n = length(data)
  if(class(data) == "paleoTS") n = length(data$mm)
  aprime = (alpha/est.probs(data, dists, true, N, N2, FALSE, ic = ic))[best.ix]
  if(correct == TRUE){
    fit =  biascorrect(data, true, best, dists, N, full = FALSE, ic)
  }
  if(correct == FALSE){
    fit = est(data, true)
  }
  data = gen.best(n, true, fit, best,dists,  N, ic)
  if(!is.matrix(data)) if(data == 0 )return(0)
  scores.data = scores2(data, dists, ic)
  difs = apply(scores.data, 2, function(x) x[best.ix] - min(x[-best.ix])) %>%  sort
  thresh = difs[ceiling(alpha*N)]
  thresh = ifelse(is.na(thresh), 0, thresh)
  thresh
}

threshpar.aprime = function(data, best, true, dists, N, N2, alpha, correct = TRUE, ic = 'AIC'){
  best.ix = which(names(dists)==best)
  true.ix = which(names(dists)==true)
  if(!class(data) == "paleoTS") n = length(data)
  if(class(data) == "paleoTS") n = length(data$mm)
  probs = est.probs(data, dists, true, N, N2, correct, ic = ic)
  aprime = (alpha/probs)[best.ix]
  thresh = threshpar.est(data, best, true, dists, N, N2, aprime, correct, ic)
  list("Threshold" = thresh, "aPrime" =aprime, "Probabilities" = probs)
}

threshpar.aprime2 = function(data, best, true, dists, N, N2, alpha, correct = TRUE, ic = 'AIC'){
  ###use the same fits i.e. only call biascorrect once
  best.ix = which(names(dists)==best)
  true.ix = which(names(dists)==true)
  if(!class(data) == "paleoTS") n = length(data)
  if(class(data) == "paleoTS") n = length(data$mm)
  probs = est.probs(data, dists, true, N, N2, correct, ic)
  aprime = (alpha/probs)[best.ix]
  thresh = threshpar.est(data, best, true, dists, N, N2, aprime, correct, ic)
  list("Threshold" = thresh, "aPrime" =aprime, "Probabilities" = probs)
}
2


####FULL PROCEDURE / TWO THRESHOLDS

threshmins = function(data, best, dists, N, N2, alpha, correct = TRUE, ic = 'AIC'){
  if(!class(data) == "paleoTS") n = length(data)
  if(class(data) == "paleoTS") n = length(data$mm)
  best.ix = which(names(dists)==best)
  best = dists[best.ix]
  threshdata = lapply(dists[-best.ix], function(x) threshpar.aprime(data, best,x, dists, N, N2, alpha, correct, ic))
  probs = sapply(threshdata, function(x) x[["Probabilities"]])
  threshes = sapply(threshdata, function(x) x[["Threshold"]])
  aprime = sapply(threshdata, function(x) x[["aPrime"]])
  list("Thresholds" = threshes, "aPrime" = aprime, "Probabilities" = probs)
}

thresh.perm = function(dists, param, n, N, alpha){
  true.thresh = unname(sapply(dists, function(x) sapply(dists[dists!=x],
                                                        function(y) threshdirect(x, y, dists, param, n, N, alpha), USE.NAMES = FALSE ))) %>% matfill(length(dists))
  dimnames(true.thresh) = list(do.call("paste" , list(dists, "best")), dists)
  true.thresh
}

thresh.est.perm = function(dists, param, n, N, alpha, correct, np){
  data = lapply(dists, function(x) lapply(dists[dists!=x], function(y) gen.best(n, x, param, y, N) %>% as.matrix(nrow = n)) )
  threshes =  sapply(names(data), function(x) sapply(names(data[[x]]), function(y) get.thresh(data[[x]][[y]],x,y, dists, N, N2, alpha, correct ,np)   )) %>% matfill(length(dists))

  dimnames(threshes) = list(do.call("paste" , list(dists, "best")), dists)
  thresh
}
###THE FULL PROCEDURE CURRENTLY TAKES THE MINIMUM OF THE CANDIDATE THRESHOLDS TO TRY TO GUARANTEE ERROR BELOW ALPHA

###SEE THE EFFECT THIS HAS IN THE CASE WHERE THE THRESHOLDS ARE KNOWN

###PERFORM THE WHOLE PROCEDURE AND TRACK ITS SUCCESS

data = gen(50000, "norm2", param) %>% matrix(nrow = 50)
ICError.newest = function(data, dists, N = 10, N2 = 10, alpha = 0.05, correct = TRUE, ic = "AIC"){
  if(!class(data) == "paleoTS") n = length(data)
  if(class(data) == "paleoTS") n = length(data$mm)
  best.ix = scores(data, dists, ic) %>% which.min
  best = names(best.ix)
  obs = scores(data,dists, ic) %>% (function(x) x[best.ix] - min(x[-best.ix]))
  threshdata = threshmins(data, best, dists, N, N2, alpha, correct, ic)
  probs = threshdata$Probabilities
  threshes = threshdata$Thresholds
  aprime = threshdata$aPrime
  decision = ifelse(obs < min(threshes), best, "Null Decision")
  list("Decision" = decision, "Observed Difference" = obs, "Thresholds" = threshes, "Min AIC" = best , "aPrime" = aprime, "Probabilities"= probs)
}



system.time(est(gen(25, "seg", param.seg), "seg"))
system.time(gen.best(25, "seg", param.seg, "smooth2",dists, 100))
gen.best(25)


###threshpar.est calls biascorrect

###optimizing


#ic.probs.true 1000 times

#threshpar.aprime - est.probs and threshpar.est in different steps
###threshpar.est and est.probs both call biascorrect
####threshpar.est: fit =  biascorrect(data, true, best, dists, N, full = FALSE)
###est.probs: fit =   biascorrect(data, true, best, dists, N, full = FALSE)


####threshmins













### Using this (mu,sigma) gives approximately 1/3 probability of each model scoring best
param = c(0.2135291, 1.1237613)

apply.vec = function(X, MARGIN, sap = TRUE, FUN) {
  if(is.null(dim(X))){
    if(sap == TRUE){
      sapply(X, FUN)
    }
    else{
      FUN = as.character(eval(ICError.sim(25, true, param, dists, 100, 100, alpha)))
      if("function(x)" %in% FUN){
        splitstr(FUN, "(")
      }
      call(as.character(quote(FUN)), X)
    }
  }
  else{
    apply(X, MARGIN, FUN)
  }
}


####
ICError.sim.range = function(ns = c(25, 50, 100), true, params, dists, N = 1000, N2 = 1000, alphas  = c(0.01, 0.05, 0.1), correct = TRUE, runs = 1000){
  names(alphas) = paste("alpha =", alphas)
  names(ns) = paste("n =", ns)
  lapply(alphas,
         function(alpha) lapply(ns,
                                function(n)ICError.sim(n, true, param, dists, N, N2, alpha, correct, runs) ))
}

### Verify model MLE, look at conditional MLEs
ic.MLE.true(true, param, dists, 50, 1000)

### Get true model probabilites
ic.probs.true(true, param2, dists, 50, 10000)
test.thresh(-1.986868, best, true, dists, param, 50, 100000)
threshdirect("norm0", "norm2", dists, param, 50, 10000, 0.05)


### Generate some data
data1 = gen(n, true, param)

### Get Estimated model probabilities, no correction
data3 = gen.best(50, "norm2", param2, "norm0", 1)

est.probs(data, dists, true, N = 1000, N2 = 1000, correct = FALSE)
est.probs(data, dists, true, N = 10000, N2 = 10000, correct = TRUE)


data = gen(25, "norm2", param)
ICError.newest(data, dists, 1000, 1000, alpha)

ICError.sim = function(n, true, param, dists, N = 1000, N2 = 1000, alpha = 0.05, correct = TRUE, runs = 1000, ic){
  data = gens(n, true, param, runs)
  results =   apply(data, 2, function(x) ICError.newest(x, dists, N, N2, alpha, correct, ic))

  list("True" = true, "Results" = results)
}

ICError.sim.reshape = function(sim){
  outs = c("Decision", "Observed Difference", "Thresholds",
           "Min AIC", "aPrime", "Probabilities")
  true = sim$True
  sim = sim$Results
  k = length(dists)
  names(outs) = outs
  sim2 = lapply(outs, function(x) sapply(sim, function(y) y[[x]]))
  spec = length(sim2$Decision[sim2$Decision == true])/length(sim2$Decision[sim2$`Min AIC`==true])

  sim2$Decision %<>% factor(levels = c(dists, "Null Decision"))
  sim.mat = lapply(sim2, function(x) matrix(x, byrow = TRUE, nrow = length(sim2$Decision)))
  colnames(sim.mat$Probabilities) = rep(dists, length(dists) - 1)

  split.ix = sapply(dists, function(x) which(sim.mat[["Min AIC"]]==x ))
  split.sim = lapply(split.ix, function(x) lapply(sim.mat, function(y) y[x,] %>% matrix(nrow = length(x))))

  psplit = split(1:(k*(k-1)), rep(1:(k-1), each = k))
  simnames = function(dist) {
    if(!is.null(split.sim[[dist]])){
      dist.ix = which(dists == dist)
      if(length(split.sim[[dist]]$Decision) == 1) {
        names(split.sim[[dist]]$Thresholds) = dists[-dist.ix]
        names(split.sim[[dist]]$aPrime) = dists[-dist.ix]
      }
      if(length(split.sim[[dist]]$Decision) > 1) {
        colnames(split.sim[[dist]]$Thresholds) = dists[-dist.ix]
        colnames(split.sim[[dist]]$aPrime) = dists[-dist.ix]
      }
      split.sim[[dist]]
    }
  }
  split.sim[[true]]$Specificity = spec
  out = lapply(dists, simnames)
  out
}

ICError.sim.compare = function(split.sim, n, true, param, dists, alpha, correct = TRUE){
  true.ix = which(dists==true)
  trueProbs = ic.probs.true(true, param, dists, n, 10, FALSE)
  trueThreshes = sapply(dists[-true.ix], function(x) threshdirect(x, true, dists, param, n, 10, alpha))
  spec = split.sim[[true]]$Specificity
  AICbests = sapply(dists, function(x) split.sim[[x]]$`Min AIC` %>% length)
  decisions = unlist(sapply(dists, function(x) split.sim[[x]]$Decision))%>%factor(levels = c(dists, "Null Decision")) %>% table
  error = sum(decisions[dists[-true.ix]])/sum(decisions)
  power = decisions[dists[true.ix]]/AICbests[dists[true.ix]]


  simProbs =  sapply(dists[-true.ix], function(x){
    if(nrow(split.sim[[x]]$Probabilities) ==1)
      split.sim[[x]]$Probabilities
    if(nrow(split.sim[[x]]$Probabilities>1))
      split.sim[[x]]$Probabilities %>% colMeans
  })
  simThreshes = sapply(dists[-true.ix], function(x){
    if(nrow(split.sim[[x]]$Thresholds) ==1)
      split.sim[[x]]$Thresholds
    if(nrow(split.sim[[x]]$Thresholds>1))
      split.sim[[x]]$Thresholds %>% colMeans
  })
  obs = as.vector(unlist(sapply(dists, function(x) split.sim[[x]]$`Observed Difference`)))
  mins = as.vector(unlist(sapply(dists, function(x) split.sim[[x]]$`Min AIC`)))
  BA.decision = sapply(1:length(obs), function(x) ifelse(obs[x] < -10, mins[x], "Null Decision")) %>%factor(levels = c(dists, "Null Decision")) %>%table
  BA.spec = length(BA.decision[BA.decision==true])/length(mins[mins==true])
  list("Decisions" = decisions, "Best AICs" = AICbests, "BA Decisions" = BA.decision, "Observed Error Rate" = error, "Specificity.EC" = spec, "Specificity.BA" = BA.spec, "Est Probabilities" = simProbs,"True Probabilities" = trueProbs,  "Est Thresholds" = simThreshes, "True Thresholds" = trueThreshes)
}



alpha.n.test = function(true, param, dists){
  results.n2.25.01  = ICError.sim(25,  true, param, dists, 1000, 1000, 0.01, TRUE, 1000)
  results.n2.50.01  = ICError.sim(50,  true, param, dists, 1000, 1000, 0.01, TRUE, 1000)
  results.n2.100.01 = ICError.sim(100, true, param, dists, 1000, 1000, 0.01, TRUE, 1000)

  results.n2.25.05 =  ICError.sim(25,  true, param, dists, 1000, 1000, 0.05, TRUE, 1000)
  results.n2.50.05 =  ICError.sim(50,  true, param, dists, 1000, 1000, 0.05, TRUE, 1000)
  results.n2.50.10 =  ICError.sim(100, true, param, dists, 1000, 1000, 0.05, TRUE, 1000)

  results.n2.25.10  = ICError.sim(25,  true, param, dists, 1000, 1000, 0.1, TRUE, 1000)
  results.n2.50.10  = ICError.sim(50,  true, param, dists, 1000, 1000, 0.1, TRUE, 1000)
  results.n2.100.10 = ICError.sim(100, true, param, dists, 1000, 1000, 0.1, TRUE, 1000)
  return(list("alpha = 0.01" = list("n = 25" = results.n2.25.01, "n = 50" = results.n2.50.01, "n = 100" = results.n2.100.01),
              "alpha = 0.05" = list("n = 25" = results.n2.25.05, "n = 50" = results.n2.50.05, "n = 100" = results.n2.100.05),
              "alpha = 0.1" =  list("n = 25" = results.n2.25.10, "n = 50" = results.n2.50.10, "n = 100" = results.n2.100.10)))
}

c.25.01  = results.n2.25.01  %>% ICError.sim.reshape %>% ICError.sim.compare(25,  true, param, dists, 0.01, TRUE)
c.50.01  = results.n2.50.01  %>% ICError.sim.reshape %>% ICError.sim.compare(50,  true, param, dists, 0.01, TRUE)
c.100.01 = results.n2.100.01 %>% ICError.sim.reshape %>% ICError.sim.compare(100, true, param, dists, 0.01, TRUE)

c.25.05  = results.n2.25.05  %>% ICError.sim.reshape %>% ICError.sim.compare(25,  true, param, dists, 0.05, TRUE)
c.50.05  = results.n2.50.05  %>% ICError.sim.reshape %>% ICError.sim.compare(50,  true, param, dists, 0.05, TRUE)
c.100.05 = results.n2.100.05 %>% ICError.sim.reshape %>% ICError.sim.compare(100, true, param, dists, 0.05, TRUE)

c.25.10  = results.n2.25.10  %>% ICError.sim.reshape %>% ICError.sim.compare(25,  true, param, dists, 0.1, TRUE)
c.50.10  = results.n2.50.10  %>% ICError.sim.reshape %>% ICError.sim.compare(50,  true, param, dists, 0.1, TRUE)
c.100.10 = results.n2.100.10 %>% ICError.sim.reshape %>% ICError.sim.compare(100, true, param, dists, 0.1, TRUE)

# return(list("alpha = 0.01" = list(c.25.01, c.50.01, c.100.01), "alpha = 0.05" = list(c.25.05, c.50.05, c.100.05), "alpha = 0.1" = list(c.25.10, c.50.10, c.100.10)))
# }

results.norms = alpha.n.test(true, param,dists)

results.norms$`alpha = 0.05`$`n = 100`
results.norms.01 = lapply(results.norms$`alpha = 0.01`, ICError.sim.reshape)
results.norms.05 = list(ICError.sim.reshape(results.norms$`alpha = 0.05`$`n = 25`), ICError.sim.reshape(results.norms$`alpha = 0.05`$`n = 50`), out)

results.norms.10 = lapply(results.norms$`alpha = 0.1`, ICError.sim.reshape)

compare.01 = lapply(1:3, function(x) ICError.sim.compare(results.norms.01[[x]], c(25, 50, 100)[x], true, param, dists, 0.01, TRUE))
compare.05 = lapply(1:3, function(x) ICError.sim.compare(results.norms.05[[x]], c(25, 50, 100)[x], true, param, dists, 0.05, TRUE))
compare.10 = lapply(1:3, function(x) ICError.sim.compare(results.norms.10[[x]], c(25, 50, 100)[x], true, param, dists, 0.10, TRUE))
compare.outs = list("Decisions", "Best AICs", "BA Decisions", "Observed Error Rate", "Specificity.EC", "Specificity.BA",
                    "Est Probabilities","True Probabilities",  "Est Thresholds", "True Thresholds")

lapply(compare.05, function(y) setNames(lapply(compare.outs, function(x) y[[x]] ), compare.outs))
threshmat = c(-1.992 ,-1.880, -1.611 , -1.990, -1.843  , -1.534 ,-1.987 , -1.972,      -1.475  ) %>% matrix(ncol = 3, byrow = TRUE)
sapply(1:3, function(col) sapply(threshmat[,col], function(x) test.thresh(x, "norm0", "norm2", dists, param, c(25, 50, 100)[col], 10000)))

sapply(c(25,50,100), function(n) threshdirect("norm0", "norm2", dists, param, n, 10000, 0.01))
sapply(c(25,50,100), function(n) threshdirect("norm0", "norm2", dists, param, n, 10000, 0.01))

sapply(c(25,50,100), function(n) threshdirect("norm0", "norm2", dists, param, n, 10000, 0.01))


p1n0 = c(0.512, 0.357, 0.154)
p1n1 = c(0.295, 0.340, 0.323)
N = 10
N2 = 10
data  = lapply(ns, function(n) gen(n*2, true, param)  %>% matrix(nrow = n))
aprime.data = lapply(setNames(1:3, ns), function(x) apply(data[[x]],2, function(data) threshdirect(best, true, dists, 100, 100, alpha = p1n1[x], correct = TRUE)$Threshold))

threshpar.aprime(data[[1]][,1], best, true, dists, 10, 10, alpha)
threshdirectdata01n1 = lapply(setNames(1:3, ns), function(x) sapply(1:1000, function(y) threshdirect(best, true, dists, param, ns[x], 10000, 0.01/p1n1[x])))
threshdirectdata05n1 = lapply(setNames(1:3, ns), function(x) sapply(1:1000, function(y) threshdirect(best, true, dists, param, ns[x], 10000, 0.05/p1n1[x])))
threshdirectdata10n1 = lapply(setNames(1:3, ns), function(x) sapply(1:1000, function(y) threshdirect(best, true, dists, param, ns[x], 10000, 0.10/p1n1[x])))

threshdirectdata01n0 = lapply(setNames(1:3, ns), function(x) sapply(1:1000, function(y) threshdirect(best, true, dists, param, ns[x], 10000, 0.01/p1n0[x])))
threshdirectdata05n0 = lapply(setNames(1:3, ns), function(x) sapply(1:1000, function(y) threshdirect(best, true, dists, param, ns[x], 10000, 0.05/p1n0[x])))
threshdirectdata10n0 = lapply(setNames(1:3, ns), function(x) sapply(1:1000, function(y) threshdirect(best, true, dists, param, ns[x], 10000, 0.10/p1n0[x])))


lapply(1:3, function(x) sapply(threshdirectdata01n1[[x]], function(y) test.thresh(y, "norm1", "norm2", dists, param, ns[x], 10000)))
lapply(1:3, function(x) sapply(threshdirectdata05n1[[x]], function(y) test.thresh(y, "norm1", "norm2", dists, param, ns[x], 10000)))
lapply(1:3, function(x) sapply(threshdirectdata10n1[[x]], function(y) test.thresh(y, "norm1", "norm2", dists, param, ns[x], 10000)))

lapply(1:3, function(x) sapply(threshdirectdata01n0[[x]], function(y) test.thresh(y, "norm0", "norm2", dists, param, ns[x], 10000)))
lapply(1:3, function(x) sapply(threshdirectdata05n0[[x]], function(y) test.thresh(y, "norm0", "norm2", dists, param, ns[x], 10000)))
lapply(1:3, function(x) sapply(threshdirectdata10n0[[x]], function(y) test.thresh(y, "norm0", "norm2", dists, param, ns[x], 10000)))


best.ix = which(names(dists)==best)
true.ix = which(names(dists)==true)
data = gen.best(n, true, param, best, N)
scores.data = apply(data, 2, function(x) scores(x, dists))
difs = apply(scores.data, 2, function(x) x[best.ix] - min(x[-best.ix])) %>%  sort

length(difs[difs < thresh])/N

ICError.sim(25, "smooth2", c(0.4, 0.1, .9),dists,100, 100, 0.1, TRUE, 20)


