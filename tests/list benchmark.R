

CombineLogliksM = function(data, models){

}

CombineLogliksL = function(data, models)


logliks(data, "norm")
flatt = sapply(models, function(x) logliks(data, x)) %>% flatten()
liks = sapply(models, function(x) logliks(data, x))
(function(x) 'colnames')
do.call(cbind, liks)

ICMultiMulti = function(data, models){
  p = length(models)
  n = nrow(data)
  N = ncol(data)
  liks = sapply(models, function(x) logliks(data, x))
  liks1 = sapply(1:p, function(ix){
    x = liks[[ix]]
    'colnames<-'(x, paste(models[ix], colnames(x), sep = '.'))
  } )
  do.call(cbind, liks1)
}



