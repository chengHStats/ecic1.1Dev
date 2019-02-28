library(microbenchmark)

for (size in c(1e3, 1e4)){
  x <- sample(1:3, size, TRUE)
  xint <-  as.integer(x)
  xfac = factor(x, levels = c("1","2","3"))
  xintchar = as.character(x)
  xnames <-  sample(models, size, TRUE)
  x = rep(1:2, 100000)
  microbenchmark(
    "tabulate" = {
      tabulate(x)
    },
    "tabulate + 1" = {
      tabulate(c(1:3,x))-1
    },
    "table int" = {
      table(x)
    },
    "table char" = {
      table(xnames)
    },
    "table fac" = {
      table(factor(x, levels = 1:3))
    }
  )%>%print
}

table(
  factor(apply(scores,1,function(x) models[which.min(x)]),
         levels = models))/N

x = runif(30000) %>% matrix(ncol = 3)
microbenchmark(
  "which.min" = {
    apply(x, 1, which.min)
  }
)

x = matrix(runif(900000), ncol = 9000)
microbenchmark(
  "anonymous" = {
    apply(x, 2, function(y) sum(y))
  },
  "named" = {
    apply(x, 2, sum)
  }
)

