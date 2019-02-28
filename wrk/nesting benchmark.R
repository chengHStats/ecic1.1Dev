x = runif(60000)%>%matrix(ncol = 3)
microbenchmark(
  "nested" = {
    mins = apply(x, 1, which.min)
    tabulate(mins)
  },
  "allocated" = {
    tabulate(apply(x,1,which.min))
  }
)


l1 <- list("a" = 1:3, "b" = 1:3)
l2 <- list("a" = 4:6, "b" = 4)

