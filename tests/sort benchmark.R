library(data.table)

x <- runif(1e5)


library(microbenchmark)

for (size in c(1e3, 1e4, 1e5, 1e6)){
  x <- runif(size)
microbenchmark(
  "data.table" = {
   data.table(x, key = "x")
  },
  "sort (quick)" = {
    sort(x, method = "quick")
  },
  "sort (radix)" = {
    sort(x, method = "radix")
  },
  "sort" = {
    sort(x)
  }
)%>%print
}



