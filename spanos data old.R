spanos.X1 = cbind(1, spanos.x)
spanos.X2 = cbind(1, spanos.x, spanos.x^2, spanos.x^3)

spanos.X3 = cbind(1, spanos.x, spanos.x^2)

bhat = solve(t(X2)%*%X2)%*%t(X2)%*%matrix(df3$pop, ncol = 1) %>%unname
data = df3$pop
sd = sqrt(t(data-X2%*%bhat)%*%(data-X2%*%bhat)/32)

system.time(est(df3$pop, "lm2"))

system.time(replicate(10000,{
  bhat = solve(t(X2)%*%X2)%*%t(X2)%*%matrix(data, ncol = 1) %>%unname
  sd = sqrt(t(data-X2%*%bhat)%*%(data-X2%*%bhat)/32)
  c(bhat, sd)
  }))
system.time(replicate(10000, est(df3$pop, "lm2")))




### predictor variable from plot in Mayo and Spanos 2004
spanos.x = c(1,2,3,5,5,6,8,
      11,11,12,13,15,16,19,19,
      21,22,22,23, 24, 25, 28, 28, 29,
      30,30,32, 34, 35, 35, 37,
      40, 41, 41, 43)


library(rvest)
theurl <- "http://www.multpl.com/united-states-population/table"
file<-read_html(theurl)
tables<-html_nodes(file, "table")
table1 <- html_table(tables[1], fill = TRUE)
df = table1[[1]]
colnames(df) = c("year", "pop")
df2 = df
df2$year = sapply(df$year, function(x) as.numeric(unname(strsplit(x, " ")[[1]][3])))
df2$pop = sapply(df$pop, function(x) as.numeric(strsplit(x, " ")[[1]][1]))
strsplit(df$year[1], " ")

df3 = df2[df2$year %in% 1955:1989,][order(35:1),]
rownames(df3) = NULL
df3$x = spanos.x

lmspanos = lm(pop ~ poly(x, 3, raw = TRUE), data = df3)

dfspanos = df3
pwr.f2.test(2, 33, 0.15, 0.05)



###BOOTSTRAP POWER ESTIMATE


library(pwr)
library(car)

#Set up the linear hypothesis
#This hypothesis matrix will evaluate
#if B2 = B3 = 0, because the
# Coefficient vector is B = c(B0, B1, B2, B3)
# And the hypothesis is of the form
# RB = rhs as below.
R = rbind(c(0,0,1,0),
          c(0,0,0,1))
sigma = summary(lmspanos)$sigma
pred = predict.lm(lmspanos, new = data.frame(x = x))
N = 1000
bootdata = sapply(1:N, function(y) pred + rnorm(35, 0, sigma))
boot.pvals = apply(bootdata, 2, function(y){
  dfboot = data.frame(pop = y, x = x)
  lmboot = lm(pop ~ poly(x, 3, raw=TRUE), data = dfboot)
  linearHypothesis(lmboot, R, rhs = c(0,0))$`Pr(>F)`[2]
})
sum(boot.pvals < 0.05)


