# Test ANOVA type I error
rm(list=ls())

#library(ez)

generate.dat <- function(N=30)
  # N is number of subjects per group
{
  dat <- data.frame(id=factor(1:(N*4)), group=factor(c(rep(c("A1", "A2"), each=N), rep(c("B1", "B2"), each=N))), value=rnorm(N*4))
  # split group into separate factors
  dat$factor1 <- factor(substring(dat$group, 1, 1))
  dat$factor2 <- factor(substring(dat$group, 2, 2))
  dat
}

#do.anova2 <- function(dat)
#{
#  results <- ezANOVA(dat, dv=value, wid=id, between=.(factor1, factor2))
#  min(results$ANOVA$p)
#}

do.anova <- function(dat)
{
  results <- aov(value ~ factor1 * factor2, data=dat)
  min( summary(results)[[1]]$`Pr(>F)`[1:3] )
}

out <- replicate(1000, do.anova(generate.dat()))
sum(out < .05) / 1000

# e.g.
# dat <- generate.dat(3)
# print(dat)