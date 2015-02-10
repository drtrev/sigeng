# Keep simple

rm(list=ls())

###
# Diagnostics:
# outliers, normality
# normality, outliers

# outliers:
# classical1, classical2, boxplot, robust
# on groups, or on residuals (between model), or on residuals (after removing b/w subj variance) <-- if using lme then lme resids?
# other possibilities: on projected within residuals

# normality:
# KS with correction or shapiro

# 2(diagnostics order) * 4(outliers) * 3(on) * 2(normality) = 48 combinations

###
# analysis:
# within-subjects anova (type 1 or 2 SS), linear mixed effect (random slope for effect)?
# if not-normal: robust, non-parametric

# 48 * 3 = 144 (normal)
# or 48 * 2 = 96 (not normal)

###
# TODO
# covariate age, gender?
# TODO regression adding params in different order?
# TODO investigate R simulate function which allows you to simulate based on an aov model.

library(ggplot2)

source("analyse.r")

generate.dat <- function(N=30)
  # kind of within subj's but may not have sphericity TODO
{
  dat <- data.frame(id=factor(1:N), group=factor(c(rep(paste0("A", 1:2), each=N), rep(paste0("B", 1:2), each=N))), value=rnorm(N*4))
  # split group into separate factors
  dat$factor1 <- factor(substring(dat$group, 1, 1))
  dat$factor2 <- factor(substring(dat$group, 2, 2))
  dat
}

generate.dat.between <- function(N=30)
  # N is per group
{
  dat <- data.frame(id=factor(1:(N*4)), group=factor(c(rep(paste0("A", 1:2), each=N), rep(paste0("B", 1:2), each=N))), value=rnorm(N*4))
  # split group into separate factors
  dat$factor1 <- factor(substring(dat$group, 1, 1))
  dat$factor2 <- factor(substring(dat$group, 2, 2))
  dat
}

generate.dat.within <- function(N=30)
  # within subjects data
{
  idmeans <- data.frame(id=factor(1:N), idmean=rnorm(N))
  dat <- data.frame(id=factor(1:N), group=factor(c(rep(paste0("A", 1:2), each=N), rep(paste0("B", 1:2), each=N))))
  dat <- join(dat, idmeans, by="id")
  #tail(dat)
  dat$value <- dat$idmean + rnorm(N*4) # add condition effect (0) to participant mean
  # split group into separate factors
  dat$factor1 <- factor(substring(dat$group, 1, 1))
  dat$factor2 <- factor(substring(dat$group, 2, 2))
  dat
}

sim <- function(analyses, N=30)
{ 
  dat <- generate.dat(N)
  
  out <- NULL
  for (i in 1:nrow(analyses))
  {
    temp <- analyse(dat, analyses[i,])
    out <- rbind(out, temp$analysis)
  }
  #out
  #head(out)

  min.p <- min(c(min(out$factor1.pval), min(out$factor2.pval), min(out$factor1.2.pval)))
  min.p
}

source("analyse.r")
# just one for testing
analyses.test <- expand.grid(diag.order=factor(c("outliers-normality")),
                        outliers=factor(c("boxplot")),# "robust")),
                        outliers.notnormal=factor(c("boxplot")), #, "robust")),
                        outliers.response=factor(c("remove")),
                        normality.func=factor(c("Shapiro-Wilk")),
                        normality.on=factor(c("resids")),
                        analysis=factor(c("anova.type3")))#, "lme")))

sim(analyses.test)
dat <- generate.dat(100); out <- analyse(dat, analyses.test[1,])
#

# TODO does SD2 make sense when not normal? could have outliers.notnormal.response
# but may be normal after outlier is transformed?

# diag=diagnostics
analyses <- expand.grid(diag.order=factor(c("outliers-normality", "normality-outliers")),
                        outliers=factor(c("classical1", "classical2", "boxplot")),# "robust")),
                        outliers.notnormal=factor(c("boxplot")), #, "robust")),
                        outliers.response=factor(c("SD2", "remove")),
                        normality.func=factor(c("KS", "Shapiro-Wilk")),
                        normality.on=factor(c("groups", "resids")),
                        analysis=factor(c("anova.type3", "lme")))
# not using anova.type2 for now (just within subjects so when there's a removed outlier we have
# to remove whole subject)
class(analyses)
analyses

#pvals <- aaply(1:100, 1, function(x) sim(analyses))
# 1min 28 (88 seconds for 5 sims)
# from analyses-timetest analyses
#save(analyses, file="analyses-timetest.RData")
#pvals <- aaply(1:5, 1, function(x) sim(analyses))
#mean(pvals < .05) # 24%

?foreach
library(foreach)
#library(doParallel)
#library(multicore)
#library(snow)
#?parallel
registerDoParallel()
#registerDoMC()
#registerDoSNOW()
# from foreach package, useful when there are no varying arguments:

# 23 % !!
# TODO what about with tiny effect when generating data?

# TODO generate data using built-in simulate() function!

# TODO sphericity

# TODO order effect

# TODO how much does each decision contribute to the false-positive rate?

# may hang without registering a parallel backend, e.g. registerDoParallel()
nreps <- 100
system.time(
out <- foreach(colcurr=names(analyses), .combine=rbind) %do%
{
  colcut <- analyses[,colcurr]
  if (length(levels(colcut))!=length(unique(colcut)))
  {
    warning("Factor levels have changed, refactoring")
    colcut <- factor(colcut)
  }
  
  foreach(l=levels(colcut), .combine=rbind) %do%
  {
    analyses.sub <- analyses[colcut==l,]
    
    # may need to use .packages or source/load.packages, seems to be needed after registerDoParallel
    pvals <- foreach(i=1:nreps, .combine=data.frame) %dopar%
    {  
      load.packages()
      pvals <- try(sim(analyses.sub))
      if (class(pvals)=="try-error")
      {
        cat("Caught error\n")
        pvals <- NA
      }
      pvals
    }
    pvals$nsig <- sum(pvals[1,] < 0.05)
    pvals$colcut <- colcurr
    pvals$collevel <- l
    print(pvals)
    pvals
  }
}
) # takes 2 hours

head(out)
#save(out, file="out-holdEachLevel.RData")
#save(analyses, file="out-holdEachLevelAnalyses.RData") # save associated analyses data frame
nrow(out)
# what range of results do we expect anyway?
out.all <- foreach(j=1:nrow(out), .combine=rbind) %do%
{
  pvals <- foreach(i=1:nreps, .combine=data.frame) %dopar%
  {  
    load.packages()
    pvals <- try(sim(analyses.sub))
    if (class(pvals)=="try-error")
    {
      cat("Caught error\n")
      pvals <- NA
    }
    pvals
  }
  pvals$nsig <- sum(pvals[1,] < 0.05)
  print(pvals)
  pvals
}
out.all
#save(out.all, file="out-all.RData")
#save(analyses, file="out-allAnalyses.RData")

# 26 +- 7 covers most of the data
#qplot(out.all$nsig)
#mean(out.all$nsig) # == 25.6
#sd(out.all$nsig) * 2 # == 6.7
#qplot(out$nsig)

# I think KS is the only one that might be above chance.
out[order(out$nsig),]

f <- aaply(.data=1:5, .margins=1, .fun=function(x) {print(x); x}, .parallel=T)
