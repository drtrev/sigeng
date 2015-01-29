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

#sim(analyses.test)
dat <- generate.dat(100); out <- analyse(dat, analyses.test[1,])
#

# TODO does SD2 make sense when not normal? could have outliers.notnormal.response

# diag=diagnostics
analyses <- expand.grid(diag.order=factor(c("outliers-normality", "normality-outliers")),
                        outliers=factor(c("classical1", "classical2", "boxplot")),# "robust")),
                        outliers.notnormal=factor(c("boxplot")), #, "robust")),
                        outliers.response=factor(c("SD2", "remove")),
                        normality.func=factor(c("KS", "Shapiro-Wilk")),
                        normality.on=factor(c("groups", "resids")),
                        analysis=factor(c("anova.type3")))#, "lme")))
# not using anova.type2 for now (just within subjects so when there's a removed outlier we have
# to remove whole subject)
class(analyses)
analyses

pvals <- aaply(1:100, 1, function(x) sim(analyses))
mean(pvals < .05) # 24%
pvals

# 23 % !!
# TODO what about with tiny effect when generating data?

# TODO sphericity

# TODO order effect

# TODO how much does each decision contribute to the false-positive rate?

pvals.all <- c(pvals, pvals2)
sum(pvals.all < .05)

