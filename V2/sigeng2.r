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

library(ggplot2)
library(lme4)
library(plyr)
library(ez)

?expand.grid
analyses <- expand.grid(diag.order=factor(c("outliers-normality", "normality-outliers")),
                        outliers=factor(c("classical1", "classical2", "boxplot", "robust")),
                        outliers.notnormal=factor(c("boxplot", "robust")),
                        outliers.response=factor(c("SD2", "remove")),
                        normality=factor(c("KS", "shapiro-wilk")),
                        analysis=factor(c("anova.type3", "lme")))
# not using anova.type2 for now (just within subjects so when there's a removed outlier we have
# to remove whole subject)
class(analyses)
analyses

source("analyse.r")

generate.dat <- function(N=30)
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

pvals <- aaply(1:100, 1, function(x) sim(analyses))
sum(pvals < .05)
pvals

# 23 % !!
# TODO what about with tiny effect when generating data?

# TODO sphericity

# TODO order effect

pvals.all <- c(pvals, pvals2)
sum(pvals.all < .05)


# What about normally:
analyses.normal <- expand.grid(diag.order=factor("outliers-normality"),
                        outliers=factor("classical1"),
                        outliers.notnormal=factor("boxplot"),
                        outliers.response=factor("SD2"),
                        normality=factor("shapiro-wilk"),
                        analysis=factor("anova.type3"))
args(sim)

pvals.normal <- aaply(1:100, 1, function(x) sim(analyses.normal))
sum(pvals.normal < .05) # 14%

just.anova <- function(dat)
{
  out <- ezANOVA(dat, dv=value, wid=id, within=.(factor1, factor2))
  min(out$ANOVA$p)
}

pvals.just.anova <- aaply(1:100, 1, function(x) just.anova(generate.dat()))
sum(pvals.just.anova < .05)

pvals.just.anova <- aaply(1:100, 1, function(x) just.anova(generate.dat.within()))
sum(pvals.just.anova < .05) # 14%

pvals.just.anova <- aaply(1:1000, 1, function(x) just.anova(generate.dat.within()))
sum(pvals.just.anova < .05) / 1000 # 13%

# why is type I error above 5% for anova?

just.anova.between <- function(dat)
{
  out <- ezANOVA(dat, dv=value, wid=id, between=.(factor1, factor2))
  min(out$ANOVA$p)
}

pvals.just.anova <- aaply(1:100, 1, function(x) just.anova.between(generate.dat.between()))
sum(pvals.just.anova < .05) # 12%


generate.dat.between2 <- function(N=30)
  # N is per group
{
  dat <- data.frame(id=factor(1:(N*2)), group=factor(c(rep("A", each=N), rep("B", each=N))), value=rnorm(N*2))
  dat
}

just.anova.between2 <- function(dat)
{
  out <- ezANOVA(dat, dv=value, wid=id, between=group)
  min(out$ANOVA$p)
}

pvals.just.anova <- aaply(1:100, 1, function(x) just.anova.between2(generate.dat.between2()))
sum(pvals.just.anova < .05) # 5%
length(pvals.just.anova)

