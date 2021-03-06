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
# TODO does SD2 make sense when not normal? could have outliers.notnormal.response
#   but may be normal after outlier is transformed?
# TODO what about with tiny effect when generating data?

# TODO generate data using built-in simulate() function!

# Sphericity: can't have (or can't detect) violations of sphericity for only two levels.
# But data can still be generated in one of two ways: first subj mean then exp effect,
# or simply exp/control effect.

# TODO order effect

# TODO how much does each decision contribute to the false-positive rate?

library(ggplot2)

source("generateData.r")
source("analyse.r")
source("investigate.r")
load.packages()

##########
# just one for testing
analyses.test <- expand.grid(diag.order=factor(c("outliers-normality")),
                        outliers=factor(c("boxplot")),# "robust")),
                        outliers.notnormal=factor(c("boxplot")), #, "robust")),
                        outliers.response=factor(c("remove")),
                        normality.func=factor(c("Shapiro-Wilk")),
                        normality.on=factor(c("resids")),
                        analysis=factor(c("anova.type3")))#, "lme")))

temp <- sim(analyses.test)

dat <- generate.dat.within(100);
out <- analyse(dat, analyses.test[1,])
out$analysis
out$dat

rm(temp)

##########

nWorkers <- 8

source("investigate.r")
#out.all <- investigateRepeatAll(remake=F)
#debugonce(investigateRepeatAll)
#out.all <- investigateRepeatAll(remake=T, nreps=2, analyses=analyses.test)
#head(out.all)

#debugonce(investigateHoldEachLevel)
#out.holdEachLevel <- investigateHoldEachLevel(loadFromCache=F, nreps=2, analyses=initAnalyses()[1:5,], nWorkers)
holdEachLevelList <- investigateHoldEachLevel(simulate=F, nWorkers=8)
interestingCols <- c("nSig", "columnName", "columnLevel", "simulationId")
head(holdEachLevelList$results[interestingCols])
tail(holdEachLevelList$results[interestingCols])

library(doBy)
library(plotrix)
resultsSummary <- summaryBy(nSig ~ columnName + columnLevel, holdEachLevelList$results, FUN=c(mean, std.error))
head(resultsSummary)
resultsSummary$level <- paste(resultsSummary$columnName, resultsSummary$columnLevel, sep="\n")
resultsSummary$ymin <- resultsSummary$nSig.mean - resultsSummary$nSig.std.error
resultsSummary$ymax <- resultsSummary$nSig.mean + resultsSummary$nSig.std.error

# Note outliers.notnormal boxplot may coincidentally be showing the expected results with all combinations
# (because there is only one level for outliers.notnormal)
ggplot(resultsSummary, aes(x=level, y=nSig.mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=ymin, ymax=ymax))

#############
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

# library(foreach)
# library(doParallel)
# # windows doesn't have fork, so have to use snow "cluster" style.
# # cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf
# cl <- makeCluster(2)
# registerDoParallel(cl)
# getDoParWorkers()
# 
# 
# 
# # I think KS is the only one that might be above chance.
# out[order(out$nsig),]
# 
# ######################
# # Check KS on its own
# 
# names(analyses)
# colcurr <- "normality.func"
# colcut <- analyses[,colcurr]
# if (length(levels(colcut))!=length(unique(colcut)))
# {
#   warning("Factor levels have changed, refactoring")
#   colcut <- factor(colcut)
# }
# 
# l <- "KS"
# analyses.sub <- analyses[colcut==l,]
# 
# # may need to use .packages or source/load.packages, seems to be needed after registerDoParallel
# pvals <- foreach(i=1:nreps, .combine=data.frame) %dopar%
# {  
#   load.packages()
#   pvals <- try(sim(analyses.sub))
#   if (class(pvals)=="try-error")
#   {
#     cat("Caught error\n")
#     pvals <- NA
#   }
#   pvals
# }
# pvals$nsig <- sum(pvals[1,] < 0.05)
# pvals$colcut <- colcurr
# pvals$collevel <- l
# print(pvals)
# pvals # nsig 29, so probably not above chance
# 
# 
# 
# f <- aaply(.data=1:5, .margins=1, .fun=function(x) {print(x); x}, .parallel=T)

