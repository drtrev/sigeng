# Keep simple

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
                        analysis=factor(c("anova.type2", "anova.type3", "lme")))
class(analyses)
analyses

source("analyse.r")

sim <- function(analyses, N=30)
{
  dat <- data.frame(id=factor(1:N), group=factor(c(rep(paste0("A", 1:2), each=N), rep(paste0("B", 1:2), each=N))), value=rnorm(N*4))
  # split group into separate factors
  dat$factor1 <- factor(substring(dat$group, 1, 1))
  dat$factor2 <- factor(substring(dat$group, 2, 2))
  
  
  #head(dat); tail(dat)
  #str(dat)
  
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

