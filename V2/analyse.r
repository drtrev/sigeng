library(lme4)
library(ez)
library(nortest)
library(plyr)

analyse <- function(dat, analysis)
  # Given analysis id, do a specific analysis/order:
{
  test <- F; if (test) analysis <- analyses[1,]
  
  outliers <- function(dat, analysis, norm=T)
  {
    outliers.sub <- function(dat, outliers.method, outliers.response)
    {
      test <- F; if (test) outliers.response <- analysis$outliers.response
      
      outliers.response.func <- function(value, MEAN, SD, outliers.response)
      {
        value.out <- NULL
        if (outliers.response=="SD2") value.out <- ifelse(value > MEAN, MEAN + 2 * SD, MEAN - 2 * SD)
        if (outliers.response=="remove") value.out <- NA
        if (is.null(value.out)) stop(paste0("Error with outliers response; not found: ", outliers.response))
        
        cat("outliers.response:", as.character(outliers.response), "\n")
        cat("Outlier converted to", as.character(value.out), "\n")
        value.out
      }
      
      classical1 <- function(datgroup, outliers.response)
      {
        classical1.sub <- function(value, MEAN, SD, outliers.response, debug.group)
        {
          if (value > MEAN + SD * 3.29 || value < MEAN - SD * 3.29)
          {
            cat("Outlier group:", debug.group, "id:", names(value), "value:", value, "\n")
            value <- outliers.response.func(value, MEAN, SD, outliers.response)
          }else{
            # cat("value:", value, " not > ", MEAN + SD * 3.29, "\n")
          }
          value
        }
        values <- datgroup$value
        names(values) <- 1:length(values)
        datgroup$value <- aaply(values, 1, classical1.sub, mean(datgroup$value), sd(datgroup$value), outliers.response, as.character(datgroup$group[1]))
        datgroup
      }
      
      classical2 <- function(datgroup, outliers.response)
      {
        classical2.sub <- function(value, MEAN, SD, outliers.response, debug.group)
        {
          if (value > MEAN + SD * 3 || value < MEAN - SD * 3)
          {
            cat("Outlier group:", debug.group, "id:", names(value), "value:", value, "\n")
            value <- outliers.response.func(value, MEAN, SD, outliers.response)
          }else{
            # cat("value:", value, " not > ", MEAN + SD * 3.29, "\n")
          }
          value
        }
        values <- datgroup$value
        names(values) <- 1:length(values)
        # TODO change to values below in mean sd and in func above
        datgroup$value <- aaply(values, 1, classical2.sub, mean(datgroup$value), sd(datgroup$value), outliers.response, as.character(datgroup$group[1]))
        datgroup
      }
      
      boxplotout <- function(datgroup, outliers.response)
      {
        # 1.5 IQR, apparantly from Tukey's work on boxplots
        boxplotout.sub <- function(value, q1, q3, MEAN, SD, outliers.response, debug.group)
        {
          iqr <- q3 - q1
          if (value > q3 + iqr * 1.5 || value < q1 - iqr * 1.5)
          {
            cat("Boxplot outlier group", debug.group, "id:", names(value), "value:", value, "\n")
            value <- outliers.response.func(value, MEAN, SD, outliers.response)
          }
          value
        }

        # test:
        #datgroup <- data.frame(group="A1", value=rnorm(100))
        #outliers.response <- "SD2"
        
        values <- datgroup$value
        names(values) <- 1:length(values)
        datgroup$value <- aaply(values, 1, boxplotout.sub, quantile(values, 1/4), quantile(values, 3/4), mean(values), sd(values), outliers.response, as.character(datgroup$group[1]))
        datgroup
      }
      
      
      dat2 <- NULL
      if (outliers.method=="classical1")
      {
        #debugonce(classical1)
        dat2 <- ddply(dat, .(group), classical1, outliers.response)
      }
      if (outliers.method=="classical2")
      {
        dat2 <- ddply(dat, .(group), classical1, outliers.response)          
      }
      if (outliers.method=="boxplot")
      {
        dat2 <- ddply(dat, .(group), boxplotout, outliers.response)
      }
      #if (outliers.method=="robust")
      #{
      #  dat2 <- dat # TODO
      #}
      if (is.null(dat2)) stop(paste0("Outliers method not found: ", outliers.method))
      
      dat2$outlier <- !(dat2$value==dat$value)
      dat2
    }
    
    if (norm)
    {
      dat <- outliers.sub(dat, analysis$outliers, analysis$outliers.response)
    }else{
      dat <- outliers.sub(dat, analysis$outliers.notnormal, analysis$outliers.response)
    }
    dat
  }
  
  
  normality <- function(dat, analysis)
  {
    normality.sub <- function(values, normality.func)
      # call for each group separately
    {
      if (normality.func=="KS")
      {
        return(lillie.test(values)$p.value)
      }
      else if (normality.func=="Shapiro-Wilk")
      {
        return(shapiro.test(values)$p.value)
      }
      else
      {
        stop(paste0("normality.func not found: ", normality.func))
      }
    }
    
    # Call normality.sub for each group
    if (analysis$normality.on=="groups")
    {
      normal.p <- ddply(dat, .(group), function(x) normality.sub(x$value, analysis$normality.func))
    
      # For now, normal if < half are sig
      if (sum(normal.p$V1 < .05) > 0) cat(paste0("Groups not-normal: ", sum(normal.p$V1 < .05), "\n"))
      if (sum(normal.p$V1 < .05) < nrow(normal.p) / 2)
      {
        normal <- T
      }else{
        normal <- F
      }
    }
    else if (analysis$normality.on=="resids")
    {
      # Remove ids that have an NA
      remove.ids <- dat[is.na(dat$value),"id"]
      dat <- dat[!(dat$id %in% remove.ids),]
      dat$id <- factor(dat$id)
      
      aov1 <- aov(value ~ factor1*factor2+Error(id/(factor1*factor2)), data=dat)
      #print(summary(aov1))
      #print(head(dat))
      pvals <-          normality.sub(resid(aov1[[3]]), analysis$normality.func) # factor1
      pvals <- c(pvals, normality.sub(resid(aov1[[4]]), analysis$normality.func)) # factor2
      pvals <- c(pvals, normality.sub(resid(aov1[[5]]), analysis$normality.func)) # factor1:2

      # From MASS:
      aov1p <- proj(aov1)
      pvals2 <-           normality.sub(aov1p[[3]][,"Residuals"], analysis$normality.func) # factor1
      pvals2 <- c(pvals2, normality.sub(aov1p[[4]][,"Residuals"], analysis$normality.func)) # factor2
      pvals2 <- c(pvals2, normality.sub(aov1p[[5]][,"Residuals"], analysis$normality.func)) # factor1:2
      
      normal <- T
      num.notnorm <- sum(pvals < .05)
      if (num.notnorm > 1)
      {
        cat(paste0(num.notnorm, " sets of resids are not normal\n"))
        normal <- F
      }
      
      num.notnorm <- sum(pvals2 < .05)
      if (num.notnorm > 1)
      {
        cat(paste0(num.notnorm, " sets of projected resids are not normal\n"))
        normal <- F
      }
      
    }else
    {
      stop(paste0("normality.on not found: ", normality.on))
    }

    normal
  }
  
  
  analysis.anova.type2 <- function(dat)
  {
    # For within-subjects, if participant has missing value then remove participant
    remove.ids <- dat[is.na(dat$value),"id"]
    dat <- dat[!(dat$id %in% remove.ids),]
    dat$id <- factor(dat$id)
    out <- ezANOVA(dat, dv=value, wid=id, within=.(factor1, factor2), type=2)
    out$ANOVA$p
  }

  analysis.anova.type3 <- function(dat)
  {
    # TODO could use aov here, faster (and type 3)
    # For within-subjects, if participant has missing value then remove participant
    remove.ids <- dat[is.na(dat$value),"id"]
    dat <- dat[!(dat$id %in% remove.ids),]
    dat$id <- factor(dat$id)
    out <- ezANOVA(dat, dv=value, wid=id, within=.(factor1, factor2), type=3)
    out$ANOVA$p
  }

  analysis.lme <- function(dat)
  {
    # TODO could also nest factors wihtin ID see field example: participants/animal or something,
    # although that seems dubious (full model?)
    lm1 <- lmer(value ~ 1 + (1|id), data=dat)
    lm2 <- lmer(value ~ factor1 + (1|id), data=dat)
    lm3 <- lmer(value ~ factor1 + factor2 + (1|id), data=dat)
    lm4 <- lmer(value ~ factor1 + factor2 + factor1:factor2 + (1|id), data=dat)
    an1 <- anova(lm1, lm2, lm3, lm4)
    #qplot(resid(lm2))
    pvals <- an1$`Pr(>Chisq)`[2:4]
    pvals
  }
  
  if (analysis$diag.order == "outliers-normality")
  {
    dat <- outliers(dat, analysis)
    norm <- normality(dat, analysis)
  }else{
    norm <- normality(dat, analysis)
    dat <- outliers(dat, analysis, norm)
    if (!norm) norm <- normality(dat, analysis) # may be normal after removing outliers
  }
  
  if (!norm)
  {
    # could transform data TODO Field
    # or just use robust test TODO Field
    cat("Would use robust test now.\n")
  }
  
  pvals <- NULL
  if (analysis$analysis=="anova.type2")
  {
    pvals <- analysis.anova.type2(dat)
  }
  if (analysis$analysis=="anova.type3")
  {
    pvals <- analysis.anova.type3(dat)
  }
  if (analysis$analysis=="lme")
  {
    pvals <- analysis.lme(dat)
  }
  if (is.null(pvals)) stop("Invalid analysis supplied")
  
  analysis$factor1.pval <- pvals[1]
  analysis$factor2.pval <- pvals[2]
  analysis$factor1.2.pval <- pvals[3] # interaction factor1:factor2
  
  analysis$star <- ""
  if (sum(pvals < .05) > 0) analysis$star <- "*"
  
  # return dat (marks outliers) and analysis used (now including p value)
  list(dat=dat, analysis=analysis)
}

test.normal.anova <- function()
{
  # Summary: It will be 14% because there are three F tests.
  
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
  
}