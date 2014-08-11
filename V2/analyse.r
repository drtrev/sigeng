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
        switch(
          outliers.response,
          SD2=ifelse(value > MEAN, MEAN + 2 * SD, MEAN - 2 * SD),
          remove=NA,
          "Error: no response found")
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
        datgroup$value <- aaply(values, 1, classical2.sub, mean(datgroup$value), sd(datgroup$value), outliers.response, as.character(datgroup$group[1]))
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
        dat2 <- dat # TODO
      }
      if (outliers.method=="robust")
      {
        dat2 <- dat # TODO
      }
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
    T # TODO
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
    # For within-subjects, if participant has missing value then remove participant
    remove.ids <- dat[is.na(dat$value),"id"]
    dat <- dat[!(dat$id %in% remove.ids),]
    dat$id <- factor(dat$id)
    out <- ezANOVA(dat, dv=value, wid=id, within=.(factor1, factor2), type=3)
    out$ANOVA$p
  }

  analysis.lme <- function(dat)
  {
    # TODO could also nest factors wihtin ID see field example: participants/animal or something
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
