############################################################
# Investigate: Run analyses to answer questions
############################################################

library(foreach)
library(doParallel)
# cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf
cl <- makeCluster(2)
registerDoParallel(cl)
#check we have 2:
#getDoParWorkers()

############################################################
# Helper functions
############################################################

getOutputFileName <- function(filePrefix=NULL)
{
  # Param:
  #   filePrefix - some kind of prefix name to give to the file,
  #     e.g. out-all, out-holdEachLevel
  if (is.null(filePrefix))
  {
    stop("Specify file prefix")
  }
  fileNumber <- 1
  fileSuffix <- ".RData"
  analysesSuffix <- "Analyses" # for the analyses data frame
  while (file.exists(paste0(filePrefix, fileNumber, fileSuffix)))
  {
    fileNumber <- fileNumber + 1
  }
  cat("fileNumber:", fileNumber, "\n")
  
  resultsFileName <- paste0(filePrefix, fileNumber, fileSuffix)
  analysesFileName <- paste0(filePrefix, analysesSuffix, fileNumber, fileSuffix)
  list(resultsFileName=resultsFileName, analysesFileName=analysesFileName)
}

############################################################
# Create a set of analyses for use in investigations
############################################################

initAnalyses <- function()
{
  # diag=diagnostics
  analyses <- expand.grid(diag.order=factor(c("outliers-normality", "normality-outliers")),
                          outliers=factor(c("classical1", "classical2", "boxplot")),# "robust")),
                          outliers.notnormal=factor(c("boxplot")), #, "robust")),
                          outliers.response=factor(c("SD2", "remove")),
                          normality.func=factor(c("KS", "Shapiro-Wilk")),
                          normality.on=factor(c("groups", "resids")),
                          analysis=factor(c("anova.type3", "lme")))
  analyses
}

############################################################
# Are some levels causing higher chance of sig?
############################################################

# If you get an error with %do% when running parallel, you
# need to include the foreach package within the loop or use
# .packages="foreach" as a parameter to foreach().
# See parallelDemo.r
holdEachLevel <- function(analyses, nreps)
{
#  system.time(
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
          source("analyse.r")
          load.packages()
          pvals <- try(sim(analyses.sub))
          if (class(pvals)=="try-error")
          {
            cat("Caught error\n")
            pvals <- NA
          }
          pvals
        }
        # pvals looks like this:
        #    result.1  result.2
        #  1   xxxxxx     xxxxx
        # Add a column showing the total number of significant values from the row:
        pvals$nsig <- sum(pvals[1,] < 0.05)
        # Add the name of the column we are currently varying:
        pvals$colcut <- colcurr
        # Add the factor level of colcurr we just used:
        pvals$collevel <- l
        print(pvals)
        pvals
      }
  }
#  ) # takes 2 hours
  
  # out looks like this:
  #         result.1     result.2 ... nsig             colcut           collevel
  #   1  0.064134139 0.0236179504 ...    1         diag.order normality-outliers
  #   2  0.284763115 0.0755375576 ...    0         diag.order outliers-normality
  #   3  0.106062614 0.1656168818 ...    0           outliers            boxplot
  #   4  ...
  out
}

############################################################

investigateHoldEachLevel <- function(remake=F, nreps=100, analyses=NULL)
{
  if (remake)
  {
    if (is.null(analyses))
    {
      cat("Using default analyses\n")
      analyses <- initAnalyses()
    }
    out.holdEachLevel <- holdEachLevel(analyses, nreps)
    fileNames <- getOutputFileName("out-holdEachLevel")
    cat("resultsFileName:", fileNames$resultsFileName, "\n")
    cat("analysesFileName:", fileNames$analysesFileName, "\n")
    
    save(out.holdEachLevel, file=fileNames$resultsFileName)
    save(analyses, file=fileNames$analysesFileName) # save associated analyses data frame
  }
  else
  {
    # TODO use file number
    load("out-holdEachLevel.RData")
    load("out-holdEachLevelAnalyses.RData")
  }
  out <- list(holdEachLevel=out.holdEachLevel, analyses=analyses)
}

############################################################
# What range of nsig results do we expect anyway?
############################################################

repeatAll <- function(analyses, out, nreps)
{
  # The following took 522 seconds (8.7mins) with 2 workers, nreps==3 and nrow(out)==14,
  #   with the outer loop as %do% and the inner loop as %dopar%.
  # 442 seconds (7.3mins) with the outer loop as %dopar% and the inner loop as %do%,
  #   but with source("analyse.r") and load.packages() in the inner loop.
  # 445 seconds same but with source and load.packages in outer loop.
  # TODO What's the SD on time measurements?
  system.time(out.all <- foreach(j=1:nrow(out), .combine=rbind) %do%
  {
    library(foreach)    
    # The following (inner) foreach loop took:
    # 37 seconds for 3 reps, parallel (cluster of 2)
    # 53 seconds for 3 reps, not parallel
    pvals <- foreach(i=1:nreps, .combine=data.frame) %do%
    {
      source("analyse.r")
      load.packages()
      
      pvals <- try(sim(analyses))
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
  })
  out.all
}

############################################################

investigateRepeatAll <- function(remake=F, nreps=100, analyses=NULL)
{
  if (remake)
  {
    # TODO make this work with fileNumber
    out <- investigateHoldEachLevel(remake=F)
    out.holdEachLevel <- out$holdEachLevel
    out.analyses <- analyses
    if (is.null(analyses))
    {
      out.analyses <- out$analyses
    }
    
    out.all <- repeatAll(out.analyses, out.holdEachLevel, nreps)
    fileNames <- getOutputFileName("out-all")
    cat("resultsFileName:", fileNames$resultsFileName, "\n")
    cat("analysesFileName:", fileNames$analysesFileName, "\n")
    
    save(out.all, file=fileNames$resultsFileName)
    save(analyses, file=fileNames$analysesFileName)
  }
  else
  {
    load("out-all.RData")    
  }
  
  # 26 +- 7 covers most of the data:
  #qplot(out.all$nsig)
  #mean(out.all$nsig) # == 25.6
  #sd(out.all$nsig) * 2 # == 6.7
  #qplot(out$nsig)
  
  out.all
}

