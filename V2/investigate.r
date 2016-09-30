############################################################
# Investigate: Run analyses to answer questions
############################################################

initializeCluster <- function()
{
  library(foreach)
  library(doParallel)
  # cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf
  cluster <- makeCluster(2)
  registerDoParallel(cluster)
  # Check we have 2:
  stopifnot(getDoParWorkers()==2)
  
  cluster
}

killCluster <- function(cluster)
{
  # If you forget to call this, the cluster processes will terminate themselves when
  # the master R session's process ends.
  stopCluster(cl)
}

############################################################
# Helper functions
############################################################

getOutputFileName <- function(filePrefix=NULL, path="results")
{
  # Param:
  #   filePrefix   some kind of prefix name to give to the file,
  #                e.g. holdEachLevel
  #   path         path to results output folder

  if (is.null(filePrefix))
  {
    stop("Specify file prefix")
  }

  getResultsFileName <- function(path, filePrefix, fileNumber, fileSuffix)
  {
    paste0(path, "/", filePrefix, fileNumber, fileSuffix)
  }
  
  fileNumber <- 1
  fileSuffix <- ".RData"

  while (
    file.exists(
      getResultsFileName(path, filePrefix, fileNumber, fileSuffix)
    )
  )
  {
    fileNumber <- fileNumber + 1
  }
  cat("fileNumber:", fileNumber, "\n")
  
  getResultsFileName(path, filePrefix, fileNumber, fileSuffix)
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
  processingTime <- system.time(
    results <- foreach(colcurr=names(analyses), .combine=rbind) %do%
    {
      colcut <- analyses[,colcurr]
      if (length(levels(colcut))!=length(unique(colcut)))
      {
        # This should only happen if analyses was cut,
        # e.g. initAnalyses()[1:2,]
        warning("Factor levels have changed, refactoring")
        colcut <- factor(colcut)
      }
      
      pvalDataFrame <- foreach(l=levels(colcut), .combine=rbind) %do%
      {
        analyses.sub <- analyses[colcut==l,]
        
        # may need to use .packages or source/load.packages, seems to be needed after registerDoParallel
        pvalRow <- foreach(i=1:nreps, .combine=data.frame) %dopar%
        {
          source("generateData.r")
          source("analyse.r")
          load.packages()
          cat("nreps: ", nreps, "\n")
          # sim returns a single p value, the minimum of the two main effects and interaction.
          pval <- try(sim(analyses.sub))
          if (class(pval)=="try-error")
          {
            cat("Caught error\n")
            warning("Caught error")
            # Do not overwrite pval because it contains error
            #pval <- NA
          }
          pval
        }
        # pvalRow looks like this (combined into a data frame via foreach):
        #    result.1  result.2
        #  1   xxxxxx     xxxxx
        stopifnot(nrow(pvalRow)==1)
        # Add a column showing the total number of significant values from the row:
        pvalRow$nsig <- sum(pvalRow[1,] < 0.05)
        # Add the name of the column we are currently varying:
        pvalRow$colcut <- colcurr
        # Add the factor level of colcurr we just used:
        pvalRow$collevel <- l
        print(pvalRow)
        pvalRow
      }
      
      pvalDataFrame
    }
  ) # takes 2 hours
  
  # results looks like this:
  #         result.1     result.2 ... nsig             colcut           collevel
  #   1  0.064134139 0.0236179504 ...    1         diag.order normality-outliers
  #   2  0.284763115 0.0755375576 ...    0         diag.order outliers-normality
  #   3  0.106062614 0.1656168818 ...    0           outliers            boxplot
  #   4  ...
  metaInfo <- list(procesingTime=processingTime)
  list(results=results, analyses=analyses, metaInfo=metaInfo)
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
    
    # Get file name ready first in case there is a problem.
    outputFileName <- getOutputFileName("holdEachLevelList")

    holdEachLevelList <- holdEachLevel(analyses, nreps)

    save(holdEachLevelList, file=outputFileName)
  }
  else
  {
    # TODO use file number
    stop("remake==T not yet implemented")
    load("out-holdEachLevel.RData")
    load("out-holdEachLevelAnalyses.RData")
  }
  holdEachLevelList
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
      source("genderateData.r")
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
    stop("check this function works")
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

