############################################################
# Investigate: Run analyses to answer questions
############################################################

initializeCluster <- function(nWorkers=2)
{
  library(foreach)
  library(doParallel)
  # cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf
  cluster <- makeCluster(nWorkers)
  registerDoParallel(cluster)
  # Check we have nWorkers:
  stopifnot(getDoParWorkers()==nWorkers)
  
  cluster
}

killCluster <- function(cluster)
{
  # If you forget to call this, the cluster processes will terminate themselves when
  # the master R session's process ends.
  stopCluster(cluster)
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
# See also .export param for function names.
# See parallelDemo.r
holdEachLevel <- function(analyses, nRepetitions, nWorkers)
{
  # Here we run through all the combinations of analyses.
  # But we do so by "holding" one level at a time.
  # This allows us to see if, e.g., the diagnostic order makes a
  # difference in the probability of finding a significant result.
  processingTime <- system.time(
    results <- foreach(currentColumnName=names(analyses), .combine=rbind) %do%
    {
      currentColumn <- analyses[,currentColumnName]
      if (length(levels(currentColumn))!=length(unique(currentColumn)))
      {
        # This should only happen if analyses was cut,
        # e.g. initAnalyses()[1:2,]
        warning("Factor levels have changed, refactoring")
        currentColumn <- factor(currentColumn)
      }
      
      pvalDataFrame <- foreach(currentLevel=levels(currentColumn), .combine=rbind) %do%
      {
        analysesWithOneLevelHeld <- analyses[currentColumn==currentLevel,]
        
        # may need to use .packages or source/load.packages, seems to be needed after
        # registerDoParallel. An alternative is to use the .export param of foreach.
        # See parallelDemo.r
        pvalRow <- foreach(i=1:nRepetitions, .combine=data.frame) %dopar%
        {
          source("generateData.r")
          source("analyse.r")
          load.packages()

          # sim returns a single p value, the minimum of the two main effects and interaction.
          pval <- try(sim(analysesWithOneLevelHeld))
          if (class(pval)=="try-error")
          {
            cat("Caught error\n")
            warning("Caught error")
            # pval contains error info.
            # Either use this information, or overwrite with NA.
            #pval <- NA
          }
          pval
        }
        # pvalRow looks like this (combined into a data frame via foreach):
        #    result.1  result.2
        #  1   xxxxxx     xxxxx
        stopifnot(nrow(pvalRow)==1)
        # Add a column showing the total number of significant values from the row:
        pvalRow$nSig <- sum(pvalRow[1,] < 0.05)
        # Add the name of the column we are currently varying:
        pvalRow$columnName <- currentColumnName
        # Add the factor level of current column we just used:
        pvalRow$columnLevel <- currentLevel
        pvalRow
      }
      # pvalDataFrame rbinds all the levels of the analyses column we are currently working on.
      
      pvalDataFrame
    }
  ) # takes 2 hours

  # results is an rbind of all the pvalDataFrames, one for each column in analyses.
  
  # i.e. results looks like this:
  #         result.1     result.2 ... nSig      currentColumn       currentLevel
  #   1  0.064134139 0.0236179504 ...    1         diag.order normality-outliers
  #   2  0.284763115 0.0755375576 ...    0         diag.order outliers-normality
  #   3  0.106062614 0.1656168818 ...    0           outliers            boxplot
  #   4  ...
  interestingCols <- c("nSig", "columnName", "columnLevel")
  cat("----------\n")
  cat("nRepetitions:", nRepetitions, "\n")
  cat("ncol(results):", ncol(results), "\n")
  print(results[interestingCols])
  cat("----------\n")
  
  metaInfo <- list(procesingTime=processingTime, nWorkers=2)
  list(results=results, analyses=analyses, metaInfo=metaInfo)
}

############################################################

investigateHoldEachLevel <- function(loadFromCache=F, nreps=100, analyses=NULL, nWorkers=2)
{
  if (!loadFromCache)
  {
    if (is.null(analyses))
    {
      cat("Using default analyses\n")
      analyses <- initAnalyses()
    }
    
    # Get file name ready first in case there is a problem.
    outputFileName <- getOutputFileName("holdEachLevelList")

    cat("Initializing cluster with", nWorkers, "workers\n")
    cluster <- initializeCluster(nWorkers)

    # Perform analyses
    holdEachLevelList <- holdEachLevel(analyses, nreps, nWorkers)
    save(holdEachLevelList, file=outputFileName)

    # Clean up
    killCluster(cluster)
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

