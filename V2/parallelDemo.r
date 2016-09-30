#####
# Testing foreach nesting operator and parallel
library(doParallel)
# windows doesn't have fork, so have to use snow "cluster" style?
# cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf
cl <- makeCluster(2)
registerDoParallel(cl)
getDoParWorkers()

# NOTE: can use .export to export global functions

# not parallel:
system.time(out <- foreach(j=1:100, .combine=cbind) %do%
{
  Sys.sleep(0.1)
})

# parallel:
system.time(out <- foreach(j=1:100, .combine=cbind) %dopar%
{
  Sys.sleep(0.1)
})

##
# not parallel:
system.time(out <- foreach(j=1:10, .combine=cbind) %do%
{
  foreach(j=1:10, .combine=cbind) %do%
{
  Sys.sleep(0.1)
}
})

# parallel:
system.time(out <- foreach(j=1:10, .combine=cbind) %do%
{
  foreach(j=1:10, .combine=cbind) %dopar%
{
  Sys.sleep(0.1)
}
})

# NOTE: when using foreach within foreach, we need to load the foreach package!
system.time(out <- foreach(j=1:10, .combine=cbind, .packages="foreach") %dopar%
{
  foreach(j=1:10, .combine=cbind) %do%
{
  Sys.sleep(0.1)
}
})

# Consider as a single stream of tasks, using the nesting operator,
# i.e. all items in the matrix are independent, and the two loops
# are turned into a "single loop" (see nested.pdf vignette)
system.time(out <- foreach(j=1:10, .combine=cbind, .packages="foreach") %:%
              foreach(j=1:10, .combine=cbind) %dopar%
{
  Sys.sleep(0.1)
}
)

# If the inner loop tasks are small, then you might want to execute the inner loop as a single task
# (to save time forking). In this case see nested.pdf vignette and look at doNWS (chunking).


#####
# trying to work out why dopar is not working in sigeng V2

source("analyse.r")
system.time(out <- foreach(j=1:10, .combine=cbind, .packages="foreach") %do% # dopar
{
  pvals <- foreach(j=1:10, .combine=data.frame) %dopar%
  {
    source("analyse.r")
    load.packages()
    
    analyses <- initAnalyses()
    analyses <- analyses[1:2,]
    
    pvals <- try(sim(analyses))
    cat("About to print pvals:----\n")
    print(pvals) 
    cat("----Printed pvals\n")
    
    if (class(pvals)=="try-error")
    {
      cat("Caught error\n")
      pvals <- NA
    }
    pvals
  }
  # returns:
  #   result.1  result.2 result.3 ...
  # 1 0.477461 0.1759253  0.17853 ...
  pvals
  pvals$nsig <- sum(pvals[1,] < 0.05)
  # Add the name of the column we are currently varying:
  #pvals$colcut <- colcurr
  # Add the factor level of colcurr we just used:
  #pvals$collevel <- l
  pvals
})

out

##
# It seems it does not work when it's in a function.
# Could be something to do with environments
test <- function()
{
  source("investigate.r")
  source("analyse.r")
  load.packages()
  a <- 1
  print(ls())

  pvalRow <- foreach(j=1:10, .combine=data.frame,
                     .export=c("initAnalyses", "load.packages", "sim", "generate.dat.within")) %dopar%
  {
    print(getwd())
    source("investigate.r")
    source("analyse.r")
    #source("generateData.r")
    load.packages()
    
    analyses <- initAnalyses()
    analyses <- analyses[1:2,]
    
    pval <- try(sim(analyses))
    cat("About to print pvals:----\n")
    print(pval)
    cat("----Printed pvals\n")
    
    if (class(pval)=="try-error")
    {
      cat("Caught error\n")
      # Don't overwrite pval because it contains the error
      # Problem was generateData.r was not loaded,
      # and it couldn't find generate.dat.within
      # Need to either source("generateData.r")
      # or add generate.dat.within to the .export param of foreach.
      #pval <- 5
    }
    pval
  }
  print("-------")
  print(pvalRow)
  print("-------")
  
}

a <- test()
?foreach

