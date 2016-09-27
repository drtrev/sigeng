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
  foreach(j=1:10, .combine=data.frame) %do%
  {
    source("analyse.r")
    load.packages()
    
    analyses <- initAnalyses()
    analyses <- analyses[1:2,]
    
    pvals <- try(sim(analyses))
    print(pvals) # <<---- pvals returns $dat and $analysis,
    # not sure where nsig comes in!!
    if (class(pvals)=="try-error")
    {
      cat("Caught error\n")
      pvals <- NA
    }
    pvals
  }
})

out
