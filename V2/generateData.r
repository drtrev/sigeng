generate.dat <- function(N=30)
  # Generate within subj's data without using idmeans.
  # Note that you cannot fail to have sphericity with only two levels of a factor.
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
  # within subjects data
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
