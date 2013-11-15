# Here I develop my own Cousineau-Morrey CIs and check they are same as Baguley's
# Along the way I discovered the purpose of difference (p. 166) and testing
# reveals that conf.level=.95 and difference=T gives error bars that can be used for inference.
# Also note that difference=T is making error bars smaller.
# I am working on copying my function into sigeng/libs.r
# and when it is finalised into sigeng itself.
# So this should be used only for testing, and the latest version will be in version control in sigeng.

#for (i in 1:1000)
#{
data <- matrix(c(rnorm(10), rnorm(10, 1)), nrow=10)
data <- data.frame(data)

# Give cm.ci wide matrix:
#cm.ci(data)

# Give my.cm.ci long data frame:
data <- melt(data)
data$id <- factor(1:10)
t.test(value ~ variable, data, paired=T)
(pval <- t.test(value ~ variable, data, paired=F)$p.value)
#if (pval > .04 && pval < .06) break;
#}

# With difference = T and conf.level=.95 the error bars are borderline sig when they are borderline overlapping
# difference=T also means error bars are smaller
test <- my.cm.ci(data, .(id), .(value), .(variable), conf.level=.95, difference=T)
#test <- my.cm.ci(data, .(id), .(value), .(variable), conf.level=pnorm(1)-pnorm(-1), difference=F)

# Calculate a real 95% between subjects confidence interval based on the t-distribution (same as error.bars in psych package):
desc <- ddply(data, .(variable), summarize, mean=mean(value), bs.std.error=std.error(value))
n <- length(levels(data$id))
alpha <- .05
desc <- within(desc, bs.ci <- bs.std.error * qt(1-alpha/2, n-1))

args(my.bs.ci)
#debugonce(my.bs.ci)
bs.ci.diff <- my.bs.ci(data, value, variable, conf.level=.95, difference=T)
bs.ci.diff

#wide <- dcast(data, id ~ variable)
#row.names(wide) <- wide$id
#wide <- wide[,2:ncol(wide)]

#bs.ci.diff <- data.frame(bs.ci(wide, difference=F)) # difference=T will make it interpretable
#bs.ci.diff$variable <- row.names(bs.ci.diff)

lm.ci.diff <- data.frame(lm.ci(wide, difference=T))
lm.ci.diff$variable <- row.names(lm.ci.diff)

#t.test(wide$X1, wide$X2, paired=T)
limits <- c(-1, 2)
ggplot(test$CI, aes(ymin=lower, ymax=upper, x=variable)) + geom_errorbar(width=.1) + scale_y_continuous(limits=limits)
ggplot(desc, aes(ymin=mean-bs.ci, ymax=mean+bs.ci, x=variable)) + geom_errorbar(width=.1) + scale_y_continuous(limits=limits)
ggplot(bs.ci.diff$CI, aes(ymin=lower, ymax=upper, x=variable)) + geom_errorbar(width=.1) + scale_y_continuous(limits=limits)
ggplot(lm.ci.diff, aes(ymin=lower, ymax=upper, x=variable)) + geom_errorbar(width=.1) + scale_y_continuous(limits=limits)
ggplot(desc, aes(ymin=mean-bs.std.error, ymax=mean+bs.std.error, x=variable)) + geom_errorbar(width=.1) + scale_y_continuous(limits=limits)

?colMeans
colMeans(data[,c("DV1", "DV2")])


my.cm.ci <- function(data, id, dv, group, conf.level = 0.95, difference=F)
  # Give me a long format data frame
  # ezANOVA has code to convert symbols to .(symbol) but for now we require dv to be .(value)
  # Difference converts to a Goldstein-Healy type plot where you can judge significance (see p. 166 of Baguley)
{
  id <- as.character(id)
  dv <- as.character(dv)
  group <- as.character(group)
  
  # test
  #id <- "id"
  #dv <- "value"
  #group <- "variable"
  #difference <- T
  #

  diff.factor <- 1
  if (difference) diff.factor <- sqrt(2)/2
  
  wide <- dcast(data, id ~ variable)
  row.names(wide) <- wide$id
  wide <- wide[,2:ncol(wide)]
  
  bag <- cm.ci(wide, conf.level=conf.level, difference=difference) # check it's the same
  
  # From Baguley's paper: normalization:
  # Y'ij = Yij - Mi + Mgrand
  Mi <- rowMeans(wide)
  Mgrand <- mean(unlist(wide))
  norm <- wide - Mi + Mgrand
  norm
    
  # Interval:
  # Mj +- t(1-alpha/2, n-1) * sqrt(J/(J-1)) * std.error(norm[,j])
  alpha <- 1-conf.level
  n <- nrow(norm)
  J <- ncol(norm)
  qtval <- qt(1-alpha/2, n-1)
  correction <- sqrt(J/(J-1))
  # it applies to columns:
  halfCI <- qtval * correction * std.error(norm) * diff.factor
  upper <- colMeans(norm) + halfCI
  lower <- colMeans(norm) - halfCI
  # Baguley calculates std.error using mean.mat/t.mat which confused me!
  # but t = m/(s/sqrt(n))  and so m/t = m/(m/(s/sqrt(n))) = (m*s/sqrt(n))/m = s/sqrt(n)
  
  CI=data.frame(cbind(lower, mean=colMeans(norm), upper))
  CI[,group] <- row.names(CI)
  list(CI=CI, bag=bag)
}


# my.cm.ci <- function(data, id, dv, group, conf.level = 0.95)
#   # Give me a long format data frame
#   # ezANOVA has code to convert symbols to .(symbol) but for now we require dv to be .(value)
# {
#   id <- as.character(id)
#   dv <- as.character(dv)
#   group <- as.character(group)
#   
#   p.means <- ddply(data, .(id), function(x)
#     {
#       data.frame(p.mean=mean(x[,dv]))
#     }
#   )
#   p.means
#   
#   n <- length(levels(data[,id]))
#   k <- length(levels(data[,group]))
#   grandmean <- mean(data[,dv])
#   
#   # subtract p.means to normalize:
#   norm.df <- ddply(data, .(id, variable), function(x)
#     {
#       pm <- with(p.means, p.mean[id==x$id])
#       data.frame(norm.dv=x[,dv] - pm + grandmean)
#     }
#   )
#   norm.df
#   
# }

cm.ci <- function(data.frame, conf.level = 0.95, difference = TRUE) {
  # test:
  #data.frame <- data
  #conf.level=.95
  #difference=T
  #
    #cousineau-morey within-subject CIs
    k = ncol(data.frame)
    if (difference == TRUE) {
        diff.factor = 2^0.5/2
    }else diff.factor = 1
    n <- nrow(data.frame)
    df.stack <- stack(data.frame)
    index <- rep(1:n, k)
    p.means <- tapply(df.stack$values, index, mean)
    norm.df <- data.frame - p.means + (sum(data.frame)/(n * k))
    t.mat <- matrix(, k, 1)
    mean.mat <- matrix(, k, 1)
    for (i in 1:k) t.mat[i, ] <- t.test(norm.df[i])$statistic[1]
    for (i in 1:k) mean.mat[i, ] <- mean(norm.df[,i])
    c.factor <- (k/(k - 1))^0.5
    moe.mat <- mean.mat/t.mat * qt(1 - (1 - conf.level)/2, n - 1) * c.factor * 
        diff.factor
    ci.mat <- matrix(, k, 2)
    dimnames(ci.mat) <- list(names(data.frame), c("lower", "upper"))
    for (i in 1:k) {
        ci.mat[i, 1] <- mean.mat[i] - moe.mat[i]
        ci.mat[i, 2] <- mean.mat[i] + moe.mat[i]
    }
    ci.mat
}

bs.ci <- function(data.frame, conf.level = 0.95, difference = FALSE) {
  # between-subject CIs
  k = ncol(data.frame)
  n <- nrow(data.frame)
  df.stack <- stack(data.frame)
  df.stack$ind <- factor(df.stack$ind, levels=names(data.frame)) # Added by TD, keep factors in same order so can be labelled below
  group.means <- colMeans(data.frame, na.rm = TRUE)
  if (difference == TRUE) 
    ci.mat <- (confint(lm(values ~ 0 + ind, df.stack), level=conf.level) - group.means) * 
    2^0.5/2 + group.means
  else ci.mat <- confint(lm(values ~ 0 + ind, df.stack), level=conf.level)
  dimnames(ci.mat) <- list(names(data.frame), c("lower", "upper"))
  ci.mat
}

my.bs.ci <- function(data, value, group, conf.level = 0.95, difference=F)
  # Give me a long data frame.
  # This borrows some code (args2symbols) from ezANOVA to allow the function to be called like this:
  # my.bs.ci(data, value=response, group=variable) where response and variable are column names in data.
{
  args <- args2symbols(match.call(), c("value", "group"))
  
  # Can use like this: data[, names(data) == args$value] 
  # Or convert to factor:
  value <- as.character(args$value)
  group <- as.character(args$group)
  
  # Make wide to test with Baguley:
  wide <- dcast(data, id ~ variable)
  row.names(wide) <- wide$id
  wide <- wide[,2:ncol(wide)]

  bag <- bs.ci(wide, conf.level=conf.level, difference=difference)
  
  
  form <- as.formula(paste(value, "~ 0 +", group))
  ci <- as.data.frame(confint(lm(form, data), conf.level=conf.level))
                      
  if (difference)
  {
    group.means <- colMeans(wide)
    diff.factor <- sqrt(2)/2
    ci <- (ci - group.means) * diff.factor + group.means
  }

  names(ci) <- c("lower", "upper")
  ci[,group] <- row.names(ci)
  # Group names returned from confint are in the form variableXX
  # Rather than just overwriting them, I think this should be safer:
  ci[,group] <- substring(ci[,group], nchar(group)+1)
  list(CI=ci, bag=bag)
}

args2symbols <- function(args, args_to_check)
  # From ezANOVA():
  # Originally args = as.list(match.call()[-1]) from caller
  # Now caller can use match.call()
  # E.g. args <- args2symbols(match.call(), c("value", "group"))
{
  args <- as.list(args[-1])
  args.out <- list()
  for (i in 1:length(args)) {
    arg_name = names(args)[i]
    if (arg_name %in% args_to_check) {
      if (is.symbol(args[[i]])) {
        code = paste("args.out <- c(args.out, list(", arg_name, "=.(", as.character(args[[i]]), 
                     ")))", sep = "")
        eval(parse(text = code))
      }
      else {
        if (is.language(args[[i]])) {
          arg_vals = as.character(args[[i]])
          arg_vals = arg_vals[2:length(arg_vals)]
          arg_vals = paste(arg_vals, collapse = ",")
          code = paste("args.out <- c(args.out, list(", arg_name, "=.(", arg_vals, ")))", 
                       sep = "")
          eval(parse(text = code))
        }
      }
    }
  }
  
  args.out
}
