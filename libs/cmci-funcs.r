library(reshape2)
library(plyr)
library(plotrix)
library(ggplot2)

my.cm.ci <- function(data, id, dv, group, conf.level = 0.95, difference=F)
  # Give me a long format data frame
  # ezANOVA has code to convert symbols to .(symbol) but for now we require dv to be .(value)
  # Difference converts to a Goldstein-Healy type plot where you can judge significance (see p. 166 of Baguley)
{
  args <- args2symbols(match.call(), c("id", "dv", "group"))
  
  # Can use like this: data[, names(data) == args$id] 
  # Or convert to character:  
  id <- as.character(args$id)
  dv <- as.character(args$dv)
  group <- as.character(args$group)
  
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

lm.ci <- function(data.frame, conf.level = 0.95, difference = FALSE) {
  #loftus-masson within-subject CIs
  k = ncol(data.frame)
  n <- nrow(data.frame)
  df.stack <- stack(data.frame)
  require(nlme)
  parts <- rep(1:n, k)
  root.ms.error <- lme(values ~ 0 + ind, random = ~1 | parts, cbind(parts, 
                                                                    df.stack))[[6]]
  detach(package:nlme)
  mean.mat <- matrix(, k, 1)
  ci.mat <- matrix(, k, 2)
  if (difference == TRUE) 
    diff.factor = 2^0.5/2
  else diff.factor = 1
  moe <- root.ms.error/n^0.5 * qt(1 - (1 - conf.level)/2, (n - 1) * (k - 
                                                                       1)) * diff.factor
  for (i in 1:k) mean.mat[i, ] <- mean(data.frame[i])
  for (i in 1:k) {
    ci.mat[i, 1] <- mean.mat[i] - moe
    ci.mat[i, 2] <- mean.mat[i] + moe
  }
  dimnames(ci.mat) <- list(names(data.frame), c("lower", "upper"))
  ci.mat
}

ml.ci <- function(data.frame, conf.level = 0.95, cov.matrix = "unstructured") {
  # CI based on multilevel model with covariance matrix unstructured or
  #   constrained to compound symmetry
  k = ncol(data.frame)
  n.parts <- nrow(data.frame)
  data.long <- reshape(data.frame, idvar = "id", direction = "long", varying = 1:k, 
                       v.names = "dv")
  require(nlme)
  if (cov.matrix == "comp.symm") 
    ml.mod <- lme(dv ~ 0 + factor(time), random = ~1 | id, na.action = na.omit, 
                  data.long)
  else ml.mod <- lme(dv ~ 0 + factor(time), random = ~0 + factor(time) | 
                       id, na.action = na.omit, data.long)
  detach(package:nlme)
  require(gmodels)
  ci.mat <- ci(ml.mod, confidence = conf.level)[, 2:3]
  detach(package:gmodels)
  ci.mat
}

my.bs.ci <- function(data, value, group, conf.level = 0.95, difference=F)
  # Give me a long data frame.
  # This borrows some code (args2symbols) from ezANOVA to allow the function to be called like this:
  # my.bs.ci(data, value=response, group=variable) where response and variable are column names in data.
{
  args <- args2symbols(match.call(), c("value", "group"))
  
  # Can use like this: data[, names(data) == args$value] 
  # Or convert to character:
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
