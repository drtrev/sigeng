# Here I develop my own Cousineau-Morrey CIs and check they are same as Baguley's
# Along the way I discovered the purpose of difference (p. 166) and testing
# reveals that conf.level=.95 and difference=T gives error bars that can be used for inference.
# Also note that difference=T is making error bars smaller.

source("~/reps/sigeng/libs/cmci-funcs.r")

for (i in 1:1000)
{
  data <- matrix(c(rnorm(10), rnorm(10, 1)), nrow=10)
  data <- data.frame(data)
  
  # Give cm.ci wide matrix:
  #cm.ci(data)
  
  # Give my.cm.ci long data frame:
  data <- melt(data)
  data$id <- factor(1:10)
  (pval <- t.test(value ~ variable, data, paired=T)$p.value)
  #(pval <- t.test(value ~ variable, data, paired=F)$p.value)
  if (pval > .04 && pval < .06) break;
}
pval

# With difference = T and conf.level=.95 the error bars are borderline sig when they are borderline overlapping
# difference=T also means error bars are smaller
debugonce(args2symbols)
test <- my.cm.ci(data, id, value, variable, conf.level=.95, difference=T)
#test <- my.cm.ci(data, .(id), .(value), .(variable), conf.level=pnorm(1)-pnorm(-1), difference=F)

# Calculate a real 95% between subjects confidence interval based on the t-distribution (same as error.bars in psych package):
desc <- ddply(data, .(variable), summarize, mean=mean(value), bs.std.error=std.error(value))
n <- length(levels(data$id))
alpha <- .05
desc <- within(desc, bs.ci.half <- bs.std.error * qt(1-alpha/2, n-1))
desc$upper <- desc$mean + desc$bs.ci.half
desc$lower <- desc$mean - desc$bs.ci.half

args(my.bs.ci)
#debugonce(my.bs.ci)
bs.ci.diff <- my.bs.ci(data, id, value, variable, conf.level=.95, difference=T)
bs.ci.diff

wide <- dcast(data, id ~ variable)
row.names(wide) <- wide$id
wide <- wide[,2:ncol(wide)]

#bs.ci.diff <- data.frame(bs.ci(wide, difference=F)) # difference=T will make it interpretable
#bs.ci.diff$variable <- row.names(bs.ci.diff)

lm.ci.diff <- data.frame(lm.ci(wide, difference=T))
lm.ci.diff$variable <- row.names(lm.ci.diff)

debugonce(ml.ci)
ml.ci.test <- data.frame(ml.ci(wide))
ml.ci.test$variable <- c("X1", "X2")
ml.ci.test$CI <- ml.ci.test$CI.upper - ml.ci.test$CI.lower
ml.ci.test.comp.symm <- data.frame(ml.ci(wide, cov.matrix="comp.symm"))
ml.ci.test.comp.symm$variable <- c("X1", "X2")
ml.ci.test.comp.symm$CI <- ml.ci.test.comp.symm$CI.upper - ml.ci.test.comp.symm$CI.lower
ml.ci.test
ml.ci.test.comp.symm # same CIs as each other when using comp.symm (i.e. one CI size is returned)
# NOTE: the ml.ci() confidence intervals with unstructured covariance are the same as qt(.975, n-1) * std.error(value) for each variable,
# i.e. the same as desc. This is actually how they are calculated from the gmodels sub function in ml.ci(). It gets the std error
# from the model summary.

#t.test(wide$X1, wide$X2, paired=T)
limits <- c(-2, 2)
ggplot(test$CI, aes(ymin=lower, ymax=upper, x=variable)) + geom_errorbar(width=.1) + scale_y_continuous(limits=limits)
ggplot(desc, aes(ymin=mean-bs.ci.half, ymax=mean+bs.ci.half, x=variable)) + geom_errorbar(width=.1) + scale_y_continuous(limits=limits)
ggplot(bs.ci.diff$CI, aes(ymin=lower, ymax=upper, x=variable)) + geom_errorbar(width=.1) + scale_y_continuous(limits=limits)
ggplot(lm.ci.diff, aes(ymin=lower, ymax=upper, x=variable)) + geom_errorbar(width=.1) + scale_y_continuous(limits=limits)
ggplot(desc, aes(ymin=mean-bs.std.error, ymax=mean+bs.std.error, x=variable)) + geom_errorbar(width=.1) + scale_y_continuous(limits=limits)
ggplot(ml.ci.test, aes(ymin=CI.lower, ymax=CI.upper, x=variable)) + geom_errorbar(width=.1) + scale_y_continuous(limits=limits)
ggplot(ml.ci.test.comp.symm, aes(ymin=CI.lower, ymax=CI.upper, x=variable)) + geom_errorbar(width=.1) + scale_y_continuous(limits=limits)

?colMeans
colMeans(data[,c("DV1", "DV2")])


