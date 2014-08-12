generate.dat <- function(N=30)
  # Generate data with 2 factors.
  # N is number of subjects per group.
{
  # id is the subject id, we use N*4 because there are four groups.
  # group is the full name of the group, A1, A2, B1, B2
  # value is the measured/observed result in our experiment, again N*4.
  # the mean value is 0 for all groups, i.e. there are no effects of group.
  dat <- data.frame(id=factor(1:(N*4)),
                    group=factor(c(
                      rep(c("A1", "A2"), each=N),
                      rep(c("B1", "B2"), each=N))),
                    value=rnorm(N*4))
  
  # split group into separate factors
  dat$factor1 <- factor(substring(dat$group, 1, 1))
  dat$factor2 <- factor(substring(dat$group, 2, 2))
  
  dat
}

# e.g.
# dat <- generate.dat(3)
# print(dat)
#   id group       value factor1 factor2
# 1   1    A1  0.42385602       A       1
# 2   2    A1 -0.34829466       A       1
# 3   3    A1 -1.40946883       A       1
# 4   4    A2 -0.09177423       A       2
# 5   5    A2 -1.26614034       A       2
# 6   6    A2 -1.18024188       A       2
# 7   7    B1 -0.86129559       B       1
# 8   8    B1 -1.30517594       B       1
# 9   9    B1  0.50221849       B       1
# 10 10    B2  0.42755864       B       2
# 11 11    B2 -0.03262990       B       2
# 12 12    B2  0.85593000       B       2

do.anova <- function(dat)
  # Perform a 2 way between-subjects anova on dat.
  # Return the lowest p value obtained.
  # (We are interested whether there is at least one significant effect
  # in this particular data set).
{
  results <- aov(value ~ factor1 * factor2, data=dat)
  min( summary(results)[[1]]$`Pr(>F)`[1:3] )
}


# Peform 1000 replications:
# each time generate some data and get the lowest p value
# from the resulting anova.
out <- replicate(1000, do.anova(generate.dat()))

# How many times out of the 1000 was there a significant effect?
sum(out < .05) / 1000 # about 14%

# 1-(1-0.05)^3 = 0.14

dat <- generate.dat()
head(dat)

model1 <- lm(value ~ 1, dat)
model2 <- lm(value ~ 1 + factor1, dat)
model3 <- lm(value ~ 1 + factor1 + factor2, dat)
model4 <- lm(value ~ 1 + factor1 + factor2 + factor1:factor2, dat)
anova(model1, model2, model3, model4)

aov1 <- aov(value ~ factor1*factor2, dat)
summary(aov1)

anova(model1, model4)
summary.lm(aov1)

# same p value:
model1 <- lm(value ~ 1, dat)
model2 <- lm(value ~ 1 + group, dat)
anova(model1, model2)

# A1 A2 B1 B2
# dummy vars:
# int  x1    x2    x3  meaning:
# 1     0     0     0   A1
# 1     1     0     0   A2
# 1     0     1     0   B1
# 1     0     0     1   B2
# no interaction (treating as levels of same factor)
