# Put all my useful functions here
# see also the corresponding R files if you want to just run them and keep the vars

# source this file (I think) to use, i.e. source("U:/R/libs/libs.R")

# to restructure from columns to rows:
# e.g.
# participant  spk_static  spk_anim  list_static  list_anim
# ...
# to
# participant  cond  anim
# ...
# note cond corresponds here to column number, which may be different to original condition numbering!
# use:
# df2 <- reshape(df1, direction="long", varying = c("spk_static", "spk_anim", "list_static", "list_anim"), v.names=c("anim"), idvar="participant", times = c(1, 2, 3, 4), timevar = "cond")
# from https://stat.ethz.ch/pipermail/r-help/2007-July/136494.html

# to extract data from a dataframe:
# df1[df1$block == 2,]

# to cut off part of axis:
# p + geom_bar(position=dodge, stat="identity") + geom_errorbar(limits, position=dodge, width = 0.25) + coord_cartesian(ylim=c(12, 17)) + scale_y_continuous(breaks=c(13, 14, 15, 16, 17))
# 

library(ggplot2)
library(plotrix) # for std.error
library(foreign) # for read.spss
library(foreach)
library(reshape)

remove.factor <- function(x, convert="character") {
  command <- paste("y <- as.", convert, "(levels(x)[x])", sep="")
  eval(parse(text=command))
  y
}

remove.levels <- function(DF, var, remove)
# e.g. DF <- remove.level(DF, "id", c(1, 2, 3)) to remove participants 1, 2 and 3
# The point of this is it recreates the factor
{
  DF <- DF[!(DF[,var] %in% remove),]
  DF[,var] <- factor(DF[,var])
  DF
}

vsubset <- function(...)
{
  subset(...)[,1]
}

calcadj <- function(df1, id="id", variable="variable", value="value", variableID="variableID")
  # This doesn't restrict participant and condition to being numeric
  # WARNING: doesn't work with two groups, i.e. grandmean will take mean from everything,
  # and adjustment will adjust mean, if there are between subject elements to the design.
  # Adjust separately, e.g. calcadj(df1[study=="1",],...)
{
  if (!(id %in% colnames(df1))) {
    stop(paste("Calcadj: Error, id column '", id, "' not found", sep=""))
  }
  if (!(variable %in% colnames(df1))) {
    stop(paste("Calcadj: Error, variable column '", variable, "' not found", sep=""))
  }

  ids <- df1[,id]
  variables <- df1[,variable]
  if (variableID %in% colnames(df1)) {
    variableIDs <- df1[,variableID]
  }else{
    variableIDs <- NULL
    #warning("No variable IDs found")
  }

  # error check
  if (!is.factor(ids)) {
    stop("Calcadj error: id must be a factor")
    #return(NULL)
  }
  if (!is.factor(variables)) {
    stop("Calcadj error: variables must be a factor")
    #return(NULL)
  }

  idmeans <- foreach(i=1:length(levels(ids)), .combine = "data.frame") %do% {
    mean( df1[ids==levels(ids)[i],value] )
  }

  # foreach doesn't return a data frame if there's only one ID
  if (!is.data.frame(idmeans)) {
    idmeans <- data.frame(idmeans)
  }

  # foreach labels the columns result.1, result.2 etc.
  # note that these won't be in a sensible order, i.e. IDs 1 to 19
  # it depends on levels
  # so we label variables here:
  colnames(idmeans) <- levels(ids)

  grandmean <- mean( df1[,value] )

  adjustment <- grandmean - idmeans

  # go through each row, add corresponding adjustment
  for (i in 1:nrow(df1)) {
    # get id
    df1[i,value] <- df1[i,value] + adjustment[1,df1[i,id]]
  }


  ##### that's adjusted the means, part 2 puts it in a df for plotting

  # calc mean and stderr for each condition (variable)
  # could just use tapply here, e.g. varmeans <- tapply(values, variables, mean)
  # varstderrs <- tapply(values, variables, std.error)
  #varmeans <- foreach(i=1:length(levels(variables)), .combine = "c") %do% {
  #  mean( df1[ variables==levels(variables)[i], value ] )
  #}
  #varstderrs <- foreach(i=1:length(levels(variables)), .combine = "c") %do% {
  #  std.error( df1[ variables==levels(variables)[i], value ] )
  #}
  values <- df1[,value] # they have changed
  varmeans <- tapply(values, variables, mean)
  varstderrs <- tapply(values, variables, std.error)

  # put into data frame
  meansdf <- data.frame(means=varmeans, stderrs=varstderrs, variable=factor(levels(variables), levels=as.character(levels(variables))))
  rownames(meansdf) <- NULL

  if (!is.null(variableIDs)) {
    # Note combining factors with c() like this gives a numeric value
    # but only after combining!
    varIDs <- foreach(i=1:length(levels(variables)), .combine = "c") %do% {
      df1[ variables==levels(variables)[i], variableID ][1]
    }
    meansdf <- cbind(meansdf, varIDs)
    meansdf <- meansdf[order(meansdf$varIDs),]
    meansdf <- rename(meansdf, c(varIDs="variableID"))
    rownames(meansdf) <- NULL
  }

  # get number of factors, guess from first variable
  #nfacs <- length(strsplit(levels(df1$variable)[1], "\\.")[[1]])
  nfacs <- length(strsplit(levels(variables)[1], "\\.")[[1]])

  if (nfacs > 1) {
    names <- paste("factor", 1, sep="")
    for (i in 2:nfacs) {
      names <- append(names, paste("factor", i, sep=""))
    }
    cols <- colsplit(meansdf$variable, "\\.", names=names)
    meansdf <- cbind(meansdf, cols)
  }
  meansdf
}

calcadj2 <- function(df1)
{
  #data <- data.frame(participant=as.numeric(df1$partn), condition=as.numeric(df1$condn), value=as.numeric(df1$value))
  dataf <- data.frame(participant=as.numeric(df1$id), condition=as.numeric(df1$condn), value=as.numeric(df1$value))
  if (is.character(df1$id) || is.character(df1$condn)) {
    cat("Error: partn and condn must be numeric\n");
    return(NULL);
  }

  # calculate the adjusted mean and stderr
  # input: 'data'
  # main output: 'meansdf' (and adds 'valueadj' column to 'data', run calcadj.R for this)

  # expects you to have a data frame called data, with a column called participant, one called value, one called condition

  # NOTE PARTICIPANT MUST BE INT NOT STRING (i.e. NOT 3A etc.)
  # NOTE ALSO DATA MUST BE SORTED BY PARTICIPANT

  # e.g. data:
  # participant     condition     value
  # 1               1             9.8
  # 1               2             0.2
  # 1               3             15.2
  # 1               4             7.5
  # 2               1             2.3
  # 2               2             10.2

  # first get the range for participants and conditions

  pmin <- min(dataf$participant)
  pmax <- max(dataf$participant)

  cmin <- min(dataf$condition)
  cmax <- max(dataf$condition)

  # get grand mean

  grandmean <- mean(dataf$value)

  # get participant means

  pmeans <- matrix()

  for (p in pmin:pmax) {
    pmeans[p - pmin + 1] <- mean(dataf$value[dataf$participant == p])
  }

  # calculate adjustment

  adjustment <- grandmean - pmeans

  # for each condition, add the adjustment for that participant

  n <- 1

  for (p in pmin:pmax) {
    for (c in cmin:cmax) {

      dataf$valueadj[n] <- dataf$value[n] + adjustment[p - pmin + 1]
      n <- n + 1

    }
  }

  means <- matrix()
  stderrs <- matrix()

  for (c in cmin:cmax) {
    means[c - cmin + 1] <- mean(dataf[dataf$condition == c,]$valueadj)
    stderrs[c - cmin + 1] <- std.error(dataf[dataf$condition == c,]$valueadj)
  }

  means <- cbind(means, stderrs)
  means <- cbind(means, c(cmin:cmax))
  colnames(means)[3] <- "cond"

  meansdf <- as.data.frame(means)
  meansdf$cond <- factor(meansdf$cond)

  # make factor variables for plotting interactions
  # max 3 factors for now
  faclen <- c(0, 0, 0, 1)
  for (facn in 3:1) {
    command <- paste("faclen[facn] <- length(levels(df1$factor", facn, "))", sep="")
    eval(parse(text=command))

    if (faclen[facn] > 0) {
      eachn <- 1
      for (i in facn:3) {
        eachn <- eachn * faclen[i+1]
      }

      command <- paste("if (faclen[facn] > 0) meansdf$factor", facn, " <- factor(rep(1:faclen[facn], each=eachn))", sep="")
      eval(parse(text=command))
    }else faclen[facn] = 1;
  }

  return(meansdf)
}

# without adjustment
# note that you can't adjust twice (no difference), so probably don't need this
calcmeans2 <- function(df1)
{
  #data <- as.data.frame(cbind(df1$partn, df1$condn, df1$value))
  data <- as.data.frame(cbind(df1$id, df1$condn, df1$value))
  colnames(data) <- c("participant", "condition", "value")

  # calculate the adjusted mean and stderr
  # input: 'data'
  # main output: 'meansdf'

  # expects you to have a data frame called data, with a column called participant, one called value, one called condition

  # NOTE PARTICIPANT MUST BE INT NOT STRING (i.e. NOT 3A etc.)
  # NOTE ALSO DATA MUST BE SORTED BY PARTICIPANT

  # e.g. data:
  # participant     condition     value
  # 1               1             9.8
  # 1               2             0.2
  # 1               3             15.2
  # 1               4             7.5
  # 2               1             2.3
  # 2               2             10.2

  cmin <- min(data$condition)
  cmax <- max(data$condition)

  means <- matrix()
  stderrs <- matrix()

  for (c in cmin:cmax) {
    means[c - cmin + 1] <- mean(data[data$condition == c,]$value)
    stderrs[c - cmin + 1] <- std.error(data[data$condition == c,]$value)
  }

  means <- cbind(means, stderrs)
  means <- cbind(means, c(cmin:cmax))
  colnames(means)[3] <- "cond"

  meansdf <- as.data.frame(means)
  meansdf$cond <- factor(meansdf$cond)

  # make factor variables for plotting interactions
  # max 3 factors for now
  faclen <- c(0, 0, 0, 1)
  for (facn in 3:1) {
    command <- paste("faclen[facn] <- length(levels(df1$factor", facn, "))", sep="")
    eval(parse(text=command))

    if (faclen[facn] > 0) {
      eachn <- 1
      for (i in facn:3) {
        eachn <- eachn * faclen[i+1]
      }

      command <- paste("if (faclen[facn] > 0) meansdf$factor", facn, " <- factor(rep(1:faclen[facn], each=eachn))", sep="")
      eval(parse(text=command))
    }else faclen[facn] = 1;
  }

  return(meansdf)
}

plotmeans <- function(my.meansdf, x="variableID", y="means", fill="factor1", dodge.width=0.9, global=FALSE)
{
  # Define in the global environment for ggplot
  # The underscore is in case meansdf is already a var
  #meansdf__ <<- data.frame(Condition=as.factor(meansdf$condx_), Means=y_, StdError=err_)

  # requires library(ggplot2)

  # NOTE: can add '+ coord_flip()' to the end of the plot command!

  # Also can change colours with: meansdf__$group <- c(1, 2, 1, 2)
  # and change to fill=factor(meansdf__$group)

  #p <- ggplot(meansdf, aes(fill=Condition, y=Means, x = Condition))
  command <- paste("p <- ggplot(my.meansdf, aes(x=", x, ", y=", y, ", fill=", fill, "))", sep="")
  eval(parse(text=command))

  limits <- aes(ymax = means + stderrs, ymin = means - stderrs)

  # dodge is for when different rows have same value for x
  dodge <- position_dodge(width=dodge.width)
  if (global) {
    p <<- p
    my.limits <<- limits
    dodge <<- dodge
    # TODO think I need to also save out my.meansdf so p works -- check this
    cat("Saved p, mylimits and dodge\n")
    cat("Use: p + geom_bar(position=dodge, stat=\"identity\") + ",
        "geom_errorbar(my.limits, position=dodge, width=0.25)\n")
  }

  p + geom_bar(position=dodge, stat="identity") + geom_errorbar(limits, position=dodge, width=0.25)
}

# plot interactions
plotints <- function(meansdf, graphfactor=1, graphvalue=1, xfactor=2)
# graphfactor is the factor number to use for the graph
# graphvalue is which graph to draw, out of all the levels in that factor
# xfactor is the factor to use as the x axis. It works out the grouping by deduction.
{
  # find number of factors
  for (i in 1:10) {
    command <- paste("if (!is.null(meansdf$factor", i, ")) factors = i", sep="")
    eval(parse(text=command))
  }

  # right now just coded for factors == 2 or 3, when I see a pattern this can be neatened

  # this deduces group/colour factor, and if necessary cuts down data frame to graphvalues only

  if (factors == 3) {
    command <- paste("meansdf <- meansdf[meansdf$factor", graphfactor, "==graphvalue,]", sep="")
    eval(parse(text=command))

    # deduce group (colour) factor
    if (xfactor == 1) {

      if (graphfactor == 2) colourfactor <- 3
      else colourfactor <- 2

    }else if (xfactor == 2) {

      if (graphfactor == 1) colourfactor <- 3
      else colourfactor <- 1

    }else if (xfactor == 3){

      if (graphfactor == 1) colourfactor <- 2
      else colourfactor <- 1

    }else cat("Error with xfactor.\n")

  }else if (factors == 2) {

    if (xfactor == 1) colourfactor <- 2
    else { if (xfactor == 2) colourfactor <- 1 }

  }else cat("Not dealing with != 3 factors yet.\n")


  # now draw the graph

  limits <- aes(ymax = means + stderrs, ymin = means - stderrs)

  command <- paste("p <- ggplot(meansdf, aes(colour=factor", colourfactor, ", y=means, x=factor", xfactor, "))", sep="")
  eval(parse(text=command))

  command <- paste("p + geom_line(aes(group=factor", colourfactor, ")) + geom_errorbar(limits, width=0.1)", sep="")
  eval(parse(text=command))
}

reshapelong <- function(df1, id="id")
{
  # could use:
  # m <- melt(videodrift2)
  # m2 <- cbind(m, colsplit(m$variable, split="\\.", names=c("factor1", "factor2")))
  # from reshape package


  cols <- colnames(df1)[colnames(df1)!=id]
  #cols <- colnames(df1)
  #df2 <- reshape(df1, direction="long", varying=cols, timevar = "cond", v.names="value", times=cols, idvar = "partn")
  # make compatible with reshape package:
  df2 <- reshape(df1, direction="long", varying=cols, timevar = "variable", v.names="value", times=cols, idvar = "id")
  #df2 <- df2[order(df2$partn),]
  df2 <- df2[order(df2$id),]
  rownames(df2) <- NULL

  #df2$condn <- c(1:length(levels(factor(df2$cond))))
  df2$variableID <- c(1:length(levels(factor(df2$variable))))

  df2$variableID <- factor(df2$variableID)
  #df2$partn <- factor(df2$partn)
  df2$id <- factor(df2$id)

  # get levels in cond
  # this gets a matrix where rows are (conceptually) factors, and columns are levels within that factor
  #facs <- sapply(strsplit(df2$cond, "\\."), unlist)
  facs <- sapply(strsplit(df2$variable, "\\."), unlist)

  # maybe there's only one factor
  if (!is.null(nrow(facs))) {

    for (facn in 1:nrow(facs)) {
      # need to initialise this for grep assign to work
      command <- paste("df2$factor", facn, " <- 1", sep="")
      eval(parse(text=command))
      
      levs <- levels(factor(facs[facn,]))
      for (leveln in 1:length(levs)) {
        # factor could appear within another factor, so check for either
        # start of string or dot, followed by factor, followed by dot or end
        # of string
        regex <- paste("(^|\\.)", levs[leveln], "(\\.|$)", sep="")
        #command <- paste("df2[grep(levs[leveln],df2$cond),]$factor", facn, " <- levs[leveln]", sep="")
        command <- paste("df2[grep(regex,df2$variable),]$factor", facn, " <- levs[leveln]", sep="")
        #cat(command, "\n", sep="")
        eval(parse(text=command))
      }

      command <- paste("df2$factor", facn, " <- factor(df2$factor", facn, ")", sep="")
      eval(parse(text=command))
    }
  }
  
  # couldn't do this before because of strsplit
  df2$variable <- factor(df2$variable)

  return(df2)
}

# repeated measures aov for a dataframe generated like the one from reshapelong
# call with e.g. summary(aovrep(DF)) or summary(aovrep(reshapelong(df1)))
aovrep <- function(df1, valuename="value")
{
  if (!is.factor(df1$id)) {
    cat("Error: id must be a factor\n");
    return(NULL);
  }

  errorfound = F
  testfactor <- paste("if (is.factor(df1$", valuename, ")) { errorfound <- T; }", sep="")
  eval(parse(text=testfactor))
  if (errorfound) {
    cat("Error: value must not be a factor\n");
    return(NULL);
  }

  # see above
  #facs <- nrow(sapply(strsplit(as.character(df1$cond), "\\."), unlist))
  facs <- nrow(sapply(strsplit(as.character(df1$variable), "\\."), unlist))

  if (is.null(df1$factor1)) {
    # only one factor, just use cond var
    #factorvar <- "cond"
    factorvar <- "variable"
    if (!is.factor(df1$variable)) {
      cat("Error: variable must be a factor\n");
      return(NULL);
    }
  }else{
    factorvar <- "factor1"
    if (!is.factor(df1$factor1)) {
      cat("Error: factor1 must be a factor\n");
      return(NULL);
    }
  }

  # aov(valuename ~ (factor1
  command <- paste("aov(", valuename, " ~ (", factorvar, sep="")

  if (!is.null(facs) && facs > 1) {
    # * factor2...
    for (i in 2:facs)  command <- paste(command, " * factor", i, sep="")
    errorfound <- F
    testfactor <- paste("if (!is.factor(df1$factor", i, ")) { errorfound <- T; cat(\"Error: factor", i, " must be a factor\\n\"); }", sep="")
    eval(parse(text=testfactor))
    if (errorfound) { return(NULL); }
  }

  # ) + Error(id/(factor1
  command <- paste(command, ") + Error(id/(", factorvar, sep="")

  if (!is.null(facs) && facs > 1) {
    # * factor2...
    for (i in 2:facs)  command <- paste(command, " * factor", i, sep="")
  }

  # )), df1)
  command <- paste(command, ")), df1)", sep="")

  cat(command, "\n", sep="")
  # last statement is returned in absence of a return statement
  eval(parse(text=command))
}


compquestion <- function(df1, df2)
  # df <- compquestion(tuq1, epflq1) ...
{
  m <- melt(df1)
  m2 <- melt(df2) # or melt(data.frame(df2))

  m3 <- rbind(m, m2)
  m4 <- cbind(m3, colsplit(m3$variable, split="\\.", names=c("factor1", "factor2", "factor3")))
  m4$id <- factor(m4$id)
  meansdf <- calcadj(m4)
  meansdf$variableID <- factor(1:nrow(meansdf))
  meansdf$order <- c(1, 2, 3, 4, 1, 2, 3, 4)
  meansdf$ordersep <- c(1, 3, 5, 7, 2, 4, 6, 8)
  list(m4, meansdf)
}

makemcnemar <- function(df1, variable="variable", value="value", firstcond=levels(df1[,variable])[1], secondcond=levels(df1[,variable])[2], responses=levels(df1[,value]))
# Take data frame structured in long format (gets factor levels properly)
# Take two columns, cond 1 and cond 2
# put into matrix for mcnemar test
# firstcond is column number of first condition
# secondcond is column number of second condition
# nresponses is number of possible responeses for the value itself,
# so 'left', 'both', 'right' would be 3
{

  # make mapping for responses, this way each can correspond to a column or row in the matrix

  mappings <- foreach(i=1:length(responses), .combine = "data.frame") %do% {
    i
  }
  colnames(mappings) <- responses

  # reshape data frame for lookups (just found this easier
  mydf <- cast(df1, id ~ variable, value=value)

  mymatrix <- matrix(0, nrow=length(responses), ncol=length(responses))

  for (r in 1:nrow(mydf)) {
    firstcondval <- mappings[1,as.character(mydf[r,firstcond])]
    secondcondval <- mappings[1,as.character(mydf[r,secondcond])]

    mymatrix[firstcondval,secondcondval] <- mymatrix[firstcondval,secondcondval] + 1

    #if (r == 11) {
    #  cat(paste("Firstcondval: ", firstcondval, ", secondcondval: ", secondcondval, "\n", sep=""))
    #  cat(paste("secondcond: ", secondcond, ", mydf[11,secondcond]: ", mydf[11,secondcond], "\n", sep=""))
    #  cat(paste("mapping: ", mappings[1,"right"], "\n", sep=""))
    #}
     
  }

  rownames(mymatrix) <- responses
  colnames(mymatrix) <- responses

  cat("Mydf:\n")
  print(mydf)
  #cat("Mappings:\n")
  #print(mappings)
  cat(paste("First cond (rows): ", firstcond, ", second cond (cols): ", secondcond, "\n", sep=""))
  mymatrix

}

