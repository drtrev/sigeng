# v0.2
# Copied back to MPI

# v0.1
# Copied from MPI to Aurelie's laptop

#source("../libs/libs.R")
library(ggplot2)
#library(hash)

checkForBetweenInteraction <- function(mysum, IDs, DVs, withinIVs, betweenIVs)
{
  # WORK IN PROGRESS

  interaction = F

  # when there are not betweens, length is two, i.e. Pr for main effect, and NULL for residual SS
  if (length(mysum[[1]][[1]]$Pr) > 2) {
    # we have betweens

    for (i in 1:length(mysum)) {
      startloop <- 1

      # special case for main effect, e.g. main effect of order
      if (i == 1 && mysum[[i]][[1]]$Pr[1] < 0.05) {
        print("+++Main effect of between IV+++")
        startloop <- 2
      }

      # Go through each P value
      for (j in startloop:(length(mysum[[i]][[1]]$Pr) - 1)) {
        print(paste("i is ", i, ", j is", j))
        if (mysum[[i]][[1]]$Pr[j] < 0.05) {
          print("+++Interaction effect for between IV+++")
          interaction=T
        }
      }
    }

  }

  interaction
}

calcSummaryDF <- function(DF, variable, value)
{
  summaries <- tapply(DF[,value], DF[,variable], summary)
  summaryDF <- data.frame(do.call(rbind, summaries))

  command <- paste("summaryDF <- cbind(", variable, "=rownames(summaryDF), summaryDF)", sep="")
  eval(parse(text=command))
  rownames(summaryDF) <- NULL

  summaryDF
}

calcSignifDF <- function(DF, sig)
# Take a list of significant data, e.g.
#sig <- list(
#     list(
#       cut    = list(name="quest", value="Q5"),
#       first  = list(body.sync="Body.Synchronous"),
#       second = list(body.sync="Body.Asynchronous"),
#       y      = 7.5,
#       height = 0.1,
#       text   = "***"
#     )
#)
# Return signifDF, i.e. x, y, xend, yend, quest=...
{
  # For each significant thing
  for (s in sig) {

    xvariable <- names(s[["first"]])[1] # e.g. body.sync

    ypos <- 7.5 # Default
    if (!is.null(s[["y"]])) ypos <- s[["y"]]

    height <- 0.1 # Default
    if (!is.null(s[["height"]])) height <- s[["height"]]

    if (!is.null(s[["cut"]])) {
      # Cut out relevant part of DF
      cutname  <- s[["cut"]][["name"]]
      cutvalue <- s[["cut"]][["value"]]
      DFcut <- DF[DF[,cutname] == cutvalue,]
    }else{
      cutname <- NULL
      cutvalue <- NULL
      DFcut <- DF
    }

    # Get start point, i.e. x pos where var is first value
    firstvalue  <- s[["first"]][[xvariable]]
    secondvalue <- s[["second"]][[xvariable]]

    # Find where this is in DF
    x <- which(levels(DFcut[,xvariable]) == firstvalue)
    xend <- which(levels(DFcut[,xvariable]) == secondvalue)

    # Make a signifDF based on this
    #  2>-------3
    #  |        V
    #  |        |
    #  ^
    #  1

    #          1  2  3                    1            2     3
    xs    <- c(x, x, xend);    ys    <- c(ypos-height, ypos, ypos)
    xends <- c(x, xend, xend); yends <- c(ypos,        ypos, ypos-height)
    

    # TODO append don't overwrite
    signifDF <- data.frame(x=xs, y=ys, xend=xends, yend=yends, tempvar=rep(cutvalue, 3))
    signifDF <- rename(signifDF, c(tempvar=cutname))

    # Store text position, i.e. x, y, cutname, label
    text <- "****"
    if (!is.null(s[["text"]])) text <- s[["text"]]
    textDF <- data.frame(x=mean(c(x, xend)), y=ypos + 0.1, tempvar=cutvalue, label=text)
    textDF <- rename(textDF, c(tempvar=cutname))
  }

  signifDFs <- list(signifDF=signifDF, textDF=textDF)

  return(signifDFs)
}

plotbar <- function(meansdf, x="variable", dvlabel, fill, filebase)
{
  command <- paste("p <- ggplot(meansdf, aes(x=", x, ", y=means", sep="")

  if (is.null(fill)) command <- paste(command, "))", sep="")
  else command <- paste(command, ", fill=", fill, "))", sep="")

  print(command)
  eval(parse(text=command))

  limits <- aes(ymin=means-stderrs, ymax=means+stderrs)
  dodge <- position_dodge(width=0.1)
  myplot <- p + geom_bar(position=dodge) + geom_errorbar(limits, width=0.1, position=dodge) + scale_y_continuous(name=dvlabel)

  print(myplot)
  #ggsave(file=paste("/tmp/", filebase, ".eps", sep=""), width=10, height=7)
}

plotit <- function(DF, meansDF, summaryDF, fit, plothash, nosave=F) #x="variable", y, fill, rows=NULL, xlabel, filebase)
{
  #print(meansDF)

  ####
  # Prepare summary DFs and signifDF

  if (length(plothash[["combine"]]) > 1) {
    # make a new variable, which is named by joining all parts of "combine", e.g. "body.sync"
    newvar <- paste(plothash[["combine"]], collapse=".")

    # Set the value to be the combination of other variables, the value of "combine", e.g. c("body", "sync")
    # would make Body.Synchronous, Body.Asynchronous etc.
    summaryDF[,newvar] <- do.call(paste, c(summaryDF[ , plothash[["combine"]] ], sep=".") )
    # Don't allow for alphabetical order, take levels from the order they appear in DF
    summaryDF[,newvar] <- factor(summaryDF[,newvar], levels=unique(summaryDF[,newvar]))

    # Same again for main DF, because that's also used for some plots
    DF[,newvar] <- do.call(paste, c(DF[ , plothash[["combine"]] ], sep=".") )
    # Don't allow for alphabetical order, take levels from the order they appear in DF
    DF[,newvar] <- factor(DF[,newvar], levels=unique(DF[,newvar]))

    # Same again for meansDF, because that's also used for some plots
    meansDF[,newvar] <- do.call(paste, c(meansDF[ , plothash[["combine"]] ], sep=".") )
    # Don't allow for alphabetical order, take levels from the order they appear in DF
    meansDF[,newvar] <- factor(meansDF[,newvar], levels=unique(meansDF[,newvar]))
  }

  signifDFs <- NULL
  if (!is.null(plothash[["sig"]])) signifDFs <- calcSignifDF(summaryDF, plothash[["sig"]])

  #
  ####

  ####
  # Start the plot

  cat("Plot:\n")

  # TODO tidy this up, most of these options are just ggplot() so can add geom_segment(data=signifDF...) later

  if (plothash[["type"]] == "bar") {

    # TODO leave ggplot() blank and do it all below so can add geom_segment(data=signifDF)
    # TODO note I've already changed fill below to not add if we're using bar here, so should
    # fill even lower down
    command <- paste("p <- ggplot()")

  }else if (plothash[["type"]] == "psycho") {

    means <- "means"
    stderrs <- "stderrs"
    if (!is.null(plothash[["means"]])) means <- plothash[["means"]]
    if (!is.null(plothash[["stderrs"]])) stderrs <- plothash[["stderrs"]]
    command <- paste("p <- ggplot() + geom_point(data=meansDF, aes(x=", plothash[["x"]], ", y=", means, sep="")

    if (is.null(plothash[["fill"]])) command <- paste(command, "))", sep="")
    else command <- paste(command, ", fill=", plothash[["fill"]], "))", sep="")

    if (!is.null(plothash[["colour"]])) command <- paste(command, ", colour=", plothash[["colour"]], sep="")

  }else if (plothash[["type"]] == "crossbar") {
    
    command <- paste("p <- ggplot()")

  }else{
    # Default: boxplot
    plothash[["type"]] <- "boxplot"
    command <- paste("p <- ggplot()")

  }

  cat(command, "\n")
  eval(parse(text=command))

  #
  ####


  ####
  # Make facet, themes and opts vars etc.

  if (!is.null(plothash[["facet"]])) facet <- paste(" + facet_", plothash[["facet"]], sep="")
  else facet <- NULL
  #else facet <- paste("facet_grid(. ~ ", x, ")", sep="")
  if (!is.null(plothash[["theme"]])) theme <- paste(" + theme_", plothash[["theme"]], sep="")
  else theme <- NULL

  if (!is.null(plothash[["scale_fill"]])) scale_fill <- paste(" + scale_fill_", plothash[["scale_fill"]], sep="")
  else scale_fill <- NULL

  # TODO opts deprecated
  myopts <- NULL
  if (!is.null(plothash[["opts"]])) myopts <- paste(" + opts(", plothash[["opts"]], ")", sep="")

  signif <- NULL
  #if (!is.null(plothash[["signif"]])) signif <- paste(" + geom_segment(", plothash[["signif"]], ")", sep="")
  if (!is.null(signifDFs)) signif <- " + geom_segment(data=signifDFs$signifDF, aes(x=x, y=y, xend=xend, yend=yend))"

  text <- NULL
  #if (!is.null(plothash[["text"]])) text <- paste(" + geom_text(", plothash[["text"]], ")", sep="")
  if (!is.null(signifDFs)) text <- " + geom_text(data=signifDFs$textDF, mapping=aes(x=x, y=y, label=label))"

  #
  ####


  ####
  # Make plot specific command, e.g. geom_boxplot

  errorbar <- NULL
  legend <- NULL

  if (plothash[["type"]] == "bar") {
    means <- "means"
    stderrs <- "stderrs"
    if (!is.null(plothash[["means"]])) means <- plothash[["means"]]
    if (!is.null(plothash[["stderrs"]])) stderrs <- plothash[["stderrs"]]
    aescommand <- paste("aes(x=", plothash[["x"]], ", y=", means, sep="")
    if (!is.null(plothash[["fill"]])) aescommand <- paste(aescommand, ", fill=", plothash[["fill"]], sep="")
    aescommand <- paste(aescommand, ")", sep="")

    dodge <- position_dodge(width=0.9)
    plottype <- paste("geom_bar(data=meansDF, mapping=", aescommand, ", position=dodge, stat=\"identity\")", sep="")
    #if (!is.null(plothash[["stderrs"]])) {
      command <- paste("limits <- aes(x=", plothash[["x"]], ", y=", means, ", ymin=", means, "-", stderrs, ", ymax=", means, "+", stderrs, ")", sep="")
      cat(command, "\n")
      eval(parse(text=command))
      errorbar <- " + geom_errorbar(data=meansDF, mapping=limits, width=0.1, position=dodge)"
    #}
  }else if (plothash[["type"]] == "point") {
    plottype <- "geom_point()"
  }else if (plothash[["type"]] == "psycho") {

    # the ../psychophysics/psychophysics.r function should give a fit with x and y
    # and another function should join together different conditions with a variable 'cond'
    # so we have already x, y, cond
    plottype <- "geom_line(data=fit, aes(x=x, y=y, group=cond, colour=cond))"

  }else if (plothash[["type"]] == "crossbar") {

    median <- "Median"
    X1stQu <- "X1st.Qu."
    X3rdQu <- "X3rd.Qu."
    if (!is.null(plothash[["median"]])) median <- plothash[["median"]]

    plottype <- "geom_crossbar("

    crossbar <- paste("data=summaryDF, mapping=aes(x=", plothash[["x"]], ", y=", median, ", ymin=", X1stQu, ", ymax=", X3rdQu, sep="")

    if (!is.null(plothash[["fill"]])) crossbar <- paste(crossbar, ", fill=", plothash[["fill"]], sep="")

    crossbar <- paste(crossbar, ")", sep="") # end aes

    if (!is.null(plothash[["width"]])) crossbar <- paste(crossbar, ", width=", plothash[["width"]], sep="")

    plottype <- paste(plottype, crossbar, ")", sep="") # end geom_crossbar

    # You can use this to draw a geom_boxplot() after the plot to replace the crossbar legend (not always useful) with a boxplot one
    if (!is.null(plothash[["legend"]])) {
      if (plothash[["legend"]] == "boxplot") legend <- paste(" + geom_boxplot(", crossbar, ")", sep="")
      if (plothash[["legend"]] == "bar") {
        legend <- " + geom_bar(data=summaryDF, mapping=aes(x=0, y=0"
        if (!is.null(plothash[["fill"]])) legend <- sprintf("%s, fill=%s", legend, plothash[["fill"]])
        legend <- paste(legend, "), stat=\"identity\", width=0)", sep="")
      }
      if (plothash[["legend"]] == "boxplot") legend <- paste(" + geom_", plothash[["legend"]], "(", crossbar, ")", sep="")
    }
  }else{
    plottype <- paste("geom_boxplot(data=DF, mapping=aes(x=", plothash[["x"]], ", y=", plothash[["y"]], sep="")
    if (!is.null(plothash[["fill"]])) plottype <- paste(plottype, ", fill=", plothash[["fill"]], sep="")
    plottype <- paste(plottype, ")", sep="") # end aes
    if (!is.null(plothash[["width"]])) plottype <- paste(plottype, ", width=", plothash[["width"]], sep="")
    plottype <- paste(plottype, ")", sep="") # end geom_boxplot
  }
  #if (!is.null(plothash[["xticlabs"]]) xticlabs <- plothash[["xticlabs"]]

  #
  ####



  ####
  # Join plot commands together, and adjust the scale
  
  command <- paste("myplot <- p + ", plottype, errorbar, facet, sep="")

  if (plothash[["type"]] == "psycho") {
    command <- paste(command, " + scale_x_continuous(name=\"", plothash[["xlab"]],
                              "\", labels=c(\"", paste(plothash[["xticlabs"]], collapse='","', sep=""), "\"),", sep="")
    # without quote marks around breaks, i.e. continuous var not factor
    command <- paste(command, " breaks=c(", paste(plothash[["xticbreaks"]], collapse=',', sep=""), "))", sep="")
  }else{
    command <- paste(command, " + scale_x_discrete(name=\"", plothash[["xlab"]],
                              "\", labels=c(\"", paste(plothash[["xticlabs"]], collapse='","', sep=""), "\"),", sep="")
    command <- paste(command, " breaks=c(\"", paste(plothash[["xticbreaks"]], collapse='","', sep=""), "\"))", sep="")
  }


  command <- paste(command, " + scale_y_continuous(name=\"", plothash[["ylab"]], "\")", sep="")

  if (is.null(plothash[["scale_fill"]])) {
    # Filllab is unique in the sense that it also gets its value from fill
    command <- paste(command, " + scale_fill_discrete(name=\"", plothash[["filllab"]], "\"", sep="")
    if (!is.null(plothash[["filllabs"]])) {
      command <- paste(command, ", labels=c(\"", paste(plothash[["filllabs"]], collapse='","', sep=""), "\")", sep="")
      command <- paste(command, ", breaks=c(\"", paste(plothash[["fillbreaks"]], collapse='","', sep=""), "\")", sep="")
    }
    command <- paste(command, ")", sep="")
  }

  command <- paste(command, scale_fill, signif, text, theme, myopts, legend, sep="")
  cat(command, "\n")
  eval(parse(text=command))

  #
  ####


  print(myplot)

  if (!nosave) {
    for (f in plothash[["filenames"]]) {
      cat("Saving to file \"", f, "\"\n", sep="")
      ggsave(file=f, width=10, height=6)
    }
    cat("\n")
  }

  # DFs are sometimes changed, e.g. combine
  outputDFs <- list(summaryDF=summaryDF, meansDF=meansDF)
  return(outputDFs)
}

plotwivbiv <- function(meansdf, dvlabel)
{
  # TODO use parameters from plothashes or just call ploteng
  if (is.null(meansdf$biv)) {
    p <- ggplot(meansdf, aes(x=wiv, y=means))
  }else{
    p <- ggplot(meansdf, aes(x=wiv, y=means, group=biv, colour=biv))
  }
  limits <- aes(ymin=means-stderrs, ymax=means+stderrs)
  dodge <- position_dodge(width=0.1)
  # TODO is it right to use limits in geom_line?
  myplot <- p + geom_line(limits, position=dodge) + geom_errorbar(limits, width=0.1, position=dodge) + scale_y_continuous(name=dvlabel)

  print(myplot)

#  # this will look for withinIVs in DF, so need to expand the vars
#  command <- paste("p <- ggplot(DF, aes(x=", iv, sep="")
#  #for (iv in withinIVs) {
#  #  if (iv == withinIVs[1]) { command <- paste(command, iv, sep="") }
#  #  else { command <- paste(command, ", ", iv, sep="") }
#  #}
#  command <- paste(command, ", y=", dv, ", group=", biv, ", colour=", biv, "))", sep="")
#
#  print(command)
#  eval(parse(text=command))
#
#  myplot <- p + geom_line()
#
#  print("About to print plot")
#  print(myplot)
#  print("Printed")
}

ploteng <- function(DF, meansDF, summaryDF, fit=NULL, plothashes, nosave=F)
{
  for (i in plothashes) {
    #if (i[["type"]] == "box")
    newsummaryDFs <- plotit(DF, meansDF, summaryDF, fit, i, nosave)
  }

  # For now return the last one
  return (newsummaryDFs)
}

makeplothashes <- function(DF, dv, withinIVs, betweenIVs)
  # make a default plothash
{
  filenames = c("/tmp/signengplot.png") # TODO increment filename

  plothashes <- NULL

  if (length(betweenIVs) == 1 && length(withinIVs) == 3) {
    # facet factor1, x factor2, fill factor3
    plothashes <- list(list(xlab=withinIVs[2], x=withinIVs[2], xticlabs=levels(DF[,withinIVs[2]]), xticbreaks=levels(DF[,withinIVs[2]]),
          ylab=dv, y=dv, fill=withinIVs[3], filllab=withinIVs[3], scale_fill=NULL, signif=NULL, text=NULL, width=NULL,
          facet=paste("grid(", betweenIVs[1], " ~ ", withinIVs[1], ")", sep=""), type="bar",
          opts="axis.title.x=theme_text(size=12, vjust=0), axis.title.y=theme_text(size=12, vjust=0.4, angle=90)", #, panel.margin=unit(0.5, \"cm\")",
          filenames=filenames))
  }

  if (length(betweenIVs) == 1 && length(withinIVs) == 2) {
    # facet factor1, x factor2, fill factor3
    plothashes <- list(list(xlab=withinIVs[1], x=withinIVs[1], xticlabs=levels(DF[,withinIVs[1]]), xticbreaks=levels(DF[,withinIVs[1]]),
          ylab=dv, y=dv, fill=withinIVs[2], filllab=withinIVs[2], filllabs=levels(DF[,withinIVs[2]]), fillbreaks=levels(DF[,withinIVs[2]]),
          facet=paste("grid(. ~ ", betweenIVs[1], ")", sep=""), type="bar", scale_fill=NULL, signif=NULL, text=NULL, width=NULL,
          opts="axis.title.x=theme_text(size=12, vjust=0), axis.title.y=theme_text(size=12, vjust=0.4, angle=90)", #, panel.margin=unit(0.5, \"cm\")",
          filenames=filenames))
  }

  if (length(betweenIVs) == 0 && length(withinIVs) == 3) {
    # facet factor1, x factor2, fill factor3
    plothashes <- list(list(xlab=withinIVs[2], x=withinIVs[2], xticlabs=levels(DF[,withinIVs[2]]), xticbreaks=levels(DF[,withinIVs[2]]),
          ylab=dv, y=dv, fill=withinIVs[3], filllab=withinIVs[3],
          facet=paste("grid(. ~ ", withinIVs[1], ")", sep=""), type="bar", scale_fill=NULL, signif=NULL, text=NULL, width=NULL,
          opts="axis.title.x=theme_text(size=12, vjust=0), axis.title.y=theme_text(size=12, vjust=0.4, angle=90)", #, panel.margin=unit(0.5, \"cm\")",
          filenames=filenames))
  }

  if (length(betweenIVs) == 0 && length(withinIVs) == 2) {
    # facet factor1, x factor2, fill factor3
    plothashes <- list(list(xlab=withinIVs[1], x=withinIVs[1], xticlabs=levels(DF[,withinIVs[1]]), xticbreaks=levels(DF[,withinIVs[1]]),
          ylab=dv, y=dv, fill=withinIVs[2], filllab=withinIVs[2], scale_fill=NULL, signif=NULL, text=NULL, width=NULL,
          facet=NULL, type="bar",
          opts=NULL,
          filenames=filenames))
  }
  if (length(betweenIVs) == 0 && length(withinIVs) == 1) {
    # x factor1
    plothashes <- list(list(xlab=withinIVs[1], x=withinIVs[1], xticlabs=levels(DF[,withinIVs[1]]), xticbreaks=levels(DF[,withinIVs[1]]),
          ylab=dv, y=dv, fill=NULL, filllab=NULL, scale_fill=NULL, signif=NULL, text=NULL, width=NULL,
          facet=NULL, type="bar",
          opts=NULL,
          filenames=filenames))
  }
  if (length(betweenIVs) == 1 && length(withinIVs) == 0) {
    # x factor1
    plothashes <- list(list(xlab=betweenIVs[1], x=betweenIVs[1], xticlabs=levels(DF[,betweenIVs[1]]), xticbreaks=levels(DF[,betweenIVs[1]]),
          ylab=dv, y=dv, fill=NULL, filllab=NULL, scale_fill=NULL, signif=NULL, text=NULL, width=NULL,
          facet=NULL, type="bar",
          opts=NULL,
          filenames=filenames))
  }

  if (is.null(plothashes)) cat("NO DEFAULT PLOT FOUND\n\n")

  return(plothashes)
}

plotWithinMult <- function(dv, params, plotPsycho)
# plot for when there are > 1 within IVs
{
  DF <- params$DF
  IDs <- params$IDs
  withinIVs <- params$withinIVs
  betweenIVs <- params$betweenIVs
  plothashes <- params$plothashes

  # concatenate vars for calcadj
  DF$conc <- getConc(DF, withinIVs)

  meansdf <- calcadj(DF, variable="conc", value=dv)

  summarydf <- calcSummaryDF(DF, variable="conc", value=dv)
  # split out our conc variables into separate ones for plotting
  summarydf <- cbind(summarydf, colsplit(summarydf$conc, "\\.", names=withinIVs))
  # make the new column factor levels in the same order as they appear in the DF
  for (i in withinIVs) { summarydf[,i] <- factor(summarydf[,i], levels=unique(summarydf[,i])) }

  cat("Summary DF\n")
  print(summarydf)

  #cat("withinIVs")
  #print(withinIVs)
  #print(DF$conc)
  #cat("Means DF before col rename\n")
  #print(meansdf)

  # means, stderrs, variable, ...
  colnames(meansdf)[3] <- "condition"
  colnames(meansdf)[4:ncol(meansdf)] <- withinIVs

  #cat("Means DF\n")
  #print(meansdf)

  # get factor levels in order they appear, not alpha order (default)
  for (idx in 4:ncol(meansdf)) {
    # if they are not strings then calcadj can leave them as not factors
    if (is.factor(meansdf[,idx])) {
      temp <- remove.factor(meansdf[,idx])
      meansdf[,idx] <- factor(meansdf[,idx], levels=unique(temp))
    }else{
      meansdf[,idx] <- factor(meansdf[,idx])
    }
  }

  # rounding error
  for (i in 1:nrow(meansdf)) {
    if (meansdf[i,"means"] > -0.001 && meansdf[i,"means"] < 0.0) {
      meansdf[i,"means"] = 0
    }
    if (meansdf[i,"means"] > 1.0 && meansdf[i,"means"] < 1.001) {
      meansdf[i,"means"] = 1
    }
  }

  cat("Means DF\n")
  print(meansdf)
  #plotbar(meansdf, "Condition", dv, withinIVs[2], name)
  
  fit <- NULL; # fit for psychometric function

  #needToFit <- F

  # an easy hack to decide whether or not to fit, is based on whether we plot it
  #for (i in 1:length(plothashes)) {
  #  if (plothashes[[i]][["type"]] == "psycho") {
  #    needToFit <- T
  #    break;
  #  }
  #}

  if (plotPsycho) {
    # fit each condition and join them together

    # assume independent variable is first within IV, i.e. column 4 of meansdf
    # assume group variable is second within IV, i.e. column 5 of meansdf
    if (ncol(meansdf) < 5) stop("Error: not enough columns in meansdf")

    # change meansdf iv (column 4) to numeric for psychophysics
    meansdf[,4] <- remove.factor(meansdf[,4], convert="numeric")

    factorlevels <- levels(meansdf[,5])

    for (i in 1:length(factorlevels)) {
      fit.data <- meansdf[meansdf[,5] == factorlevels[i],]
      #fit.data.response <- fit.data$mean
      #fit.data.iv <- remove.factor(fit.data[,4], convert="numeric")
      #fit.cond <- psycho(fit.data.response, fit.data.iv)
      cat("about to fit\n")
      fit.cond <- psycho(response = fit.data$mean, iv = fit.data[,4])
      cat("fit\n")
      fit.cond$cond <- factorlevels[i]
      fit <- rbind(fit, fit.cond)
    }

  }

  # assume we have a variable with concatenated conditions:
  #plotbox(DF, withinIVs[1], dv, withinIVs[2], withinIVs[3], xlabel, filebase)
  # summarydf may be updated, e.g. combine
  summaryDFs <- ploteng(DF, meansdf, summarydf, fit, plothashes=plothashes, nosave=params$settings$noplotsave)
  #summaryDFs <- list(meansDF=meansdf, summaryDF=summarydf)

  summaryDFs
}

addPlotOpt <- function(plothash, str)
{
  if (is.null(plothash[["opts"]]) || plothash[["opts"]] == "") {
    plothash[["opts"]] = str
  }else{
    plothash[["opts"]] = paste(plothash[["opts"]], ", ", str, sep="")
  }

  return(plothash)
}

getConc <- function(DF, withinIVs)
  # concatenate vars for calcadj
  # TODO callers of this function are clobbering conc by using DF$conc
  # What if user has a variable called conc?
{

  # This function basically does this
  #conc <- do.call(paste, c(DF[,withinIVs], sep="."))
  #DF$conc <- factor(conc, levels=unique(conc))

  # but we replace existing dots first

  if (length(withinIVs) > 1) {
    # if withinIVs contain the separator char, replace it
    tempDF <- DF
    for (i in withinIVs) {
      tempDF[,i] <- gsub('\\.', '_', tempDF[,i])
    }
    conc <- do.call(paste, c(tempDF[,withinIVs], sep="."))
    conc <- factor(conc, levels=unique(conc))
  }else{
    conc <- DF[,withinIVs]
    if (!is.factor(conc)) conc <- factor(conc)
  }

  return(conc)
}

diagnosticplot <- function(dv, resids)
{
  p <- qplot( resids ) + opts(title=paste(dv, "residuals"))
  print(p)
}

doplot <- function(dv, params)
  # this is the main wrapper for all plotting
  # ploteng plots based on the hash information, and accepts a list of hashes (one for each plot)
  # the other plot functions are more experimental
{
  # get summary, means and stderrs data frame

  # output
  meansdf <- NULL
  summarydf <- NULL
  summaryDFs <- NULL

  # make params easier
  DF <- params$DF
  IDs <- params$IDs
  withinIVs <- params$withinIVs
  betweenIVs <- params$betweenIVs
  plothashes <- params$plothashes

  if (is.null(plothashes)) {
    # make a default plot: guess at what will be useful

    plothashes <- makeplothashes(DF, dv, withinIVs, betweenIVs)
  }else{
    plothashesdef <- makeplothashes(DF, dv, withinIVs, betweenIVs)
    # fix filllab - if not given then use the value from 'fill', otherwise default from makeplothashes will not match
    for (i in 1:length(plothashes)) {
      if (is.null(plothashes[[i]][["filllab"]])) plothashes[[i]][["filllab"]] <- plothashes[[i]][["fill"]]
    }

    # modify the defaults
    for (i in 1:length(plothashes)) plothashes[[i]] <- modifyList(plothashesdef[[1]], plothashes[[i]])

    # set colour to fill by default if plottype is psycho
    for (i in 1:length(plothashes)) {
      if (plothashes[[i]][["type"]] == "psycho") {
        if (is.null(plothashes[[i]][["colour"]])) plothashes[[i]][["colour"]] <- plothashes[[i]][["fill"]]
      }
    }
  }

  # store changes for subfunctions
  params$plothashes <- plothashes

  if (!is.null(betweenIVs)) {
    # TODO loop for all between IVs and within IVs
    #wiv <- withinIVs[1]
    biv <- betweenIVs[1]
    meansdf <- NULL
    #print(paste("wiv: ", wiv, ", biv: ", biv))
    # Could probably use summaryBy here:

    if (length(betweenIVs) > 0) cat("Only using one biv right now\n")

    if (!is.null(withinIVs)) {

      DF$conc <- getConc(DF, withinIVs)

      for (i in levels(DF[,biv])) {
        DFsub <- DF[DF[,biv]==i,]
        #print(head(DFsub))

        # correct within error bars for subset only
        temp <- calcadj(DFsub, variable="conc", value=dv)
        temp[,biv] <- unique(DFsub[,biv])

        if (is.null(meansdf)) meansdf <- temp
        else meansdf <- rbind(meansdf, temp)
        
        #for (wiv in withinIVs) {
        #  temp <- data.frame(means=tapply(DFsub[,dv], DFsub[,wiv], mean), stderrs=tapply(DFsub[,dv], DFsub[,wiv], std.error), wiv=levels(DFsub[,wiv]), biv=i)
        #  if (is.null(meansdf)) meansdf <- temp
        #  else meansdf <- rbind(meansdf, temp)
        #  print(meansdf)
        #}
      }

      rownames(meansdf) <- NULL
      colnames(meansdf)[3] <- "condition"
      colnames(meansdf)[4:(ncol(meansdf)-1)] <- withinIVs # -1 for biv
      for (tempcol in withinIVs) meansdf[,tempcol] <- factor(meansdf[,tempcol], levels=levels(DF[,tempcol]))
      meansdf[,biv] <- factor(meansdf[,biv], levels=levels(DF[,biv]))
      row.names(meansdf) <- NULL
      #attr(meansdf$means, "names") <- NULL # calcadj should do this really
      #attr(meansdf$stderrs, "names") <- NULL
      print(meansdf)

      summarydf <- calcSummaryDF(DF, biv, dv)

      if (length(withinIVs) == 1) {
        colnames(meansdf)[4] <- "wiv"
        plotwivbiv(meansdf, dv)
      }else{
        ploteng(DF, meansdf, summarydf, plothashes=plothashes, nosave=params$settings$noplotsave)
      }

    }else{

       # for one biv no within ivs
       means <- unname(tapply(DF[,dv], DF[,biv], mean))
       stderrs <- unname(tapply(DF[,dv], DF[,biv], std.error))
       meansdf <- data.frame(means=means, stderrs=stderrs, condition=levels(DF[,biv]))
       colnames(meansdf)[3] <- biv

       summarydf <- calcSummaryDF(DF, biv, dv)

       ploteng(DF, meansdf, summarydf, plothashes=plothashes, nosave=params$settings$noplotsave)

    }

  }else{ # If we don't have between IVs

    print("Correcting error bars with old Loftus and Masson method (no between IVs)...")
    if (length(withinIVs)==1) {
      meansdf <- calcadj(DF, variable=withinIVs[1], value=dv)
      colnames(meansdf)[3] <- "wiv"
      #print(meansdf)
      #summarydf <- calcSummaryDF(DF, withinIVs[1], dv)
      plotwivbiv(meansdf, dv)
    }else{
      # if psychophysic do this for each participant first then overall

      plotPsycho <- F

      for (i in 1:length(plothashes)) {
        if (plothashes[[i]][["type"]] == "psycho") {
          plotPsycho <- T
          break;
        }
      }

      if (plotPsycho) {
        # for each participant
        for (p in levels(DF$id)) {
          # get id
          params$DF <- DF[DF$id==p,] # danger overwriting params here, but it is put back below
          for (i in 1:length(params$plothashes)) {
            # TODO should change existing title if it's already an opt
            params$plothashes[[i]] <- addPlotOpt(params$plothashes[[i]], paste("title=\"Participant ", p, "\"", sep=""))
          }
          summaryDFs <- plotWithinMult(dv, params, plotPsycho)
          params$plothashes <- plothashes
        }
        # put params back
        params$DF <- DF
      }

      # overall
      summaryDFs <- plotWithinMult(dv, params, plotPsycho)

    }

  }

  # summaryDFs is either generated in a sub function or here
  if (is.null(summaryDFs)) summaryDFs <- list(meansDF=meansdf, summaryDF=summarydf)

  return(summaryDFs)
} # end doplot

# When there's something significant found we can split data and redo sigeng
redo <- function(sigs, params)
{
  retval <- NULL

  sink(params$logfile, split=T, append=T)
  #print("Starting redo")
  for (i in sigs) {
    # check for interaction
    if (length(grep(":", i)) > 0) {
      cat(params$redolevelsp, "## Redoing for interaction: ", i, "\n\n", sep="")
      covariates <- strsplit(i, ":")

      for (my.cov in covariates[[1]]) {

        cat(params$redolevelsp, "### Covariate: ", my.cov, "\n\n", params$redolevelsp, "- ", sep="")
        redoparams <- params
        redoparams$withinIVs <- params$withinIVs[params$withinIVs!=my.cov]
        redoparams$betweenIVs <- params$betweenIVs[params$betweenIVs!=my.cov]

        do.call(cat, list(levels(params$DF[,my.cov]), sep=", "))
        cat("\n\n")

        for (my.level in levels(params$DF[,my.cov])) {
          redoparams$redolevel <- params$redolevel + 1
          redoparams$redolevelsp <- makeredolevelsp(redoparams$redolevel)
          redoparams$DF <- params$DF # get DF back
          cat(redoparams$redolevelsp, "#### Level: ", my.level, "\n", sep="")
          cat(redoparams$redolevelsp, "\n", sep="")
          redoparams$DF <- redoparams$DF[redoparams$DF[,my.cov]==my.level,]
          redoparams$DF[,my.cov] <- factor(redoparams$DF[,my.cov]) # 'condition' variable won't have reliable factors
          sink()
            redoparams$settings$noplot=T # plots aren't labelled and get confusing when redoing
            retval <- sigengdv(redoparams)

          sink(params$logfile, split=T, append=T)

        }

      }

    }

  }

  #print("Exiting redo")
  sink()

  return(retval)
}

errorcheck <- function(params)
{
  if (length(params$IDs) > 1) {
    cat("Too many ID vars")
    return("Error")
  }

  for (i in params$DVs) {
    if (!(i %in% colnames(params$DF))) {
      cat("Sigeng error, DV: \"", i, "\" not found\n", sep="")
      return("Error")
    }
  }
  for (i in params$withinIVs) {
    if (!(i %in% colnames(params$DF))) {
      cat("Sigeng error, withinIV: \"", i, "\" not found\n", sep="")
      return("Error")
    }
  }
  for (i in params$betweenIVs) {
    if (!(i %in% colnames(params$DF))) {
      cat("Sigeng error, betweenIV: \"", i, "\" not found\n", sep="")
      return("Error")
    }
  }

  if (is.null(params$blocks)) {
    cat("No blocks\n")
  }else{
    if (params$blocks %in% colnames(params$DF)) {
      cat("Found blocks\n")
    }
  }

  return("OK")
}

error <- function()
{
  # clean up
  #detach(params)
  return("Error")
}

robustanova <- function(dv, params)
{
  DF <- params$DF
  IDs <- params$IDs
  withinIVs <- params$withinIVs
  betweenIVs <- params$betweenIVs

  robust <- NULL

  # robust

  # TODO could also use ranked based methods - this here is trimmed means
  if (length(betweenIVs) == 0 && length(withinIVs) == 1) {
    cat("Detected no bivs and 1 wiv\n")
    print(withinIVs)
    # Just within, can use rmanova from WRS package
    # Get groups as columns, pass it as matrix
    m <- melt(DF)
    m <- m[m$variable==dv,]
    # If we're redoing then this withinIV might not cover all the data, so first get the columns we're interested in
    m <- m[,c(IDs, "value", "variable", withinIVs[1])]
    command <- paste("DFwide <- as.matrix(cast(m, ... ~ variable + ", withinIVs[1], ", fun.aggregate=mean))", sep="")
    eval(parse(text=command))
    print(head(DFwide))

    robust <- rmanovab(DFwide, tr=.2)

    if (!params$settings$noanovaprint) {
      cat("Robust analysis trimmed means bootstrap:\n")
      print(robust)
    }

    robust <- friedman.test(DFwide)

    if (!params$settings$noanovaprint) {
      cat("Robust analysis friedman:\n")
      print(robust)
    }

    robust <- bprm(DFwide)

    if (!params$settings$noanovaprint) {
      cat("Robust analysis bprm:\n")
      print(robust)
    }

    robust <- rmanova(DFwide, tr=.2)

    if (!params$settings$noanovaprint) {
      cat("Robust analysis trimed means:\n")
      print(robust)
    }

    if (!params$settings$noanovasave) sink(params$logfile, split=T, append=T)
      if (robust$siglevel <= 0.05) 
        cat(params$redolevelsp, "+ **Robust significant**\n\n", sep="")
    if (!params$settings$noanovasave) sink()
  }else if (length(betweenIVs) == 1 && length(withinIVs) == 0) {
    
    # Get groups as columns, pass it as matrix
    m <- melt(DF)
    m <- m[m$variable==dv,]
    cat("Robust note: taking mean to aggregate trials\n")
    command <- paste("DFwide <- as.matrix(cast(m, ... ~ variable + ", betweenIVs[1], ", fun.aggregate=mean))", sep="")
    eval(parse(text=command))
    print(head(DFwide))
    robust <- rananova(DFwide, tr=.2)

    if (!params$settings$noanovaprint) {
      cat("Robust analysis:\n")
      print(robust)
    }

    if (!params$settings$noanovasave) sink(params$logfile, split=T, append=T)
      if (robust$siglevel <= 0.05) 
        cat(params$redolevelsp, "+ **Robust significant**\n\n", sep="")
    if (!params$settings$noanovasave) sink()

  }else{
    cat("Robust method not implemented for this design\n")
  }

  return(robust)
}

doanova <- function(dv, params)
{
  DF <- params$DF
  IDs <- params$IDs
  withinIVs <- params$withinIVs
  betweenIVs <- params$betweenIVs
  
  sigs <- NULL
  trends <- NULL

  command <- paste("my.aov <- aov(", dv, " ~ ", sep="")
  for (iv in withinIVs) {
    if (iv == withinIVs[1]) { command <- paste(command, iv, sep="") }
    else { command <- paste(command, " * ", iv, sep="") }
  }
  for (iv in betweenIVs) {
    if (iv == betweenIVs[1] && length(withinIVs)==0) { command <- paste(command, iv, sep="") }
    else { command <- paste(command, " * ", iv, sep="") }
  }
  if (length(withinIVs)>0) {
    command <- paste(command, " + Error(", IDs, "/(", sep="")
    for (iv in withinIVs) {
      if (iv == withinIVs[1]) { command <- paste(command, iv, sep="") }
      else { command <- paste(command, " * ", iv, sep="") }
    }

    command <- paste(command, "))", sep="")
  }
  # TODO add on betweenIVs here, see
  # http://wiki.stdout.org/rcookbook/Statistical%20analysis/ANOVA/
  command <- paste(command, ", data=DF)", sep="")
  # TODO noanovasave was meant for not saving anova file
  # Should add a nolog and change all sinks to log() and unlog()
  if (!params$settings$noanovasave) sink(params$logfile, split=T, append=T)
    cat(params$redolevelsp, "+ ", command, "\n", sep="")
  if (!params$settings$noanovasave) sink()

  eval(parse(text=command))

  if (!params$settings$noanovaprint) {
    #sink(params$logfile, split=T, append=T)
      print(summary(my.aov)) # TODO indent / make code block
      cat("\n")
    #sink()
  }

  # try to be clever and work out what the significant effects are
  if (!params$settings$noanovasave) {
    anovafilename <- params$anovafile

    sink(anovafilename)
      print(summary(my.aov))
    sink() # reset to stdout

    temp <- scan(anovafilename, "character", sep="\n")
    temp <- temp[grep("Signif", temp, invert=T)] # remove Signif legend because we're gonna search for *'s and .'s
    tempsig <- temp[grep("\\*", temp)] # search for *'s and .'s
    temptrend <- temp[grep("\\.$", temp)]

    sink(params$logfile, split=T, append=T)

      if (length(tempsig) > 0) {

        tempsig <- strsplit(tempsig, "[[:space:]]+")
        # now we've got a list of the significant lines, we want the first element of each - that's the main effect
        for (i in tempsig) { sigs <- c(sigs, i[1]) }

          cat(params$redolevelsp, "+ **Found significant effects**: ", sep="")
          do.call(cat, list(sigs, sep=", "))
          cat("\n")

      }else cat(params$redolevelsp, "+ Found no significant effects.\n", sep="")

      if (length(temptrend) > 0) {
        temptrend <- strsplit(temptrend, "[[:space:]]+")
        # now we've got a list of the significant lines, we want the first element of each - that's the main effect
        for (i in temptrend) { trends <- c(trends, i[1]) }

        cat(params$redolevelsp, "+ *Found trends*: ", sep="")
        do.call(cat, list(trends, sep=", "))
        cat("\n")

      }else cat(params$redolevelsp, "+ Found no trends.\n", sep="")

      cat("\n")
    sink()
  }

  # Store residuals for diagnostic plots
  my.resids <- NULL
  if (length(withinIVs) > 0) {
    # model for residuals
    command <- paste("my.lm <- lm(", dv, " ~ ", sep="")
    for (iv in withinIVs) {
      if (iv == withinIVs[1]) { command <- paste(command, iv, sep="") }
      else { command <- paste(command, " * ", iv, sep="") }
    }
    for (iv in betweenIVs) {
      if (iv == betweenIVs[1] && length(withinIVs)==0) { command <- paste(command, iv, sep="") }
      else { command <- paste(command, " * ", iv, sep="") }
    }
    command <- paste(command, " + ", IDs, ", data=DF)", sep="")
    if (!params$settings$noanovasave) sink(params$logfile, split=T, append=T)
      cat(params$redolevelsp, "+ ", command, "\n", sep="")
      eval(parse(text=command))
      cat(params$redolevelsp, "+ my.resids <- residuals(my.lm)\n", sep="")
      my.resids <- residuals(my.lm)
    if (!params$settings$noanovasave) sink()
  }else{
    if (!params$settings$noanovasave) sink(params$logfile, split=T, append=T)
      cat(params$redolevelsp, "+ my.resids <- residuals(my.aov)\n", sep="")
      my.resids <- residuals(my.aov)
    if (!params$settings$noanovasave) sink()
  }

  if (!params$settings$noanovaprint) {
    if (!params$settings$noanovasave) sink(params$logfile, split=T, append=T)
      cat(params$redolevelsp, "+ shapiro.test(my.resids)\n", sep="")
      my.shapiro <- shapiro.test(my.resids)
      print(my.shapiro)
      cat("\n")

      if (my.shapiro$p.value <= 0.05) {
        cat(params$redolevelsp, "+ **Found significant shapiro**\n", sep="")
      }

      cat("\n")
       
    if (!params$settings$noanovasave) sink()
  }
  # end residuals

  robust <- robustanova(dv, params)
  
  return(list(my.aov=my.aov, sigs=sigs, trends=trends, resids=my.resids, robust=robust))
}

makeredolevelsp <- function(redolevel)
{
  redolevelsp <- NULL

  # quote output based on which redo level we're on (0 is top level)
  if (redolevel > 0) {
    for (i in 1:redolevel) redolevelsp <- paste(redolevelsp, "> ", sep="")
  }

  return(redolevelsp)
}

sigengdv <- function(params)
# sigeng for one dependent variable, sub function of sigeng, used by redo
# no point in trying to keep params$DVs, it will be clobbered by redo
{
  # let's be clear:
  dv <- params$DVs # should only be one DV now (chosen by sigeng)

  # TODO ordinal analysis / convert to numeric

  # TODO sphericity test!
  # TODO check residuals (can plot with plot(residuals(myaov)))

  anova.results <- NULL
  if (!params$settings$noanova) {
    anova.results <- doanova(dv, params)
  }

  if (!params$settings$noplot) {
    if (!params$settings$noanova) diagnosticplot(dv, anova.results$resids)
    summaryList <- doplot(dv, params)
  }

  #checkForBetweenInteraction(summary(my.aov), IDs=IDs, DVs=DVs, withinIVs=withinIVs, betweenIVs=betweenIVs)

  # when we've got an interaction, split the data, recalculate factors and look at each factor level independently
  if (!params$settings$noredo && !is.null(anova.results$sigs)) {
    redo(anova.results$sigs, params)
  }

  sigengout <- list(anova.results=anova.results, summaryDFs=summaryList)
  return(sigengout)
}

# Here blocks means the experiment was repeated in more than one block
# noanova noplot and nosave (figs) are to speed things up when testing
sigeng <- function(DF, IDs="id", DVs=c("value"), withinIVs=c("cond"), betweenIVs=NULL, blocks=NULL, plothashes=NULL, logfile="/tmp/sigenglog.txt", clobberlog=T, anovafile="/tmp/anova.txt", ...)
{
  sigeng.version=0.9
  
  

  settingsdefault <- list(noanova=F, noanovasave=F, noanovaprint=F, noredo=F, noplot=F, noplotsave=F, noonlysig=F, noperl=F, nolatex=F)
  if (length(list(...)) > 0) {
    if (!is.list(list(...)[[1]])) settings <- list(...)
    else settings <- list(...)[[1]]

    settings <- modifyList(settingsdefault, settings)
  }else settings <- settingsdefault

  # make plothashes a list if it's not already
  if (!is.null(plothashes) && !is.null(names(plothashes))) plothashes <- list(plothashes)

  # copy them for easy passing to functions
  # NOTE We include settings in params
  params <- list(DF=DF, IDs=IDs, DVs=DVs, withinIVs=withinIVs, betweenIVs=betweenIVs, blocks=blocks, plothashes=plothashes, logfile=logfile, clobberlog=clobberlog, anovafile=anovafile, redolevel=0, redolevelsp=NULL, settings=settings)
  #rm(DF, IDs, DVs, withinIVs, betweenIVs, blocks, plothashes, logfile, anovafile)

  # Error check
  if (errorcheck(params)=="Error") return(error())

  # Check if we can actually do anova
  # If we have only 1 participant then we can't
  if (length(levels(DF[,params$IDs])) == 1 && !settings$noanova) {
    params$settings$noanova = T
    warning("Only 1 ID; setting noanova")
  }

  if (!settings$noanova && !settings$noanovasave) {
    myappend=T
    if (params$clobberlog) myappend=F
    sink(params$logfile, split=T, append=myappend)
    cat("% Sigeng (v", sigeng.version, ")\n% Trevor Dodds\n% September 2012\n\n", sep="")
    cat("# Overall analysis\n\n", sep="")
    sink()
  }

  #attach(params) -- problematic if we crash, leaves it attached

  # TODO outlier check for removing participants based on other vars

  # TODO MANOVA

  # TODO test with between subjects

  # TODO can't we do a MANOVA for within?
  #if (!is.null(withinIVs)) {
  #  cat("No MANOVA for within\n")
  #}

  for (dv in DVs) {
    if (is.numeric(DF[,dv])) {
      cat("+ Using dv: ", dv, "\n", sep="")
      params$DVs <- dv # clobber param DVs (would have happened anyway from redo)
      params$redolevelsp <- makeredolevelsp(params$redolevel)
      #anova.results <- sigengdv(params)
      sigengout <- sigengdv(params)
    }
  }

  #detach(params)

  if (!settings$noperl && !settings$noanova) {
    onlysig <- "--onlysig "
    nolatex <- NULL
    if (settings$noonlysig) onlysig <- NULL
    if (settings$nolatex) nolatex <- "--nolatex "
# TODO set base path and path to log.pl otherwise this will not work
    # quick attempt:
    basepath <- "~/reps/sigeng/extras/"
    command <- paste("perl ", basepath, "log.pl ", onlysig, nolatex, params$logfile, sep="")
    cat(command, "\n")
    system(command)
  }

  #summary(anova.results$my.aov)
  sigengout

  #if (is.null(myredo)) { return(DF) } else { return(myredo) }
}

