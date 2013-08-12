# v0.2
# Copied back to MPI

# v0.1
# Copied from MPI to Aurelie's laptop

#source("../libs/libs.R")
library(ggplot2)
library(WRS)
#library(hash)

remove.factor <- function(x, convert="character") {
  command <- paste("y <- as.", convert, "(levels(x)[x])", sep="")
  eval(parse(text=command))
  y
}

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

calcSignifDF <- function(DF, dv, sig)
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
  signifDF <- NULL
  textDF <- NULL
  xposStore <- NULL
  yposStore <- -Inf
  lastcut <- list(name="", value="")

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

    if (nrow(DFcut) > 0) {
      # Get start point, i.e. x pos where var is first value
      firstvalue  <- s[["first"]][[xvariable]]
      secondvalue <- s[["second"]][[xvariable]]

      if (!(xvariable %in% colnames(DFcut))) stop(paste0("calcSignifDF error: first variable '", xvariable,"' not found in DF."))

      # Find where this is in DF
      x <- which(levels(DFcut[,xvariable]) == firstvalue)
      xend <- which(levels(DFcut[,xvariable]) == secondvalue)

      ####
      # Dodge vertically
      # Go through from x to xend in DFcut and push it above the max value
      for (i in x:xend) {
        current <- DFcut[DFcut[,xvariable] == levels(DFcut[,xvariable])[i] , dv]
        cat("Cut value: ", s[["cut"]][["value"]], ", x: ", i, ", max dv: ", max(current), "\n")
        if (max(current) > ypos - height - 0.3) ypos <- max(current) + height + 0.3
      }

      # Now dodge prev ypos's

      if (s[["cut"]][["name"]] != lastcut[["name"]] || s[["cut"]][["value"]] != lastcut[["value"]]) { xposStore <- NULL; yposStore <- -Inf } # reset
      lastcut <- list(name=s[["cut"]][["name"]], value=s[["cut"]][["value"]]) # store last cut

      # pull out yposStore values where the corresponding xposStore values are within the range x:xend
      yposDodge <- NULL
      if (length(xposStore) > 0) {
        for (i in 1:length(xposStore)) {
          if (xposStore[i] %in% x:xend) yposDodge <- c(yposDodge, yposStore[i])
        }
        if (length(yposDodge) > 0) {
          if (max(yposDodge) > ypos - height - 0.3) ypos <- max(yposDodge) + height + 0.3 # dodge text also
        }
      }
      xposStore <- c(xposStore, x:xend)
      yposStore <- c(yposStore, rep(ypos, length(x:xend)))

      #
      ####


      # Make a signifDF based on this
      #  2>-------3
      #  |        V
      #  |        |
      #  ^
      #  1

      #          1  2  3                    1            2     3
      xs    <- c(x, x, xend);    ys    <- c(ypos-height, ypos, ypos)
      xends <- c(x, xend, xend); yends <- c(ypos,        ypos, ypos-height)
      

      temp <- data.frame(x=xs, y=ys, xend=xends, yend=yends, tempvar=rep(cutvalue, 3))
      temp <- rename(temp, c(tempvar=cutname))
      signifDF <- rbind(signifDF, temp)

      # Store text position, i.e. x, y, cutname, label
      text <- "****"
      if (!is.null(s[["text"]])) text <- s[["text"]]
      if (!is.null(s[["p"]])) {
        if (s[["p"]] <= 0.05) text <- "*"
        if (s[["p"]] <= 0.01) text <- "**"
        if (s[["p"]] <= 0.001) text <- "***"
      }
      temp <- data.frame(x=mean(c(x, xend)), y=ypos + 0.1, tempvar=cutvalue, label=text)
      temp <- rename(temp, c(tempvar=cutname))
      textDF <- rbind(textDF, temp)
    }
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

makeOptions <- function(plothash, optnames, addname=F)
  # Take optnames e.g. facet, theme, scale_fill
  # And produce a list, e.g. list(facet="facet_grid(...)", ...)
{

  # grab relevant options
  options <- plothash[optnames]

  # remove na's
  options <- options[!is.na(names(options))]

  # remove null's
  options <- options[!sapply(options, is.null)]

  if (addname)
  {
    # add name to value, so it becomes facet_grid(...) etc.
    options2 <- paste(names(options), options, sep="_")
    names(options2) <- names(options)
    options <- as.list(options2)
  }

  return(options)
}

# Inspired by http://stackoverflow.com/questions/3737619/how-can-i-pass-a-ggplot2-aesthetic-from-a-variable
# and aes_string
aes_list <- function(mapping)
{
  mapping[sapply(mapping, is.null)] <- "NULL"
  parsed <- lapply(mapping, function(x) parse(text = x)[[1]])
  structure(parsed, class = "uneval")
}

plotit <- function(DF, meansDF, summaryDF, fit, plothash, nosave=F) #x="variable", y, fill, rows=NULL, xlabel, filebase)
{
  #print(meansDF)

  # This function should be called with an individual plothash and it plots as follows:

  # 1. Combine columns if requested, updating summaryDFs as appropriate
  # 2. Make layers depending on plottype, e.g. bar makes geom_bar and geom_errorbar; Calc signifDFs with relevant DVs
  # 3. Add on user defined layers, e.g. facet_grid(~ ques)
  # 4. Evaluate and plot

  # NB: 'layers' stores a list of evaluated layers, e.g. geom_boxplot()
  #     'plotOptions' stores a list of unevaluated strings, e.g. "scale_x_discrete(...)" which could come from the user

  ####
  # 1. Prepare summary DFs

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
    DF[,newvar] <- factor(DF[,newvar], levels=unique(DF[,newvar]))

    # Same again for meansDF, because that's also used for some plots
    meansDF[,newvar] <- do.call(paste, c(meansDF[ , plothash[["combine"]] ], sep=".") )
    meansDF[,newvar] <- factor(meansDF[,newvar], levels=unique(meansDF[,newvar]))
  }


  ####
  # 3. Make layers

  # First start with plottype specific layers
  layers <- list()
  plotOptions <- list()
  signifDFs <- NULL
  xContinuous <- F

  # Default: boxplot
  if (is.null(plothash[["type"]])) plothash[["type"]] <- "boxplot"

  cat("Plot type:", plothash[["type"]], "\n")

  # Default: width
  if (is.null(plothash[["width"]])) plothash[["width"]] <- "0.9" # Confirmed default for boxplot

  # If ylimits does not cover area used by geom_segment (i.e. signifDFs) then the whole graph won't display.
  # So we take into account the users specification, but take min/max also including the signif values.
  # After initialising mins and maxs here, we append to them below depending on plot type.
  ylimits <- NULL
  if (!is.null(plothash[["ylimits"]])) {
    ylimits <- plothash[["ylimits"]] # This is technically the min/max value, but now we put all values together and use min() and max() below
  }


  # TODO could make these separate functions
  if (plothash[["type"]] == "bar") {
    means <- "means" # default variable names, can be overriden with plothash options
    stderrs <- "stderrs"
    if (!is.null(plothash[["means"]])) means <- plothash[["means"]]
    if (!is.null(plothash[["stderrs"]])) stderrs <- plothash[["stderrs"]]
    mapping <- list(x=plothash[["x"]], y=means)
    if (!is.null(plothash[["fill"]])) mapping <- c(mapping, fill=plothash[["fill"]])

    dodge <- position_dodge(width=0.9)

    # You can still plot without converting to a factor and it works, but I left this in just in case
    if (!is.null(plothash[["xfactor"]]) && plothash[["xfactor"]]) {
      cat("Converting x to factor for plotting...\n")
      meansDF[,plothash[["x"]]] <- factor(meansDF[,plothash[["x"]]], levels=unique(meansDF[,plothash[["x"]]]))
      str(meansDF)
    }

    # Here is how you can do it by converting mapping to a str and evaluating the code
    #mappingstr <- paste(names(mapping), mapping, collapse=", ", sep="=")
    #plottype <- list(paste("geom_bar(data=meansDF, mapping=aes(", mappingstr, "), position=dodge, stat=\"identity\")", sep=""))

    #layers <- geom_bar(data=meansDF, mapping=aes_list(mapping), stat="identity", position=dodge)
    mapping <- aes_string(x=plothash[["x"]], y=means, fill=plothash[["fill"]])
    layers <- geom_bar(data=meansDF, mapping=mapping, stat="identity", position=dodge)

    if (!is.factor(meansDF[,plothash[["x"]]])) xContinuous <- T
    cat("xContinuous:\n")
    print(xContinuous)
    cat("plothash[[x]]:\n")
    print(plothash[["x"]])

    ##command <- paste("limits <- aes(x=", plothash[["x"]], ", y=", means, ", ymin=", means, "-", stderrs, ", ymax=", means, "+", stderrs, ")", sep="")
    ##cat(command, "\n")
    ##eval(parse(text=command))
    #mappingstr <- paste(names(mapping), mapping, collapse=", ", sep="=")
    #plottype <- c(plottype, paste0("geom_errorbar(data=meansDF, mapping=aes(", mappingstr, "), width=0.1, position=dodge)"))


    #mapping <- list(x=plothash[["x"]], y=means, ymin=paste0(means, "-", stderrs), ymax=paste0(means, "+", stderrs))
    mapping <- aes_string(x=plothash[["x"]], y=means, ymin=paste0(means, "-", stderrs), ymax=paste0(means, "+", stderrs), fill=plothash[["fill"]])
    layers <- c(layers, geom_errorbar(data=meansDF, position=dodge, mapping=mapping, width=0.1))

    ylimits <- c(ylimits, meansDF[,means] - meansDF[,stderrs] - 0.1, meansDF[,means] + meansDF[,stderrs] + 0.1)

    if (!is.null(plothash[["sig"]])) {
      meansDF$tops <- meansDF[,means] + meansDF[,stderrs] # highest values
      signifDFs <- calcSignifDF(meansDF, "tops", plothash[["sig"]])
    }

  }else if (plothash[["type"]] == "point") {
    #plottype <- "geom_point()"
    mapping <- aes_string(x=plothash[["x"]], y=plothash[["y"]], colour=plothash[["colour"]], fill=plothash[["fill"]])
    layers <- geom_point(data=DF, mapping=mapping)
    ylimits <- c(ylimits, DF[,plothash[["y"]]])
  }else if (plothash[["type"]] == "psycho") {

    means <- "means"
    stderrs <- "stderrs"
    if (!is.null(plothash[["means"]])) means <- plothash[["means"]]
    if (!is.null(plothash[["stderrs"]])) stderrs <- plothash[["stderrs"]]

    #cat("debug psycho:\n")
    #cat("plothash[[\"x\"]]: ", plothash[["x"]], ", fill: ", plothash[["fill"]], ", colour: ", plothash[["colour"]], "\n", sep="")
    #print(unname(meansDF[,means]))
    #print(meansDF)
    # testing
    #meansDF[,means] <- unname(meansDF[,means])
    mapping <- aes_string(x=plothash[["x"]], y=means, fill=plothash[["fill"]], colour=plothash[["colour"]])
    layers <- geom_point(data=meansDF, mapping=mapping)
    ylimits <- c(ylimits, meansDF[,means])
    #temp.plot <- ggplot() + geom_point(data=meansDF, mapping=mapping) + geom_line(data=fit, aes(x=x, y=y, group=cond, colour=cond))
    #cat("Testing temp.plot\n")
    #print(temp.plot)
    #stop("debug stop")
    layers <- c(layers, geom_line(data=fit, aes(x=x, y=y, group=cond, colour=cond)))

    #print(fit)

    #mapping <- list(x=plothash[["x"]], y=means)

    #if (!is.null(plothash[["fill"]])) mapping <- c(mapping, list(fill=plothash[["fill"]]))
    #if (!is.null(plothash[["colour"]])) mapping <- c(mapping, list(colour=plothash[["colour"]]))

    ##plottype <- paste("geom_point(data=meansDF, aes(x=", plothash[["x"]], ", y=", means, sep="")
    #mappingstr <- paste(names(mapping), mapping, collapse=", ", sep="=")
    #plottype <- paste0("geom_point(data=meansDF, mapping=aes(", mappingstr, "))")

    # the ../psychophysics/psychophysics.r function should give a fit with x and y
    # and another function should join together different conditions with a variable 'cond'
    # so we have already x, y, cond
    #plottype <- c(plottype, "geom_line(data=fit, aes(x=x, y=y, group=cond, colour=cond))")

  }else if (plothash[["type"]] == "crossbar") {

    # hack levels for body illusion
    #summaryDF$body.sync <- factor(summaryDF$body.sync, levels=c("Body.Synchronous", "Body.Asynchronous", "Object.Synchronous", "Object.Asynchronous"))

    median <- "Median"
    X1stQu <- "X1st.Qu."
    X3rdQu <- "X3rd.Qu."
    if (!is.null(plothash[["median"]])) median <- plothash[["median"]]

    #plottype <- "geom_crossbar("

    mapping <- aes_string(x=plothash[["x"]], y=median, ymin=X1stQu, ymax=X3rdQu, fill=plothash[["fill"]])
    #cat("Width:\n")
    #print(as.numeric(plothash[["width"]]))
    layers <- geom_crossbar(data=summaryDF, mapping=mapping, width=as.numeric(plothash[["width"]]))

    #crossbar <- paste("data=summaryDF, mapping=aes(x=", plothash[["x"]], ", y=", median, ", ymin=", X1stQu, ", ymax=", X3rdQu, sep="")

    #if (!is.null(plothash[["fill"]])) crossbar <- paste(crossbar, ", fill=", plothash[["fill"]], sep="")

    #crossbar <- paste(crossbar, ")", sep="") # end aes

    #if (!is.null(plothash[["width"]])) crossbar <- paste(crossbar, ", width=", plothash[["width"]], sep="")

    #plottype <- paste(plottype, crossbar, ")", sep="") # end geom_crossbar

    ylimits <- c(ylimits, summaryDF[,X1stQu], summaryDF[,X3rdQu])
    if (!is.null(plothash[["sig"]])) signifDFs <- calcSignifDF(summaryDF, X3rdQu, plothash[["sig"]])

    # You can use this to draw a geom_boxplot() after the plot to replace the crossbar legend (not always useful) with a boxplot one
    if (!is.null(plothash[["legend"]])) {
      if (plothash[["legend"]] == "boxplot") layers <- c(layers, geom_boxplot(data=summaryDF, mapping=mapping, width=0))
      if (plothash[["legend"]] == "bar") {
        #layers <- c(layers, geom_bar(data=summaryDF, mapping=aes_string(x=0, y=0, fill=plothash[["fill"]]), position="dodge", stat="identity", width=0))

        plotOptions <- c(plotOptions, 'geom_bar(data=summaryDF, mapping=aes_string(x=1, y=0, fill=plothash[["fill"]]), position="dodge", stat="identity", width=0)')
        #legend <- "geom_bar(data=summaryDF, mapping=aes(x=0, y=0"
        #if (!is.null(plothash[["fill"]])) legend <- sprintf("%s, fill=%s", legend, plothash[["fill"]])
        #legend <- paste(legend, "), stat=\"identity\", width=0)", sep="")
      }
      #if (plothash[["legend"]] == "boxplot") legend <- paste("geom_", plothash[["legend"]], "(", crossbar, ")", sep="")
      #options <- c(options, legend)
    }
  }else{

    mapping <- aes_string(x=plothash[["x"]], y=plothash[["y"]], fill=plothash[["fill"]])
    layers <- geom_boxplot(data=DF, mapping=mapping, width=as.numeric(plothash[["width"]]))
    
    if (!is.null(plothash[["sig"]])) signifDFs <- calcSignifDF(DF, plothash[["y"]], plothash[["sig"]])
    ylimits <- c(ylimits, DF[,plothash[["y"]]])

    #plottype <- paste("geom_boxplot(data=DF, mapping=aes(x=", plothash[["x"]], ", y=", plothash[["y"]], sep="")
    #if (!is.null(plothash[["fill"]])) plottype <- paste(plottype, ", fill=", plothash[["fill"]], sep="")
    #plottype <- paste(plottype, ")", sep="") # end aes
    #if (!is.null(plothash[["width"]])) plottype <- paste(plottype, ", width=", plothash[["width"]], sep="")
    #plottype <- paste(plottype, ")", sep="") # end geom_boxplot
  }
  #if (!is.null(plothash[["xticlabs"]]) xticlabs <- plothash[["xticlabs"]]

  # Then add layers that are usually the same, with some checks on plottype where necessary
  
  if (plothash[["type"]] == "psycho" || xContinuous) {
    ## TODO just add to options
    # Note this was just for psycho:
    ##command <- paste(command, " + scale_x_continuous(name=\"", plothash[["xlab"]],
    ##                          "\", labels=c(\"", paste(plothash[["xticlabs"]], collapse='","', sep=""), "\"),", sep="")
    ## without quote marks around breaks, i.e. continuous var not factor
    ##command <- paste(command, " breaks=c(", paste(plothash[["xticbreaks"]], collapse=',', sep=""), "))", sep="")

    # Continuous
    plotOptions <- c(plotOptions, "scale_x_continuous(name=plothash[[\"xlab\"]], labels=plothash[[\"xticlabs\"]], breaks=as.numeric(plothash[[\"xticbreaks\"]]))")
  }else{
    # Discrete
    #command <- paste(command, " + scale_x_discrete(name=\"", plothash[["xlab"]],
    #                          "\", labels=c(\"", paste(plothash[["xticlabs"]], collapse='","', sep=""), "\"),", sep="")
    #command <- paste(command, " breaks=c(\"", paste(plothash[["xticbreaks"]], collapse='","', sep=""), "\"))", sep="")

    plotOptions <- c(plotOptions, "scale_x_discrete(name=plothash[[\"xlab\"]], labels=plothash[[\"xticlabs\"]], breaks=plothash[[\"xticbreaks\"]])")
  }

  if (is.null(plothash[["scale_fill"]])) {

    # Filllab is unique in the sense that it also gets its value from fill
    #fillparams <- list(name=plothash[["filllab"]], labels=plothash[["filllabs"]], breaks=plothash[["fillbreaks"]])
    #fillparamsstr <- paste(names(fillparams), fillparams, collapse=", ", sep="=")

    #plotOptions <- c(plotOptions, paste0("scale_fill_discrete(", fillparamsstr, ")")

    cat("Filllabs:\n")
    print(plothash[["filllabs"]])
    plotOptions <- c(plotOptions, 'scale_fill_discrete(name=plothash[["filllab"]], labels=plothash[["filllabs"]], breaks=plothash[["fillbreaks"]])')

   #   fillstr <- paste(fillstr, ", 
    #command <- paste(command, " + scale_fill_discrete(name=\"", plothash[["filllab"]], "\"", sep="")
    #if (!is.null(plothash[["filllabs"]])) {
    #  command <- paste(command, ", labels=c(\"", paste(plothash[["filllabs"]], collapse='","', sep=""), "\")", sep="")
    #  command <- paste(command, ", breaks=c(\"", paste(plothash[["fillbreaks"]], collapse='","', sep=""), "\")", sep="")
    #}
    #command <- paste(command, ")", sep="")
  }

  # Output signifDFs, change ylimits

  cat("signifDF\n")
  print(signifDFs$signifDF)
  print(signifDFs$textDF)

  if (!is.null(signifDFs)) ylimits <- c(ylimits, signifDFs$signifDF$y+0.1) # This 0.1 is to take into account text position - could also take max from text y's
  plothash[["ylimits"]]=c(min(ylimits), max(ylimits))
  #cat("debug psycho2:\n")
#plothash[["ylimits"]]=c(0, 1)
  # the null ones will be removed on the next line... clever!
  scaley <- list(name='plothash[["ylab"]]', limits=plothash[["ylimits"]], breaks=plothash[["yticbreaks"]], labels=plothash[["yticlabs"]])
  scaley <- scaley[!sapply(scaley, function (x) is.null(x))]
  scaleystr <- paste(names(scaley), scaley, collapse=", ", sep="=")
  plotOptions <- c(plotOptions, paste0("scale_y_continuous(", scaleystr, ")"))

  ####
  # 4. Add user defined layers

  plotOptions <- c(plotOptions, makeOptions(plothash, c("facet", "theme", "scale_fill"), addname=T))
  #layers <- c(layers, lapply(plotOptions, function (x) eval(parse(text=x))))

  # TODO opts deprecated
  if (!is.null(plothash[["opts"]]))
  {
    myopts <- paste("opts(", plothash[["opts"]], ")", sep="")
    plotOptions <- c(plotOptions, myopts)
  }

  # This is set with the sig option
  if (!is.null(signifDFs))
  {
    signif <- "geom_segment(data=signifDFs$signifDF, aes(x=x, y=y, xend=xend, yend=yend))"
    plotOptions <- c(plotOptions, signif)
  }

  # Text also comes from the sig option
  if (!is.null(signifDFs))
  {
    text <- "geom_text(data=signifDFs$textDF, mapping=aes(x=x, y=y, label=label))"
    plotOptions <- c(plotOptions, text)
  }

  #
  ####




  ####
  # Join plot commands together, and adjust the scale
  
  cat("Plot:\n")

  #command <- "myplot <- ggplot()"

  ##command <- paste("myplot <- p + ", plottype, errorbar, facet, sep="")
  #command <- paste(command, paste(plottype, collapse=" + "), sep=" + ")
  #command <- paste(command, paste(layers, collapse=" + "))

  print(plotOptions)

  #cat(command, "\n")
  cat("Making plot\n")
  myplot <- ggplot() + layers + lapply(plotOptions, function (x) eval(parse(text=x)))
  #myplot <- ggplot() + layers
  #cat("Layers:\n")
  #print(layers)
  #myplot <- ggplot() + lapply(layers, function (x) eval(parse(text=x))) + lapply(plotOptions, function (x) eval(parse(text=x)))

  #
  ####

  cat("Printing plot\n")
  print(myplot)
  cat("Printed plot\n")

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
  cat("Plot wivbiv\n")
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

  if (length(betweenIVs) == 2 && length(withinIVs) == 1) {
    # facet biv1, fill biv2, x wiv
    plothashes <- list(list(xlab=withinIVs[1], x=withinIVs[1], xticlabs=levels(DF[,withinIVs[1]]), xticbreaks=levels(DF[,withinIVs[1]]),
          ylab=dv, y=dv, fill=betweenIVs[2], filllab=betweenIVs[2], filllabs=levels(DF[,betweenIVs[2]]), fillbreaks=levels(DF[,betweenIVs[2]]), scale_fill=NULL, signif=NULL, text=NULL, width=NULL,
          facet=paste("grid(", betweenIVs[1], " ~ .)", sep=""), type="bar",
          opts="axis.title.x=theme_text(size=12, vjust=0), axis.title.y=theme_text(size=12, vjust=0.4, angle=90)", #, panel.margin=unit(0.5, \"cm\")",
          filenames=filenames))
  }

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
  if (length(betweenIVs) == 1 && length(withinIVs) == 1) {
    if (!is.factor(withinIVs)) {
      # put it on x axis
      # x within, fill between
      plothashes <- list(list(xlab=withinIVs[1], x=withinIVs[1], xticlabs=unique(DF[,withinIVs[1]]), xticbreaks=unique(DF[,withinIVs[1]]),
            ylab=dv, y=dv, fill=betweenIVs[1], filllab=betweenIVs[1], filllabs=levels(DF[,betweenIVs[1]]), fillbreaks=levels(DF[,betweenIVs[1]]), scale_fill=NULL, signif=NULL, text=NULL, width=NULL,
            facet=NULL, type="point",
            opts=NULL,
            filenames=filenames))
    }else{
      # x between, fill within
      plothashes <- list(list(xlab=betweenIVs[1], x=betweenIVs[1], xticlabs=levels(DF[,betweenIVs[1]]), xticbreaks=levels(DF[,betweenIVs[1]]),
            ylab=dv, y=dv, fill=withinIVs[1], filllab=withinIVs[1], scale_fill=NULL, signif=NULL, text=NULL, width=NULL,
            facet=NULL, type="bar",
            opts=NULL,
            filenames=filenames))
    }
  }
  if (length(betweenIVs) == 1 && length(withinIVs) == 0) {
    # x factor1
    plothashes <- list(list(xlab=betweenIVs[1], x=betweenIVs[1], xticlabs=levels(DF[,betweenIVs[1]]), xticbreaks=levels(DF[,betweenIVs[1]]),
          ylab=dv, y=dv, fill=NULL, filllab=NULL, scale_fill=NULL, signif=NULL, text=NULL, width=NULL,
          facet=NULL, type="bar",
          opts=NULL,
          filenames=filenames))
  }

  if (length(betweenIVs) == 0 && length(withinIVs) == 3) {
    # facet factor1, x factor2, fill factor3
    plothashes <- list(list(xlab=withinIVs[2], x=withinIVs[2], xticlabs=levels(DF[,withinIVs[2]]), xticbreaks=levels(DF[,withinIVs[2]]),
          ylab=dv, y=dv, fill=withinIVs[3], filllab=withinIVs[3], filllabs=levels(DF[,withinIVs[3]]), fillbreaks=levels(DF[,withinIVs[3]]),
          facet=paste("grid(. ~ ", withinIVs[1], ")", sep=""), type="bar", scale_fill=NULL, signif=NULL, text=NULL, width=NULL,
          opts="axis.title.x=theme_text(size=12, vjust=0), axis.title.y=theme_text(size=12, vjust=0.4, angle=90)", #, panel.margin=unit(0.5, \"cm\")",
          filenames=filenames))
  }

  if (length(betweenIVs) == 0 && length(withinIVs) == 2) {
    # facet factor1, x factor2, fill factor3
    plothashes <- list(list(xlab=withinIVs[1], x=withinIVs[1], xticlabs=levels(DF[,withinIVs[1]]), xticbreaks=levels(DF[,withinIVs[1]]),
          ylab=dv, y=dv, fill=withinIVs[2], filllab=withinIVs[2], filllabs=levels(DF[,withinIVs[2]]), fillbreaks=levels(DF[,withinIVs[2]]), scale_fill=NULL, signif=NULL, text=NULL, width=NULL,
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

  cat("Correcting stderr using old Loftus and Masson method by concatenating all withinIVs\n")
  if (length(betweenIVs) > 0) cat("TODO not sure what this correction is doing with between IVs\n")
  if (length(IDs) > 1) stop("sigeng/plotWithinMult: not sure how to deal with > 1 ID yet")
  meansdf <- calcadj(DF, id=IDs, variable="conc", value=dv)

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
  cat("Means DF before col rename\n")
  print(meansdf)

  # means, stderrs, variable, ...
  colnames(meansdf)[3] <- "condition"
  # when there's > 1 wiv
  cat("Got here", ncol(meansdf), "\n")
  if (ncol(meansdf) < 4) meansdf <- cbind(meansdf, unique(DF[,withinIVs[1]]))

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
  cat("Correcting values near zero and 1\n")
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

      if (!is.factor(DF[,i])) {
        warning("getConc: conc is converting a withinIV to a factor; concatenation not implemented for numeric.")
        tempDF[,i] <- factor(tempDF[,i])
      }

      if (is.factor(DF[,i]))
        my.levels.temp <- levels(DF[,i])
      else
        my.levels.temp <- levels(tempDF[,i]) # preserve levels if we can

      tempDF[,i] <- factor(gsub('\\.', '_', tempDF[,i]), levels=my.levels.temp) # need factor here because the levels are used below
    }
    #cat("head tempDF:\n")
    #print(head(tempDF))
    conc <- do.call(paste, c(tempDF[,withinIVs], sep="."))
    #cat("conc:\n")
    #print(conc)
    # get levels in same order as the previous levels (hopefully what the user wants)
    my.levels.temp <- lapply(tempDF[,withinIVs], levels)
    #print(my.levels.temp)
    # vary the last one first, so set others to length of one to the right
    for (i in (length(withinIVs)-1):1) {
      my.levels.temp[withinIVs[i]][[1]] <- rep(my.levels.temp[withinIVs[i]][[1]], each=length(my.levels.temp[withinIVs[i+1]][[1]]))
    }
    my.levels <- do.call(paste, c(my.levels.temp, sep="."))
    #print(my.levels)
    #conc <- factor(conc, levels=unique(conc))
    conc <- factor(conc, levels=my.levels)
  }else{
    conc <- DF[,withinIVs]
    if (!is.factor(conc)) conc <- factor(conc)
  }

  return(conc)
}

diagnosticplot <- function(dv, resids)
{
  p <- qplot( resids ) + ggtitle(paste(dv, "residuals"))
  print(p)
}

calcadj.biv <- function(DF, IDs, dv, betweenIVs)
# Adjustment for within subject error is done for each group of subjects (e.g. each level of the between subject variable)
{
  if (length(betweenIVs)==1) {
    # base case
    #cat("Base case\n")
    biv <- betweenIVs[1]
    #print(biv)
    #print(levels(DF[,biv]))
    meansdf <- NULL
    # TODO optimise with tapply
    for (i in levels(DF[,biv])) {
      DFsub <- DF[DF[,biv]==i,]
      #cat("DFsub:\n")
      #print(DFsub)

      # correct within error bars for subset only
      temp <- calcadj(DFsub, id=IDs, variable="conc", value=dv)
      temp[,biv] <- unique(DFsub[,biv])

      meansdf <- rbind(meansdf, temp)
      
      #print(meansdf)
      #for (wiv in withinIVs) {
      #  temp <- data.frame(means=tapply(DFsub[,dv], DFsub[,wiv], mean), stderrs=tapply(DFsub[,dv], DFsub[,wiv], std.error), wiv=levels(DFsub[,wiv]), biv=i)
      #  if (is.null(meansdf)) meansdf <- temp
      #  else meansdf <- rbind(meansdf, temp)
      #  print(meansdf)
      #}
    }

    row.names(meansdf) <- NULL
    return(meansdf)
  }else{
    # cut down DF and recurse for each level
    # start with outermost (last) biv
    cat("Not base case\n")
    biv <- betweenIVs[length(betweenIVs)]
    #print(biv)
    #print(levels(DF[,biv]))
    meansdf <- NULL
    betweenIVs.temp <- betweenIVs[1:(length(betweenIVs)-1)]
    # TODO optimise with tapply
    for (i in levels(DF[,biv])) {
      DFsub <- DF[DF[,biv]==i,]
      DFsub$biv <- NULL
      #print(DFsub)
      #print(betweenIVs.temp)
      meansdf.temp <- calcadj.biv(DFsub, IDs, dv, betweenIVs.temp)
      meansdf.temp <- cbind(meansdf.temp, i)
      colnames(meansdf.temp)[length(colnames(meansdf.temp))] <- biv
      meansdf <- rbind(meansdf, meansdf.temp)
      #print(meansdf)
    }
    return(meansdf)
  }

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
    # same for filllabs and fillbreaks
    for (i in 1:length(plothashes)) {
      if (is.null(plothashes[[i]][["filllab"]])) plothashes[[i]][["filllab"]] <- plothashes[[i]][["fill"]]
      if (is.null(plothashes[[i]][["filllabs"]])) plothashes[[i]][["filllabs"]] <- levels(DF[,plothashes[[i]][["fill"]]])
      if (is.null(plothashes[[i]][["filllabbreaks"]])) plothashes[[i]][["fillbreaks"]] <- levels(DF[,plothashes[[i]][["fill"]]])
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

    #if (length(betweenIVs) > 2) cat("MAX 2 BIVs RIGHT NOW\n")

    if (!is.null(withinIVs)) {

      DF$conc <- getConc(DF, withinIVs)

      print(betweenIVs)
      meansdf <- calcadj.biv(DF, params$IDs, dv, betweenIVs)
      cat("MeansDF:\n")
      print(meansdf)
      if (length(withinIVs) > 1) # we have conc
        colnames(meansdf)[4:(ncol(meansdf)-length(betweenIVs))] <- withinIVs
      else
        colnames(meansdf)[3:(ncol(meansdf)-length(betweenIVs))] <- withinIVs
      for (i in withinIVs) {
        if (is.factor(DF[,i])) meansdf[,i] <- factor(meansdf[,i], levels=levels(DF[,i]))
      }
      cat("MeansDF:\n")
      print(meansdf)
      str(meansdf)

      ##rownames(meansdf) <- NULL
      ##colnames(meansdf)[3] <- "condition"
      #for (tempcol in withinIVs) {
      #  # If it was a factor, then make it a factor in meansdf also.
      #  # It will actually always be a factor, because that was needed to adjust the error bars (calcadj) and to generate conc
      #  # but here we sort out the levels, and remove the factor if it originally wasn't one
      #  if (is.factor(DF[,tempcol])) {
      #    meansdf[,tempcol] <- factor(meansdf[,tempcol], levels=levels(DF[,tempcol]))
      #  }else{
      #    meansdf[,tempcol] <- remove.factor(meansdf[,tempcol], convert=class(DF[,tempcol]))
      #  }
      #}
      #if (is.factor(DF[,biv])) meansdf[,biv] <- factor(meansdf[,biv], levels=levels(DF[,biv]))
      ##print(meansdf)

      summarydf <- calcSummaryDF(DF, biv, dv)

      #if (length(withinIVs) == 1) {
      #  colnames(meansdf)[4] <- "wiv"
      #  plotwivbiv(meansdf, dv)
      #}else{
        ploteng(DF, meansdf, summarydf, plothashes=plothashes, nosave=params$settings$noplotsave)
      #}

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

    #cat("Correcting error bars with old Loftus and Masson method (no between IVs)...\n")
    ##if (length(withinIVs)==1) {
    #  meansdf <- calcadj(DF, variable=withinIVs[1], value=dv)
    #  colnames(meansdf)[3] <- "wiv"
    #  #print(meansdf)
    #  #summarydf <- calcSummaryDF(DF, withinIVs[1], dv)
    #  plotwivbiv(meansdf, dv)
    #}else{
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

    #}

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
    print(withinIVs[1])
    # Just within, can use rmanova from WRS package
    # Get groups as columns, pass it as matrix
    # Get rid of variables we don't need
    DFtemp <- DF[,c(IDs, withinIVs[1], dv)]
    m <- melt(DFtemp, id.vars=c(IDs, withinIVs[1]))
    # We've already got rid of variables we don't need, so no need for this now
    #m <- m[m$variable==dv,]
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

  # Do linear model, useful e.g. in case IVs are numeric (not factor)
  command <- paste("my.lm <- lm(", dv, " ~ ", sep="")
  for (iv in withinIVs) {
    if (iv == withinIVs[1]) { command <- paste(command, iv, sep="") }
    else { command <- paste(command, " * ", iv, sep="") }
  }
  for (iv in betweenIVs) {
    if (iv == betweenIVs[1] && length(withinIVs)==0) { command <- paste(command, iv, sep="") }
    else { command <- paste(command, " * ", iv, sep="") }
  }
  command <- paste(command, ", data=DF)", sep="")
  eval(parse(text=command))
  if (!params$settings$noanovasave) sink(params$logfile, split=T, append=T)
    cat(params$redolevelsp, "+ ", command, "\n", sep="")
    print(summary(my.lm))
  if (!params$settings$noanovasave) sink()
  

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
      #print(summary(my.lm))
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

  robust <- NULL
  if (!params$settings$norobust) robust <- robustanova(dv, params)
  
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

create.dir <- function(dirname)
{
  errorcheck <- dir.create(dirname)
  if (!errorcheck[1]) stop(paste0("Error, could not create directory: ", dirname))
}

cache.save <- function(filename, obj, path)
{
  if (!file.exists(path)) {
    # If we use the default cache path, we can help out and create it all
    # Otherwise it will create max one dir
    if (!file.exists("~/.sigeng")) create.dir("~/.sigeng")
    create.dir(path)
  }

  save(obj, file=paste0(path, "/", filename))
}

cache.load <- function(filename, path)
{
  cat("Attempting cache load\n")
  filenameall <- paste0(path, "/", filename)
  print(filenameall)
  if (!file.exists(filenameall)) {
    cat(paste0("Could not find cache file: ", filenameall, "\n"))
    return (NULL) # Error
  }else{
    temp.space <- new.env()
    foo <- load(filenameall, temp.space) # this will probably be called "obj"
    the.object <- get(foo, temp.space)
    rm(temp.space)
  }

  return(the.object)
}

settings.fast <- function(norobust=T)
{
  return(list(noanova=F, noanovasave=T, noanovaprint=T, noredo=T, noplot=T, noplotsave=T, noonlysig=T, noperl=T, nolatex=T, norobust=norobust))
}

brute.outs <- function(params)
{
  
  cat("Brute: testing outs\n")
  if (params$brute$outs$test > 3)
  {
    warning("brute.outs: Not dealing with test > 3 yet.\n")
    params$brute$outs$test <- 3
  }

  finished = F
  myouts <- NULL

  # speed up local params
  params$settings <- settings.fast(norobust=T)

  cat(paste0("ID var: ", params$IDs, "\n"))
  IDs <- params$DF[, params$IDs]
  print(IDs)
  #if (is.factor(IDs)) IDs <- as.integer(levels(IDs)) # TODO deal with levels that are string, i.e. change this function a bit
  #else IDs <- unique(IDs)

  IDs.numeric <- IDs
  if (is.factor(IDs)) {
    # generate numeric IDs, because algorithm needs to start from i to endID
    IDs.numeric <- rep(-1, times=length(levels(IDs))) # make array, default -1
    for (i in 1:length(levels(IDs))) { # set elements of the array to increase
      IDs.numeric[i] <- i
    }
    print(IDs.numeric)
  }
  startID <- min(IDs.numeric)
  endID <- max(IDs.numeric)

  cat(paste0("startID: ", startID, "\n"))
  cat(paste0("endID: ", endID, "\n"))

  for (i in startID:endID) {
    out1 <- i
    for (j in i:endID) {
      #if (i == j) next # allow duplicates because then it will just remove one element
      out2 <- c(out1, j)
      
      for (k in j:endID) {
        outs <- c(out2, k)
        if (params$brute$outs$test == 1) outs <- c(k)
        if (params$brute$outs$test == 2) outs <- c(j, k)
        cat("Leaving out: ")
        do.call(cat, list(outs, sep=","))
        cat("\n")
        DFbak <- params$DF
          if (is.factor(IDs)) {
            # lookup outlier number, i.e. levels(IDs)[outs]
            cat("Removing the following outliers:\n")
            print(levels(IDs)[outs])
            params$DF <- params$DF[!(params$DF$id %in% levels(IDs)[outs]),]
          }else{
            # if numeric, this is fine
            params$DF <- params$DF[!(params$DF$id %in% outs),]
          }
          myeng <- sigengdv(params)
        params$DF <- DFbak
        #cat("Anova in brute:\n")
        #print(myeng$anova.results$my.aov)
        pval <- summary(myeng$anova.results$my.aov)[[2]][[1]][["Pr(>F)"]][1]
        if (is.null(myeng$anova.results$robust)) pval2 <- 1
        else pval2 <- myeng$anova.results$robust$siglevel

        matchedToIgnore <- F
        if (is.factor(IDs)) matchedToIgnore <- identical(levels(IDs)[outs], params$brute$outs$ignore)
        else matchedToIgnore <- identical(outs, as.integer(params$brute$outs$ignore))

        if ((pval <= 0.05 || pval2 <= 0.05) && !matchedToIgnore) {
          #print(outs)
          sigtest <- NULL
          if (pval <= 0.05) sigtest <- "normal"
          if (pval2 <= 0.05) sigtest <- c(sigtest, "robust")
          #print(sigtest)
          #myouts <- c(myouts, list(myeng=myeng, outs=outs, sigtest=sigtest))
          if (is.factor(IDs)) {
            myouts <- c(myouts, list(outs=levels(IDs)[outs], sigtest=sigtest))
          }else{
            myouts <- c(myouts, list(outs=outs, sigtest=sigtest))
          }
          #finished=T
          #break;
        }
      }

      if (params$brute$outs$test==1) finished=T
      if (finished) break;
    }

    if (params$brute$outs$test==2) finished=T
    if (finished) break;
  }

  if (finished) cat("Brute.outs: Finished with test < 3.\n")
  
  for (i in myouts) {
    cat("outs[i]:")
    print(i)
  }

  if (is.null(myouts)) {
    myouts <- "No brute outliers" # code for no outliers
    cat("No brute outliers\n")
  }

  cache.save(cache.mkfilename("brute.outs", params), myouts, params$cachepath)
  # Store settings when outs were generated
  outs.settings <- list(ignore=params$brute$outs$ignore, norobust=params$settings$norobust)
  cache.save(cache.mkfilename("brute.outs.settings", params), outs.settings, params$cachepath)
  return(myouts)
}

cache.mkfilename <- function(name, params)
# You can optionally use this to make a standard file name
{
  return(paste(name, params$analysisID, params$DVs, "Rdata", sep="."))
}

sigengbrutedv <- function(params)
# First check if we need to do any brute stuff
# Then call sigengdv, e.g. after removing outliers
# Return the same as sigengdv
{
  outs <- NULL
  if (!is.null(params$brute$outs$test) && params$brute$outs$test > 0) {
# test for outliers and cache
# note we have to store the DV along with the cache
    outs <- cache.load(cache.mkfilename("brute.outs", params), params$cachepath)
    if (is.null(outs)) outs <- brute.outs(params) # do for real
    else cat("Loaded brute.outs from cache\n")
    if (length(outs) == 1 && outs == "No brute outliers") cat(outs, "\n")
    else {
      for (i in outs) {
        cat("outs[i]:")
          print(i)
      }
    }
  }
  if (!is.null(params$brute$outs$retest) && params$brute$outs$retest > 0) {
    if (!is.null(outs)) warning("brute$outs$test and brute$outs$retest both set")
    params$brute$outs$test <- params$brute$outs$retest
    outs <- brute.outs(params)
  }

  # Now run sigengdv 'for real' with outs removed
  if (!is.null(outs) && outs != "No brute outliers") {
    params$DF <- params$DF[!(params$DF$id %in% outs),]
  }
  sigengdv(params)
}

sigengdv <- function(params)
# sigeng for one dependent variable, sub function of sigeng, used by redo
# no point in trying to keep params$DVs, it will be clobbered by redo
{
  # let's be clear:
  dv <- params$DVs # should only be one DV now (chosen by sigeng)

  # TODO ordinal analysis / convert to numeric for DV

  # TODO convert independent variable to numeric, e.g. xzratio factor in ANOVA or numeric in regression

  # TODO sphericity test!

  anova.results <- NULL
  if (!params$settings$noanova) {
    anova.results <- doanova(dv, params)
  }

  summaryList <- NULL
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

# AnalysisID is used for caching/reloading cache
# Here blocks means the experiment was repeated in more than one block
# noanova noplot and nosave (figs) are to speed things up when testing
# brute is for brute force checking, where outs$test checks the cache for outliers (=no. of outliers to remove), outs$retest forces a rebuild of the cache, pitac is for problematic collaborators, generates loads of figures TODO
# outs$ignore will not remove those IDs, used with outs$test and outs$retest
# TODO check IVs continuous or factor
sigeng <- function(DF, analysisID="default", IDs="id", DVs=c("value"), withinIVs=c("cond"), betweenIVs=NULL, blocks=NULL, plothashes=NULL,
                   brute=list(outs=list(test=0, retest=0, ignore=NULL), pitac=F), cachepath="~/.sigeng/cache",
                   logfile="/tmp/sigenglog.txt", clobberlog=T, anovafile="/tmp/anova.txt", ...)
{
  sigeng.version=0.9
  
  

  settingsdefault <- list(noanova=F, noanovasave=F, noanovaprint=F, noredo=F, noplot=F, noplotsave=F, noonlysig=F, noperl=F, nolatex=F, norobust=F)
  if (length(list(...)) > 0) {
    if (!is.list(list(...)[[1]])) settings <- list(...)
    else settings <- list(...)[[1]]

    settings <- modifyList(settingsdefault, settings)
  }else settings <- settingsdefault

  # make plothashes a list if it's not already
  if (!is.null(plothashes) && !is.null(names(plothashes))) plothashes <- list(plothashes)

  # copy them for easy passing to functions
  # NOTE We include settings in params
  params <- list(analysisID=analysisID, DF=DF, IDs=IDs, DVs=DVs, withinIVs=withinIVs, betweenIVs=betweenIVs, blocks=blocks, plothashes=plothashes, brute=brute, cachepath=cachepath, logfile=logfile, clobberlog=clobberlog, anovafile=anovafile, redolevel=0, redolevelsp=NULL, settings=settings)
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
  # This should have a global outlier check here based on all DVs, e.g. using out() from WRS
  # Note: can use a function for this, currently working on outliers()

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
      
      sigengout <- sigengbrutedv(params)
      #sigengout <- sigengdv(params)
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



# useful utils
file.overwrite.prompt <- function(filename=stop("file.overwrite.prompt: no filename specified"))
{
  overwrite <- T # default, i.e. also if file does not exist

  if (file.exists(filename)) {
    cat("File exists:", filename, "\n")
    cat("Overwrite?\n")
    response <- readLines(stdin(), 1)
    if (response=="y") {
      cat("Overwriting.\n")
    }else{
      cat("Not overwriting.\n")
      overwrite <- F
    }
  }else{
    cat("Creating new file:", filename, "\n")
  }

  return(overwrite)
}

remake.check <- function(fun.call, cache.file, clobber=F)
  # e.g. remake <- remake | remake.check(match.call(), "cache/analyse.fun.call.Rdata", clobber=T) # Need to call remake.check to clobber the cache. Single | does not short-circuit.

  # cache.file <- "cache/analyse.Rdata"
  # if (remake) {
  #   ... # reload/process/format DF
  #   save(DF, file=cache.file)
  # }else{
  #   load(cache.file)
  # }
  #
  # If clobber is true and we need a remake then the new function call will be saved (and file clobbered if it exists).
  #
  # Otherwise you'll need to do this instead:
  #
  # cache.fun.call <- "cache/analyse.fun.call.Rdata"
  # fun.call <- match.call()
  # if (!remake) remake <- remake.check(fun.call, cache.fun.call)
  # ...
  # cache.file <- "cache/analyse.Rdata"
  # if (remake) {
  #   save(fun.call, file=cache.fun.call)
  #   ...
  #   ... # reload/process/format DF
  #   save(DF, file=cache.file)
  # }else{
  #   load(cache.file)
  # }
  #
  # General workflow:
  # Raw data (e.g. *.xml, potentially appended to in case of a crash)
  # -> Files ready to be read and processed (e.g. *-r.xml with only the last <?xml> cut out. Remove these files to remake them.)
  # -> "Intermediary files": Relevant data extracted from *-r.xml and converted into R data frame, with some initial calculations
  # (e.g. thresholds); saved to a cache.
  # The intermediary files only need to be remade (remake=T) if the arguments change (e.g. number of reversals
  # to ignore for calculating the threshold) or if the raw data change.
  # The actual arguments are cached separately from the intermediary files (e.g. cache/analyse-fun.call.Rdata).
  # Here, we check if the actual arguments from the function call changed by comparing this call to the one in the cache.
  # We could also check if raw data changes, e.g. by having an md5sum or sha1 for each subject, or just let the user specify
  # remake=T in the analyse() function when they update the raw data. This is probably quicker and easier than checking all the md5sums
  # each time we run the analysis.
{
  remake <- F

  # First, check if the cache even exists.
  # If it doesn't, then we can presume this is a first run and we need to make ('remake') our intermediary files.
  if (file.exists(cache.file)) {
    # Load cached arguments into a new environment, so as not to mask anything.
    newenv <- new.env()
    load(cache.file, newenv)

    # Check if loaded call matches current call.
    # Actually cannot just check calls like this:
    # !identical(newenv$fun.call, fun.call)
    # because then if the user adds remake=F it would actually remake (one time only).


    # So we check the actual args individually
    mismatch <- F
    fun.call.list <- as.list(fun.call)
    cached.call.list <- as.list(newenv$fun.call)

    # Go through the longer one and check they match (if length is different and the new arg is just remake this is fine)
    if (length(fun.call.list) > length(cached.call.list)) longest.names <- names(fun.call.list)
    else longest.names <- names(cached.call.list)

    for (i in longest.names) {
      if (i != "remake" && !identical(fun.call.list[[i]], cached.call.list[[i]])) { mismatch <- T; break }
    }

    if (mismatch) {
      cat("Detected function call change, will remake.\n")
      remake <- T
    }

  }else{
    cat("No cache of function call, will remake.\n")
    remake <- T
  }

  if (remake && clobber)
  {
    # Create/update cache.file
    save(fun.call, file=cache.file)
  }
  
  remake
}

outliers.getidx <- function(dv, method, method.param=NULL)
# Simply gets the indexes of outliers. See outliers() for details of method and method.param
{
  idx <- NULL
  highlow <- NULL
  if (method=="classic")
  {
    cat("Using classic method for detecting outliers.\n")
    if (is.null(method.param)) stop("sigeng: outliers.sub() method.param not defined")
    if (!is.numeric(method.param)) stop("sigeng: outliers.sub() method.param is not numeric")
    out1 <- dv[ dv > mean(dv) + sd(dv) * method.param ]
    out2 <- dv[ dv < mean(dv) - sd(dv) * method.param ]
    if (length(out1) > 0 || length(out2) > 0)
    {
      cat("Outlier values:\n")
      print(out1)
      print(out2)
    }

    # Get idx's of outliers
    idx     <- which(dv > mean(dv) + sd(dv) * method.param)
    if (length(idx) > 0) highlow <- rep("high", length(idx))
    idx2    <- which(dv < mean(dv) - sd(dv) * method.param)
    if (length(idx2) > 0) highlow <- c(highlow, rep("low", length(idx2)))
    idx <- c(idx, idx2)
    #cat("outliers.getidx() highlow:\n")
    #print(highlow)
  }
  if (method=="boxplot")
  {
    idx <- which(dv %in% boxplot(dv)$out)
    # TODO need to set highlow for boxplot()$out
  }
  if (method=="out")
  {
    # TODO note could do this in matrix form, col1 is var1 and col2 is var2 (then it plots it for you!)
    idx <- out(dv)$out.id
    # TODO need to set highlow for out()
  }
  if (method=="ggboxplot") stop("sigeng: outliers.getidx() method ggboxplot not yet implemented.")
  if (is.null(idx)) stop(paste("sigeng: outliers.getidx() method not found:", method))

  if (length(idx) > 0)
  {
    cat("Outlier indexes:\n")
    print(idx)
  }

  if (length(highlow) != length(idx)) stop("outliers.getidx(): highlow and idx are of different lengths.")

  highlow <- factor(highlow, levels=c("high", "low"))

  return(list(idx=idx, highlow=highlow))
}

find.row <- function(DF, searchlist, return.idx=F)
# e.g. find.row(thresholds, find.row(thresholds, list(staircaseID=1, sid="01", xzratio=2))
# return.idx means return the row numbers (integer) instead of the rows themselves
# Note: R will warn if some objects in searchlist are not a multiple of shorter object lengths.
{
  #for (i in 1:nrow(df)) if (all(df[i,] == row)) return(i)
  #return(NULL)
  mat <- DF[,names(searchlist)]==searchlist
  logv <- aaply(mat, 1, all)
  if (!return.idx) return(DF[logv,])
  else return(which(unname(logv)))
}

outliers.respond <- function(params, split.blocks, method, method.param, response, response.param)
# See outliers() for details
{
  if (length(params$DVs) > 1) stop("sigeng: outliers.respond() only deals with 1 DV for now")

  DF <- params$DF # easier

  if (!split.blocks) {
    temp <- outliers.getidx(DF[,params$DVs[1]], method, method.param)
    idx <- temp$idx
    highlow <- temp$highlow
  }else{
    # Get index for each block, but update idx based on position in original DF
    # It actually uses the row names of the DF for the new idx, which could be character strings but must be unique
    #cat("Check row numbers:\n")
    #print(head(DF))
    origidx <- NULL
    highlow <- factor(NA, levels=c("high", "low")) # need to start with NA to append with unlist(list())
    for (b in levels(DF[,params$blocks]))
    {
      DFcut <- DF[DF[,params$blocks]==b,]
      cat("Block:", b, "\n")
      temp <- outliers.getidx(DFcut[,params$DVs[1]], method, method.param)
      idx <- temp$idx
      if (length(temp$highlow) > 0) highlow <- unlist(list(highlow, temp$highlow))
      #cat("highlow5:\n")
      #print(highlow)
      #print(temp$highlow)
      for (i in idx)
      {
        #cat("idx:", i, "\n")
        #cat("row:\n")
        #print(DFcut[i,])
        #cat("row number:\n")
        rown <- row.names(DFcut[i,]) 
        #print(rown)
        #cat("original row:\n")
        #print(DF[rown,])
        origidx <- c(origidx, rown)
      }
      #cat("Done split blocks outlier check. origidx:", origidx, "\n")
    }
    idx <- origidx
    #cat("idx is now:", idx, "\n")
    if (length(highlow) > 1) highlow <- highlow[2:length(highlow)] # remove first NA
    else highlow <- NULL
  }

  #cat("######\n")
  #cat("idx:\n")
  #print(idx)
  #cat("highlow:\n")
  #print(highlow)
  #cat("######\n")

  if (!is.null(idx))
  {
    responded <- F

    #cat("Choosing response...\n")

    if (response=="meanotherblocks" || response=="meanblocks")
    {
      if (is.null(params$blocks)) stop("sigeng: outliers.respond() no blocks defined")

      responded <- T

      cat("Responding with mean of other blocks\n")
      #print(idx)

      # Go through idx. Note idx could be character string or numeric, and can be used like this: DF[i,]
      for (i in idx)
      {
        # TODO can I just use find.row (above function)?
        cat("Outlier row:\n")
        print(DF[i,])
        # lookup other blocks for same ID and IVs
        # So first get current block, ID, and IVs
        currentBlock <- DF[i,params$blocks]
        #cat("currentBlock:\n")
        #print(currentBlock)
        currentId <- DF[i,params$IDs]
        #cat("CurrentId:\n")
        #print(currentId)
        #cat("params$withinIVs:\n")
        #print(params$withinIVs)
        currentWithinIVs <- DF[i,params$withinIVs, drop=F] # don't drop, so if there's only one column it will stay as a DF
        #cat("CurrentWithinIVs:\n")
        #print(currentWithinIVs)
        #print(class(currentWithinIVs))
        currentBetweenIVs <- DF[i,params$betweenIVs, drop=F]
        lookupBlocks <- levels(DF[,params$blocks])[levels(DF[,params$blocks])!=currentBlock]
        #cat("lookupBlocks:\n")
        #print(lookupBlocks)

        # Now look them all up
        # Make a DF of currentWithinIVs for comparison
        currentWithinIVsDF  <- ldply(1:nrow(DF), function(x) currentWithinIVs)
        currentBetweenIVsDF <- ldply(1:nrow(DF), function(x) currentBetweenIVs)
        # Set matching cells to True
        #print(class(currentWithinIVsDF))
        #print(class(DF[,params$withinIVs, drop=F]))

        if (length(params$withinIVs) > 0) {
          #print(currentWithinIVsDF)
          logicalWithinIVMatch <- DF[,params$withinIVs, drop=F] == currentWithinIVsDF
          #print(logicalWithinIVMatch)
          rowsWithin <- aaply(logicalWithinIVMatch, 1, all)
          #print(rowsWithin)
        } else rowsWithin <- rep(T, nrow(DF))

        if (length(params$betweenIVs) > 0) {
          logicalBetweenIVMatch <- DF[,params$betweenIVs, drop=F] == currentBetweenIVsDF
          # When whole row matches set it to True
          rowsBetween <- aaply(logicalBetweenIVMatch, 1, all)
        } else rowsBetween <- rep(T, nrow(DF))

        # When both within and between rows match:
        rows <- rowsWithin & rowsBetween
        values <- DF[DF[,params$IDs] == currentId & DF[,params$blocks] %in% lookupBlocks & rows,]

        # Check the values being used for replacement are not themselves outliers!
        for (replacementidx in row.names(values))
        {
          if (replacementidx %in% idx) warning("Replacement is also an outlier! That also means one of the replacements will be incorrect")
        }

        #cat("Replacement values found:\n")
        #print(values)
        #cat("Those were the values.\n")
        replace.value <- mean(values[,params$DVs[1]]) # response=="meanotherblocks"
        if (response=="meanblocks") replace.value <- mean(c(values[,params$DVs[1]], DF[i,params$DVs[1]]))

        # Check the values being used for replacement are not higher if the outlier is "high" or lower if the outlier is "low"
        #cat("highlow:\n")
        #print(highlow)
        if (!is.null(highlow)) {
          currenthighlow <- highlow[which(i==idx)]
          warn <- F
          if (currenthighlow=="high" && replace.value > DF[i, params$DVs[1]]) warn <- T
          if (currenthighlow=="low" && replace.value < DF[i, params$DVs[1]]) warn <- T
          if (warn) warning(paste0("Replacing outlier with further outlying value! idx: ", idx, ", original value: ", DF[i, params$DVs[1]], ", replacement: ", replace.value))
        }

        #DF[DF[,params$IDs] == currentId & DF[,params$blocks] == currentBlock & rows,params$DVs[1]] <- replace.value
        cat(paste("Original value:", DF[i, params$DVs[1]], "replaced with:", replace.value, "\n"))
        DF[i, params$DVs[1]] <- replace.value
      }
    }

    if (response=="NA")
    {
      DF[idx,params$DVs[1]] <- NA
      responded <- T
    }

    if (response=="convert")
    {
      # Check we've got a response.param
      if (is.null(response.param) || !is.numeric(response.param))
        stop("sigeng: outliers.respond() response is convert but missing response param")

      responded <- T

      # Get mean and sd
      m <- mean(DF[,params$DVs[1]])
      s <- sd(DF[,params$DVs[1]])
      for (i in idx)
      {
        cat("Outlier row:\n")
        print(DF[i,])
        cat(paste("Original value:", DF[i, params$DVs[1]]))
        # first is it a high or low outlier?
        if (DF[i, params$DVs[1]] < m)
          DF[i, params$DVs[1]] <- m - response.param * s
        else
          DF[i, params$DVs[1]] <- m + response.param * s

        cat(paste(" replaced with:", DF[i, params$DVs[1]], "\n"))
      }
    }

    if (!responded) stop(paste("sigeng: outliers.respond() response not found:", response))
  }else{
    cat("No outliers found.\n")
  }

  return(DF)
}

outliers <- function(params, method="classic", method.param=2, split.blocks=F, response="meanotherblocks", response.param=NULL)
# requires DF, DVs, blocks, IDs, withinIVs, betweenIVs
# method is either "classic", "boxplot", or TODO "ggboxplot", "out" from WRS
# method.param is currently just for classic, and corresponds to the number of SDs to consider it an outlier
# split.blocks: either check for outliers with all blocks together (pretend they are all independent observations)
#             or set split.blocks=T and check for each block separately
# response is either "meanotherblocks" (use the mean score from all other blocks)
#                    "meanblocks" (use the mean score of all blocks including the outlier)
#                    "removeblock" (remove that row)
#                    "removeparticipant" (remove whole participant i.e. maybe > 1 row)
#                    "NA" (set to NA)
#                    (from Field, 5.7.1, p. 153):
#                    "convert" with response.param number of SDs (= Z score) e.g. convert to mean + 2SDs
#                  
{
  # Error check (useful if calling from outside sigeng() )
  if (is.null(params$DF) || is.null(params$DVs) || is.null(params$IDs))
    stop("sigeng: outliers() missing DF, DVs or IDs")

  # this is done in outliers.respond: if (length(params$DVs) > 1) stop("sigeng: outliers() only deals with 1 DV for now")
  #if (!split.blocks)
  #{
    params$DF <- outliers.respond(params, split.blocks, method, method.param, response, response.param)
    # TODO after responding with mean of other blocks or of all blocks, check new point is not an outlier!
  #}

  #if (response=="meanotherblocks" || response=="meanblocks") cat("NOTE: Look at the output. Check the response did not replace an outlier by an outlier from another block!\n")

  return(params$DF)
}

psignifit.summariseFor <- function(DF, idvar="id", intensityvar="x", responsevar="y", condvar=NULL, outfile="summary.txt", outfile.overall="summary-overall.txt", auto.write=T)
  # Usage e.g.:
  # mysum <- psignifit.summariseFor(...)
  # write.table(mysum, file="summary.txt", row.names=F)
  # Write auto.write==T, we will automatically write the summary to outfile, but first we check it exists and ask whether to overwrite it.
{
  ## We want x, y, n for each participant
  ## intensity, number of positive response (for a yes/no task) or no. correct responses (nAFC), number of trials
  
  # cover possibilities of a numeric or factor id
  if (is.factor(DF[,idvar])) ids <- levels(DF[,idvar])
  else ids <- unique(DF[,idvar])
  
  xs <- levels(DF[,intensityvar])

  #nocond <- F
  if (is.null(condvar)) {
    #nocond <- T
    condvar <- "cond"
    DF[,condvar] <- factor(1) # set up a dummy cond var, so it can be easily processed by psignifit script
  }else{
    if (!is.factor(DF[,condvar])) stop("Error, condvar must be a factor")
  }
  
  res <- NULL
  
  for (id in ids) {
    DFid <- DF[DF$id==id,]
    
    for (cond in levels(DF[,condvar])) {
      DFidcond <- DFid[DFid[,condvar]==cond,]
      ys <- NULL
      ns <- NULL
      for (x in xs) {
        DFidcondx <- DFidcond[DFidcond[,intensityvar]==x,]
        ys <- c(ys, sum(DFidcondx[,responsevar]))
        #ys <- c(ys, sum(DFidcondx$correct))
        ns <- c(ns, nrow(DFidcondx))
      }
      residcond <- data.frame(id=id, cond=cond, x=xs, y=ys, n=ns)
      residcond <- rename(residcond, c("cond"=condvar))
      res <- rbind(res, residcond)
    }
  }

  # leave the cond, so it can be easily processed by the psignifit script
  #if (nocond) res[,condvar] <- NULL

  # We have the mean for each participant. Now calculate the variance for each condition and intensity,
  # and the standard error
  # We could here just create a new variable in the data frame called value, which is y/n, rather than
  # calculate it each loop. That's what we end up doing below anyway to call calcadj.
  # So this works but could be neater, but I am far too exhausted/disinterested to work on this right now.
  ressum <- NULL
  for (cond in levels(DF[,condvar])) {
    res.cond <- res[res[,condvar]==cond,]
    for (x in xs) {
      res.condx <- res.cond[res.cond$x==x,]

      # we have the data. get the mean, sd, std.error
      my.mean <- mean(res.condx$y/res.condx$n)
      my.sd <- sd(res.condx$y/res.condx$n)
      my.std.error <- std.error(res.condx$y/res.condx$n)
      
      temp <- data.frame(cond=cond, x=x, mean=my.mean, sd=my.sd, std.error=my.std.error)
      ressum <- rbind(ressum, temp)
      #res[res[,condvar]==cond & res$x==x,]$all.mean <- my.mean
      #res[res[,condvar]==cond & res$x==x,]$all.sd <- my.sd
      #res[res[,condvar]==cond & res$x==x,]$all.std.error <- my.std.error
    }
  }

  # now with adjustment
  # Note I've left the mean.adj in the adjusted values even though it should be the same
  # This allows for a quick check (e.g. it could be different if there are between subject elements)
  #print(res[,condvar])
  #print(head(res))
  res.adj <- res
  res.adj$variable <- paste(res.adj[,condvar], res.adj$x, sep=".")
  res.adj$variable <- factor(res.adj$variable)
  res.adj$value <- res.adj$y/res.adj$n
  #print(res.adj$id)
  res.adj <- calcadj(res.adj)
  #print(res.adj)
  ressum$mean.adj <- -10000
  ressum$std.error.adj <- -10000
  for (cond in levels(DF[,condvar])) {
    for (x in xs) {
      ressum[ressum$cond==cond & ressum$x==x,]$mean.adj <- res.adj[res.adj$variable==paste(cond, x, sep="."),]$means
      ressum[ressum$cond==cond & ressum$x==x,]$std.error.adj <- res.adj[res.adj$variable==paste(cond, x, sep="."),]$stderrs
    }
  }

  #ressum <- rename(ressum, c("cond"=condvar))

  # Calculate for all participants in main psignifit table, allow for fitting/plotting overall data

  for (cond in levels(DF[,condvar])) {
    DFcond <- DF[DF[,condvar]==cond,]
    ys <- NULL
    ns <- NULL
    for (x in xs) {
      DFcondx <- DFcond[DFcond[,intensityvar]==x,]
      ys <- c(ys, sum(DFcondx[,responsevar]))
      ns <- c(ns, nrow(DFcondx))
    }
    residcond <- data.frame(id="all", cond=cond, x=xs, y=ys, n=ns)
    residcond <- rename(residcond, c("cond"=condvar))
    res <- rbind(res, residcond)
  }
  
  if (auto.write) {
    if (file.overwrite.prompt(outfile)) write.table(res, outfile, row.names=F)
    if (file.overwrite.prompt(outfile.overall)) write.table(ressum, outfile.overall, row.names=F)
  }

  return (list(psignifit=res, overall=ressum))
}

psignifit.run <- function(psignifit="~/ollik-home/utils/psignifit_3.0_beta.20120611.1/tests/dodds/analysesummary.py",
                          summaryfilein="summary.txt",
                          threshslopesfileout="rectwidth-data/threshslopes.txt",
                          nafc=1,
                          bootstrapbayes="bootstrap",
                          gammaislambda=T,
                          constraints="'unconstrained'/'unconstrained'/'Uniform(0,0.1)'") # constrains separated by '/' cos comma used in prior functions
  # Get psignifit to analyse summary file and output thresholds, slopes, and other params
  # Returns false when not run (e.g. output file exists, chose not to overwrite) or true otherwise
{
  # First run e.g. rectwidth-analysis.r, and call mysum <- summariseForPsignifit(...), then
  # write.table(mysum, file="summary.txt", row.names=F); i.e. make summary.txt

  # This will warn before overwriting threshslopesfileout
  if (file.exists(threshslopesfileout)) {
    cat("File exists:", threshslopesfileout, "\n")
    cat("Overwrite?\n")
    response <- readLines(stdin(), 1)
    if (response=="y") {
      cat("Overwriting.\n")
    }else{
      cat("Not overwriting.\n")
      return (list(run=F, retval=0)) # Not run
      #stop("Stopping due to threshslopes file existing.")
    }
  }
  
  
  cat("Analysing with psignifit...\n")
  
  # Then analyse it with psignifit
  if (gammaislambda) gammaislambda <- 'T'
  else gammaislambda <- 'F'
  retval <- system(paste(psignifit, summaryfilein, threshslopesfileout, nafc, bootstrapbayes, gammaislambda, constraints, sep=" "))
  cat("Retval:")
  print(retval)
  
  cat("Finished psignifit\n")

  return (list(run=T, retval=retval)) # run, value returned
}

psignifit.readthreshslope <- function(threshslopesfile="rectwidth-data/threshslopes.txt", idvar="id", condvar=NULL)
  # Get data frame of threshslopes file (the file output from psignifit analysis script), see psignifit.run
  # Set idvar or condvar to null to keep original values
{
  
  #######
  # Read it in
  #rm(list=ls())
  DF <- data.frame(read.table(threshslopesfile, header=T))

  if (is.null(idvar)) idvar <- "pid"
  if (is.null(condvar)) condvar <- "cond"

  # This actually works with null values but I had to set idvar and condvar anyway for the code below
  DF <- rename(DF, c("pid"=idvar, "cond"=condvar))

  DF[,idvar] <- factor(DF[,idvar])
  DF[,condvar] <- factor(DF[,condvar])
  #str(DF)

  return(DF)
}

psignifit.predict <- function(x, alpha, beta, gamma, lambda)
  # Predict the resulting function based on parameter estimates
  # This function corresponds to the logistic function (Hill 2001)
{
  y <- 1/(1+exp(-(x-alpha)/beta))
  y <- gamma + (1 - gamma - lambda) * y
  cat("gamma: ", gamma, ", lambda: ", lambda, "\n")

  return(y)
}

psignifit.threshold <- function(f=0.5, alpha, beta)
# This function calculates Tf, and corresponds to the function for calculating the threshold, see Appendix of Hill (2001)
# f is the threshold requested e.g. 0.5
{
  return (alpha - beta * log(1/f-1))
}

psignifit.plot <- function(DFthresh, DFpsi, DFoverall, ids=NULL, conds=NULL, idvar="id", condvar=NULL, nafc=1, xmin=NULL,
                           xlabel="Intensity", xticbreaks="auto", xticlabs="auto", colourlabel=NULL, titleprefix=NULL, title="default", threshold=T, nothreshleft=F, nothreshdown=F)
  # Parameters:
  # DFthresh from psignifit.readthreshslope
  # DFpsi from DFsum$psignifit from DFsum <- psignifit.summariseFor
  # DFoverall from DFsum$overall
  # Plot, with only ids or conds.
  # When null, use them all.
  # threshold is T/F, display threshold line
  # nothreshleft is T/F, hide horizontal part of line
  # nothreshdown is T/F, hide vertical part of line
{

  if (is.null(ids))   { if (is.factor(DFthresh[,idvar]))   ids   <- levels(DFthresh[,idvar])   else ids   <- DFthresh[,idvar] }
  if (is.null(conds)) { if (is.factor(DFthresh[,condvar])) conds <- levels(DFthresh[,condvar]) else conds <- DFthresh[,condvar] }

  if (is.null(colourlabel)) colourlabel <- condvar

  for (i in ids) {

    p <- ggplot()
    if (xticbreaks != "auto" && xticlabs != "auto") p <- p + scale_x_continuous(breaks=xticbreaks, labels=xticlabs)

    for (cond in conds) {
      #DFidcond <- DF[DF[,idvar]==i & DF[,condvar]==cond,]
      DFthreshidcond <- DFthresh[DFthresh[,idvar]==i & DFthresh[,condvar]==cond,]
      DFpsiidcond <- DFpsi[DFpsi[,idvar]==i & DFpsi[,condvar]==cond,]

      temp.x <- remove.factor(DFpsiidcond$x, convert="numeric")
      temp.xmin <- min(temp.x)
      if (!is.null(xmin)) temp.xmin <- xmin # allow user to override, e.g. set to 0
      my.model <- data.frame(x=temp.xmin:max(temp.x))
      my.model$y <- psignifit.predict(my.model$x, DFthreshidcond$alpha, DFthreshidcond$beta, DFthreshidcond$gamma, DFthreshidcond$lambda)

     

      colnames(DFpsiidcond)[colnames(DFpsiidcond)==condvar] <- "cond"
      my.model$cond <- cond
      #print(head(DFpsiidcond))
      #print(head(my.model))

      

      if (i == "all") {
        # use now the DFoverall for stderror etc.

        p <- p + geom_point(data=DFpsiidcond, mapping=aes(x=remove.factor(x, "numeric"), y=y/n, colour=cond)) +
             geom_errorbar(data=DFoverall, mapping=aes(x=remove.factor(x, "numeric"), y=mean, ymin=mean-std.error.adj, ymax=mean+std.error.adj, colour=cond, width=0.9)) +
             geom_line(data=my.model, mapping=aes(x=x, y=y, colour=cond))
      }else{
        p <- p + geom_point(data=DFpsiidcond, mapping=aes(x=remove.factor(x, "numeric"), y=y/n, colour=cond)) +
             geom_line(data=my.model, mapping=aes(x=x, y=y, colour=cond))
      }

      if (nothreshleft && nothreshdown) threshold <- F
      if (threshold) {
        # predict the y value of the threshold
        # I thought psignifit is returning the threshold based on the resulting y value, i.e. including the guesses and lapses!
        # So I calculated it myself based on the fit from psignifit
        # Then I found it's the same - but it makes sense that it results as 0.5 because gammaislambda!!
        threshold.x <- psignifit.threshold(0.5, DFthreshidcond$alpha, DFthreshidcond$beta) # same as DFthreshidcond$thresh_0.5
        print(threshold.x)
        threshold.y <- psignifit.predict(threshold.x, DFthreshidcond$alpha, DFthreshidcond$beta, DFthreshidcond$gamma, DFthreshidcond$lambda)
        print(threshold.y)
        # (min(xs), ty) to end (tx, ty)
        # (tx, ty) to end (tx, lowbound)
        lowbound <- 0
        if (nafc > 1) lowbound <- 1 / nafc
        minxs <- min(remove.factor(DFpsi$x, "numeric"))
        DFthresh.temp <- data.frame(x=c(minxs, threshold.x), xend=rep(threshold.x, 2), y=rep(threshold.y, 2), yend=c(threshold.y, lowbound), cond=cond)
        if (nothreshleft) DFthresh.temp <- DFthresh.temp[2,]
        if (nothreshdown) DFthresh.temp <- DFthresh.temp[1,]

        p <- p + geom_segment(data=DFthresh.temp, mapping=aes(x=x, xend=xend, y=y, yend=yend, colour=cond), linetype=2)
      }

    }
    if (is.null(title)) {
      p <- p + labs(x=xlabel, colour=colourlabel)
    }else{
      if (title == "default") {
        if (is.null(titleprefix)) title.current <- paste("Participant", i, sep=" ")
        else title.current <- paste(titleprefix, "Participant", i, sep=" ")
      }else title.current <- title
      p <- p + ggtitle(title.current) + labs(x=xlabel, colour=colourlabel)
    }
    print(p)
  }

}


# Adapted from: http://www.r-bloggers.com/visual-debugging-with-rstudio/
debugonce.fun <- function(f){ 
  fname = deparse(substitute(f))
  debugfile <- tempfile("debug-", ".", fileext=".r")
  dump(fname, file = debugfile) 
  source(debugfile)
  cat("Debug file:", debugfile, "\n")
  do.call("debugonce", args = list(fname), envir = globalenv()) 
  invisible(NULL) 
}

