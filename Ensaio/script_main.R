### I: load libraries: #########################################################
rm(list = ls())

install <- F
load.pkg <- T
runLongFunctions <- F

# libs <- c("readxl", "dplyr", "plyr", "zoo", "svMisc", "data.table", "progress",
#           "bigmemory", "ff")
# 
# libs <- c("readxl", "dplyr", "plyr", "zoo", "svMisc", "vars", "Jmisc",
#           "dlm", "ggplot2", "ggfortify", "parallel", "progress", "xts", "FarmSelect",
#           "POET", "nFactors", "devtools", "copulaedas")


libs <- c("readxl", "dplyr", "plyr", "zoo", "svMisc", "vars", "Jmisc",
          "dlm", "ggplot2", "ggfortify", "progress", "xts", "FarmSelect",
          "POET", "nFactors", "devtools", "copulaedas", "forecast", "astsa",
          "glmnet", "plotly", "lubridate", "caret")

# invert order of libs, to avoid masking functions of more important 
# libraries (here: plotly's select was messing around with dplyr's) 
libs <- libs[length(libs):1]

if (install == T) {
  new.pkgs <- c()
  for (i in 1:length(libs)) {
    if (libs[i] %in% installed.packages() == F) {
      new.pkgs[i] <- libs[i]
      if (is.null(new.pkgs) == T) {
        cat("no packages needed to be installed", "\n")
      } else if (!is.null(new.pkgs) == T){
        cat("some packages need to be installed:", length(new.pkgs), "\n")
      }
    }
  }
  install.packages(new.pkgs, dependencies = T)
  if (load.pkg == T) {
    for (i in 1:length(libs)) {
      library(libs[i]) 
    }
  }
} else if (install == F) {
  if (load.pkg == T) {
    for (i in 1:length(libs)) {
      sapply(libs, library, character.only = T) 
    }
  }
}
### A: user selection: #########################################################
### A.1: set working directory: ################################################

# save standard par() parameters:
oldParSettings <- par()

memLimit.old <- memory.limit()
memLimit.old

memSize.old <- memory.size()
memSize.old

### A.2: define fundamental functions ##############################################

# create 'not in' function:
'%!in%' <- function(x,y)!('%in%'(x,y))

# function for selecting the number of factors (Ahn and Horenstein 2013, Econometrica)
ahnhorenstein <- function(data = "", K.factors = NULL, extractFactors=FALSE) {
  n <- ncol(data) #dimension
  p <- nrow(data) #sample size
  
  if (extractFactors == FALSE) {
    values <- extract(data, K = ncol(data))$eval
  
    values = pmax(values, 0) #values
  } else {
    values <- data
  }
  
  ratio = c()  
  K.factors <- if (is.null(K.factors)) {
    for(i in 1:(floor(min(n,p)/2))){
      ratio = append(ratio, values[i+1]/values[i])
    }
    ratio = ratio[is.finite(ratio)]
    K.factors = which.min(ratio)} else {
      for (i in 1:K.factors) {
        ratio = append(ratio, values[i+1]/values[i])
      }
      ratio = ratio[is.finite(ratio)]
      K.factors = which.min(ratio)
    }
  
  return(K.factors)
  cat("Number of factors with eigenvalue-ratio test (Ahn and Hohrenstein, 2013, Econometrica): ", K.factors, "\n")
}

# lag data:
lagdata <-function(y, lags, intercept = FALSE){
  T <- nrow(y)
  K <- ncol(y)
  obs <- T-lags
  x  <- embed(y, dimension = lags + 1)[,-(1:K)]
  if(intercept == TRUE){
    x <- cbind(1, x)
  }
  yi <- y[(lags+1):T,]
  return(list(y = yi, x = x, obs = obs, T = T, K = K));
}

# exctract factors (Stock and Watson 2002, JASA):
extract <- function(x, K, compMarR2 = T){
  f_call <- match.call()
  n <- ncol(x)
  x <- as.matrix(x)
  x.x <- t(x)%*%x
  
  evectors <- eigen(x.x)$vectors
  evalues <- eigen(x.x)$values
  
  varexplained <- evalues/sum(evalues)*100
  cumvarexplained <- cumsum(varexplained)
  
  ret.evectors <- sqrt(n)*evectors[, 1:K] #loadings
  fac <- x%*%ret.evectors/n #factors
  error = x - (fac %*% t(ret.evectors)) #error
  
  if (compMarR2 == T) {
    print("Compute marginal R-squared: TRUE")
    
    dataReg <- as.data.frame(cbind(fac, x))
    fhat <- c(paste0("fhat", seq(1, K, by = 1)))
    colnames(dataReg) <- c(fhat, colnames(x))
    dataFhat <- matrix(NA, nrow = K, ncol = ncol(dataReg))
    colnames(dataFhat) <- colnames(dataReg)
    
    pbar <- progress_bar$new(format = " :spin estimation process [:bar] :percent estimated: :eta",
                             total = ncol(dataFhat),
                             clear = FALSE,
                             width = 60)
    
    for (i in 1:ncol(dataFhat)) {
      pbar$tick()
      for (j in 1:K) {
        linmod.fhat <- lm(data = dataReg,
                          dataReg[, which(colnames(dataReg) == fhat[j])] ~ dataReg[, i])
        dataFhat[j,i] <- summary(linmod.fhat)$r.squared*100
      }
    }
    
    
  }
  
  return(list(lam = ret.evectors,
              fac = fac,
              err = error,
              cumvar = cumvarexplained,
              eval = evalues,
              marR2 = dataFhat,
              K = K,
              f_call = f_call))
}

transData <- function(x, standardize = T){
  
  data <- x
  for (i in seq(2, ncol(data), by = 1)) {
    if (tcodeuse[1,i] == 1) { #no transformation
      cat("tcode: 1 -->", colnames(data[i]), "\n")
      next()
    } else if (tcodeuse[1,i] == 2) { #first-difference
      cat("tcode: 2 -->", colnames(data[i]), "\n")
      data[,i] <- c(NA, diff(data[,i], lag = 1, differences = 1))
    } else if (tcodeuse[1,i] == 3) { #second-difference
      cat("tcode: 3 -->", colnames(data[i]), "\n")
      data[,i] <- c(NA, diff(data[,i], lag = 1, differences = 2))
    } else if (tcodeuse[1,i] == 4) { #log-transform
      cat("tcode: 4 -->", colnames(data[i]), "\n")
      data[,i] <- log10(data[,i])
    } else if (tcodeuse[1,i] == 5) { #first log-differences
      cat("tcode: 5 -->", colnames(data[i]), "\n")
      data[,i] <- c(NA, diff(log10(data[,i]), lag = 1, differences = 1))
    } else if (tcodeuse[1,i] == 6) { #second log-differences
      cat("tcode: 6 -->", colnames(data[i]), "\n")
      data[,i] <- c(NA, NA, diff(log10(data[,i]), lag = 1, differences = 2))
    } else if (tcodeuse[1,i] == 7) { #(x_t / x_t-1 - 1.0)
      cat("tcode: 7 -->", colnames(data[i]), "\n")
      data[,i] <- (c(diff(data[,i], lag = 1, differences = 1), NA) / lag(data[,i], lag = 1)) - 1
      data[,i] <- c(NA, data[-length(data[,i]), i])
    }
  }
  
  if (standardize == T) {
    print("standardize data: TRUE")
    # standardize stationary data:
    for (i in seq(2, ncol(data), by = 1)) {
      data[,i] <- scale(data[,i], center = T, scale = T)
    }
  } else if (standardize == F) {
    print("standardize data: FALSE")
  }
  return(data)
}

estFac3d <- function(data = ..., from = forecastStart, to = forecastEnd, K = 20, standardize = T) { #start function
  
  R <- array(NA, dim = c(length(seq.Date(from = as.Date(from), to = as.Date(to), by = "month")), ncol(data)-1, K))
  
  if (standardize == T) { #start if I
    dateTime2 <- c(seq.Date(from = as.Date(from), to = as.Date(to), by = "month"))
    for (t in dateTime2) { #start for-loop I
      print(as.Date(t))
      ndata <- data[which(data[, "sasdate"] <= t),]
      ndata <- scale(ndata[, which(colnames(data) != "sasdate")], center = T, scale = T)
      fac <- extract(as.matrix(ndata), K = K)
      marR2 <- fac$marR2[, -(1:K)]
      
      tt <- which(dateTime2 == t)
      for (kk in 1:K) { #start for-loop I.1
        R[tt, 1:ncol(ndata), kk] <- marR2[kk,]
      } #end for-loop I.1
      
    } #end for-loop I
    return(R)
  } # end if I
} #end function

plotFac3d <- function(data = factors3d[,,], K = 1){
  plot_ly(z = ~data[,,K],
          colorbar = list( title = ""),
          colors = grey(level = c(0.2, 0.5, 0.8), alpha = 1), linetypes = 5) %>% 
    add_surface() %>%
    layout(title = paste("Marginal R-squares for Factor", K),
           scene = list(
             xaxis = list(title = paste0("Variables (1 - ", ncol(data[,,K]), ")"),
                          gridcolor = "rgb(255, 255, 255)",
                          zerolinecolor = "rgb(255, 255, 255)",
                          showbackground = TRUE,
                          backgroundcolor = "rgb(240, 240, 240)",
                          autorange = "reversed",
                          showticklabels = F),
             yaxis = list(title = paste0("Time (",
                                         substr(forecastStart,1,4), "/M", substr(forecastStart, 9,10), " - ",
                                         substr(forecastEnd,1,4), "/M", substr(forecastEnd, 9,10), ")"),
                          gridcolor = "rgb(255, 255, 255)",
                          zerolinecolor = "rgb(255, 255, 255)",
                          showbackground = TRUE,
                          backgroundcolor = "rgb(230, 230, 230)",
                          showticklabels = F),
             zaxis = list(title = "R-squared",
                          range = c(0,100),
                          gridcolor = "rgb(255, 255, 255)",
                          zerolinecolor = "rgb(255, 255, 255)",
                          showbackground = TRUE,
                          backgroundcolor = "rgb(220, 220, 220)")
           )
    ) 
}
# date to split = 2000,
# horizon = 2002
splitTrainTestData <- function(data, dateTimeTravel, forecastHorizon=60,
                               verbose=TRUE, .yvar=yvar) {
  # Returns:
  #   A list with 4 data.frames. x.train, y.train, x.test, y.test. 
  #
  # Args:
  #   dateTimeTravel
  #       cut the data with information at time t,
  #   forecastHorizon
  #       how long should be the test horizon (in months) looking back from
  #       time t (dateTimeTravel)
  #   .yvar
  #       which variable to be modelled 
  
  # make sure dateTimetravel is in date format
  if (!is.Date(dateTimeTravel)) {
    dateTimeTravel <- as.Date(dateTimeTravel)
  }
  
  # first cut data up to "time travel" period. 
  data <- data %>% 
    filter(sasdate <= dateTimeTravel)
  
  # split remaining into 2 sections. Begin til split; split til dateTimeTravel
  splitDate <- dateTimeTravel %m-% months(forecastHorizon)
  startDate <- head(data$sasdate, 1)
  endDate <- tail(data$sasdate, 1)
  
  ## split train und test (read %>% as in "then")
  x.train <- data %>% 
    filter(sasdate <= splitDate) %>% select(-.yvar)
  y.train <- data %>% 
    filter(sasdate <= splitDate) %>% select(yvar) %>% pull()  # pull: extract
  x.test <- data %>%                                            # as numeric 
    filter(sasdate > splitDate) %>% select(-.yvar)
  y.test <- data %>% 
    filter(sasdate > splitDate) %>% select(.yvar) %>% pull()  # pull: extract
                                                                # as numeric 
  if (nrow(x.train) != length(y.train)) {
    warning("check dimensions of training and test")
  } 
  if (isTRUE(verbose)) {
    message(paste0("Splitting Dates...\n",
                   "Start: ", startDate,
                   " || Split: ", splitDate,
                   " || End: ", endDate, 
                   " || Horizon: ", forecastHorizon, "Ms"))
    }
  return(list(x.train = x.train,
              x.test = x.test, 
              y.train = y.train,
              y.test = y.test, 
              metaSplits = c(startDate = startDate,
                             splitDate = splitDate,
                             endDate   = endDate)
              )
         )
}

alphaSearchFun <- function(splittedData, alphaSearchMin = 0.1, alphaSearchMax = 0.9, 
                           alphaSearchBy = 0.01, lambdaseq = NULL, verbose=FALSE) {

  ## alpha search
  if (verbose) {
    message(paste("Alpha Search: from ", alphaSearchMin, "to ",
                  alphaSearchMax , "by ", alphaSearchBy))  
  }

  list.of.fits <- list()
  results <- data.frame()
  
  for(i in seq(alphaSearchMin, alphaSearchMax, by = alphaSearchBy)) {
    fit.name <- paste0(("alpha"), i)
    list.of.fits[[fit.name]] <- glmnet(data.matrix(splittedData$x.train[,-1]),
                                       splittedData$y.train,
                                       alpha = i,
                                       lambda = lambdaseq,
                                       family = "gaussian")

    # predicted <- predict(list.of.fits[[fit.name]],
    #                      s = list.of.fits[[fit.name]]$lambda,
    #                      newx = as.matrix(splittedData$x.test[,-1]))
    
    # mse <- mean((as.numeric(splittedData$y.test) - predicted)^2)
    temp <- data.frame(alpha = i,  fit.name = fit.name)
    results <- rbind(results, temp)
  }
  list.of.fits$results <- results
  return(list.of.fits)
}

EN_SelectVars <- function(data, dateTimeTravel = forecastStart, forecastHorizon,
                          yvar = yvar, lambdaseq = NULL,
                          alphaSearchMin = 0.4, alphaSearchMax = 0.6, 
                          alphaSearchBy = 0.01, maxVars = 60,
                          returnOnlyDfr = NULL, verbose=FALSE) {
  call <- sys.call()
  whichData <- call$data
  
  
  
  # split data
  splittedData <- splitTrainTestData(data, dateTimeTravel,
                                    forecastHorizon=forecastHorizon,
                                    verbose=verbose, .yvar=yvar)
  
  list.of.fits <- alphaSearchFun(splittedData, alphaSearchMin, alphaSearchMax, 
                                 alphaSearchBy, lambdaseq, verbose)
  
  # grid search alha and lambda
  l <- list()
  for (i in seq(alphaSearchMin, alphaSearchMax, by = alphaSearchBy)) {
    fit.name <- paste0(("alpha"), i)
    assign(paste0("predtest_mse", i), 
           predict(list.of.fits[[fit.name]], newx = as.matrix(splittedData$x.test[,-1])))
    xx <- eval(parse(text = paste0("predtest_mse", i)))
    sqerrors <- ((splittedData$y.test) - xx)^2
    mse <- colMeans(sqerrors)
    mse <- c(mse, rep(NA, 100 - length(mse))) #just to make all same dimension
    #assign(paste0("predtest_mse", i), mse) 
    l[[fit.name]] <- mse
  }
  
  # transform to list to data frame
  asDF <- data.frame(l)
  # get location of min MSE among all alphas and lambdas
  loc<-which(asDF == min(asDF, na.rm = TRUE), arr.ind = TRUE)
  # extract name of best fit
  bestfitname <- colnames(asDF[loc[2]])
  # extract best fit of list of fits
  bestfit <- list.of.fits[[bestfitname]]
  # get MSE of best fit, lambda and alpha
  bestMSE <- asDF[loc]
  bestLambda <- bestfit$lambda[loc[1]]
  bestAlpha <- list.of.fits$results$alpha[loc[2]]
  
  # extract the Coefs based on the above selected best lambda
  coefsOfBestLambda <- bestfit$beta[, bestfit$lambda == bestLambda]
  nrOfNonNullCoefs <- sum(abs(coefsOfBestLambda) > 0)
  
  # get colnames of all variables with non-zero coeficients estimates
  resVarsSortedAbs <- sort(abs(coefsOfBestLambda[abs(coefsOfBestLambda) > 0 ]),
                           decreasing = T)
  resVarsNames <- names(resVarsSortedAbs)
  resVarsNamesClean <- gsub("\\..*$","", resVarsNames)
  resVarsNamesUniq <- unique(resVarsNamesClean)
  numberOfUniqVars <- length(resVarsNamesUniq)
  # generate vector of NON duplicatse (true,false)
  noDups <- !duplicated(resVarsNamesClean)
  resVarsSortedAbsUniq <- resVarsSortedAbs[noDups]
  names(resVarsSortedAbsUniq) <- gsub("\\..*$","", names(resVarsSortedAbsUniq))
  
  if (verbose) {
    message("Resulting in...")
    message(paste("alpha: ", bestAlpha, "\tmse: ", bestMSE, "\tlambda: ", bestLambda))
    message("Nr of selected Vars: ",length(resVarsNames),
            "\tNr of unique Vars: ", numberOfUniqVars)
    message(rep("=", 80))
    
  }
  if (length(resVarsNamesUniq) > maxVars) {
    # if more than maxVars, restrict to only the 
    #resVarsSortedAbs  <- head(resVarsSortedAbs, maxVars)
    resVarsNames      <- head(resVarsNames, maxVars)
    resVarsNamesClean <- head(resVarsNamesClean, maxVars)
    resVarsNamesUniq  <- head(resVarsNamesUniq, maxVars)
    resVarsSortedAbsUniq <- head(resVarsSortedAbsUniq, maxVars)
  }
  
  outData <- data[, names(resVarsSortedAbs)]
  # if only interested in final DF with relevant vars ... 
  if (isTRUE(returnOnlyDfr)) {
    return(outData)
  }
  # ... otherwise return extra info as well in a list
  meta <- list(whichData = whichData,
               splitDates = splittedData$metaSplits,
               dimOutData = dim(outData),
               resVarsSortedAbs = resVarsSortedAbs,
               resVarsSortedAbsUniq = resVarsSortedAbsUniq,
               resVarsNames = resVarsNames,
               resVarsNamesUniq = resVarsNamesUniq,
               resVarsNamesClean = resVarsNamesClean,
               numberOfUniqVars = numberOfUniqVars,
               maxVars = maxVars,
               coefsOfBestLambda = coefsOfBestLambda,  # coefsOfBestLambda inclusive 0 valued coefs!
               nrOfNonNullCoefs = nrOfNonNullCoefs,
               bestMSE = bestMSE, 
               bestLambda = bestLambda,
               bestAlpha = bestAlpha)
            
  out <- list(outData = outData, meta = meta, bestfit = bestfit)
  class(out) <- "EN_Select"
  return(out)
}


# create generic barplot function for EN_SElect object
barplot.EN_Select <- function(x, unique=FALSE, ...) {
  leg = paste0(
    "start date\t", "split date\t\t", " end date\n", 
    x$meta$splitDates[1],"\t", x$meta$splitDates[2], "\t",x$meta$splitDates[3],
    "\n# of selected vars:  ", x$meta$nrOfNonNullCoefs, 
    "\n# of unique vars:  ", x$meta$numberOfUniqVars,
    "\n max vars:  ", x$meta$maxVars,
    "\n data:  ", x$meta$whichData,
    "\nlambda:  ", format(x$meta$bestLambda, digits=3),
    "\nmse:  ", format(x$meta$bestMSE, digits=3),
    "\nalpha:  ", format(x$meta$bestAlpha, digits=3)
  )
  
  varsToPlotClean <- gsub("\\..*$","", x$meta$resVarsNames)
  groups <- metaDataFrame[6, ][varsToPlotClean]
  xtoPlot <- x$meta$resVarsSortedAbs
  xAxis <- paste0(names(xtoPlot), " (", groups, ")")
  color = as.character(metaDataFrame[metaDataFrame$sasdate == "color", ][varsToPlotClean])
  
  #uniqueVarstoPlot <- x$meta$resVarsNamesUniq %in% names(xtoPlot) ### PROBLEM!!!
  
  
  if (isTRUE(unique)) {
    warning("unique is in beta still!!! Check for repeated vars!!!")
    xtoPlotUniq <- x$meta$resVarsSortedAbsUniq
    noDups <- !duplicated(x$meta$resVarsNamesClean)
    colorUniq = color[noDups]
    
    varsToPlotCleanUniq <- gsub("\\..*$","", names(x$meta$resVarsSortedAbsUniq))
    groupsUniq <- metaDataFrame[6, ][varsToPlotCleanUniq]
    xAxisUniq <- paste0(names(x$meta$coefsOfBestLambda[x$meta$resVarsNamesUniq]), " (", groupsUniq, ")")
    
    if (any(duplicated(names(xtoPlotUniq)))) {
      stop("duplicated in unique, there is something wrong.")
    }
    
    main = "Selected Variables (Unique)"
    invisible(barplot(main = main, xtoPlotUniq, names.arg = xAxisUniq,
                             col = colorUniq, las=2, cex.axis = .60, cex.names = .6,
                             cex.main= .7, ...))
    return(mtext(leg, adj = 1, side=4, cex = .7, las=1, padj = -.9, las=2, line=-1.6))
  }
  
  main = "Selected Variables (Duplicates)"
  par(mar=c(8, 4, 4, 2) + 0.1 )
  invisible(barplot(xtoPlot, main = main, names.arg = xAxis, col = color, las=2,
                    cex.axis = .60, cex.names = .6, cex.main= .7, ...))
  mtext(leg, adj = 1, side=4, cex = .7, las=1, padj = -.8, las=2, line=-1.6)
}

genLaggedData <- function(data, lags=12) {
  tmplist <- as.matrix(data[,-1]) %>% lagdata(lags = lags)
  x <- as.data.frame(tmplist$x)
  y <- as.data.frame(tmplist$y)
  varnames <- colnames(tmplist$y)
  egrid <- expand.grid(varnames, paste0(".", 1:lags)) 
  expVarNames <- do.call(paste0, as.list(egrid))
  x <- as.data.frame(x)
  names(x) <- expVarNames
  dates <- data$sasdate[-(1:lags)]
  x$sasdate <- dates
  y$sasdate <- dates
  
  dfMerge <- merge(y, x)
}

## Sanity checks!
sanityCheckLagData <- function(ldata) {
  # basically this functions checks if the "var" in t is equal
  # to "var.lag", accounting for the shift in the data 
  # The variable and obs is chosen randomly
  i <- sample(1:125, 1)
  rowNbr <- sample(1:(nrow(ldata)-24),1)
  p <- sample(1:12, 1)
  toCheck <- names(ldata[,-1])[i]
  
  c1 <- ldata %>% select(toCheck) %>% 
    filter(row_number() == rowNbr) %>% as.numeric()
  c2 <- ldata %>% select( paste0(toCheck,".", p) ) %>% 
    filter(row_number() == rowNbr+p) %>% as.numeric()
  
  assertthat::are_equal(c1, c2)
}


### A.3: select parameters: ####################################################

# select transformation codes to be used:
tcodeselec <- "tcodebbe" # 'tcodebbe' for Bernanke, Boivon and Eliasz (2005, QJE) and 'tcode' for McCracken and Ng (2014, JBES)

# select sample size (SW: 1970 - 1998):
customTimeStart <- as.Date("1970-01-01")
customTimeEnd <- as.Date("2018-01-01")

forecastStart <- as.Date("2005-01-01")
forecastEnd <- as.Date("2017-01-01") #end of dataset!

forecastHorizon = 12 * 5

# select variable to be forecasted:
yvar <- "CPIAUCSL" # or: CPIULFSL (alternative)

# select sheet:
sheetSelect <- "data"

### B: load data: ##############################################################

# read csv:
rawdata <- read_xlsx("2019-05.xlsx", sheet = sheetSelect)

data <- rawdata[, -1] # delete first column from import, which has no data!

# drop variables:
dropvars <- matrix(0, nrow = 1, ncol = ncol(data))
colnames(dropvars) <- colnames(data)

for (i in seq(2, ncol(data), by = 1)) { #find variables that are not usable!
  if (data[which(data[,1] == "usage"),i] == 0) {
    cat("drop variable (because timedate) -->", colnames(data[i]), "\n")
    
    dropvars[1, colnames(data[i])] <- 1
  }
}

dropvars <- t(dropvars[, !(colSums(dropvars == 0))])
dropvars <- c(colnames(dropvars))
data <- data[, -(which(colnames(data) %in% dropvars))]

varlabels   <- data[1,]
usedate     <- data[2,]
useavail    <- data[3,]
usage       <- data[4,]
slow        <- data[5,]
block       <- data[6,]
tcode       <- data[7,]
tcodebbe    <- data[8,]

# save metadata in a object to be able to open from other scripts
# generate misc object with plotting arguments
source("./R/genMiscObjForPlots.R")
metaDataFrame <- data[1:8, ]
metaDataFrame <- rbind(metaDataFrame, c("color", misc$colorsCode[misc$blockAsNum]))
saveRDS(metaDataFrame, "./data/metaDataFrameFromXLS.RDS")

data        <- data[-(1:8),]

data <- as.data.frame(sapply(data, as.numeric))
data[, 1] <- seq.Date(from = as.Date(as.Date(as.numeric(data[1, "sasdate"]), origin = "1899-12-30")),
                      length.out = nrow(data), #or: to = as.Date(as.numeric(data[nrow(data), "sasdate"]), origin = "1899-12-30")
                      by = "month")

# drop variables that are not available:
dropvars <- matrix(0, nrow = 1, ncol = ncol(data))
colnames(dropvars) <- colnames(data)

# transform data to be approx. stationary and standardize data (data, odata, sdata):
tcodeuse <- tcodebbe

# check data!
odata <- transData(data, standardize = F) #odata: original stationary data
sdata <- transData(data, standardize = T) #sdata: standardized stationary data

# split sample (datetime, only split after transform due to NA!):
data <- data[c(seq(which(data[, "sasdate"] == customTimeStart),
                   which(data[, "sasdate"] == customTimeEnd), by = 1)), ]

odata <- odata[c(seq(which(odata[, "sasdate"] == customTimeStart),
                     which(odata[, "sasdate"] == customTimeEnd), by = 1)), ]

sdata <- sdata[c(seq(which(sdata[, "sasdate"] == customTimeStart),
                     which(sdata[, "sasdate"] == customTimeEnd), by = 1)), ]


### B: 2 - Generate lagged data ################################################

l_data <- genLaggedData(data)
lodata <- genLaggedData(odata)
lsdata <- genLaggedData(sdata)

checks <- data.frame()
for (i in seq(10)) {
  checks[i,1] <- sanityCheckLagData(l_data)
  checks[i,2] <- sanityCheckLagData(lodata)
  checks[i,3] <- sanityCheckLagData(lsdata)
}

if (all(unlist(checks))) {
  message("All sanity checks ran fine")
  rm(checks)
} else {
  warning("Some sanity check went wrong")
}




### C: extract and visualize factors ###########################################
if (isTRUE(runLongFunctions)) {
  # extract factors (full sample, for exact dates check above!):
  factors <- extract(sdata[, which(colnames(sdata) %!in% c("sasdate", yvar))], K = 5)
  saveRDS(factors, paste0("./data/factors_", Sys.Date() ,".RDS"))
  # bar plot 
  if (!exists("misc")) {
    source("./R/genMiscObjForPlots.R")
  }
  
  for (i in factors$K:1) {
    barplot(factors$marR2[i, -(1:factors$K)],
            main = paste("Marginal R-squares for Factor", i),
            names.arg = misc$tableNames2,
            col = misc$colorsCode[misc$blockAsNum],
            ylim = c(0, 80),
            cex.axis = .7,
            cex.names = .7,
            las = 2)
    box()
  }
}

if (isTRUE(runLongFunctions)) {
  # plot factor stability over time (use odata as data is standardized new at each point in time!):
  factors3d <- estFac3d(data = odata[, which(colnames(odata) != yvar)], K = 5) #estimate factors recursively for each point in time
  for (i in dim(factors3d)[3]:1) {
    persp(z = factors3d[,,i], main = paste("Marginal R-squares for Factor", i),
          zlab = "R-squared",
          xlab = paste0("Time (",
                        substr(forecastStart,1,4), "/M", substr(forecastStart, 9,10), " - ",
                        substr(forecastEnd,1,4), "/M", substr(forecastEnd, 9,10), ")"),
          ylab = paste0("Variables (1 - ", ncol(factors3d[,,i]), ")"),
          theta = 40, phi = 30, expand = 0.5, zlim = c(0,100),
          col = "lightgrey")
  } #3dplot of factors over time (open in new window to see anything due to grid size!)
  
  plotObj <- plotFac3d(data = factors3d[,,], K = 1) #plotly interactive 3dplot (here K chooses the single factor to be plotted!)
  saveRDS(plotObj, file = paste0("./data/plot3D_", Sys.Date(), ".RDS")) # save plot object to be able to
  # check it without running code again
  plotObj
} 
### D. Elastic Net (Simple  Examples) ##########################################

# debugonce(alphaSearchFun)
# debugonce(EN_SelectVars)
lambdaseq <- 10^seq(-7, 1, length.out = 60)
ab <- EN_SelectVars(data = lodata, 
                    dateTimeTravel = forecastStart,
                    forecastHorizon = 24,
                    lambdaseq = lambdaseq,
                    yvar = yvar,
                    alphaSearchMin =  .4,
                    alphaSearchMax = .6,
                    alphaSearchBy = .01,
                    verbose = T)
barplot(ab, unique = F, ylim = c(0, .02))

## test ab heir (simple training example)
if (FALSE) {
  splittedData <- splitTrainTestData(odata, dateTimeTravel = forecastStart)
  xTrain <- ts(splittedData$x.train[,-1],
               start = as.yearmon(customTimeStart), frequency = 12)
  xTest <- ts(splittedData$x.test[,-1],
              start = as.yearmon(forecastStart %m+% months(1)), frequency = 12)
  yTrain <- ts(splittedData$y.train, 
               start = as.yearmon(customTimeStart), frequency = 12)
  yTest <- ts(splittedData$y.test,
              start = as.yearmon(forecastStart %m+% months(1)), frequency = 12)
  
  fit <- glmnet(x = xTrain, y = yTrain)
  predd <- predict(fit, newx = xTest)
  
  predd1 <- ts(predd[,55], start = as.yearmon(forecastStart %m+% months(1)), frequency = 12)
  plot.ts(cbind(predd1, yTest), plot.type = "s", col=c("blue", "red"))
  plot.ts(predd1 - yTest, plot.type = "s")
  abline(h=0)
}


### E: Elastic Net (Looping) ###################################################

check <- TRUE
# generate data.frame for storing the binary data. 
dfBinary <- data %>%
  # TODO: check if > or >= is right!!
  filter(sasdate >  forecastStart, 
         sasdate <= forecastEnd) %>% 
  mutate_if(is.numeric, function(x) {x<-NA})

if (check) {
  message("start: ", head(dfBinary$sasdate, 1))
  message("end:   ", tail(dfBinary$sasdate, 1))
} 

dfBinary <- data.frame(matrix(NA, nrow = ncol(data), ncol = 0))
dfBinary$variables <- names(data)
diffStartEndMonths <- (as.yearmon(forecastEnd) - as.yearmon(forecastStart)) * 12
for (i in seq_len(diffStartEndMonths)) {
  dateTimeTravel <- forecastStart %m+% months(i)
  forecastHorizon <- 12 * 5
  res <- EN_SelectVars(data = lodata,
                       dateTimeTravel = dateTimeTravel,
                       forecastHorizon = forecastHorizon,
                       lambdaseq = lambdaseq,
                       yvar = yvar,
                       alphaSearchMin =  .4,
                       alphaSearchMax = .6,
                       alphaSearchBy = .1,
                       verbose=T, 
                       maxVars = 100)
  barplot(res, unique = T)
  dateAsNum <- as.numeric(dateTimeTravel)
  row <- c(sasdate=dateAsNum,res$meta$resVarsSortedAbsUniq)
  dfBinary[i+1] <- row[dfBinary$variables]
}
dfBinary <- replace(dfBinary, is.na(dfBinary), 0)
# saveRDS(dfBinary, "./data/dfBinary.RDS")
dff <- dfBinary
colnames(dff) <- c("sasdate", as.character(as.Date(as.numeric(dfBinary[1,-1]))))
dff <- dff[-1,]
rownames(dff) <- dff$sasdate
dff <- dff[,-1]

dfrank <- apply(dff, MARGIN = 2, rank, ties.method = "max" )
dfrank <- t(floor(dfrank))
dfrank_ts <- ts(dfrank, start = as.yearmon(head(rownames(dfrank),1)), frequency = 12)
plot(dfrank_ts, plot.type = "single", col = misc$colorsCode[misc$blockAsNum], 
     type="l")
legend("bottomright", legend = unique(misc$blockAsNum),  col = misc$colorsCode,
       lty=1, cex=.5, horiz = T)


dfrank_ts %>% heatmap(Rowv = NA)
dim(dff)
dffBin <- dff > 0
dffBin <- apply(dffBin, as.numeric, MARGIN = 2)
rownames(dffBin) <- rownames(dff)
heatmap(as.matrix(dff),Colv = NA, Rowv = NA,
        keep.dendro=FALSE, cexRow = .5, scale = "c")

pheatmap::pheatmap(dffBin, cluster_cols = F, cluster_rows = F,show_rownames = T,
                   show_colnames = T, fontsize = 5, legend = T)
####
pheatmap::pheatmap(dff, cluster_cols = F, cluster_rows = F,show_rownames = T,
                   show_colnames = T, fontsize = 5, legend = T)


