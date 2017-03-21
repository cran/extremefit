
######################################################################################
# Copyright Ion Grama (2008)
# This package was used in the paper Grama and Spokoiny Annal of Statist (2008).
#
#
# we introduce two functions hill and  hill.adapt
# and some related functions to compute quantiles
# rm(list = ls())
######################################################################################

######################################################################################
## compute the Hill estimator as a process

#' Hill estimator
#'
#' @description Compute the weighted Hill estimator.
#'
#' @param X a vector of data.
#' @param weights a vector of weights assiociated to \eqn{x}.
#' @param grid a vector of values for which the Hill estimator is computed.
#'
#' @details Compute the weighted Hill estimator for vectors \eqn{grid}, data and weights (see references below).
#' @return
#' \item{xsort}{the sorted data.}
#' \item{wsort}{the weights assiociated to \eqn{xsort}.}
#' \item{grid}{the grid for which the Hill estimator is computed.}
#' \item{hill}{the Hill estimators.}
#'
#' @references
#' Grama, I. and Spokoiny, V. (2008). Statistics of extremes by oracle estimation. Ann. of Statist., 36, 1619-1648.
#'
#' Durrieu, G. and Grama, I. and Pham, Q. and Tricot, J.- M (2015). Nonparametric adaptive estimator of extreme conditional tail probabilities quantiles. Extremes, 18, 437-478.
#'
#' Hill, B.M. (1975). A simple general approach to inference about the tail of a distribution. Annals of Statistics, 3, 1163-1174.
#'
#' @author Ion Grama
#' @export
#'
#' @examples
#' X <- abs(rcauchy(100))
#' weights <- rep(1, length(X))
#' wh <- hill(X, w = weights)
hill <- function(X, weights = rep(1, length(X)), grid = X){
  # by Ion Grama (2006)
  # compute weighted Hill estimator
  Xsort <- rev(sort(X))
  n <- length(X)
  k <- 2:n
  Xord <- rev(order(X))
  word <- weights[Xord]
  #cat('word[k] = ', word[k], '\n')
  hillpoint <- function(tau){
    sum(log(Xsort[which(Xsort>tau)])*word[which(Xsort>tau)])/sum(word[which(Xsort>tau)])-log(tau)
  }
  #if(identical(grid, x)){
    #Y <- cumsum(word[k-1] * log(Xsort[k-1]))/cumsum(word[k-1]) - log(Xsort[k])
  #}else{
    Y <- as.numeric(lapply(rev(sort(grid)), hillpoint))
    #}
  #Y <- hillpoint(grid)
  result <- list(xsort = Xsort, wsort = word, grid = rev(sort(grid)), hill = c(Y))
  class(result) <- "hill"
  result
}



#########################################################################
# the following function is the main tool for the adaptive choice (by Ion Grama)
#' Compute the extreme quantile procedure
#'
#' @param X a numeric vector of data values.
#' @param weights a numeric vector of weigths associated to the vector \eqn{X}.
#' @param initprop the initial proportion at which we begin to test the model.
#' @param gridlen the length of the grid for which the test is done.
#' @param r1 a proportion value of the data from the right that we skip in the test statistic.
#' @param r2 a proportion value of the data from the left that we skip in the test statistic.
#' @param CritVal the critical value assiociated to the weights.
#' @param plot If \code{TRUE}, the results are plotted.
#'
#' @return
#' \item{Xsort}{the sorted vector of the data.}
#' \item{sortweights}{the weights associated to Xsort.}
#' \item{wh}{the weighted Hill estimator associated to X (output of the function hill).}
#' \item{TestingGrid}{the grid used for the statistic test.}
#' \item{TS,TS1,TS.max,TS1.max}{respectively the test statistic, the likelihood ratio test, the maximum of the test statistic and the maximum likelihood ratio test.}
#' \item{Paretodata}{logical: if TRUE the distribution of the data is a Pareto distribution.}
#' \item{Paretotail}{logical: if TRUE a Pareto tail was detected.}
#' \item{madapt}{the first indice of the TestingGrid for which the test statistic exceeds the critical value.}
#' \item{kadapt}{the adaptive indice of the threshold.}
#' \item{kadapt.maxlik}{the maximum likelihood corresponding to the adaptive threshold in the selected testing grid.}
#' \item{hadapt}{the adaptive weighted parameter of the Pareto distribution after the threshold.}
#' \item{Xadapt}{the adaptive threshold.}
#'
#' @details Given a vector of data and assiociated weights, this function compute the adaptive procedure described in Grama and Spokoiny (2008) and Durrieu et al. (2015).
#'
#' We suppose that the data are in the domain of attraction of the Frechet-Pareto type. Otherwise, the procedure will not work.
#'
#' @references
#' Grama, I. and Spokoiny, V. (2008). Statistics of extremes by oracle estimation. Ann. of Statist., 36, 1619-1648.
#'
#' Durrieu, G. and Grama, I. and Pham, Q. and Tricot, J.- M (2015). Nonparametric adaptive estimator of extreme conditional tail probabilities quantiles. Extremes, 18, 437-478.
#'
#' @author Ion Grama
#'
#' @export
#' @examples
#'
#' x <- abs(rcauchy(100))
#' HH <- hill.adapt(x, weights=rep(1, length(x)), initprop = 0.1,
#'                gridlen = 100 , r1 = 0.25, r2 = 0.05, CritVal=10,plot=TRUE)
#' #the critical value 10 is assiociated to the rectangular kernel.
#' HH$Xadapt # is the adaptive threshold
#' HH$hadapt # is the adaptive parameter of the Pareto distribution
#'
hill.adapt <- function(X,  weights = rep(1, length(X)) ,  initprop = 1/10,  gridlen = 100,  r1 = 1/4,  r2 = 1/20,  CritVal = 10 ,  plot = F)
{
  # CPTMLE version 1.001 by Ion Grama (2006)
  # Input:
  # X the observations
  # weights
  # initprop
  # gridlen
  # r1 is the proportion to skip from the right in the test statistic
  # r2 is the proportion to skip from the left in the test statistic
  # CritVal is the critical value
  # plot  alows to visualise the procedure
  # Output:
  # the output is the list:
  #  list(
  #	initprop = initprop,
  #	r1 = r1,
  #	r2 = r2,
  #	CritVal = CritVal,
  #	Xsort = Xsort,  # sorted data X
  #	sortweights = TTSS$sortweights,  # sorted weights
  #	n = n,
  #	wh = wh,  # hill estimator
  #	gridlen = gridlen,
  #	window1 = window1,
  #	window2 = window2,
  #		TS = TS,
  #		TS1 = TS1,
  #    	TS.max = TS.max,
  #    	TS1.max = TS1.max,
  #	Paretodata = Paretodata,  # logical: if  = F the data has no Pareto tail
  #	Paretotail = Paretotail,  # logical: if  = T a Pareto tail was detected
  #	madapt = mad,          # the adaptive m
  #	kadapt = kadapt.maxlik.Final,   # the adaptive threshold index k
  #	kadapt.maxlik = kadapt.maxlik,  # maximum point for TS1
  #	hadapt =  ifelse( is.na(kadapt.maxlik.Final) , NA,  wh[kadapt.maxlik.Final] ),   # the adaptive value of the weighted Hill estimator
  #	Xadapt =  ifelse( is.na(kadapt.maxlik.Final) , NA,  Xsort[kadapt.maxlik.Final])  # the adaptive threshold \tau
  #	)


  #########################################################################
  # functions used in the main body
  #########################################################################
  hill = function(x, weights = rep(1, length(x))){
    # by Ion Grama (2006)
    # compute weighted Hill estimator
    n <- length(x)
    k <- 2:n
    Xsort <- rev(sort(x))
    Xord <- rev(order(x))
    word <- weights[Xord]
    #cat('word[k] = ', word[k], '\n')
    #Y <- hillpoint(grid)
    Y <- cumsum(word[k-1] * log(Xsort[k-1]))/cumsum(word[k-1]) - log(Xsort[k])
    c(NA, Y)
  }
  #########################################################################
  # repeat columns and rows
  #repcol = function(x, nrep){x%*%t(rep(1, nrep))}
  #reprow = function(x, nrep){t(x%*%t(rep(1, nrep)))}

  #########################################################################

  WeightedTestStatist = function(X, weights, m){
    # by Ion Grama (2006)
    # Input:
    # X the observations
    # weights
    # m is the length of the tested interval (given as the number of the order statistic)
    # Output:
    #this function compute the weighted test statistics WTS1,  WTS2 and the estimators wh1 and wh2
    # the output is a list containing WTS1,  WTS2,  wh1,  wh2

    G.KL = function(x){x - log(1 + x)}		# define a function to compute Kullback-Leibler information
    #######################################

    # order log(X) and weights
    n <- length(X)
    order <- rev(order(X))
    Xsort <- rev(sort(X))
    lnX <- log(Xsort)
    weights <- weights[order]

    Sweights <- cumsum(weights[1:n])	  # compute the cumulative weights
    SX <- pmax(cumsum( weights[1:(n-1)] * lnX[1:(n-1)]) - lnX[2:n] * Sweights[1:(n-1)], 0)   # compute SX at points X_2, X_3, ..., X_n (attention! X_1 is excluded )
    wh1 <- SX[1:(n-1)]/Sweights[1:(n-1)]        # compute wighted Hill's estimator at X_2, X_3, ..., X_n (attention! X_1 is excluded )
    wh2 <- abs(SX[m-1]-SX[1:(m-1)])/(Sweights[m-1]-Sweights[1:(m-1)]); wh2[m-1] = NA  # compute weighted wh2  at X_2, X_3, ..., X_m (attention! X_1 is excluded )
    # do not delete the following comment
    # sometimes the difference SX[m-1]-SX[1:(m-1)] can be negative
    # because of the computing error (of order -10^-17)
    # that is why in the above formula we used abs(...)

    ## the following lines are introduced to check the correctness
    ##if (  is.number(SX[m-1]) == F  ){ #  this line is to check correctness
    ##cat("m-1 = ", m-1, "\n") #  this line is to check correctness
    ##cat("SX = ", SX, "\n") #  this line is to check correctness
    ##cat("abs(SX[m-1]-SX[1:(m-1)]) = ", abs(SX[m-1]-SX[1:(m-1)]), "\n") #  this line is to check correctness
    ##cat("Sweights[m-1]-Sweights[1:(m-1)] = ", Sweights[m-1]-Sweights[1:(m-1)], "\n") #  this line is to check correctness
    ##cat("wh2[(m-10):(m-1)] = ", round(wh2[(m-10):(m-1)], 5), "\n") #  this line is to check correctness
    ##cat("Xsort[(m-10):(m-1)] = ", Xsort[(m-10):(m-1)], "\n") #  this line is to check correctness
    ##}  #  this line is to check correctness

    #############################################
    # The previous computation computes wh2 only for a given m and requires a vector of size m
    # In the following lines here we propose a new vector computation of wh2:
    # this computation requires introduction of a matrix of size m x n
    # and compute all the values wh2 at once (this is not used in this version)
    #	Dh = repcol(Sh[1:m], nrep = m)-reprow(Sh[1:m], nrep = n)    # the same computation:  E = expand.grid(Sh[1:m], Sh[1:m]); matrix(E[, 1]-E[, 2], nrow = m)
    #	Dw = repcol(Sw[1:m], nrep = m)-reprow(Sw[1:m], nrep = n)
    #	wh2matr = Dh/Dw
    #	wh2matr[, m]  #  this line is to check correctness
    #	wh2           #  this line is to check correctness
    ############################################

    TS1 <- Sweights[1:(m-1)] * G.KL(  wh1[1:(m-1)]/wh1[m-1] - 1  )                  # compute weighted TS1  at X_2, X_3, ..., X_m (attention! X_1 is excluded )
    TS2 <- (Sweights[m-1]-Sweights[1:(m-1)]) * G.KL(  wh2[1:(m-1)]/wh1[m-1] - 1  )  # compute weighted TS2  at X_2, X_3, ..., X_m (attention! X_1 is excluded )

    # we can compute normalized TS1 and TS2 if necessary (not used in the sequel)
    #varW = sqrt(mean(weights[1:m]^2))   # the normalisation is not used here
    #TS1 = TS1/varW   # the normalisation is not used here
    #TS2 = TS2/varW # the normalisation is not used here

    ###############
    # the output

    list(WTS1 = TS1,  WTS2 = TS2,  wh1 = wh1[1:(m-1)],  wh2 = wh2[1:(m-1)],  sortweights = weights)

  } # end of WeightedTestStatist

  ############################################################################################################################################################
  #########################################################################


  argmaxw <- function(x,  window){ window[rev(order(x[window]))]}



  ############################################################################################################################################################
  ############################################################################################################################################################
  ############################################################################################################################################################

  #if(plot == T) pdf("hill.adapt.pdf")

  ############################################################################################################################################################
  # program body of the function hill.adapt

  n <- length(X)
  idx <- split(1:n, X)
  test<-data.frame(X=sapply(idx, function(i) X[i[1]]),
                   w=sapply(idx, function(i) sum(weights[i])))
  Xunique<-test$X
  Wunique<-test$w
  n <- length(Xunique)
  wh <- hill(Xunique, w = Wunique)
  Xsort <- sort(Xunique, decreasing = TRUE)
  order <- rev(order(Xunique))
  sortweights <- Wunique[order]
  cumsumw <- cumsum(sortweights)
  initprop <- max(initprop,3/n)  # = ifelse(floor(n*initprop) > 3, initprop, 3/n)
  gridlen <- min(gridlen,n-1)
  #cat("this line is to check correctness in hill.adapt: gridlen = ", gridlen, "\n")

  Paretodata <- TRUE   # Paretodata = TRUE if the whole data are from Pareto law
  Paretotail <- FALSE	# this variable will indicate if the tail is
  FirstTime <- TRUE	# FirstTime is a control valiable to verify if the test statistic exceeds CritVal from the very beginning
  # i.e. FirstTime = TRUE if it is the first passage of the loop in  m in TestingGrid

  mad <- n	# the adaptive value mad (from "m adaptive") is first initialized as equal to n
  # if the the Pareto model is accepted and
  # a change from Pareto tail is not detected the final result will be  m = n

  Grid <- floor(seq(1,  n-1,  length = gridlen))				# define the grid of values for the index k
  #cat("this line is to check correctness in hill.adapt: n = ", n, "\n")
  #cat("this line is to check correctness in hill.adapt: Grid = ", Grid, "\n")

  if(n > 3)   {  # begin if n>3
    #############################################################################################
    #####  logSpacings = (1:(n-1))*log( Xsort[1:(n-1)]/Xsort[2:n] );   # compute the log spacings
    #############################################################################################

    #cat("this line is to check correctness in hill.adapt: floor(n * initprop) = ", floor(n * initprop), "\n")


    while(floor(n * initprop) >= 3) {		# at the end of the while if  m0 = n * initprop is too large we resume with initprop = 0.8 * initprop
      # n * initprop cannot take value 3
      #cat("this line is to check correctness in hill.adapt: n * initprop = ", n * initprop, "\n") #  this line is to print the init propotion
      #cat("this line is to check correctness in hill.adapt: Grid = ", Grid, "\n")

      count_k0_too_smal = 0
      TS.max <- rep(NA,  n)  # to store TS.max  for each m in TesingGrid
      TS1.max <- rep(NA,  n) # to store TS1.max for each m in TesingGrid
      kadapt.maxlik <- rep(NA,  n) # to store kadapt for each m in TesingGrid

      FirstTime <- TRUE
      TestingGrid <- Grid[Grid >=   max(floor(n * initprop), 3)   ] 	# select testing points > =  m0 = n * initprop on the grid (attention! initprop can change if  m0 is too large at the next repeat loop)

      #cat(" TestingGrid = ", TestingGrid, "\n")
      #cat(" max(floor(n * initprop), 3) = ", max(floor(n * initprop), 3), "\n")


      #print(paste(      "floor(n * initprop) = ", floor(n * initprop)            ))
      #print(paste("n = ", n * initprop))

      #for(m in (TestingGrid)) {   # this loop compute iteratively the test statistics TS.max and TS1.max
      for(indm in 1:length(TestingGrid)) {   # this loop compute iteratively the test statistics TS.max and TS1.max
        # the previous line was introduces by Pham ?????? apparently to have acces to the index \hat m -1 the last accepted

        m <- TestingGrid[indm]
        #cat("this line is to check correctness in hill.adapt: m = ", m, "\n") #  this line is to print m
        #cat("this line is to check correctness in hill.adapt: indm = ", indm, "\n") #  this line is to print m

        TTSS <- WeightedTestStatist(Xunique, Wunique, m)  # compute the test statistic
        TS <- TTSS$WTS1 + TTSS$WTS2		# computed at     X_2, X_3, ..., X_m
        TS1 <- TTSS$WTS1					# computed at     X_2, X_3, ..., X_m
        TS <- c(NA, TS)						# computed at X_1, X_2, X_3, ..., X_m
        TS1 <- c(NA, TS1)				    # computed at X_1, X_2, X_3, ..., X_m

        #cat("TS = ", TS, "\n") #  this line is to check correctness

        restrictedw <- which(cumsumw >=  r1*cumsumw[m] & cumsumw <=  (1-r2)*cumsumw[m] )
        w1 <- min(restrictedw)
        w2 <- max(restrictedw)


        window1 <- max(2,    w1  ) 	#define left  end of the testing window in the set X_1, X_2, X_3, ..., X_m
        window2 <- min(m-2,  w2  ) 	#define right end of the testing window in the set X_1, X_2, X_3, ..., X_m

        #cat("window1 = ", window1, "\n") 	#  this line is to check correctness
        #cat("window2 = ", window2, "\n") 	#  this line is to check correctness
        #print(paste("m = ", m))			#  this line is to check correctness
        #print(window1:window2)			#  this line is to check correctness

        kadapt.maxlik[m] <- argmaxw(TS1,  window = window1:window2)[1]
        ### kadapt.maxlik[m] is the adaptive value of k
        ### computed by penalized maxlikelihood (PML) for each m of the TesingGrid

        # the next 2 lines are included to plot the outputs (they are not necessary to compute the adaptive threshold  kadapt)
        TS.max[m] <- max(TS[window1:window2])     # store TS.max for each m
        TS1.max[m] <- max(TS1[window1:window2])	 # store TS1.max for each m


        #############################################################
        #cat("TS.max[m] = ", round(TS.max[m], 5), "\n") #  this line is to check correctness
        #if(is.number(TS.max[m]) == F){cat("TS = ", round(TS, 5), "\n")} #  this line is to check correctness
        #############################################################

        #############
        # kadapt.changepoint[m] = (window1:window2)[rev(order(TS[window1:window2]))[1]]
        #	### this is the alternative adaptive value of  k
        #	### computed as the change point (CP)
        #	### we found by simulations that the CP choice is not better than the PML choice
        #############

        #############################################################
        # visualise the procedure  #######
        if(plot == T){

          par(mfrow = c(2, 2))
          ylim0 <- c(0,  1.2*CritVal )

          # plot TS = TS1 + TS2
          #plot(log(Xsort[1:m]), xlim = xrange, TS.max[1:m], ylim = ylim0,  ylab = "Test Statist");
          plot(TS.max[1:m], ylim = ylim0,  ylab = "Test Statist: TS.max")
          #abline(v = c(window1, window2), col = "red");
          abline(v = m, col = "blue")
          #abline(h = sqrt(mean(weights[1:m]^2)), col = 8);
          abline(h = CritVal, col = 8)
          title(paste( "start.pt", round(n * initprop, 0),   ",  m = ",  m))

          # plot the hill estimator
          #plot(log(Xsort[1:m]), wh[1:m], xlim = xrange, ylim = c(0, max(wh[2:m])),  ylab = "hill estimator");
          plot(wh[1:m], ylim = c(0, max(wh[which(!is.na(wh[1:m]))])),  ylab = "hill estimator: hill")
          #abline(v = c(window1, window2));
          abline(v = kadapt.maxlik[m], col = "red")
          abline(h = wh[kadapt.maxlik[m]], col = "red")
          title( paste("k.ad = ",  kadapt.maxlik[m], ",  hill.ad = ",  round(wh[ kadapt.maxlik[m] ], 4))   )

          # plot TS1
          #plot(log(Xsort[1:m]), TS1[1:m], xlim = xrange, ylim = ylim0,  ylab = "TS1");
          plot(TS1[1:m], ylim = ylim0,  ylab = "Pen Lik: TS1")
          abline(v = c(window1, window2), col = "red"); abline(v = m, col = "blue")
          #abline(v = log(Xsort[kadapt.maxlik[m]]), col = 8);
          abline(v = kadapt.maxlik[m], col = "red", lwd = 2)
          title(paste("k.ad = ",  kadapt.maxlik[m]))

          # plot TS2
          #plot(TTSS$WTS2[1:m], ylim = ylim0, ylab = "TS2"); abline(v = c(window1, window2));

          par(mfrow = c(1, 1))
        }
        #############################################################
        # end: visualise the procedure  #######
        #############################################################

        #############################################################
        # this line is to check correctness
        #plot(TS1); points(window1:window2, TS1[window1:window2], col = 5, ylim = c(0, 10));
        #############################################################
        #print(paste(" ==  ==  ==  ==  ==  ==  ==  =  TSMAX = ", TS.max[m])) # this line is to check correctness
        #############################################################


        if( (TS.max[m] >=  CritVal) & (!is.na(TS.max[m])) ) {
          if(FirstTime & (count_k0_too_smal<10) ) {

            #cat("k0 is too small: resume with the next element on the grid\n");
            count_k0_too_smal <- count_k0_too_smal+1
            next;
          } else if(FirstTime == FALSE){
            mad <- TestingGrid[ifelse(indm == 1,  1,  indm-1)]
            Paretodata <- FALSE
            Paretotail <- TRUE
            break
          } else {
            Paretodata <- FALSE
            Paretotail <- FALSE
            break
          }
        }

        FirstTime <- FALSE
      }#end for(m in TestingGrid) # this is the end of loop on the testing grid


      ### the final part of the program

      ### in the following line we verify if the test statistic exceeds CritVal
      ### from the very beginning
      ### and if TS.max[m0] is > CritVal at m0 we resume with a diminished initprop

      if(FirstTime  ==  FALSE | initprop < 5/n) break

      initprop <- initprop * 0.8
      #cat("initprop is too large: resume with initprop = 0.8*initprop\n")
    } #end while(n*initprop > 3)

    if(Paretodata) {
      if(FirstTime) {
        #cat("No Pareto tail detected\n")
        Paretotail <- FALSE
        mad <- NA
        kadapt.maxlik.Final <- NA
        hdapt.maxlik <- NA
        Xdapt.maxlik <- NA
      }
      else {
        mad <- n
        kadapt.maxlik.Final <- n
      }
    }
    else {
      kadapt.maxlik.Final <- kadapt.maxlik[mad]		# the adaptive  k
    }

    ### end of the final part of the program

    ### the following line was included to choose between MLE and CP approach: do not use in this version
    ### if(MLE == T){kad = kml; had = hml; Xad = Xml} else {kad = kls; had = hls; Xad = Xls}
  }  # end if n>3

  else{  # begin if n<=3
    TestingGrid <- NA
    TS <- NA
    TS1 <- NA
    TS.max <- NA
    TS1.max <- NA
    Paretodata <- NA  # logical: if  = F the data has no Pareto tail
    Paretotail <- NA  # logical: if  = T a Pareto tail was detected
    mad <- NA # madapt = mad,          # the adaptive m
    kadapt.maxlik.Final <- n  #kadapt = kadapt.maxlik.Final,   # the adaptive k
    kadapt.maxlik <- NA  #kadapt.maxlik = kadapt.maxlik,  # maximum points for TS1
  }  # end else if n>3

  ############################
  ### compute the output values
  res <- list(
    n = n,
    initprop = initprop,
    r1 = r1,
    r2 = r2,
    CritVal = CritVal,
    Xsort = Xsort,  # sorted data X
    sortweights = sortweights,  # sorted weights
    wh = wh,  # hill estimator
    TestingGrid = TestingGrid,
    TS = TS,
    TS1 = TS1,
    TS.max = TS.max,
    TS1.max = TS1.max,
    Paretodata = Paretodata,  # logical: if  = F the data has no Pareto tail
    Paretotail = Paretotail,  # logical: if  = T a Pareto tail was detected
    madapt = mad,          # the adaptive m
    kadapt = kadapt.maxlik.Final,   # the adaptive k
    kadapt.maxlik = kadapt.maxlik,  # maximum points for TS1
    hadapt =  ifelse( is.na(kadapt.maxlik.Final) , NA,  wh[kadapt.maxlik.Final] ),   # the adaptive weighted h
    Xadapt =  ifelse( is.na(kadapt.maxlik.Final) , NA,  Xsort[kadapt.maxlik.Final])  # the adaptive X
  )

  class(res) <- "hill.adapt"
  #attr(res,  "call")  <-  sys.call()
  res
  #if(plot == T) dev.off()
} # end of the function  hill.adapt

##############################################################################################
# the following function visualise the selection procedure


#' Hill.adapt Plot
#'
#' @description Graphical representation of the hill.adapt function last iteration
#'
#' @param x output object of the function hill.adapt.
#' @param ... further arguments passed to or from other methods.
#'
#' @details The weighted hill estimator, the test statistic, the penalized likelihood graphs of the last iteration and the survival function are given. The blue line corresponds to the threshold (indice or value). The magenta lines correspond to the window (r1, r2) where the estimation is computed. The red lines corresponds to the initial proportion (initprop) and the last non rejected point of the statistic test (madapt).
#'
#' @seealso
#' \code{\link{hill.adapt}}, \code{\link{plot}}
#'
#' @export
#' @examples
#'
#' x <- abs(rcauchy(100))
#' HH <- hill.adapt(x, weights=rep(1, length(x)), initprop = 0.1,
#'                gridlen = 50 , r1 = 0.25, r2 = 0.05, CritVal=10)
#' plot(HH)
#'
#'
plot.hill.adapt <- function(x, ...){
  k0 <- x$n*x$initprop
  m <- x$madapt
  madapt <- x$madapt
  kadapt <- x$kadapt
  #Xsort = hill.adapt$Xsort
  m1 <- x$madapt*x$r1
  m2 <- x$madapt*(1-x$r2)

  par(mfrow = c(2, 2))
  plot(x$wh[1:m],  ..., xlab = "k", ylab = "Hill estimator")
  abline(h =  x$hadapt, col = "blue", lwd = 2)
  abline(v =  kadapt, col = "blue", lwd = 2, lty = 2)
  abline(v =  madapt, col = "red", lty = 2, lwd = 2)
  abline(v =   k0, col = "red", lty = 2, lwd = 2)
  title( paste("Hill = ",  round(x$hadapt, 4),  "kadapt = ",  round(kadapt, 0))  );
  plot(x$TS.max[1:m],  ..., xlab = "k", ylim = c(0, x$CritVal), ylab = "Test Statistic")
  abline(h = x$CritVal, lty = 2, lwd = 2)
  abline(v =   madapt, col = "red", lty = 2, lwd = 2)
  abline(v =   k0, col = "red", lty = 2, lwd = 2)
  title(paste("Test Statistic,  start pt = ",  k0, " break pt = ",  madapt));
  plot(x$TS1[1:m],  ..., xlab = "k", ylab = "Penalized MaxLikelihood")
  abline(v =  kadapt, col = "blue", lwd = 2, lty = 2)
  abline(v =  madapt, col = "red", lty = 2, lwd = 2)
  #abline(v =  k0, col = "red", lty = 2, lwd = 2)
  abline(v =  c(m1, m2), col = "magenta", lty = 2)
  title( paste("Pen Lik,  kadapt = ",  round(kadapt, 0))  );
  plot(rev(log(x$Xsort)), 1-wecdf(x$Xsort, rev(x$Xsort), x$sortweights), type = "s", xlab = "log(X)", ylab = "survival function",  ...)
  abline(v =  log(x$Xadapt), col = "blue", lty = 2, lwd = 2)
  title(paste("Survival function,  threshold = ", round(x$Xadapt, 4)));
  par(mfrow = c(1, 1))
}



#' Weighted quantile
#'
#' @description Compute the weighted quantile of order p.
#'
#' @param X a vector of data.
#' @param p a vector of probabilities.
#' @param weights the weights assiociated to the vector \eqn{X}.
#'
#' @details Give the weighted quantile for a given \eqn{p}
#'
#' @return A vector of quantile assiociated to the probabilities vector given in input.
#'
#' @export
#' @examples
#'
#' X <- rpareto(10)
#' p <- seq(0.01, 0.99, 0.01)
#' plot(p, wquantile(X, p, rep(1,length(X))), type = "s")
#'
wquantile <- function(X, p, weights = rep(1, length(X))){
  #Input : X is the data
  # weights is the vector of weights assiociated to X
  # p is the vector of grid for which we want the weighted quantile
  # output: vector of weighted ecdf assiociated to x
  ord <- (order(X))
  wdf <- cumsum(weights[ord])/sum(weights)
  rank.p <- function(prob){
    k <- which(wdf >=  prob)[1]
    x <- c(sort(X))[k]
  }
  px <- sapply(p, rank.p)
  return(px)
}

#' Weighted Empirical Cumulative Distribution Function
#'
#' @description Calculate the values of the weighted empirical cumulative distribution function for a given vector of data
#'
#' @param X the vector of data to create the wecdf.
#' @param x the vector of data that you want the corresponding wecdf values.
#' @param weights the weights applicated to the vector \eqn{X}.
#'
#' @details Give the value of the wecdf.
#' If the weights are 1 (the default value), the wecdf become the ecdf of \eqn{X}.
#'
#' @return Return a vector of the wecdf values corresponding to \eqn{x} given a reference vector \eqn{X} with weights \eqn{weights}.
#'
#' @export
#' @examples
#'
#' X <- rpareto(10)
#' x <- seq(0.8, 50, 0.01)
#' plot(x, wecdf(X, x, rep(1,length(X))))
#'
#' #to compare with the ecdf function
#' f <- ecdf(X)
#' lines(x, f(x), col = "red", type = "s")
#'
wecdf <- function(X, x, weights = rep(1, length(X))){
  #by Kevin Jaunatre
  #Input : X is the data
  # weights is the vector of weights assiociated to X
  # x is the vector of grid for which we want the weighted ecdf
  # output: vector of weighted ecdf assiociated to x
  ord <- (order(X))
  wdf <- cumsum(weights[ord])/sum(weights)
  rank.x <- function(x){
    k <- ceiling(rank(c(x, X)))[1]
    p <- c(0, wdf)[k]
  }
  px <- sapply(x, rank.x)
  return(px)
}


#' Compute the extreme quantile procedure on a time dependent data
#'
#' @description Compute the function hill.adapt on time dependent data.
#'
#' @param X a vector of the observed values.
#' @param t a vector of time covariates which should have the same length as X.
#' @param Tgrid a grid of time (can be any sequence in the interval \code{[min(t) , max(t)]} ).
#' @param h a bandwidth value (vector values are not admitted).
#' @param kernel a kernel function used to compute the weights in the time domain, with default the truncated Gaussian kernel.
#' @param kpar a value for the kernel function parameter, with no default value.
#' @param CritVal a critical value associated to the kernel function given by \code{\link{CriticalValue}}. The default value is 3.6 corresponding to the truncated Gaussian kernel.
#' @param gridlen the gridlen parameter used in the function hill.adapt. The length of the grid for which the test will be done.
#' @param initprop the initprop parameter used in the function hill.adapt. The initial proportion at which we will begin to test the model.
#' @param r1 the r1 parameter used in the function hill.adapt. The proportion from the right that we will skip in the test statistic.
#' @param r2 the r2 parameter used in the function hill.adapt. The proportion from the left that we will skip in the test statistic.
#' @param x the result of the hill.ts function
#' @param ... further arguments to be passed from or to other methods.
#'
#' @details For a given time serie and kernel function, the function hill.ts will give the results of the adaptive procedure for each \eqn{t}. The adaptive procedure is described in Durrieu et al. (2005).
#'
#' The kernel implemented in this packages are : Biweight kernel, Epanechnikov kernel, Rectangular kernel, Triangular kernel and the truncated Gaussian kernel.
#'
#' @return
#' \item{Tgrid}{the given vector \eqn{Tgrid}.}
#' \item{h}{the given value \eqn{h}.}
#' \item{Threshold}{the adaptive threshold \eqn{\tau} for each \eqn{t} in \eqn{Tgrid}.}
#' \item{Theta}{the adaptive estimator of \eqn{\theta} for each \eqn{t} in \eqn{Tgrid}.}
#'
#'
#' @references Durrieu, G. and Grama, I. and Pham, Q. and Tricot, J.- M (2015). Nonparametric adaptive estimator of extreme conditional tail probabilities quantiles. Extremes, 18, 437-478.
#'
#' @export
#' @seealso \code{\link{hill.adapt}}, \code{\link{Biweight.kernel}}, \code{\link{Epa.kernel}}, \code{\link{Rectangular.kernel}}, \code{\link{Triang.kernel}}, \code{\link{TruncGauss.kernel}}
#'
#' @examples
#'
#' theta <- function(t){
#'    0.5+0.25*sin(2*pi*t)
#'  }
#' n <- 5000
#' t <- 1:n/n
#' Theta <- theta(t)
#' Data <- NULL
#' Tgrid <- seq(0.01, 0.99, 0.01)
#' #example with fixed bandwidth
#' \dontrun{ #For computing time purpose
#'   for(i in 1:n){
#'      Data[i] <- rparetomix(1, a = 1/Theta[i], b = 5/Theta[i]+5, c = 0.75, precision = 10^(-5))
#'    }
#'
#'   #example
#'   hgrid <- bandwidth.grid(0.009, 0.2, 20, type = "geometric")
#'   TgridCV <- seq(0.01, 0.99, 0.1)
#'   hcv <- bandwidth.CV(Data, t, TgridCV, hgrid, pcv = 0.99, TruncGauss.kernel,
#'                      kpar = c(sigma = 1), CritVal = 3.6, plot = TRUE)
#'
#'   Tgrid <- seq(0.01, 0.99, 0.01)
#'   hillTs <- hill.ts(Data, t, Tgrid, h = hcv$h.cv, kernel = TruncGauss.kernel,
#'              kpar = c(sigma = 1), CritVal = 3.6,gridlen = 100, initprop = 1/10, r1 = 1/4, r2 = 1/20)
#'   plot(hillTs$Tgrid, hillTs$Theta, xlab = "t", ylab = "Estimator of theta")
#'   lines(t, Theta, col = "red")
#'
#' }
#'
#'
hill.ts <- function(X, t, Tgrid, h, kernel = TruncGauss.kernel, kpar = NULL, CritVal = 3.6, gridlen = 100, initprop = 1/10,  r1 = 1/4,  r2 = 1/20){
  # By Kevin Jaunatre 2015
  # Input : X.ts is a time serie
  # kernel is a kernel function for the weights
  # h is the bandwidth for the kernel function
  # t is a vector of time for which we want to calculate the adaptive model
  # Tmax is the
  # Output : list of t,  Theta and threshold assiocated to each t
  #
  if( identical(kernel,TruncGauss.kernel) & is.null(kpar) ){
    kpar = c(sigma = 1)
  }

  timeSeriehill.adapt <- function(Tgrid){
    indices <- which((t>= Tgrid-h)&(t<= Tgrid+h))
    seq.X <- X[indices]
    tgrid <- t[indices]
    tx <- (tgrid-Tgrid) / h
    if( is.null(kpar) ){
      par <- list(tx)
    }else{
      par <- list(tx, kpar)
    }
    weights <- do.call(kernel, par)
    hillAdapt <- hill.adapt(seq.X,  weights, CritVal = CritVal, gridlen = gridlen, initprop = initprop,  r1 = r1,  r2 = r2)
    tau.adapt <- hillAdapt$Xadapt
    theta.adapt <- hillAdapt$hadapt
    k.adapt <- hillAdapt$kadapt
    return(list(tau = tau.adapt, theta = theta.adapt, kadapt = k.adapt))
  }
  result <- sapply(Tgrid,  timeSeriehill.adapt)
  res <- list(X = X, t = t, Tgrid = Tgrid, h = h, kernel = kernel,
              kpar = kpar, CritVal = CritVal, r1 = r1, r2 = r2, threshold = as.numeric(result[1, ]),
              Theta = as.numeric(result[2, ]), kadapt = as.numeric(result[3, ]))
  class(res) <- "hill.ts"
  #attr(res,  "call")  <-  sys.call()
  res
}


#' Prints a hill.ts object
#'
#'
#' @rdname hill.ts
#' @export print hill.ts
print.hill.ts<-function(x,...){
  result <- list(Tgrid = x$Tgrid, threshold = x$threshold, Theta = x$Theta)
  result
}




#' Gaussian kernel function
#'
#' @description Gaussian kernel function
#'
#' @param x a vector.
#'
#' @details Gaussian Kernel with the value of standard deviation equal to 1/3.
#' \deqn{
#'   K(x) = (1/{(1/3)*sqrt(2 \pi)}  exp(-(3*x)^2/2)) (abs(x) <= 1)
#' }
#' We recommend a critical value of 8.3 for this kernel.
#'
#'
#' @export
#'
#' @examples
#' plot(function(x) Gaussian.kernel(x), -2, 2,
#' main = " Gaussian kernel")
#'
Gaussian.kernel <- function(x){
  return((1/((1/3)*sqrt(2*pi))*exp(-(3*x)^2/2))*(abs(x)<= 1))
}




#' Truncated Gaussian kernel function
#'
#' @description Truncated Gaussian kernel function
#'
#' @param x a vector.
#' @param sigma the standard deviation of the truncated gaussian kernel.
#'
#' @details Truncated Gaussian Kernel with \eqn{sigma} the standard deviation parameter with default value \eqn{1}.
#' \deqn{
#'   K(x) = (1/{sigma*sqrt(2 \pi)}  exp(-(x/sigma)^2/2)) (abs(x) <= 1)
#' }
#' We recommend a critical value of 3.6 for this kernel with sigma=1.
#'
#'
#' @export
#'
#' @examples
#' plot(function(x) TruncGauss.kernel(x), -2, 2,
#' main = " Truncated Gaussian kernel")
#'
TruncGauss.kernel <- function(x,sigma=1){
  return((1/(sigma*sqrt(2*pi))*exp(-(x/sigma)^2/2))*(abs(x)<= 1))
}


#' Epanechnikov kernel function
#'
#' @description Epanechnikov kernel function.
#'
#' @param x vector.
#'
#' @details Epanechnikov kernel:
#' \deqn{
#'   K(x) = 3/4 ( 1 - x^2 ) (abs(x)<=1)
#' }
#' We recommend a critical value of 6.1 for this kernel function.
#' @export
#'
#' @examples
#' plot(function(x) Epa.kernel(x), -2, 2, ylab = "Epanechnikov kernel function")
#'
Epa.kernel <- function(x){
  return((3/4*(1-x^2))*(abs(x)<= 1))
}


#' Triangular kernel function
#'
#' @description Triangular kernel function
#'
#' @param x a vector.
#'
#' @details Triangular Kernel
#' \deqn{
#'   K(x) = ( 1 - abs(x) )  (abs(x) <= 1)
#' }
#' We recommend a critical value of 6.9 for this kernel.
#'
#' @export
#'
#' @examples
#' plot(function(x) Triang.kernel(x), -2, 2,
#' main = " Triangular kernel")
#'
Triang.kernel <- function(x){
  return((1-abs(x))*(abs(x)<= 1))
}

#' Biweight kernel function
#'
#' @description Biweight kernel function.
#'
#' @param x a vector.
#'
#' @details Biweight kernel:
#' \deqn{
#'   K(x) = 15/16 ( 1 - x^2 )^2  (abs(x)<=1)
#' }
#' We recommend a critical value of 7 for this kernel function.
#'
#' @export
#'
#' @examples
#' plot(function(x) Biweight.kernel(x),-2, 2,
#' main = " Biweight kernel ")
#'
Biweight.kernel <- function(x){
  return(15/16*(1-x^2)^2*(abs(x)<= 1))
}

#' Rectangular kernel function
#'
#' @details Rectangular kernel function
#'
#' @param x a vector.
#'
#' @details Rectangular Kernel
#' \deqn{
#'   K(x) = 1 (abs(x) <= 1)
#' }
#' We recommend a critical value of 10 for this kernel.
#'
#' @export
#'
#' @examples
#' plot(function(x) Rectangular.kernel(x), -2, 2,
#' main = " Rectangular kernel ")
#'
Rectangular.kernel <- function(x){
  return(1*(abs(x)<= 1))
}

#' Goodness of fit test statistics
#'
#' @description goftest is a generic function whose application depends on the class of its argument.
#'
#' @param object model object.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' The form of the value returned by goftest depends on the class of its argument. See the documentation of the particular methods for details of what is produced by that method.
#'
#' @references Grama, I. and Spokoiny, V. (2008). Statistics of extremes by oracle estimation. Ann. of Statist., 36, 1619-1648.
#'
#' Durrieu, G. and Grama, I. and Pham, Q. and Tricot, J.- M (2015). Nonparametric adaptive estimator of extreme conditional tail probabilities quantiles. Extremes, 18, 437-478.
#'
#' @export
#'
#' @seealso \code{\link{goftest.hill.adapt}}, \code{\link{goftest.hill.ts}}
#'
goftest <- function(object, ...) UseMethod("goftest",  object)


#' Goodness of fit test statistics
#'
#' @description Give the results of the goodness of fit tests for testing the null hypothesis that the tail is fitted by a Pareto distribution, starting from the adaptive threshold, against the Pareto change point distribution for all possible change points (for more details see pages 447 and 448 of Durrieu et al. (2015)).
#'
#' @param object output of the function hill.adapt.
#' @param plot If \code{TRUE}, the test statistics are plotted.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' \item{TS.window}{the test statistic inside the window. (pages 447 and 448 of Durrieu et al.(2015))}
#' \item{TS}{the test statistic.}
#' \item{CritVal}{the critical value of the test.}
#'
#' @export
#'
#' @references Grama, I. and Spokoiny, V. (2008). Statistics of extremes by oracle estimation. Ann. of Statist., 36, 1619-1648.
#'
#' Durrieu, G. and Grama, I. and Pham, Q. and Tricot, J.- M (2015). Nonparametric adaptive estimator of extreme conditional tail probabilities quantiles. Extremes, 18, 437-478.
#'
#' @seealso \code{\link{hill.adapt}}, \code{\link{goftest}}
#'
#' @examples
#' x <- abs(rcauchy(100))
#' HH <- hill.adapt(x, weights=rep(1, length(x)), initprop = 0.1,
#'                gridlen = 100 , r1 = 0.25, r2 = 0.05, CritVal=10)
#'
#' #the critical value 10 is assiociated to the rectangular kernel.
#'
#' goftest(HH, plot = TRUE)
#'
#' # we observe that for this data, the null hypothesis that the tail
#' # is fitted by a Pareto distribution is not rejected as the maximal
#' # value in the graph does not exceed the critical value.
#'
#'
goftest.hill.adapt <- function(object, plot = FALSE, ...){
  #Function returning the test statistic corresponding to the threshold for one T in Tgrid (output of hill.adapt)
  #hill is the output of the function hill.adapt
  #Output :
  #TS.window : the test statistic within the window (r1 and r2 parameters)
  #TS : the test statistic without window
  WeightedTestStatist <- function(X, weights, m){
    # by Ion Grama (2006)
    # Input:
    # X the observations
    # weights
    # m is the change point (given as the number of the order statistic)
    # Output:
    #this function compute the weighted test statistics WTS1,  WTS2 and the estimators wh1 and wh2
    # the output is a list containing WTS1,  WTS2,  wh1,  wh2

    G.KL <- function(x){x - log(1 + x)}		# define a function to compute Kullback-Leibler information
    #######################################

    # order log(X) and weights
    n <- length(X)
    order <- rev(order(X))
    Xsort <- rev(sort(X))
    lnX <- log(Xsort)
    weights <- weights[order]

    Sweights <- cumsum(weights[1:n])	  # compute the cumulative weights
    SX <- cumsum( weights[1:(n-1)] * lnX[1:(n-1)]) - lnX[2:n] * Sweights[1:(n-1)]   # compute SX at points X_2, X_3, ..., X_n (attention! X_1 is excluded )
    wh1 <- SX[1:(n-1)]/Sweights[1:(n-1)]        # compute wighted Hill's estimator at X_2, X_3, ..., X_n (attention! X_1 is excluded )
    wh2 <- abs(SX[m-1]-SX[1:(m-1)])/(Sweights[m-1]-Sweights[1:(m-1)]); wh2[m-1] = NA;  # compute weighted wh2  at X_2, X_3, ..., X_m (attention! X_1 is excluded )
    # do not delete the following comments
    # sometimes the difference SX[m-1]-SX[1:(m-1)] can be negative
    # because of the computing error (of order -10^-17)
    # that is why in the above formula we used abs(...)

    ## the following lines are introduced to check the correctness
    ##if (  is.number(SX[m-1]) == F  ){ #  this line is to check correctness
    ##cat("m-1 = ", m-1, "\n") #  this line is to check correctness
    ##cat("SX = ", SX, "\n") #  this line is to check correctness
    ##cat("abs(SX[m-1]-SX[1:(m-1)]) = ", abs(SX[m-1]-SX[1:(m-1)]), "\n") #  this line is to check correctness
    ##cat("Sweights[m-1]-Sweights[1:(m-1)] = ", Sweights[m-1]-Sweights[1:(m-1)], "\n") #  this line is to check correctness
    ##cat("wh2[(m-10):(m-1)] = ", round(wh2[(m-10):(m-1)], 5), "\n") #  this line is to check correctness
    ##cat("Xsort[(m-10):(m-1)] = ", Xsort[(m-10):(m-1)], "\n") #  this line is to check correctness
    ##}  #  this line is to check correctness

    #############################################
    # The previous computation computes wh2 only for a given m and requires a vector of size m
    # In the following lines here we propose a new vector computation of wh2:
    # this computation requires introduction of a matrix of size m x n
    # and compute all the values wh2 at once (this is not used in this version)
    #	Dh = repcol(Sh[1:m], nrep = m)-reprow(Sh[1:m], nrep = n)    # the same computation:  E = expand.grid(Sh[1:m], Sh[1:m]); matrix(E[, 1]-E[, 2], nrow = m)
    #	Dw = repcol(Sw[1:m], nrep = m)-reprow(Sw[1:m], nrep = n)
    #	wh2matr = Dh/Dw
    #	wh2matr[, m]  # to check
    #	wh2          # to check
    ############################################

    TS1 <- Sweights[1:(m-1)] * G.KL(  wh1[1:(m-1)]/wh1[m-1] - 1  )                  # compute weighted TS1  at X_2, X_3, ..., X_m (attention! X_1 is excluded )
    TS2 <- (Sweights[m-1]-Sweights[1:(m-1)]) * G.KL(  wh2[1:(m-1)]/wh1[m-1] - 1  )  # compute weighted TS2  at X_2, X_3, ..., X_m (attention! X_1 is excluded )

    # we can compute normalized TS1 and TS2 if necessary
    #varW = sqrt(mean(weights[1:m]^2))   # the normalisation is not used here
    #TS1 = TS1/varW   # the normalisation is not used here
    #TS2 = TS2/varW # the normalisation is not used here

    ###############
    # the output

    list(WTS1 = TS1,  WTS2 = TS2,  wh1 = wh1[1:(m-1)],  wh2 = wh2[1:(m-1)],  sortweights = weights)

  }
  argmaxw <- function(x,  window){ window[rev(order(x[window]))]}
  X <- object$Xsort
  weights <- object$sortweights
  wh.kadapt <- object$kadapt
  TTSS <- WeightedTestStatist(X, weights, wh.kadapt)
  TS <- TTSS$WTS1 + TTSS$WTS2		# computed at     X_2, X_3, ..., X_m
  TS1 <- TTSS$WTS1					# computed at     X_2, X_3, ..., X_m
  TS <- c(NA, TS)						# computed at X_1, X_2, X_3, ..., X_m
  TS1 <- c(NA, TS1)
  cumsumw <- cumsum(TTSS$sortweights)
  restrictedw <- which(cumsumw >=  object$r1*cumsumw[wh.kadapt] & cumsumw <=  (1-object$r2)*cumsumw[wh.kadapt] )
  w1 <- min(restrictedw)
  w2 <- max(restrictedw)
  window1 <- max(2,    w1  ) 	#define left  end of the testing window in the set X_1, X_2, X_3, ..., X_m
  window2 <- min(wh.kadapt-2,  w2  )
  kadapt.maxlik <- argmaxw(TS1,  window = window1:window2)[1]
  TS.window <- TS[window1:window2]
  TS <- TS[which(!is.na(TS))]
  if(plot == TRUE){
    plot(TS, ylim = c(0, object$CritVal+5))
    if(object$r1*length(TS)-floor(object$r1*length(TS))<= 0.5){
      ind <- floor(object$r1*length(TS)):(floor(object$r1*length(TS)-1)+length(TS.window))
    }else{ind <- round(object$r1*length(TS)):(round(object$r1*length(TS)-1)+length(TS.window))}
    points(ind, TS.window, col = "red")
    abline(h = object$CritVal, col = "blue")
    legend("topleft",c("Critical Value","test statistic inside the window","test statistic"),col=c("blue","red","black"),pch=c("-","o","o"))
    }
  return(list(TS.window = TS.window, TS = TS, CritVal = object$CritVal))
}

#' Goodness of fit test statistics for time series
#'
#' @description Give the results of the goodness of fit test for testing the null hypothesis that the tail is fitted by a Pareto distribution starting from the adaptive threshold (for more details see pages 447 and 448 of Durrieu et al. (2015)).
#'
#' @param object output of the hill.ts function.
#' @param X a vector of the observed values.
#' @param t a vector of time covariates which should have the same length as X.
#' @param plot If \code{TRUE}, the test statistic are plotted.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' \item{TS.window}{the maximum value of test statistics inside the window for each t in Tgrid (see help(hill.ts) ).}
#' \item{TS.max}{the maximum value of test statistics for each t in Tgrid (see help(hill.ts) ).}
#' \item{CritVal}{the critical value of the test.}
#'
#' @references Grama, I. and Spokoiny, V. (2008). Statistics of extremes by oracle estimation. Ann. of Statist., 36, 1619-1648.
#'
#' Durrieu, G. and Grama, I. and Pham, Q. and Tricot, J.- M (2015). Nonparametric adaptive estimator of extreme conditional tail probabilities quantiles. Extremes, 18, 437-478.
#'
#' @seealso \code{\link{hill.ts}}, \code{\link{goftest}}
#'
#' @export
#'
#' @examples
#' theta<-function(t){0.5+0.25*sin(2*pi*t)}
#' n<-5000
#' t<-1:n/n
#' Theta<-theta(t)
#' Data<-NULL
#' Tgrid<-seq(0.01,0.99,0.01)
#' #example with fixed bandwidth
#' for(i in 1:n){Data[i]<-rparetomix(1,a=1/Theta[i],b=5/Theta[i]+5,c=0.75,precision=10^(-5))}
#' \dontrun{ #For computing time purpose
#'   #example
#'   hgrid <- bandwidth.grid(0.009, 0.2, 20, type = "geometric")
#'   TgridCV <- seq(0.01, 0.99, 0.1)
#'   hcv <- bandwidth.CV(Data, t, TgridCV, hgrid, pcv = 0.99,
#'          TruncGauss.kernel, kpar = c(sigma = 1), CritVal = 3.6, plot = TRUE)
#'
#'   Tgrid <- seq(0.01,0.99,0.01)
#'   hillTs <- hill.ts(Data, t, Tgrid, h = hcv$h.cv, TruncGauss.kernel, kpar = c(sigma = 1),
#'                    CritVal = 3.6, gridlen = 100, initprop = 1/10, r1 = 1/4, r2 = 1/20)
#'   goftest(hillTs, Data, t, plot = TRUE)
#'
#'   # we observe that for this data, the null hypothesis that the tail
#'   # is fitted by a Pareto distribution is not rejected
#'   # for all points on the Tgrid
#'
#' }
#'
#'
goftest.hill.ts <- function(object, X, t, plot = FALSE, ...){
  #Function returning the test statistic corresponding to the threshold for each T in Tgrid (output of hillts)
  #hillts is the output of the function hill.ts
  #X.ts is the data
  #t is the time of the data
  #Output :
  #TS.window : the test statistic maximum within the window (r1 and r2 parameters)
  #TS.max : the test statistic maximum without window
  TS.max <- NULL
  WeightedTestStatist <- function(X, weights, m){
    # by Ion Grama (2006)
    # Input:
    # X the observations
    # weights
    # m is the change point (given as the number of the order statistic)
    # Output:
    #this function compute the weighted test statistics WTS1,  WTS2 and the estimators wh1 and wh2
    # the output is a list containing WTS1,  WTS2,  wh1,  wh2

    G.KL <- function(x){x - log(1 + x)}		# define a function to compute Kullback-Leibler information
    #######################################

    # order log(X) and weights
    n <- length(X)
    order <- rev(order(X))
    Xsort <- rev(sort(X))
    lnX <- log(Xsort)
    weights <- weights[order]

    Sweights <- cumsum(weights[1:n])	  # compute the cumulative weights
    SX <- cumsum( weights[1:(n-1)] * lnX[1:(n-1)]) - lnX[2:n] * Sweights[1:(n-1)]   # compute SX at points X_2, X_3, ..., X_n (attention! X_1 is excluded )
    wh1 <- SX[1:(n-1)]/Sweights[1:(n-1)]        # compute wighted Hill's estimator at X_2, X_3, ..., X_n (attention! X_1 is excluded )
    wh2 <- abs(SX[m-1]-SX[1:(m-1)])/(Sweights[m-1]-Sweights[1:(m-1)])
    wh2[m-1] <- NA  # compute weighted wh2  at X_2, X_3, ..., X_m (attention! X_1 is excluded )
    # do not delete the following comments
    # sometimes the difference SX[m-1]-SX[1:(m-1)] can be negative
    # because of the computing error (of order -10^-17)
    # that is why in the above formula we used abs(...)

    ## the following lines are introduced to check the correctness
    ##if (  is.number(SX[m-1]) == F  ){ #  this line is to check correctness
    ##cat("m-1 = ", m-1, "\n") #  this line is to check correctness
    ##cat("SX = ", SX, "\n") #  this line is to check correctness
    ##cat("abs(SX[m-1]-SX[1:(m-1)]) = ", abs(SX[m-1]-SX[1:(m-1)]), "\n") #  this line is to check correctness
    ##cat("Sweights[m-1]-Sweights[1:(m-1)] = ", Sweights[m-1]-Sweights[1:(m-1)], "\n") #  this line is to check correctness
    ##cat("wh2[(m-10):(m-1)] = ", round(wh2[(m-10):(m-1)], 5), "\n") #  this line is to check correctness
    ##cat("Xsort[(m-10):(m-1)] = ", Xsort[(m-10):(m-1)], "\n") #  this line is to check correctness
    ##}  #  this line is to check correctness

    #############################################
    # The previous computation computes wh2 only for a given m and requires a vector of size m
    # In the following lines here we propose a new vector computation of wh2:
    # this computation requires introduction of a matrix of size m x n
    # and compute all the values wh2 at once (this is not used in this version)
    #	Dh = repcol(Sh[1:m], nrep = m)-reprow(Sh[1:m], nrep = n)    # the same computation:  E = expand.grid(Sh[1:m], Sh[1:m]); matrix(E[, 1]-E[, 2], nrow = m)
    #	Dw = repcol(Sw[1:m], nrep = m)-reprow(Sw[1:m], nrep = n)
    #	wh2matr = Dh/Dw
    #	wh2matr[, m]  # to check
    #	wh2          # to check
    ############################################

    TS1 <- Sweights[1:(m-1)] * G.KL(  wh1[1:(m-1)]/wh1[m-1] - 1  )                  # compute weighted TS1  at X_2, X_3, ..., X_m (attention! X_1 is excluded )
    TS2 <- (Sweights[m-1]-Sweights[1:(m-1)]) * G.KL(  wh2[1:(m-1)]/wh1[m-1] - 1  )  # compute weighted TS2  at X_2, X_3, ..., X_m (attention! X_1 is excluded )

    # we can compute normalized TS1 and TS2 if necessary
    #varW = sqrt(mean(weights[1:m]^2))   # the normalisation is not used here
    #TS1 = TS1/varW   # the normalisation is not used here
    #TS2 = TS2/varW # the normalisation is not used here

    ###############
    # the output

    list(WTS1 = TS1,  WTS2 = TS2,  wh1 = wh1[1:(m-1)],  wh2 = wh2[1:(m-1)],  sortweights = weights)

  }
  argmaxw <- function(x,  window){ window[rev(order(x[window]))]}
  tstau <- function(i){
    indices <- which((t>= object$Tgrid[i]-object$h)&(t<= object$Tgrid[i]+object$h))
    seq.X <- X[indices]
    tgrid <- t[indices]
    tx <- (tgrid-object$Tgrid[i])/object$h
    if( is.null(object$kpar) ){
      par <- list(tx)
    }else{
      par <- list(tx, object$kpar)
    }
    weights <- do.call(object$kernel,par)
    wh.kadapt <- object$kadapt[i]
    TTSS <- WeightedTestStatist(seq.X, weights, wh.kadapt)
    TS <- TTSS$WTS1 + TTSS$WTS2		# computed at     X_2, X_3, ..., X_m
    TS1 <- TTSS$WTS1					# computed at     X_2, X_3, ..., X_m
    TS <- c(NA, TS)						# computed at X_1, X_2, X_3, ..., X_m
    TS1 <- c(NA, TS1)
    cumsumw <- cumsum(TTSS$sortweights)
    restrictedw <- which(cumsumw >=  object$r1*cumsumw[wh.kadapt] & cumsumw <=  (1-object$r2)*cumsumw[wh.kadapt] )
    w1 <- min(restrictedw)
    w2 <- max(restrictedw)
    window1 <- max(2,    w1  ) 	#define left  end of the testing window in the set X_1, X_2, X_3, ..., X_m
    window2 <- min(wh.kadapt-2,  w2  )
    kadapt.maxlik <- argmaxw(TS1,  window = window1:window2)[1]
    TS.maxwindow <- max(TS[window1:window2])
    TS.max <- max(TS[which(!is.na(TS))])
    return(list(TS.maxwindow = TS.maxwindow, TS.max = TS.max))
  }
  i <- 1:length(object$Tgrid)
  result <- sapply(i,  tstau)
  if(plot == TRUE){
    plot(object$Tgrid, as.numeric(result[1, ]), ylim = c(0, object$CritVal+2), ylab = "Test statistics", xlab="t")
    abline(h = object$CritVal, col = "blue")
    legend("topleft",c("Critical Value"),col=c("blue"),pch=c("-"))
  }
  return(list(TS.window = as.numeric(result[1, ]), TS.max = as.numeric(result[2, ], CritVal = object$CritVal)))
}

##############################################################################################
# the following functions compute the weighed adaptive quantiles and survival functions

#' Predict the adaptive Survival or Quantile function
#'
#' @description Give the adaptive survival function or quantile function
#'
#' @param object output object of the function hill.adapt.
#' @param newdata optionally, a vector with which to predict. If omitted, the original data points are used.
#' @param type either "quantile" or "survival".
#' @param ... further arguments passed to or from other methods.
#'
#' @details If type = "quantile", \eqn{newdata} must be between 0 and 1. If type = "survival", \eqn{newdata} must be in the domain of the data from the \code{hill.adapt} function.
#'
#' @return The function provide the quantile assiociated to the adaptive model for the probability grid (transformed to -log(1-p) in the output) if type = "quantile". And the survival function assiociated to the adaptive model for the quantile grid if type = "survival".
#'
#' @seealso \code{\link{hill.adapt}}
#'
#' @export
#'
#' @examples
#' x <- rparetoCP(1000)
#'
#' HH <- hill.adapt(x, weights=rep(1, length(x)), initprop = 0.1,
#'                gridlen = 100 , r1 = 0.25, r2 = 0.05, CritVal=10)
#'
#' newdata <- probgrid(p1 = 0.01, p2 = 0.999, length = 100)
#' pred.quantile <- predict(HH, newdata, type = "quantile")
#' newdata <- seq(0, 50, 0.1)
#' pred.survival <- predict(HH, newdata, type = "survival")#survival function
#'
#' #compare the theorical quantile and the adaptive one.
#' predict(HH, 0.9999, type = "quantile")
#' qparetoCP(0.9999)
#'
#' #compare the theorical probability and the adaptive one assiociated to a quantile.
#' predict(HH, 20, type = "survival")
#' 1 - pparetoCP(20)
#'
predict.hill.adapt <- function(object, newdata = NULL, type = "quantile",  ...){
  # Input: hill.adapt is the output of the function hill.adapt
  # grid is a probabilities grid if type = "quantile" or a quantile grid if type = "survival"
  # type is either quantile to estimate the quantiles or survival to estimate the survival function
  # Output: matrix of grid and their quantile/probabilities assiociated

  if(type == "quantile"){
    n <- length(object$Xsort)
    X <- object$Xsort
    Theta <- object$hadapt
    pgrid <- newdata
    if(is.null(pgrid)){pgrid <- 1:(length(object$Xsort)-1)/length(object$Xsort)}
    if(!is.na(Theta)){
      weights <- object$sortweights
      Xsort <- sort(X)
      Xord <- (order(X))
      word <- weights[Xord]
      p0 <- wecdf(Xsort, object$Xadapt, word)
      x1 <- wquantile(Xsort, pgrid[which(pgrid<p0)], word)
      x2 <- object$Xadapt*((1-p0)/(1-pgrid[which(pgrid>= p0)]))^(Theta)
      y <- as.numeric(c(x1, x2))
    }else{
      weights <- object$sortweights
      Xsort <- sort(X)
      x1 <- wquantile(Xsort, pgrid, weights)
      y <- x1
    }
    res <- list(p = pgrid, y = y[!is.na(y)])
    class(res) <- "predict.adapt"
    #attr(res,  "call")  <-  sys.call()
    res
  }else if(type == "survival"){
    n <- length(object$Xsort)
    Xsort <- sort(object$Xsort)
    xgrid <- newdata
    if(is.null(xgrid)){xgrid <- object$Xsort}
    NX <- sort(unique(c(object$Xadapt, xgrid)))
    proGrid <- wecdf(object$Xsort, NX, object$sortweights)
    tau <- object$Xadapt
    Theta <- object$hadapt
    if(!is.na(Theta)){
      x1 <- 1-proGrid[which(NX<tau)]
      x2 <- (1-ppareto(NX[which(NX>= tau)], a = 1/Theta, scale = tau))*(1-proGrid[which(NX == tau)])
      y <- c(x1, x2)
      if(length(xgrid[which(xgrid == tau)]) == 0){
        y <- sort(y[-which(NX == tau)], decreasing = TRUE)
        NX <- NX[-which(NX == tau)]
      }
    }else{
      proGrid <- wecdf(Xsort, NX, weights)
      x1 <- 1-proGrid
      x <- c(x1)
      y <- sort(x, decreasing = TRUE)
    }
    res <- list(x = NX, p = y)
    class(res) <- "predict.adapt"
    #attr(res,  "call")  <-  sys.call()
    res
  }else{
    cat("please choose a type between quantile and survival")
  }
}#end of function predict.hill.adapt

#' Predict the adaptive Survival or Quantile function
#'
#' @description Give the adaptive survival function or quantile function
#'
#' @param object output  object of the function hill.
#' @param newdata optionally, a vector with which to predict. If omitted, the original data points are used.
#' @param type either "quantile" or "survival".
#' @param threshold.rank the rank value for the hill output of the threshold, with default value 0.
#' @param threshold the value of threshold, with default value 0.
#' @param ... further arguments passed to or from other methods.
#'
#' @details If type = "quantile", \eqn{newdata} must be between 0 and 1. If type = "survival", \eqn{newdata} must be in the domain of the data from the \code{hill} function.
#'
#' @return The function provide the quantile assiociated to the adaptive model for the probability grid (transformed to -log(1-p) in the output) if type = "quantile". And the survival function assiociated to the adaptive model for the quantile grid if type = "survival".
#' @export
#'
#' @seealso \code{\link{hill}}
#'
#' @examples
#' x <- abs(rcauchy(100))
#' hh <- hill(x)
#' #example for a fixed value of threshold
#' predict(hh, threshold = 3)
#' #example for a fixed rank value of threshold
#' predict(hh, threshold.rank = 30)
#'
predict.hill <- function(object, newdata = NULL, type = "quantile", threshold.rank = 0, threshold = 0, ...){
  # Input: object is the output of the function hill
  # threshold.rank is the rank we use as a threshold
  # type is either quantile to estimate the quantiles or survival to estimate the survival function
  # Output: matrix of grid and their quantile/probabilities assiociated
  n <- length(object$xsort)
  Xsort <- object$xsort
  Xord <- (order((object$xsort)))
  word <- object$wsort
  if(threshold.rank == 0){
    if(threshold == 0){
      cat("please enter a value for either threshold.rank or threshold. \n")
    }else{
      Theta <- hill(object$xsort, object$wsort, threshold)$hill}
  }else{
    Theta <- object$hill[threshold.rank]
  }
  if(threshold.rank == 0){
    tau <- threshold
  }else{
    tau <- Xsort[threshold.rank]
  }
  if(type == "quantile"){
    p0 <- wecdf(Xsort, tau, word)
    pgrid <- newdata
    if(is.null(pgrid)){pgrid <- 1:(length(object$xsort)-1)/length(object$xsort)}
    #pgrid <- sort(c(p0, pgrid))
    x1 <- wquantile(Xsort, pgrid[which(pgrid<p0)], word)
    x2 <- tau*((1-p0)/(1-pgrid[which(pgrid>= p0)]))^(Theta)
    y <- c(x1, x2)
    return(list(p = pgrid, y = y[!is.na(y)]))
  }
  if(type == "survival"){
    xgrid <- newdata
    if(is.null(xgrid)){xgrid <- object$xsort}
    NX <- sort(unique(c(tau, xgrid)))
    proGrid <- wecdf(Xsort, NX, word)
    x1 <- 1-proGrid[which(NX<= tau)]
    x2 <- (1-ppareto(NX[which(NX>tau)], a = 1/Theta, scale = tau))*(1-proGrid[which(NX == tau)])
    y <- c(x1, x2)
    if(length(xgrid[which(xgrid == tau)]) == 0){
      y <- sort(y[-which(NX == tau)], decreasing = TRUE)
      NX <- NX[-which(NX == tau)]
    }
    return(list(x = NX, p = y))
  }else{cat("please choose a type between quantile and survival")}
}#end of function predict.hill

#' Predict the adaptive Survival or Quantile function for a time serie
#'
#' @description  Give the adaptive survival function or quantile function of a time serie
#'
#' @param object output object of the function hill.ts.
#' @param newdata optionally, a vector with which to predict. If omitted, the original data points are used.
#' @param type either "quantile" or "survival".
#' @param ... further arguments passed to or from other methods.
#'
#' @details If type = "quantile", \eqn{newdata} must be between 0 and 1. If type = "survival", \eqn{newdata} must be in the domain of the data from the function \code{hill.ts}.
#'
#' @return
#'   \item{p}{the input vector of probabilities.}
#'   \item{x}{the input vector of values.}
#'   \item{Tgrid}{Tgrid output of the function hill.ts.}
#'   \item{quantiles}{the estimted quantiles assiociated to newdata.}
#'   \item{survival}{the estimated survival function assiociated to newdata.}
#'
#' @seealso \code{\link{hill.ts}}
#'
#' @export
#'
#' @examples
#' #Generate a pareto mixture sample of size n with a time varying parameter
#' theta <- function(t){
#'    0.5+0.25*sin(2*pi*t)
#'  }
#' n <- 4000
#' t <- 1:n/n
#' Theta <- theta(t)
#' Data <- NULL
#' set.seed(1240)
#' for(i in 1:n){
#'    Data[i] <- rparetomix(1, a = 1/Theta[i], b = 1/Theta[i]+5, c = 0.75, precision = 10^(-5))
#'  }
#' \dontrun{ #For computing time purpose
#'   #choose the bandwidth by cross validation
#'   Tgrid <- seq(0, 1, 0.1)#few points to improve the computing time
#'   hgrid <- bandwidth.grid(0.01, 0.2, 20, type = "geometric")
#'   hcv <- bandwidth.CV(Data, t, Tgrid, hgrid, TruncGauss.kernel,
#'          kpar = c(sigma = 1), pcv = 0.99, CritVal = 3.6, plot = TRUE)
#'   h.cv <- hcv$h.cv
#'
#'   #we modify the Tgrid to cover the data set
#'   Tgrid <- seq(0, 1, 0.02)
#'   hillTs <- hill.ts(Data, t, Tgrid, h = h.cv, TruncGauss.kernel, kpar = c(sigma = 1),
#'            CritVal = 3.6, gridlen = 100, initprop = 1/10, r1 = 1/4, r2 = 1/20)
#'   p <- c(0.999)
#'   pred.quantile.ts <- predict(hillTs, newdata = p, type = "quantile")
#'   true.quantile <- NULL
#'   for(i in 1:n){
#'      true.quantile[i] <- qparetomix(p, a = 1/Theta[i], b = 1/Theta[i]+5, c = 0.75)
#'    }
#'   plot(Tgrid, log(as.numeric(pred.quantile.ts$y)),
#'        ylim = c(0, max(log(as.numeric(pred.quantile.ts$y)))), ylab = "log(0.999-quantiles)")
#'   lines(t, log(true.quantile), col = "red")
#'   lines(t, log(Data), col = "blue")
#'
#'
#'   #comparison with other fixed bandwidths
#'
#'   plot(Tgrid, log(as.numeric(pred.quantile.ts$y)),
#'        ylim = c(0, max(log(as.numeric(pred.quantile.ts$y)))), ylab = "log(0.999-quantiles)")
#'   lines(t, log(true.quantile), col = "red")
#'
#'   hillTs <- hill.ts(Data, t, Tgrid, h = 0.1, TruncGauss.kernel, kpar = c(sigma = 1),
#'                     CritVal = 3.6, gridlen = 100,initprop = 1/10, r1 = 1/4, r2 = 1/20)
#'   pred.quantile.ts <- predict(hillTs, p, type = "quantile")
#'   lines(Tgrid, log(as.numeric(pred.quantile.ts$y)), col = "green")
#'
#'
#'   hillTs <- hill.ts(Data, t, Tgrid, h = 0.3, TruncGauss.kernel, kpar = c(sigma = 1),
#'                CritVal = 3.6, gridlen = 100, initprop = 1/10, r1 = 1/4, r2 = 1/20)
#'   pred.quantile.ts <- predict(hillTs, p, type = "quantile")
#'   lines(Tgrid, log(as.numeric(pred.quantile.ts$y)), col = "blue")
#'
#'
#'   hillTs <- hill.ts(Data, t, Tgrid, h = 0.04, TruncGauss.kernel, kpar = c(sigma = 1),
#'              CritVal = 3.6, gridlen = 100, initprop = 1/10, r1 = 1/4, r2 = 1/20)
#'   pred.quantile.ts <- predict(hillTs ,p, type = "quantile")
#'   lines(Tgrid, log(as.numeric(pred.quantile.ts$y)), col = "magenta")
#' }
#'
#'
#'
predict.hill.ts <- function(object, newdata = NULL, type = "quantile", ...){
  # By Kevin Jaunatre 2015
  # Input : X is a vector of data
  # hill.ts is the output of the function hill.ts
  # grid is a probabilities grid if type = "quantile" or a quantile grid if type = "survival"
  # type is either quantile to estimate the quantiles or survival to estimate the survival function
  # Output : list with grid and the survival or quantile function assiociated for each t
  X <- object$X
  t <- object$t
  if(type == "quantile"){
    n <- length(X)
    Theta <- object$Theta
    tau <- object$threshold
    Tgrid <- object$Tgrid
    kernel <- object$kernel
    pgrid <- newdata
    if(is.null(pgrid)){pgrid <- 1:(length(X)-1)/length(X)}
    if(length(object$h) == 1){
      h <- rep(object$h, length(Tgrid))
    }else{h <- object$h}
    y <- NULL
    for(i in 1:length(Tgrid)){
      if(!is.na(Theta[i])){
        indices <- which((t>= Tgrid[i]-h[i])&(t<= Tgrid[i]+h[i]))
        t.ind <- t[indices]
        tx <- (t.ind-Tgrid[i])/h[i]
        if( is.null(object$kpar) ){
          par <- list(tx)
        }else{
          par <- list(tx, object$kpar)
        }
        weights <- do.call(kernel,par)
        Xt <- X[indices]
        p0 <- wecdf(Xt, tau[i], weights)
        #pgrid <- sort(c(p0, pgrid))
        x1 <- wquantile(Xt, pgrid[which(pgrid<= p0)], weights)
        x2 <- tau[i]*((1-p0)/(1-pgrid[which(pgrid>p0)]))^(Theta[i])
        y <- cbind(y, c(x1, x2))
      }else{
        indices <- which((t>= Tgrid[i]-h[i])&(t<= Tgrid[i]+h[i]))
        t.ind <- t[indices]
        tx <- (t.ind-Tgrid[i])/h[i]
        if( is.null(object$kpar) ){
          par <- list(tx)
        }else{
          par <- list(tx, object$kpar)
        }
        weights <- do.call(kernel,par)
        Xt <- X[indices]
        x1 <- wquantile(Xt, pgrid, weights)
        y <- cbind(y, x1)
      }
    }
    return(list(p = pgrid, Tgrid = Tgrid, y = y))
  }
  if(type == "survival"){
    n <- length(X)
    tau <- object$threshold
    Theta <- object$Theta
    Tgrid <- object$Tgrid
    kernel <- object$kernel
    xgrid <- newdata
    if(is.null(xgrid)){xgrid <- X}
    if(length(object$h) == 1){
      h <- rep(object$h, length(t))
    }else{h <- object$h}
    y <- NULL
    for(i in 1:length(Tgrid)){
      if(!is.na(Theta[i])){
        NX <- sort(unique(c(tau[i], xgrid)))
        indices <- which((t>= Tgrid[i]-h[i])&(t<= Tgrid[i]+h[i]))
        t.ind <- t[indices]
        tx <- (t.ind-Tgrid[i])/h[i]
        if( is.null(object$kpar) ){
          par <- list(tx)
        }else{
          par <- list(tx, object$kpar)
        }
        weights <- do.call(kernel,par)
        proGrid <- wecdf(X[indices], NX, weights)
        x1 <- 1-proGrid[which(NX<= tau[i])]
        x2 <- (1-ppareto(NX[which(NX>tau[i])], a = 1/Theta[i], scale = tau[i]))*(1-proGrid[which(NX == tau[i])])
        x <- c(x1, x2)
        if((length(xgrid[which(xgrid == tau[i])]) == 0)){
          y <- cbind(y, sort(x[-which(NX == tau[i])], decreasing = TRUE))
        }else{y <- cbind(y, sort(x, decreasing = TRUE))}
      }else{
        NX <- sort(xgrid)
        indices <- which((t>= Tgrid[i]-h[i])&(t<= Tgrid[i]+h[i]))
        t.ind <- t[indices]
        tx <- (t.ind-Tgrid[i])/h[i]
        if( is.null(object$kpar) ){
          par <- list(tx)
        }else{
          par <- list(tx, object$kpar)
        }
        weights <- do.call(kernel,par)
        proGrid <- wecdf(X[indices], NX, weights)
        x1 <- 1-proGrid
        x <- c(x1)
        y <- cbind(y, sort(x, decreasing = TRUE))
      }
    }
    return(list(x = sort(xgrid), Tgrid = Tgrid, p = y))
  }else{cat("please choose a type between quantile and survival")}
}#end of the function predict.hill.ts

#############################################################################################
# the following functions are plotting results

#' Hill Plot
#'
#' @description Graphical representation of the hill estimator.
#'
#' @param x output object of the function hill.
#' @param xaxis either "ranks" or "xsort".
#' @param ... further arguments passed to or from other methods.
#'
#' @details If xaxis="ranks", the function draws the Hill estimators for each ranks of the grid output of the function hill.
#' If xaxis="xsort", the function draws the Hill estimators for each data of the grid output of the function hill.
#'
#' @seealso \code{\link{hill}}
#'
#'
#' @export
#'
#' @examples
#' x <- abs(rcauchy(100))
#' hh <- hill(x)
#' par(mfrow = c(2, 1))
#' plot(hh, xaxis = "ranks")
#' plot(hh, xaxis = "xsort")
#'
plot.hill <- function(x, xaxis = "ranks", ...){
  if(xaxis == "ranks"){
    plot(x$hill, ...)
  }else if(xaxis == "xsort"){
    plot(log(x$grid), x$hill, ...)
  }
}


##############################################################################################
# The following functions are helping us on the choice of the critical value and the bandwidth

#' Computation of the critical value in the hill.adapt function
#'
#' @description For a given kernel function, compute the critical value (CritVal) of the test statistic in the hill.adapt function by Monte-Carlo simulations.
#'
#' @param NMC the number of Monte-Carlo simulations.
#' @param n the sample size.
#' @param kernel a kernel function for which the critical value is computed. The available kernel functions are Epanechnikov, Triangular, Truncated Gaussian, Biweight and Rectangular. The truncated gaussian kernel is by default.
#' @param kpar a value for the kernel function parameter, with no default value.
#' @param prob a vector of type 1 errors.
#' @param gridlen,initprop,r1,r2 parameters used in the function hill.adapt (see \code{\link{hill.adapt}}).
#' @param plot If \code{TRUE}, the empirical cummulative distribution function and the critical values are plotted.
#'
#' @return For the type 1 errors \eqn{prob}, this function returns the critical values.
#' @export
#'
#' @seealso \code{\link{hill.adapt}}
#'
#' @references Durrieu, G. and Grama, I. and Pham, Q. and Tricot, J.- M (2015). Nonparametric adaptive estimator of extreme conditional tail probabilities quantiles. Extremes, 18, 437-478.
#'
#' @examples
#' n <- 1000
#' NMC <- 500
#' prob <- c(0.99)
#' \dontrun{ #For computing time purpose
#'   CriticalValue(NMC, n, TruncGauss.kernel, kpar = c(sigma = 1), prob, gridlen = 100 ,
#'                 initprop = 1/10, r1 = 1/4, r2 = 1/20, plot = TRUE)
#' }
#'
CriticalValue <- function(NMC, n, kernel = TruncGauss.kernel,kpar = NULL, prob = 0.95, gridlen = 100, initprop = 0.1, r1 = 0.25, r2 = 0.05, plot = FALSE){
  # By Kevin Jaunatre 2015
  # Input : NMC is the number of Monte-Carlo iterations we want to do
  # n is the length of the pareto sample
  # kernel is the kernel function on which we want to calculate the critical value
  # prob is a vector of p for each we want to calculate the critical value
  # warning : each value must be between max(hl)+step (ex: 0.01) and 1-max(hl)-step
  # initprop
  # gridlen
  # r1 is the proportion to skip from the right in the test statistic
  # r1 is the proportion to skip from the left in the test statistic
  # Output : critical value of the kernel function for each p
  if( identical(kernel,TruncGauss.kernel) & is.null(kpar) ){
    kpar = c(sigma = 1)
  }

  t <- n/2
  h <- n/2
  TS <- NULL
  MCfunction <- function(Nmc){
    x <- rpareto(n, 1)
    Tgrid <- (t-h):(t+h)
    tx <- (Tgrid-t)/h
    if( is.null(kpar) ){
      par <- list(tx)
    }else{
      par <- list(tx, kpar)
    }
    weights <- do.call(kernel, par)
    xgrid <- x[Tgrid]
    hill <- hill.adapt(xgrid,  weights, CritVal = Inf, gridlen = gridlen, initprop = initprop,  r1 = r1,  r2 = r2)
    TS <- max(hill$TS.max[which(!is.na(hill$TS.max))])
    return(TS)
  }
  Nmc <- 1:NMC
  TS <- sapply(Nmc,  MCfunction)
  p <- 1:NMC/NMC
  if(plot == TRUE){
    plot(sort(TS), p, type = "l", xlab = "Test Statistic", ylab = "Empirical distribution function")
    for(i in 1:length(prob)){abline(h = prob[i], col = i,lty=2)}
    for(i in 1:length(prob)){abline(v = wquantile(TS, prob[i]), col = i,lty=2)}
    #legend("topleft",  paste("p = ", prob), col = 1:length(prob), pch = "-")
  }
  return(wquantile(TS, prob))
}#end of CriticalValue



#' Pointwise confidence intervals by bootstrap
#'
#' @description Pointwise quantiles and survival probabilities confidence intervals using bootstrap.
#'
#' @param X a numeric vector of data values.
#' @param weights a numeric vector of weights.
#' @param probs used if type = "quantile", a numeric vector of probabilities with values in \eqn{[0,1]}.
#' @param xgrid used if type = "survival", a numeric vector with values in the domain of X.
#' @param B an integer giving the number of bootstrap iterations.
#' @param alpha the type 1 error of the bootstrap (1-\eqn{alpha})-confidence interval.
#' @param type type is either "quantile" or "survival".
#' @param CritVal a critical value associated to the kernel function given by \code{\link{CriticalValue}}. The default value is 10 corresponding to the rectangular kernel.
#' @param gridlen,initprop,r1,r2 parameters used in the function hill.adapt (see \code{\link{hill.adapt}}).
#' @param plot If \code{TRUE}, the bootstrap confidence interval is plotted.
#'
#' @details Generate B samples of \eqn{X} with replacement to estimate the quantiles of orders \eqn{probs} or the survival probability corresponding to \eqn{xgrid}. Determine the bootstrap pointwise (1-\eqn{alpha})-confidence interval for the quantiles or the survival probabilities.
#'
#' @return
#'   \item{LowBound}{the lower bound of the bootstrap (1-\eqn{alpha})-confidence interval.}
#' \item{UppBound}{the upper bound of the bootstrap (1-\eqn{alpha})-confidence interval of level.}
#'
#' @seealso \code{\link{hill.adapt}},\code{\link{CriticalValue}},\code{\link{predict.hill.adapt}}
#'
#' @export
#'
#' @examples
#' X <- abs(rcauchy(400))
#' hh <- hill.adapt(X)
#' probs <- probgrid(0.1, 0.999999, length = 100)
#' B <- 200
#' \dontrun{ #For computing time purpose
#'   bootCI(X, weights = rep(1, length(X)), probs = probs, B = B, plot = TRUE)
#'   xgrid <- sort(sample(X, 100))
#'   bootCI(X, weights = rep(1, length(X)), xgrid = xgrid, type = "survival", B = B, plot = TRUE)
#' }
#'
bootCI <- function(X, weights = rep(1, length(X)), probs = 1:(length(X)-1)/length(X), xgrid = sort(X), B = 100, alpha = 0.05, type = "quantile", CritVal = 10, initprop = 1/10,  gridlen = 100, r1 = 1/4,  r2 = 1/20, plot = F){
  hx <- hill.adapt(X, weights, initprop = initprop, gridlen = gridlen, r1 = r1,  r2 = r2,  CritVal = CritVal)
  if(type == "quantile"){
    pred <- matrix(0, ncol = B, nrow = length(probs))
    for(b in 1:B){
      ind <- sample(1:length(X), length(X), replace = TRUE)
      xb <- X[ind]
      wb <- weights[ind]
      hm <- hill.adapt(xb, wb, initprop = initprop, gridlen = gridlen, r1 = r1,  r2 = r2,  CritVal = CritVal)
      pred[, b] <- as.numeric(predict(hm, newdata = probs, type = "quantile")$y)
    }
    quantInf <- NULL
    quantSup <- NULL
    for(j in 1:length(probs)){
      quantInf[j] <- quantile(pred[j, ], alpha/2)
      quantSup[j] <- quantile(pred[j, ], 1-alpha/2)
    }
    if(plot == T){
      pred.true <- predict(hx, newdata = probs, type = "quantile")
      plot(-log(1-pred.true$p), log(pred.true$y), type = "l", xlab = "-log(1-p)", ylab = "log(quantiles)")
      lines(-log(1-probs), log(quantInf), col = "green")
      lines(-log(1-probs), log(quantSup), col = "green")
    }
    return(list(LowBound = as.numeric(quantInf), UppBound = as.numeric(quantSup)))
  }else if(type == "survival"){
    pred <- matrix(0, ncol = B, nrow = length(xgrid))
    for(b in 1:B){
      ind <- sample(1:length(X), length(X), replace = TRUE)
      xb <- X[ind]
      wb <- weights[ind]
      hm <- hill.adapt(xb, wb, initprop = initprop, gridlen = gridlen, r1 = r1,  r2 = r2,  CritVal = CritVal)
      pred[, b] <- as.numeric(predict(hm, newdata = xgrid, type = "survival")$p)
    }
    survInf <- NULL
    survSup <- NULL
    for(j in 1:length(xgrid)){
      survInf[j] <- quantile(pred[j, ], alpha/2)
      survSup[j] <- quantile(pred[j, ], 1-alpha/2)
    }
    if(plot == T){
      pred.true <- predict(hx, newdata = xgrid, type = "survival")
      plot(log(pred.true$x), pred.true$p, type = "l", xlab = "log(x)", ylab = "Survival probability")
      lines(log(xgrid), survInf, col = "green")
      lines(log(xgrid), survSup, col = "green")
    }
    return(list(LowBound = as.numeric(survInf), UppBound = as.numeric(survSup)))
  }
}


#' Pointwise confidence intervals by bootstrap
#'
#' @description Pointwise quantiles and survival probabilities confidence intervals using bootstrap.
#'
#' @param X a vector of the observed values.
#' @param t a vector of time covariates which should have the same length as X.
#' @param Tgrid a sequence of times used to perform the cross validation (can be any sequence in the interval \code{[min(t) , max(t)]} ).
#' @param h a bandwidth value (vector values are not admitted).
#' @param kernel a kernel function used to compute the weights in the time domain, with default the truncated gaussian kernel.
#' @param kpar a value for the kernel function parameter, with no default value.
#' @param prob used if type = "quantile", a scalar value in \eqn{[0,1]} which determines the quantile order (vector values are not admitted).
#' @param threshold used if type = "survival", a scalar value in the domain of X.
#' @param B an integer giving the number of bootstrap iterations.
#' @param alpha the type 1 error of the bootstrap (1-\eqn{alpha})-confidence interval.
#' @param type type is either "quantile" or "survival".
#' @param CritVal a critical value associated to the kernel function given by \code{\link{CriticalValue}}. The default value is 3.6 corresponding to the truncated Gaussian kernel.
#' @param gridlen,initprop,r1,r2 parameters used in the function hill.adapt (see \code{\link{hill.adapt}}).
#' @param plot If \code{TRUE}, the bootstrap confidence interval is plotted.
#'
#' @details For each point in \eqn{Tgrid}, generate B samples of \eqn{X} with replacement to estimate the quantile of order \eqn{prob} or the survival probability beyond \eqn{threshold}. Determine the bootstrap pointwise (1-\eqn{alpha})-confidence interval for the quantiles or the survival probabilities.
#'
#' The kernel implemented in this packages are : Biweight kernel, Epanechnikov kernel, Rectangular kernel, Triangular kernel and the truncated Gaussian kernel.
#'
#' @return
#'   \item{LowBound}{the lower bound of the bootstrap (1-\eqn{alpha})-confidence interval.}
#' \item{UppBound}{the upper bound of the bootstrap (1-\eqn{alpha})-confidence interval of level.}
#'
#' @seealso \code{\link{hill.ts}},\code{\link{predict.hill.ts}}, \code{\link{Biweight.kernel}}, \code{\link{Epa.kernel}}, \code{\link{Rectangular.kernel}}, \code{\link{Triang.kernel}}, \code{\link{TruncGauss.kernel}}
#'
#' @section Warning:
#' The executing time of the function can be time consuming if the B parameter or the sample size are high (B=100 and the sample size = 5000 for example) .
#'
#' @export
#'
#' @examples
#' theta <- function(t){
#'    0.5+0.25*sin(2*pi*t)
#'  }
#' n <- 5000
#' t <- 1:n/n
#' Theta <- theta(t)
#' set.seed(123)
#' Data <- NULL
#' for(i in 1:n){
#'    Data[i] <- rparetomix(1, a = 1/Theta[i], b = 1/Theta[i]+5, c = 0.75)
#'  }
#' Tgrid <- seq(1, length(Data)-1, length = 20)/n
#' h <- 0.1
#' \dontrun{ #For computing time purpose
#'   bootCI.ts(Data, t, Tgrid, h, kernel = TruncGauss.kernel, kpar = c(sigma = 1),
#'             CritVal = 3.6, threshold = 2, type = "survival", B = 100, plot = TRUE)
#'   true.p <- NULL
#'   for(i in 1:n){
#'      true.p[i] <- 1-pparetomix(2, a = 1/Theta[i], b = 1/Theta[i]+5, c = 0.75)
#'    }
#'   lines(t, true.p, col = "red")
#'   bootCI.ts(Data, t, Tgrid, h, kernel = TruncGauss.kernel, kpar = c(sigma = 1),
#'  prob = 0.999, type = "quantile", B = 100, plot = TRUE)
#'   true.quantile <- NULL
#'   for(i in 1:n){
#'      true.quantile[i] <- qparetomix(0.999, a = 1/Theta[i], b = 1/Theta[i]+5, c = 0.75)
#'    }
#'   lines(t, log(true.quantile), col = "red")
#' }
#'
#'
bootCI.ts <- function(X, t, Tgrid, h, kernel =  TruncGauss.kernel, kpar = NULL, prob = 0.99, threshold = quantile(X, 0.99), B = 100, alpha = 0.05, type = "quantile", CritVal = 3.6, initprop = 1/10,  gridlen = 100, r1 = 1/4,  r2 = 1/20, plot = F){


  hh <- hill.ts(X, t, Tgrid, h, kernel = kernel, kpar = kpar, initprop = initprop, gridlen = gridlen, r1 = r1,  r2 = r2,  CritVal = CritVal)

  if( identical(kernel,TruncGauss.kernel) & is.null(kpar) ){
    kpar = c(sigma = 1)
  }

  if(type == "quantile"){
    Qpred <- matrix(0, ncol = B, nrow = length(Tgrid))
    for(j in 1:length(Tgrid)){
      x <- X[which((t>= Tgrid[j]-h)&(t<= Tgrid[j]+h))]
      tx <- ( t[which((t >= Tgrid[j] - h)&(t <= Tgrid[j] + h))] - Tgrid[j] ) / h
      if( is.null(kpar) ){
        par <- list(tx)
      }else{
        par <- list(tx, kpar)
      }
      weights <- do.call(kernel, par)
      for(b in 1:B){
        ind <- sample(1:length(x), length(x), replace = TRUE)
        xb <- x[ind]
        wb <- weights[ind]
        hm <- hill.adapt(xb, wb, initprop = initprop, gridlen = gridlen, r1 = r1,  r2 = r2,  CritVal = CritVal)
        Qpred[j, b] <- as.numeric(predict(hm, newdata = prob, type = "quantile")$y)
      }
    }
    quantInf <- NULL
    quantSup <- NULL
    #quantEmpInf <- NULL
    #quantEmpSup <- NULL
    for(jj in 1:length(Tgrid)){
      quantInf[jj] <- quantile(Qpred[jj, ], alpha/2)
      quantSup[jj] <- quantile(Qpred[jj, ], 1-alpha/2)
      #quantEmpInf[j] <- quantile(Qemp[j, ], alpha)
      #quantEmpSup[j] <- quantile(Qemp[j, ], 1-alpha)
    }
    if(plot == T){
      pred.true <- predict(hh, newdata = prob, type = "quantile")
      plot(Tgrid, log(as.numeric(pred.true$y)), type = "l", ylim = c(0, round(max(log(quantSup)))), xlab = "Tgrid", ylab = "log(quantiles)")
      lines(Tgrid, log(quantInf), col = "green")
      lines(Tgrid, log(quantSup), col = "green")
      #lines(-log(1-pg), log(quantEmpInf), col = "blue", lty = 2)
      #lines(-log(1-pg), log(quantEmpSup), col = "blue", lty = 2)
    }
    return(list(LowBound = as.numeric(quantInf), UppBound = as.numeric(quantSup)))
  }else if(type == "survival"){
    Ppred <- matrix(0, ncol = B, nrow = length(Tgrid))
    for(j in 1:length(Tgrid)){
      x <- X[which((t>= Tgrid[j]-h)&(t<= Tgrid[j]+h))]
      tx <- ( t[which((t >= Tgrid[j] - h)&(t <= Tgrid[j] + h))] - Tgrid[j] ) / h
      if( is.null(kpar) ){
        par <- list(tx)
      }else{
        par <- list(tx, kpar)
      }
      weights <- do.call(kernel, par)
      for(b in 1:B){
        ind <- sample(1:length(x), length(x), replace = TRUE)
        xb <- x[ind]
        wb <- weights[ind]
        hm <- hill.adapt(xb, wb, initprop = initprop, gridlen = gridlen, r1 = r1,  r2 = r2,  CritVal = CritVal)
        Ppred[j, b] <- as.numeric(predict(hm, newdata = threshold, type = "survival")$p)
      }
   }
    SurvInf <- NULL
    SurvSup <- NULL
    #quantEmpInf <- NULL
    #quantEmpSup <- NULL
    for(jj in 1:length(Tgrid)){
      SurvInf[jj] <- quantile(Ppred[jj, ], alpha/2)
      SurvSup[jj] <- quantile(Ppred[jj, ], 1-alpha/2)
      #quantEmpInf[j] <- quantile(Qemp[j, ], alpha)
      #quantEmpSup[j] <- quantile(Qemp[j, ], 1-alpha)
    }
    if(plot == T){
      pred.true <- predict(hh, newdata = threshold, type = "survival")
      plot(Tgrid, as.numeric(pred.true$p), type = "l", ylim = c(min(SurvInf)-0.1, max(SurvSup)+0.1), xlab = "Tgrid", ylab = "Survival probability")
      lines(Tgrid, SurvInf, col = "green")
      lines(Tgrid, SurvSup, col = "green")
      #lines(-log(1-pg), log(quantEmpInf), col = "blue", lty = 2)
      #lines(-log(1-pg), log(quantEmpSup), col = "blue", lty = 2)
     }
    return(list(LowBound = as.numeric(SurvInf), UppBound = as.numeric(SurvSup)))
  }
}



#' Probability grid
#'
#' @description Create a geometric grid of probabilities
#'
#' @param p1 the first element of the grid.
#' @param p2 the last element of the grid.
#' @param length the length of the grid.
#'
#' @details Create a geometric grid of length \eqn{length} between \eqn{p1} and \eqn{p2}.The default value of \eqn{length} is \eqn{50}.
#'
#' @return A vector of probabilities between \eqn{p1} and \eqn{p2} and length \eqn{length}.
#' @export
#'
#' @examples
#' p1 <- 0.01
#' p2 <- 0.99
#' length <- 500
#' pgrid <- probgrid(p1, p2, length)
#'
probgrid = function(p1, p2, length = 50){
  #p1 starting point of the grid
  #p2 end point of the grod
  u1 <- -log(1-p1)
  u2 <- -log(1-p2)
  ugrid <- seq(u1, u2, length = length)
  probs <- 1-exp(-ugrid)
  probs
}

#' Choice of the bandwidth by cross validation.
#'
#' @description Choose a bandwidth by minimizing the cross validation function.
#'
#' @param X a vector of the observed values.
#' @param t a vector of time covariates which should have the same length as X.
#' @param Tgrid a sequence of times used to perform the cross validation (can be any sequence in the interval \code{[min(t) , max(t)]} ).
#' @param hgrid a sequence of values from which the bandwidth is selected.
#' @param pcv a probability value which determines the level of quantiles used to perform the cross validation, with default 0.99.
#' @param kernel a kernel function used to compute the weights in the time domain, with default the truncated gaussian kernel.
#' @param kpar a value for the kernel function parameter, with no default value.
#' @param CritVal a critical value associated to the kernel function computed from the function \code{CriticalValue}, with default 3.6 corresponding to the truncated Gaussian kernel.
#' @param plot If \code{TRUE}, the cross validation function is plotted.
#'
#' @details The sequence \eqn{hgrid} must be geometric. (see \code{\link{bandwidth.grid}} to generate a geometric grid of bandwidths).
#'
#' The value \eqn{pcv} should be scalar (vector values are not admitted).
#'
#'
#'
#' @return
#' \item{hgrid}{the sequence of bandwidth given in input.}
#' \item{CV}{the values of the cross validation function for \code{hgrid}.}
#' \item{h.cv}{the bandwidth that minimizes the cross-validation function.}
#'
#' @author Durrieu, G., Grama, I., Jaunatre, K., Pham, Q. and Tricot, J.- M
#'
#' @references Durrieu, G. and Grama, I. and Pham, Q. and Tricot, J.- M (2015). Nonparametric adaptive estimator of extreme conditional tail probabilities quantiles. Extremes, 18, 437-478.
#'
#' @seealso \code{\link{bandwidth.grid}} , \code{\link{CriticalValue}}
#'
#' @export
#'
#' @examples
#' #Generate the data
#' theta <- function(t){
#'    0.5+0.25*sin(2*pi*t)
#'  }
#' n <- 5000
#' t <- 1:n/n
#' Theta <- theta(t)
#' Data <- NULL
#' for(i in 1:n){
#'    Data[i] <- rparetomix(1, a = 1/Theta[i], b = 1/Theta[i]+5, c = 0.75)
#'  }
#'
#' #compute the cross validation bandwidth
#' Tgrid <- seq(0, 1, 0.02) #define a grid to perform the cross validation
#' hgrid <- bandwidth.grid(0.1, 0.3, 20) #define a grid of bandwidths
#' \dontrun{ #For computation time purpose
#'   Hcv <- bandwidth.CV(Data, t, Tgrid, hgrid, pcv = 0.99, plot = TRUE)
#'   #The computing time can be long
#'   Hcv
#' }
#'
#'
bandwidth.CV <- function(X, t, Tgrid, hgrid, pcv = 0.99, kernel = TruncGauss.kernel, kpar = NULL, CritVal = 3.6, plot = FALSE){

  # Input : X is a vector of data
  # hgrid is the vector of bandwidth that we want to compare
  # pcv is the vector of p that we want to estimate the quantile
  # Tgrid is the grid of time t we want to calculate the cross validation function
  # Tmax is the maximum t of the data
  # CritVal: Critical value assiociated to the kernel function
  # plot : give a graphic output of the cross validation
  # Output : matrix of value for the cross validation fucntion for each hgrid and each pcv
  if( identical(kernel,TruncGauss.kernel) & is.null(kpar) ){
    kpar = c(sigma = 1)
  }

  repcol <- function(x, nrep){
    # by Ion Grama (2006)
    # Input:
    # x a vector
    # nrep the number of repetitions
    # Ouput:
    # a matrix matrix ontained by repeating nrep times a column
    x%*%t(rep(1, nrep))
  }
  n <- length(X)
  EmpQuant <- function(Tgrid){
    quant <- NULL
    for(l in 1:length(hgrid)){
      Xhgrid <- X[which((t>= Tgrid-hgrid[l])&(t<= Tgrid+hgrid[l]))]
      tx <- (Tgrid-t[which((t>= Tgrid-hgrid[l])&(t<= Tgrid+hgrid[l]))])/hgrid[l]
      if( is.null(kpar) ){
        par <- list(tx)
      }else{
        par <- list(tx, kpar)
      }
      weights <- do.call(kernel, par)
      quant[l] <- wquantile(Xhgrid[!is.na(Xhgrid)], pcv, weights[!is.na(Xhgrid)])
    }
    return(quant)
  }
  quant <- sapply(Tgrid,  EmpQuant)
  quant <- t(quant)
  CVCalcul <- function(m){
    hm <- hgrid[m]
    qres <- NULL
    for(i in 1:length(Tgrid)){
      Xtgrid <- X[which((t>= round(Tgrid[i]-hm, 10))&(t<= round(Tgrid[i]+hm, 10)))]
      tx <- (Tgrid[i]-t[which((t>= round(Tgrid[i]-hm, 10))&(t<= round(Tgrid[i]+hm, 10)))])/hm
      if( is.null(kpar) ){
        par <- list(tx)
      }else{
        par <- list(tx, kpar)
      }
      Wtgrid <- do.call(kernel, par)
      QQ <- wquantile(Xtgrid[!is.na(Xtgrid)], pcv, Wtgrid[!is.na(Xtgrid)])
      ind <- which(Xtgrid == QQ)
      Wtgrid <- Wtgrid[-ind]
      Xtgrid <- Xtgrid[-ind]
      HH <- hill.adapt(Xtgrid,  weights = Wtgrid,  initprop = 0.1,  gridlen = min(length(Xtgrid), 100) ,  r1 = 0.25,  r2 = 0.05,  CritVal = CritVal)
      qres[i] <- predict(HH, pcv, type = "quantile")$y
    }
    Qres <- repcol(qres, length(hgrid))
    CV <- sum(abs(log(Qres/quant)))/length(hgrid)/length(Tgrid)
    return(CV)
  }
  m  <- 1:length(hgrid)
  CV <- sapply(m,  CVCalcul)
  if(plot == TRUE){
    plot(hgrid, CV, type = "l")
  }
  return(list(hgrid = hgrid, CV = CV, h.cv = hgrid[which(CV == min(CV))]))
}#end of bandwidth.CV


#' Bandwidth Grid
#'
#' @description Create either a geometric or an uniform grid of bandwidths between two values
#'
#' @param hmin the minimum value of the grid.
#' @param hmax the maximum value of the grid.
#' @param length the length of the grid.
#' @param type the type of grid, either \eqn{geometric} or \eqn{uniform}.
#'
#' @details The geometric grid is defined by:
#' \deqn{
#'   grid(l) = hmin * exp( ( log(hmax)-log(hmin) ) / (length -1) ) ^ l  ,  l = 0 , ... , (length -1)
#' }
#'
#' @return Return a geometric or uniform grid of size \eqn{length} between \eqn{hmin} and \eqn{hmax}.
#' @export
#'
#' @examples
#' hmin <- 0.05
#' hmax <- 0.2
#' length <- 20
#' (h.geometric <- bandwidth.grid(hmin, hmax, length, type = "geometric"))
#' (h.uniform <- bandwidth.grid(hmin, hmax, length, type = "uniform"))
#'
bandwidth.grid <- function(hmin, hmax, length = 20, type = "geometric"){
  #generate a gri of bandwidth between hmin and hmax
  if(type == "geometric"){
    l <- 0:(length-1)
    q <- exp((log(hmax)-log(hmin))/(length-1))
    hl <- hmin*q^l
    return(hl)
  }
  if(type == "uniform"){
    hl <- seq(hmin, hmax, length = length)
    return(hl)
  }else{cat("please choose a type between geometric and uniform")}
}

# ########################################################################
# ## Simulations ###
# ########################################################################


#' The Pareto mix Distribution
#'
#' @name Pareto mix
#'
#' @description Density, distribution function, quantile function and random generation for the Pareto mixture distribution with \eqn{a} equal to the shape of the first Pareto Distribution, \eqn{b} equal to the shape of the second Pareto Distribution and \eqn{c} is the mixture proportion. The locations and the scales parameters are equals to \eqn{0} and \eqn{1}.
#'
#' @param x a vector of quantiles.
#' @param q a vector of quantiles.
#' @param p a vector of probabilities.
#' @param n the number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param a the shape parameter of the first Pareto Distribution.
#' @param b the shape parameter of the second Pareto Distribution.
#' @param c the value of the mixture proportion.
#' @param precision the precision of the Newton method.
#' @param initvalue the initial value of the Newton method.
#' @param Nmax the maximum of iteration done for the Newton method.
#'
#' @details If the \eqn{a}, \eqn{b} and \eqn{c} are not specified, they respectively take the default values \eqn{1}, \eqn{2} and \eqn{0.75}.
#'
#' The cumulative Pareto mixture distribution is
#' \deqn{
#'   F(x) = c  (1- x ^ {-a}) + ( 1 - c ) (1 - x ^ {-b}),  x \ge 1,  a >0,  b > 0, 0 \le c \le 1
#' }
#' where \eqn{a} and \eqn{b} are the shapes of the distribution and \eqn{c} is the mixture proportion.
#'
#' @return dparetomix gives the density, pparetomix gives the distribution function, qparetomix gives the quantile function, and rparetomix generates random deviates.
#'
#' The length of the result is determined by n for rparetomix, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.
#'
#' @export
#'
#' @examples
#' par(mfrow = c(3,1))
#' plot(function(x) dparetomix(x), 0, 5,ylab="density",
#'      main = " Pareto mixture density ")
#' mtext("dparetomix(x)", adj = 0)
#'
#' plot(function(x) pparetomix(x), 0, 5,ylab="distribution function",
#'      main = " Pareto mixture Cumulative ")
#' mtext("pparetomix(x)", adj = 0)
#'
#' plot(function(x) qparetomix(x), 0, 1,ylim=c(0,5),ylab="quantiles",
#'      main = " Pareto mixture Quantile ")
#' mtext("qparetomix(x)", adj = 0)
#'
#' #generate a sample of the Pareto mix distribution of size n
#' n <- 100
#' x <- rparetomix(n)
#'
#'
pparetomix <- function(q, a = 1, b = 2, c = 0.75){
  q <- pmax(q, 1)
  1-(c*q^(-a)+(1-c)*q^(-b))
}

#' @rdname Pareto-mix
dparetomix <- function(x, a = 1, b = 2, c = 0.75){
  #x <- pmax(x, 0)
  (a*c*x^(-a-1)+b*(1-c)*x^(-b-1))*(x>= 1)
}

#' @rdname Pareto-mix
qparetomix <- function(p, a = 1, b = 2, c = 0.75, precision = 10^(-10), initvalue = 0.5, Nmax = 1000){
  p <- 1-p
  t0 <- initvalue
  for(k in 1:Nmax){
    t1 <- t0 * (1+ (c*t0^(-a)+(1-c)*t0^(-b)-p) / (c*a*t0^(-a)+(1-c)*b*t0^(-b)) )
    if(abs(max(t1-t0))<precision) {break}
    t0 <- t1
    #cat("maximal number of iterations exeeded")
  }
  t1
}

#' @rdname Pareto-mix
rparetomix <- function(n, a = 1, b = 2, c = 0.75, precision = 10^(-10)){
  x <- runif(n)
  x <- pmin(1, pmax(0, x))
  qparetomix(x, a, b, c, precision)
}



#' Pareto Change Point distribution
#'
#' @description Distribution function, quantile function and random generation for the Pareto change point distribution with \eqn{a0} equal to the shape of the first pareto distribution, \eqn{a1} equal to the shape of the second pareto distribution, \eqn{x0} equal to the scale and \eqn{x1} equal to the change point.
#'
#' @param x a vector of quantiles.
#' @param p a vector of probabilities.
#' @param n a number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param a0 a vector of shape parameter of the Pareto distribution before \eqn{x1}.
#' @param a1 a vector of shape parameter of the Pareto distribution after \eqn{x1}.
#' @param x0 a vector of scale parameter of the function.
#' @param x1 a vector of change point value.
#'
#' @details If not specified, \eqn{a0, a1, x0} and \eqn{x1} are taking respectively the values \eqn{1, 2, 1} and \eqn{6}
#'
#' The cumulative Pareto change point distribution is given by :
#'   \deqn{
#'     F(x) = (x <= x1)* (1 - x^{-a0}) + (x > x1) * ( 1 - x^{-a1} * x1^{-a0 + a1})
#'   }
#'
#' @return pparetoCP gives the distribution function, qparetoCP gives the quantile function, and rparetoCP generates random deviates.
#'
#' The length of the result is determined by n for rparetoCP, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.
#'
#' @export
#'
#' @examples
#' par(mfrow = c(2,1))
#'
#' plot(function(x) pparetoCP(x), 0, 5,ylab="distribution function",
#'      main = " Pareto change point Cumulative ")
#' mtext("pparetoCP(x)", adj = 0)
#'
#' plot(function(x) qparetoCP(x), 0, 1,ylab="quantiles",
#'      main = " Pareto change point Quantile")
#' mtext("qparetoCP(x)", adj = 0)
#'
#' #generate a sample of pareto distribution of size n
#' n <- 100
#' x <- rparetoCP(n)
#'
pparetoCP <- function(x, a0 = 1, a1 = 2, x0 = 1, x1 = 6){
  x <- pmax(x, x0)
  x <- x/x0
  x1 <- x1/x0
  (x <= x1)*(1- x^(-a0))+ (x > x1)*(1- x^(-a1)*x1^(-a0+a1) )
}


#' @rdname pparetoCP
qparetoCP <- function(p, a0 = 1, a1 = 2, x0 = 1, x1 = 6){
  p <- pmin(1, pmax(0, p))
  y1  <-  1-(x1/x0)^(-a0)
  x0*((p <= y1)* ( (1-p)^(-1/a0) ) + (p > y1)* ((1-p)^(-1/a1) * (x1/x0)^(1-a0/a1)   ))
}

#' @rdname pparetoCP
rparetoCP <- function(n, a0 = 1, a1 = 2, x0 = 1, x1 = 6){
  x <- runif(n)
  x <- pmin(1, pmax(0, x))
  y1  <-  1-(x1/x0)^(-a0)
  x0*((x <= y1)*  ((1-x)^(-1/a0)) +  (x > y1)* ((1-x)^(-1/a1) * (x1/x0)^(1-a0/a1) )  )
}

#' The Pareto Distribution
#'
#' @name Pareto Distribution
#'
#' @description Density, distribution function, quantile function and random generation for the Pareto distribution where \eqn{a}, \eqn{loc} and \eqn{scale} are respectively the shape, the location and the scale parameters.
#'
#' @param x a vector of quantiles.
#' @param q a vector of quantiles.
#' @param p a vector of probabilities.
#' @param n a number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param a a vector of shape parameter of the Pareto distribution.
#' @param loc a vector of location parameter of the Pareto distribution.
#' @param scale a vector of scale parameter of the Pareto distribution.
#'
#' @details If \eqn{shape}, \eqn{loc} or \eqn{scale} parameters are not specified, the respective default values are \eqn{1}, \eqn{0} and \eqn{1}.
#'
#' The cumulative Pareto distribution is
#' \deqn{
#'   F(x) = 1- ((x-loc)/scale) ^ {-1/a},  x > 0,  a > 0, scale > 0
#' }
#' where \eqn{a} is the shape of the distribution.
#'
#' @return dpareto gives the density, ppareto gives the distribution function, qpareto gives the quantile function, and rpareto generates random deviates.
#'
#' The length of the result is determined by n for rpareto, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.
#' @export
#'
#' @examples
#' par(mfrow = c(3,1))
#' plot(function(x) dpareto(x), 1, 5,ylab="density",
#'      main = " Pareto density ")
#'
#' plot(function(x) ppareto(x), 1, 5,ylab="distribution function",
#'      main = " Pareto Cumulative ")
#'
#' plot(function(x) qpareto(x), 0, 1,ylab="quantile",
#'      main = " Pareto Quantile ")
#'
#' #generate a sample of pareto distribution of size n
#' n <- 100
#' x <- rpareto(n)
#'
ppareto  <-  function(q,  a = 1,  loc = 0,  scale = 1)
{
  #a = invpow = 1/pow

  1-((q-loc)/scale)^( - a)
}

#' @rdname Pareto-Distribution
dpareto  <-  function(x,  a = 1,  loc = 0,  scale = 1)
{
  #a = invpow = 1/pow

  (((x-loc)/scale)^( - a - 1)*a/scale)*(x-loc >= 1)
}

#' @rdname Pareto-Distribution
qpareto  <-  function(p,  a = 1,  loc = 0,  scale = 1)
{
  #a = invpow = 1/pow

  Q <- (1 - pmax(0, pmin(1, p)))^(-1/a)
  Q*scale+loc
}

#' @rdname Pareto-Distribution
rpareto  <-  function(n,  a = 1,  loc = 0,  scale = 1)
{
  #a = invpow = 1/pow

  Q <- (1 - runif(n,  0,  1))^(-1/a)
  Q*scale+loc
}

####################################################################
# Burr distribution

#' The Burr Distribution
#'
#' @name Burr Distribution
#'
#' @description Density, distribution function, quantile function and random generation for the Burr distribution with \eqn{a} and \eqn{k} two parameters.
#'
#' @param x a vector of quantiles.
#' @param q a vector of quantiles.
#' @param p a vector of probabilities.
#' @param n a number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param a a parameter of the burr distribution
#' @param k a parameter of the burr distribution
#'
#' @details The cumulative Burr distribution is
#' \deqn{
#'   F(x) = 1-( 1 + (x ^ a) ) ^{- k },  x >0,  a >0,  k > 0
#' }
#'
#'
#' @return dburr gives the density, pburr gives the distribution function, qburr gives the quantile function, and rburr generates random deviates.
#'
#' The length of the result is determined by n for rburr, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.
#' @export
#'
#' @examples
#' plot(function(x) dburr(x,3,1), 0, 5,ylab="density",
#' main = " burr density ")
#'
#' plot(function(x) pburr(x,3,1), 0, 5,ylab="distribution function",
#'      main = " burr Cumulative ")
#'
#' plot(function(x) qburr(x,3,1), 0, 1,ylab="quantile",
#'      main = " burr Quantile ")
#'
#' #generate a sample of burr distribution of size n
#' n <- 100
#' x <- rburr(n, 1, 1)
#'
#'
rburr <- function(n, a, k){
  ((1-runif(n))^(-1/k)-1)^(1/a)
}

#' @rdname Burr-Distribution
dburr <- function(x, a, k){
  k*a*x^(a-1)*(1+x^(a))^(-k-1)
}

#' @rdname Burr-Distribution
pburr <- function(q, a, k){
  1-(1+q^(a))^(-k)
}

#' @rdname Burr-Distribution
qburr <- function(p, a, k){
  ((1-p)^(-1/k)-1)^(1/a)
}


######################################################################
# dependent data of Burr distribution

#' Generate Burr dependent data
#'
#' @description Random generation function for the dependent Burr with a, b two shapes parameters and alpha the dependence parameter.
#'
#' @param n the number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param a a parameter of the function.
#' @param b a parameter of the function.
#' @param alpha the dependence parameter. It is defined by a single value between 0 and 1. The value 1 corresponds to the full independence. The closer to 0 the value of alpha is, the stronger is the dependence. \eqn{alpha} cannot take the value 0.
#'
#' @details The description of the dependence is described in \emph{Fawcett and Walshaw (2007)}. The Burr distribution is :
#' \eqn{
#'   F(x) = 1 - ( 1 + (x ^ a) ) ^ { - b }, x > 0, a > 0, b > 0
#' }
#' where a and b are shapes of the distribution.
#'
#' @return Generates a vector of random deviates. The length of the result is determined by n.
#'
#' @references Fawcett, D. and Walshaw, D. (2007). Improved estimation for temporally clustered extremes. Environmetrics, 18.2, 173-188.
#'
#' @export
#'
#' @examples
#' theta <- function(t){
#'    1/2*(1/10+sin(pi*t))*(11/10-1/2*exp(-64*(t-1/2)^2))
#'  }
#' n <- 200
#' t <- 1:n/n
#' Theta <- theta(t)
#' plot(theta)
#' alpha <- 0.6
#' Burr.dependent <- rburr.dependent(n, 1/Theta, 1, alpha)
#'
#'
rburr.dependent <- function(n, a, b, alpha){
  #generate n dependent burr distribution data with :
  # alpha is the dependence parameter: if  = 1 the data are independent and if close to 0,  the dependence increases
  # a and b are the parameter of the burr distribution
  rfrechet = function(m, a){
    U <- runif(m,  0, 1)
    (log(1/U))^(-a)
  }
  GyX <- function(x, y, alpha){
    Z <- x^(-1/alpha)+y^(-1/alpha)
    return( exp(1/y) * y^(-1/alpha+1) * Z^(alpha-1) * exp(-Z^alpha) )
  }
  derGyX <- function(x, y, alpha){
    Z <- x^(-1/alpha)+y^(-1/alpha)
    return( exp(1/y) * y^(-1/alpha+1) * exp(-Z^alpha) * x^(-1/alpha-1) * (Z^(2*alpha-2) - (((alpha-1)/alpha)*Z^(alpha-2))) )
  }

  y0 <- rfrechet(1, 1)
  if(n>1){
    y <- y0
    for(i in 2:n){
      p <- runif(1)
      x <- seq(0.1, 100, 0.01)
      #x0 <- 0.2
      #x0 <- x[which(abs(p-GyX(x, y, alpha)) == min(abs(p-GyX(x, y, alpha))))]
      x0 <- x[which(derGyX(x, y, alpha) == max(derGyX(x, y, alpha)))]
      m <- 1
      for (k in 1:1000){
        x1 = x0-(GyX(x0, y, alpha)-p)/derGyX(x0, y, alpha)
        #if(x1<0){x0 <- 0.2*m;m <- m+1;next}
        if (abs(x1-x0)< 10^-5)  {break}
        x0 <- x1
      }
      y <- x1
      y0[i] <- x1
    }
  }
  invY <- exp(-1/y0)
  if(length(a)>1){
    res <- NULL
    for(i in 1:n){res[i] <- qburr(invY[i], a[i], b)}
  }else{res <- qburr(invY, a, b)}
  return(res)
}
