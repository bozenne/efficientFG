### FCT-robustE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 26 2020 (16:57) 
## Version: 
## Last-Updated: mar 11 2020 (11:37) 
##           By: Brice Ozenne
##     Update #: 103
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * examples
##' @examples
##' ## load mets functions
##' library(mets)
##' butils.base:::sourcePackage("efficientFG", c.code = TRUE, trace = TRUE)
##'
##' ## prepare data
##' data(bmt)
##' bmt$id <- 1:nrow(bmt)
##' bmt$age2 <- bmt$age
##' bmt$platelet2 <- as.character(bmt$platelet)
##'
##' ## initial estimator of the cif
##' e.cif <- doubleFGR(Event(time,cause)~age+platelet,data=bmt)
##' 
##' ## initial estimator of the censoring 
##' e.censoring <- phreg(Surv(time,cause==0)~age+platelet,data=bmt)
##'
##' ## construct double robust estimator
##' e.eFG <- efficientFineGray(event = e.cif, censoring = e.censoring)
##' 


## * efficientFineGray
efficientFineGray <- function(event, censoring, type.Delta = "Delta"){

    type.Delta <- match.arg(type.Delta, c("Delta","Deltat"))
    class(event) <- append("doubleFGR",class(event))
    
    ## ** extract information from input
    ## *** order of the observations
    order.12 <- event$ord[,1]
    order.0 <- censoring$ord[,1]
    
    ## *** event times
    time <- event$time[,1] ## ordered according to event time, and then cause (first 1,2 and then 0)
    time <- time[order(order.12)] ## restaure original ordering
    time.ties <- unique(time[duplicated(time)])
    Utime <- unique(time)
    n.Utime <- length(Utime)

    ## *** type of event
    varepsilon <- event$cause ## ordered according to event time, and then cause (first 1,2 and then 0)
    varepsilon <- varepsilon[order(order.12)] ## restaure original ordering
    if(any(varepsilon %in% 0:2 == FALSE)){
        stop("Event type indicator must take its values in 0,1, or 2. \n")
    }
    jumptime.12 <- event$jumptimes
    jumptime.0 <- censoring$jumptimes

    if(!is.null(event$strata.entry) || !is.null(censoring$strata.entry)){
        stop("Cannot handle delayed entry \n")
    }

    ## *** strata
    if(event$nstrata>1 || censoring$nstrata>1){
        stop("Can only handle one strata \n")
        ## event$strata[,1]
        ## censoring$strata[,1]
    }
        
    ## *** regression parameters
    p.0 <- censoring$p
    indexp.0 <- 1:p.0
    p.1 <- event$p
    indexp.1 <- 1:p.1
    p.2 <- event$p
    indexp.2 <- p.1 + (1:p.2)

    beta.0 <- coef(censoring)[indexp.0]
    beta.1 <- event$coef[indexp.1]
    beta.2 <- event$coef[indexp.2]

    ## *** design matrix
    X.1 <- event$X ## original ordering
    X.2 <- event$X ## original ordering
    X.0 <- censoring$X  ## original ordering

    n.obs <- NROW(X.1)
    if(NROW(X.0)!=n.obs){
        stop("Unequal number of rows in the design matrix between the event model and the censoring model \n")
    }
    ## ## *** baseline hazards
    Lambda0.0 <- rep(NA, n.Utime)
    Lambda0.1 <- rep(NA, n.Utime)
    Lambda0.2 <- rep(NA, n.Utime)
    Lambda0.0[match(censoring$cumhaz[,1],Utime)] <- censoring$cumhaz[,2]
    Lambda0.1[match(event$cumhaz[,1],Utime)] <- event$cumhaz[,2]
    Lambda0.2[match(event$cumhaz[,1],Utime)] <- event$cumhaz[,3]

    
    Lambda0.0 <- NA2previousValue(Lambda0.0, start = 0)
    Lambda0.1 <- NA2previousValue(Lambda0.1, start = 0)
    Lambda0.2 <- NA2previousValue(Lambda0.2, start = 0)
    
    ## ** prepare relevant quantities

    ## *** dN: counting process
    tdN.0 <- matrix(0, nrow = n.Utime, ncol = n.obs, 
                    dimnames = list(Utime, NULL))
    tdN.1 <- matrix(0, nrow = n.Utime, ncol = n.obs, 
                    dimnames = list(Utime, NULL))
    tdN.2 <- matrix(0, nrow = n.Utime, ncol = n.obs, 
                    dimnames = list(Utime, NULL))

    indexJump.0 <- which(varepsilon == 0)
    indexJump.1 <- which(varepsilon == 1)
    indexJump.2 <- which(varepsilon == 2)
    Utime.0 <- unique(time[indexJump.0])
    n.Utime.0 <- length(Utime.0)
    Utime.1 <- unique(time[indexJump.1])
    n.Utime.1 <- length(Utime.1)
    Utime.2 <- unique(time[indexJump.2])
    n.Utime.2 <- length(Utime.2)
    
    match_timeUtime <- match(time,Utime)

    tdN.0[(indexJump.0-1) * n.Utime + match_timeUtime[indexJump.0]] <- 1
    tdN.1[(indexJump.1-1) * n.Utime + match_timeUtime[indexJump.1]] <- 1
    tdN.2[(indexJump.2-1) * n.Utime + match_timeUtime[indexJump.2]] <- 1

    dN.0 <- base::t(tdN.0)
    dN.1 <- base::t(tdN.1)
    dN.2 <- base::t(tdN.2)

    ## *** Y: at risk indicator
    ## table(rowSums(dN.1))
    Y.0 <- (1-riskRegression::rowCumSum(dN.0)) + dN.0
    Y.1 <- (1-riskRegression::rowCumSum(dN.1)) + dN.1
    Y.2 <- (1-riskRegression::rowCumSum(dN.2)) + dN.2

    ## *** Delta: censoring indicator
    time.censoring <- time
    time.censoring[varepsilon!=0] <- Inf

    Delta <- varepsilon!=0
    if(type.Delta=="Deltat"){
        Deltat <- do.call(cbind,lapply(Utime, function(t){
            return((Delta + (t < time.censoring))>0)
        }))
        colnames(Deltat) <- Utime
    }else if(type.Delta == "Delta"){
        Deltat <- matrix(Delta, nrow = n.obs, ncol = n.Utime, byrow = FALSE,
                         dimnames = list(NULL, Utime))
    }
    
    ## *** eXB: linear predictor
    eXb.0 <- exp(X.0 %*% beta.0)
    eXb.1 <- exp(X.1 %*% beta.1)
    eXb.2 <- exp(X.2 %*% beta.2)
    
    Lambda.0 <- tcrossprod(eXb.0,Lambda0.0)
    Lambda.1 <- tcrossprod(eXb.1,Lambda0.1)
    Lambda.2 <- tcrossprod(eXb.2,Lambda0.2)

    ## *** G: survival relative to the censoring mechanism
    G.0 <- exp(-Lambda.0)
    ## range(G.0 - predict(censoring, time = Utime, se = FALSE)$surv[order.12,])
    G.0atT <- G.0[which(dN.0+dN.1+dN.2==1)]
    ## matrix(1:n.obs,nrow=n.obs,ncol=n.Utime)[which(dN.0+dN.1+dN.2==1)]
    ## fields::image.plot(G.0)
    Gt.0 <- matrix(G.0atT, nrow = n.obs, ncol = n.Utime, byrow = FALSE,
                   dimnames = list(NULL, Utime))
    ## fields::image.plot(Y.0+Y.1+Y.2)
    if(type.Delta == "Deltat"){
        for(iObs in 1:n.obs){ ## iObs <- 390
            Gt.0[iObs, 1:match_timeUtime[iObs]] <- G.0[iObs, 1:match_timeUtime[iObs]]
            ## Gt.0[iObs, match_timeUtime[iObs]] - G.0[iObs, match_timeUtime[iObs]]
        }
    }else if(type.Delta == "Delta"){
        ## nothing to be done
    }

    ## *** F1/F2: cumulative incidence function relative to each event
    F1 <- 1 - exp(-Lambda.1)
    F2 <- 1 - exp(-Lambda.2)
    ## xx <- predict(event, newdata = bmt, times = Utime)
    ## range(xx[[1]] - F1)
    ## range(xx[[2]] - F2)
    ## range(F1 + F2)

    ## *** dMc: martingale relative to the censoring mechanism
    lambda0.0 <- diff(c(0,Lambda0.0))
    lambda.0 <- tcrossprod(eXb.0, lambda0.0)
    ## range(cumsum(lambda0.0)-Lambda0.0)
    dM.0 <- dN.0 - Y.0 * lambda0.0
    
    ## ** double robust estimate of the e
    drS_0 <- .S0DR(Y = Y.1, eXb = eXb.1, F1 = F1, F2 = F2, Delta = Deltat, dMc = dM.0, Gc = G.0, Gc_T = Gt.0)
    drS_1 <- .S1DR(X = X.1, Y = Y.1, eXb = eXb.1, F1 = F1, F2 = F2, Delta = Deltat, dMc = dM.0, Gc = G.0, Gc_T = Gt.0)
    
    ## ** double robust estimate of the baseline hazard
    drLamba0.1 <- .LambdaDR(Y = Y.1, eXb = eXb.1, dN = dN.1, S0 = drS_0, Delta = Deltat, Gc_T = Gt.0)

    browser()
    ## ** double robust estimate of the \beta
}

##
.S0DR <- function(Y, eXb, F1, F2, Delta, dMc, Gc, Gc_T){

    int_dMcGc <- matrix(rowSums(dMc/Gc), nrow = NROW(F1), ncol = NCOL(F1), byrow = FALSE)
    int_F1frac2dMcGc <- F1*riskRegression::rowCumSum(dMc/(Gc*(1-F1-F2)))
    int_frac1dMcGc <- riskRegression::rowCumSum(F1*dMc/(Gc*(1-F1-F2)))

    S0iDR <- riskRegression::colMultiply_cpp(Y*Delta/Gc_T + int_dMcGc - int_F1frac2dMcGc + int_frac1dMcGc, eXb)
    ## S0iDR <- riskRegression::colMultiply_cpp(Y*Delta/Gc_T, eXb)
    S0DR <- colSums(S0iDR)
    return(S0DR)
}

.S1DR <- function(X, Y, eXb, dMc, F1, F2, Delta, Gc, Gc_T){

    S1DR <- apply(X, 2, function(iX){ ## iX <- X[,1]
        .S0DR(Y = Y, eXb = eXb*iX, dMc = dMc, F1 = F1, F2 = F2, Delta = Delta, Gc = Gc, Gc_T = Gc_T)
    })
    return(S1DR)
}

.LambdaDR <- function(Y, eXb, dN, S0, Delta, Gc_T){
    browser()    
    psi <- riskRegression::rowCumSum(dN / S0)
    I <- riskRegression::colMultiply_cpp(psi + I, eXb)
}

.drBeta <- function(){
}


## * rev_rowCumSum [not used]
##' @examples
##' rowCumSum(matrix(1,3,3))
##' rev_rowCumSum(matrix(1,3,3))
rev_rowCumSum <- function(x){
    riskRegression::rowCumSum(x[,NCOL(x):1,drop=FALSE])[,NCOL(x):1,drop=FALSE]
}


######################################################################
### FCT-robustE.R ends here
