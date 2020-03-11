### predict.doubleFGR.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 11 2020 (09:49) 
## Version: 
## Last-Updated: mar 11 2020 (10:54) 
##           By: Brice Ozenne
##     Update #: 47
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * predict.doubleFG (examples)
##' @examples
##' ## load mets functions
##' library(mets)
##' library(data.table)
##' library(ggplot2)
##' butils.base:::sourcePackage("efficientFG", c.code = TRUE, trace = TRUE)
##'
##' ## prepare data
##' data(bmt)
##' bmt$id <- 1:nrow(bmt)
##' bmt$age2 <- bmt$age
##' bmt$platelet2 <- bmt$platelet
##'
##' bmt.sorted <- bmt[order(bmt$time,bmt$cause),]
##' ## initial estimator of the cif
##' e.cif <- doubleFGR(Event(time,cause)~age+platelet,data=bmt)
##' ls.pred <- predict(e.cif, time = c(1:100), newdata  = bmt[1:2,])
##' dt.pred <- predict(e.cif, time = c(1:100), newdata  = bmt[1:2,,drop=FALSE], format = "data.table")
##' range(ls.pred[[1]] + ls.pred[[2]])
##'
##' gg <- ggplot(dt.pred, aes(x = time, y = cif, group = cause, color = cause))
##' gg <- gg + geom_point() + geom_line()
##' gg <- gg + facet_wrap(~id)
##' gg

## * predict.doubleFG
predict.doubleFG <- function(object, newdata, times, format = "matrix"){

    format <- match.arg(format, c("matrix","data.table"))
    
    ## ** event times
    obstime <- object$time[,1]
    Uobstime <- unique(obstime)
    n.Uobstime <- length(Uobstime)

    ## ** type of event
    varepsilon <- object$cause ## ordered according to event time (i.e. first 1,2 and then 0)
    Uvarepsilon <- sort(unique(varepsilon))
    
    ## ** regression parameters
    p.1 <- object$p
    indexp.1 <- 1:p.1
    p.2 <- object$p
    indexp.2 <- p.1 + (1:p.2)

    beta.1 <- object$coef[indexp.1]
    beta.2 <- object$coef[indexp.2]

    ## ** design matrix
    X <- model.matrix(object$formula, newdata, intercept = FALSE)
    ## head(X)

    ## ** baseline hazards
    Lambda0.1 <- rep(NA, n.Uobstime)
    Lambda0.2 <- rep(NA, n.Uobstime)
    Lambda0.1[match(object$cumhaz[,1],Uobstime)] <- object$cumhaz[,2]
    Lambda0.2[match(object$cumhaz[,1],Uobstime)] <- object$cumhaz[,3]

    ## for timepoint where there is no jump, replace NA by the previous cumhazard value 
    Lambda0.1 <- NA2previousValue(Lambda0.1, start = 0)
    Lambda0.2 <- NA2previousValue(Lambda0.2, start = 0)

    ## ** eXb: linear predictor
    eXb.1 <- exp(X[,names(beta.1),drop=FALSE] %*% beta.1)
    eXb.2 <- exp(X[,names(beta.2),drop=FALSE] %*% beta.2)
    
    indexJump <- prodlim::sindex(jump.times = Uobstime, eval.times = times)
    Lambda.1 <- eXb.1 %*% rbind(c(0,Lambda0.1)[indexJump+1])
    Lambda.2 <- eXb.2 %*% rbind(c(0,Lambda0.2)[indexJump+1])
    if(any(times > max(obstime))){
        Lambda.1[,times > max(obstime)] <- NA
        Lambda.2[,times > max(obstime)] <- NA
    }
    
    ## ** cumulative incidence function    
    out <- setNames(vector(mode = "list", length = 3), c(Uvarepsilon[Uvarepsilon!=0], "time"))
    out[[1]] <- 1 - exp(-Lambda.1)
    out[[2]] <- 1 - exp(-Lambda.2)
    out[[3]] <- times
    if(format == "data.table"){
        ## cause 1
        dt1 <- data.table::melt(data.table::data.table(id = 1:NROW(newdata), cause = names(out)[1], out[[1]]), id.vars = c("id","cause"),
                                variable.name = "time", value.name = "cif")
        dt1[, c("time") := times[as.numeric(as.factor(.SD$time))]]

        ##  cause 2
        dt2 <- data.table::melt(data.table::data.table(id = 1:NROW(newdata), cause = names(out)[2], out[[2]]), id.vars = c("id","cause"),
                                variable.name = "time", value.name = "cif")
        dt2[, c("time") := times[as.numeric(as.factor(.SD$time))]]

        ## all causes
        dt3 <- data.table::copy(dt2)
        dt3$cif <- dt1$cif + dt2$cif
        dt3$cause <- paste(dt1$cause,"+",dt2$cause)

        ## assemble
        out <- rbind(dt1,dt2,dt3)
    }
    
    return(out)
}

## * NA2previousValue
##' @examples
##' NA2previousValue(c(1,NA,NA,2,NA,3,NA))
NA2previousValue <- function(x, start){
    indexNA <- which(is.na(x))
    if(length(indexNA)>0){
        for(iX in indexNA){
            if(iX==1){
                x[iX] <- start
            }else{
                x[iX] <- x[iX-1]
            }
        }
    }
    return(x)
}

######################################################################
### predict.doubleFGR.R ends here
