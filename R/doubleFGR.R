
##' @examples
##' library(mets)
##' data(bmt)
##' bmt$time <- bmt$time+runif(408)*0.01
##' ocif1 <- cifreg(Event(time,cause)~platelet+age,bmt,cause=1,propodds=NULL)
##' ocif2 <- cifreg(Event(time,cause)~platelet+age,bmt,cause=2,propodds=NULL)
##' ###
##' ocr<- phreg(Surv(time,cause==0)~platelet+age,bmt)
##' ## change coefficients to be same as ocif1
##' ###ocr$coef <- ocif1$coef[-2]
##' ocr$coef <- c(-1,1)*ocr$coef
##' n <- 500000
##'
##' dats <- sim.cifs(list(ocif1,ocif2),n,bmt,cens=NULL)
##' dsummary(dats,time+age+platelet~status)
##' dats <- sim.cifsRestrict(list(ocif1,ocif2),n,bmt)
##' dsummary(dats,time+age+platelet~status)
##'
##' mm <- model.matrix(~platelet+age,dats)[,-1]
##' datc <- rchaz(rbind(c(0,0),ocr$cumhaz),exp(mm %*% c(2*coef(ocr))))
##' dats$time <- pmin(dats$time,datc$time)
##' dats$status <- ifelse(dats$time<datc$time,dats$status,0)
##' dtable(dats,~status+cause)
##'
##' bplot(ocif1)
##' bplot(ocif2,add=TRUE,col=1)
##' ### on simulated data FG models
##' cif1 <- cifreg(Event(time,status)~platelet+age,data=dats,cause=1,propodds=NULL)
##' cif2 <- cifreg(Event(time,status)~platelet+age,data=dats,cause=2,propodds=NULL)
##' bplot(cif1,add=TRUE,col=2)
##' bplot(cif2,add=TRUE,col=2)
##'
##' cbind(coef(ocif1),coef(ocif2))
##' cbind(coef(cif1),coef(cif2))
##'
##' ## loop for restrict model (mix with recursive, algorithm)
##' pr2 <- doubleFGR(Event(time,status)~platelet+age,data=dats,restrict=3)
##' cbind(pr2$coef[1:2],pr2$coef[3:4])
##' matlines(pr2$cumhaz[,1],pr2$cumhaz[,2:3],type="s",lty=2,col=3)
##' ###
##'
##' pr0 <- doubleFGR(Event(time,status)~platelet+age,data=dats,restrict=0)
##' cbind(pr0$coef[1:2],pr0$coef[3:4])
##' matlines(pr0$cumhaz[,1],pr0$cumhaz[,2:3],type="s",lty=2,col=2)
##' ###


## * doubleFGR
doubleFGR <- function(formula,data,offset=NULL,weights=NULL,...) {# {{{
    cl <- match.call()
    m <- match.call(expand.dots = TRUE)[1:3]
    special <- c("strata", "cluster","offset")
    Terms <- terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Y <- model.extract(m, "response")
### if (!is.Surv(Y)) stop("Expected a 'Surv'-object")
    if (ncol(Y)==2) {
        exit <- Y[,1]
        entry <- NULL ## rep(0,nrow(Y))
        status <- Y[,2]
    } else {
        entry <- Y[,1]
        exit <- Y[,2]
        status <- Y[,3]
    }
    id <- strata <- NULL
    if (!is.null(attributes(Terms)$specials$cluster)) {
        ts <- survival::untangle.specials(Terms, "cluster")
        pos.cluster <- ts$terms
        Terms <- Terms[-ts$terms]
        id <- m[[ts$vars]]
    } else pos.cluster <- NULL
    if (!is.null(attributes(Terms)$specials$strata)) {
        ts <- survival::untangle.specials(Terms, "strata")
        pos.strata <- ts$terms
        Terms <- Terms[-ts$terms]
        strata <- m[[ts$vars]]
        strata.name <- ts$vars
    } else { strata.name <- NULL; pos.strata <- NULL}
### if (!is.null(attributes(Terms)$specials$offset)) {
### ts <- survival::untangle.specials(Terms, "offset")
### pos.offset <- ts$terms
### Terms <- Terms[-ts$terms]
### offset <- m[[ts$vars]]
### } else pos.offset <- NULL
    X <- model.matrix(Terms, m)
    if (!is.null(intpos <- attributes(Terms)$intercept))
        X <- X[,-intpos,drop=FALSE]
    if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
    res <- c(doubleFG01R(X,entry,exit,status,id,strata,offset,weights,strata.name,...),
             list(call=cl,model.frame=m,formula=formula,strata.pos=pos.strata,cluster.pos=pos.cluster))
    class(res) <- "doubleFG"

res
}

## * doubleFG01R
doubleFG01R <- function(X,entry,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,
                        strata.name=NULL,cumhaz=TRUE,
                        beta,stderr=TRUE,method="NR",no.opt=FALSE,Z=NULL,propodds=NULL,AddGam=NULL,
                        restrict=0,case.weights=NULL,...) {
    p <- ncol(X)
    if (missing(beta)) beta <- rep(0,2*p)
    if (p==0) X <- cbind(rep(0,length(exit)))
    if (is.null(strata)) { strata <- rep(0,length(exit)); nstrata <- 1; strata.level <- NULL; } else {
                                                                                                  strata.level <- levels(strata)
                                                                                                  ustrata <- sort(unique(strata))
                                                                                                  nstrata <- length(ustrata)
                                                                                                  strata.values <- ustrata
                                                                                                  if (is.numeric(strata)) strata <- fast.approx(ustrata,strata)-1 else {
                                                                                                                                                                      strata <- as.integer(factor(strata,labels=seq(nstrata)))-1
                                                                                                                                                                  }
                                                                                              }
    if (is.null(offset)) offset <- rep(0,length(exit))
    if (is.null(weights)) weights <- rep(1,length(exit))
    strata.call <- strata
    Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given
    if (!is.null(Z)) Zcall <- Z

## possible casewights to use for bootstrapping and other things
if (is.null(case.weights)) case.weights <- rep(1,length(exit))

trunc <- (!is.null(entry))
if (!trunc) entry <- rep(0,length(exit))

if (!is.null(id)) {
ids <- unique(id)
nid <- length(ids)
if (is.numeric(id)) id <- fast.approx(ids,id)-1 else {
id <- as.integer(factor(id,labels=seq(nid)))-1
}
} else id <- as.integer(seq_along(entry))-1;
## orginal id coding into integers
id.orig <- id+1;

dd <- .Call("FastCoxPrepStrata",entry,exit,status,X,id,trunc,strata,weights,offset,Zcall,case.weights,PACKAGE="mets")

cause <- status[dd$ord+1]
dd$cause <- cause
dd$nstrata <- nstrata

obj <- function(pp,U=FALSE,all=FALSE) {# {{{
val <- with(dd, doubleFGstrataR(pp,X,XX,sign,cause,jumps,strata,nstrata,weights,offset,ZX,caseweights,restrict))

if (all) {
val$time <- dd$time
val$cause <- dd$cause
val$ord <- dd$ord+1
val$jumps <- dd$jumps+1
val$jumptimes <- val$time[val$jumps]
val$weightsJ <- dd$weights[val$jumps]
val$case.weights <- dd$case.weights[val$jumps]
val$strata.jumps <- val$strata[val$jumps]
val$nevent <- length(val$S0)
val$nstrata <- dd$nstrata
val$strata <- dd$strata
return(val)
}
with(val,structure(-ploglik,gradient=-gradient,hessian=-hessian))
}# }}}

opt <- NULL
if (p>0) {
if (no.opt==FALSE) {
if (tolower(method)=="nr") {
tim <- system.time(opt <- lava::NR(beta,obj,...))
opt$timing <- tim
opt$estimate <- opt$par
} else {
opt <- nlm(obj,beta,...)
opt$method <- "nlm"
}
cc <- opt$estimate; names(cc) <- rep(colnames(X),2)
if (!stderr) return(cc)
val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
} else val <- c(list(coef=beta),obj(beta,all=TRUE))
} else {
val <- obj(0,all=TRUE)
}

se.cumhaz <- lcumhaz <- lse.cumhaz <- NULL
II <- NULL
### computes Breslow estimator
if (cumhaz==TRUE) { # {{{
if (no.opt==FALSE & p!=0) {
II <- - tryCatch(solve(val$hessian),error=
function(e) matrix(0,nrow(val$hessian),ncol(val$hessian)) )
} else II <- matrix(0,p,p)
strata <- val$strata[val$jumps]
nstrata <- val$nstrata
jumptimes <- val$jumptimes

## Brewslow estimator
### cumhaz <- cbind(jumptimes,cumsumstrata(1/val$S0,strata,nstrata))
cumhaz <- cbind(jumptimes,val$base12$base)

### if ((no.opt==FALSE & p!=0)) {
### DLambeta.t <- apply(val$E/c(val$S0),2,cumsumstrata,strata,nstrata)
### varbetat <- rowSums((DLambeta.t %*% II)*DLambeta.t)
### ### covv <- apply(covv*DLambeta.t,1,sum) Covariance is "0" by construction
### } else varbetat <- 0
### var.cumhaz <- cumsumstrata(1/val$S0^2,strata,nstrata)+varbetat
### se.cumhaz <- cbind(jumptimes,(var.cumhaz)^.5)
se.cumhaz <- cumhaz

### colnames(cumhaz1) <- c("time","cumhaz")
### colnames(cumhaz2) <- c("time","cumhaz")
### colnames(se.cumhaz) <- c("time","se.cumhaz")
} # }}}
else {cumhaz <- se.cumhaz <- lcumhaz <- lse.cumhaz <- NULL}

res <- c(val,
list(cox.prep=dd,
strata.call=strata.call, strata.level=strata.level,
entry=entry,
exit=exit,
status=status,
p=p,
X=X,
offsets=offset,
weights=weights,
id=id.orig,
opt=opt,
cumhaz=cumhaz, se.cumhaz=se.cumhaz,
lcumhaz=lcumhaz, lse.cumhaz=lse.cumhaz,
II=II,strata.name=strata.name,propodds=propodds))
class(res) <- "phreg"
res
}

## * doubleFGstrataR
doubleFGstrataR <- function(beta,
                            X,
                            XX,
                            Sign,
                            cause,
                            Jumps,
                            strata,
                            nstrata,
                            weights,
                            offsets,
                            ZX,
                            caseweights,
                            restrict){# {{{
    p=length(beta)/2
    strata=c(strata)
### X1 <- X[,1:p,drop=FALSE]
### X2 <- X[,(p+1):(2*p),drop=FALSE]
    Xb1 = c(X %*% beta[1:p]+offsets)
    Xb2 = c(X %*% beta[(p+1):(2*p)]+offsets)
    eXb1 = c(exp(Xb1)*weights);
    eXb2 = c(exp(Xb2)*weights);
    if (nrow((Sign))==length(eXb1)) { ## Truncation
        eXb1 = c(Sign)*eXb1;
        eXb2 = c(Sign)*eXb2;
    }

    S01 = c(revcumsumstrata(eXb1,strata,nstrata))
    S02 = c(revcumsumstrata(eXb2,strata,nstrata))
    E1=apply(eXb1*as.matrix(X),2,revcumsumstrata,strata,nstrata)/S01;
    E2=apply(eXb2*as.matrix(X),2,revcumsumstrata,strata,nstrata)/S02;
    Jumps=Jumps+1
    causeJ <- cause[Jumps]
    Jumps1 <- Jumps[causeJ==1]
    Jumps2 <- Jumps[causeJ==2]

    ## both E's, SO's, X same
    S0J <- cbind(S01,S02)[Jumps,]
    eXbJ <- cbind(eXb1,eXb2)[Jumps,]
    EJ <- cbind(E1,E2)[Jumps,]
    XJ <- X[Jumps,]

    E1 = E1[Jumps1,,drop=FALSE];
    E2 = E2[Jumps2,,drop=FALSE];
    E21=.Call("vecMatMat",E1,E1)$vXZ;
    E22=.Call("vecMatMat",E2,E2)$vXZ;

    XX21=apply(XX*eXb1,2,revcumsumstrata,strata,nstrata)/S01;
    XX22=apply(XX*eXb2,2,revcumsumstrata,strata,nstrata)/S02;
    ## mat XX2=revcumsumstrataMatCols(XX,eXb,S0,strata,nstrata);
    XX21 = XX21[Jumps1,,drop=FALSE];
    XX22 = XX22[Jumps2,,drop=FALSE];

    weightsJ=weights[Jumps];
    caseweightsJ=caseweights[Jumps];

    ## compute recursive weights given beta1, beta2
    base12 <- .Call("cumsumstrataDFGR",weightsJ,S0J,causeJ,strata[Jumps],nstrata,eXbJ)
    if (restrict>0) {
        for (i in 1:restrict) {
            Lam1inf <- tail(base12$base[,1],1);
            base12 <- .Call("cumsumstrataDFGRestrictR",weightsJ,S0J,causeJ,strata[Jumps],nstrata,eXbJ,Lam1inf);
        }
    }

### DLam <- .Call("DLambetaDFGR",weightsJ,S0J,causeJ,EJ,XJ,strata[Jumps],nstrata,eXbJ);

    grad1 = (X[Jumps1,,drop=FALSE]-E1); ## Score
    grad2 = (X[Jumps2,,drop=FALSE]-E2); ## Score

### weights for Fg1 an Fg2 models
    pow1 <- base12$pow1[causeJ==1]
    pow2 <- base12$pow2[causeJ==2]
    cw1 <- (caseweightsJ*weightsJ)[causeJ==1]
    cw2 <- (caseweightsJ*weightsJ)[causeJ==2]
    S012 = S01[Jumps1]/(pow1*cw1); ## S0 with weights to estimate baseline
    S022 = S02[Jumps2]/(pow2*cw2); ## S0 with weights to estimate baseline

    val = sum(pow1*cw1*(Xb1[Jumps1]-log(S012)))+
        sum(pow2*cw2*(Xb2[Jumps2]-log(S022))); ## Partial log-likelihood

    grad21= grad1*(pow1*cw1); ## score with weights
    grad22= grad2*(pow2*cw2); ## score with weights

    grad <- c(apply(grad21,2,sum),apply(grad22,2,sum))
    gradient <- grad

    ## no weights
### val2 = caseweightsJ*weightsJ*val; ## Partial log-likelihood with weights
    val2 = val; ## Partial log-likelihood with weights

    hesst1 = -(XX21-E21); ## hessian contributions in jump times
    hesst2 = -(XX22-E22); ## hessian contributions in jump times
### hess = matrix(apply(hesst,2,sum),p,p);
    hesst12 = hesst1*(cw1)*pow1; ## hessian over time with weights
    hesst22 = hesst2*(cw2)*pow2; ## hessian over time with weights

    ## missing some derivative terms for hessian (due to pow=w(beta,\Lam(beta,t-))

    ## setup hessian matrix
    hess12 = matrix(apply(hesst12,2,sum),p,p); ## hessian with weights
    hess22 = matrix(apply(hesst22,2,sum),p,p); ## hessian with weights
    hess2 <- matrix(0,2*p,2*p)
    pd2 <- p
    hess2[1:pd2,1:pd2] <- hess12
    hess2[(pd2+1):(2*p),(pd2+1):(2*p)] <- hess22

    out=list(jumps=Jumps, ploglik=sum(val2),U=grad2,base12=base12,
             gradient=matrix(gradient,1,2*p), hessian=hess2, ##hessianttime=hesst2,
### S2S0=XX2,
             E=EJ, S0=S0J
             )
    return(out)
}# }}}


