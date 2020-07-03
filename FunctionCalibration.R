# function for initialization of the range parameters in the discrepancy
# rep to account for possible replications, in this case XF has to be unique location and yF must me averaged,
# we assume that the replications are organized by concatenation of batches of data corresponding to unique location
RangeEstim = function(XF,yF,GP,u,s2,rep=1)
{
  XF = unique.matrix(XF)
  yFmat = matrix(yF,nrow=nrow(XF),ncol=rep,byrow =FALSE)
  yF = rowMeans(yFmat)
  XFU <- cbind(XF, matrix(rep(u, nrow(XF)), ncol=length(u), byrow=TRUE)) 
  p <- predict(GP, XFU, xprime=XFU)
  C <- s2*diag(nrow(p$cov)) + (p$cov + t(p$cov))/2 #+ diag(p$nugs)
  res = yF - p$mean
  GPdisc = mleHomGP(X=XF,Z=res,covtype = "Gaussian",known = list(g=s2/rep,beta0=0))
  Sigdisc = cov_gen(XF,theta=u,type="Gaussian")
  return(list(psi=GPdisc$theta,sigb=GPdisc$nu_hat,Sigdisc=Sigdisc))
}


# generate realization of a GP
simdiscrepancy = function(X,s2d=.02,psi = c(.1,.2))
{
  SIG = s2d * cov_gen(X,theta=psi,type="Gaussian")
  as.vector(rmvnorm(1,rep(0,nrow(X)),SIG))
}

# function for calibration
#priorUpBounds2b for uniform prior on [0,priorUpBounds2b] for s2b variance discrepancy parameter 
postcalibrationwithdisc = function(theta,XF,yF,GP,Sigdisc=NULL,priorUpBounds2b=NULL,logvar=FALSE)
{
  
  # bound on variances 
  s2upbound = .1
  
  if (is.null(Sigdisc)) bRange = 1
  else bRange = 0
  
  ## input processing and checking
  if(length(theta) != (ncol(GP$X0) - ncol(XF))*(bRange+1) + 2) 
    stop("length(theta), ncol(XF), ncol(GP$X0) mismatch")
  u <- theta[1:ncol(XF)]
  
 
  
  if (logvar==TRUE)  
  {
    s2f <- exp(theta[ncol(XF)+1])
    s2b = exp(theta[ncol(XF)+2])
  }
    else {
      s2f <- theta[ncol(XF)+1]
      s2b = theta[ncol(XF)+2]
    }
    
  if (bRange==1) thdisc = theta[(ncol(XF)+3):length(theta)]
  

  
  
  ## prior checking  
  if(any(u < 0 | u > 1)) return (-Inf)
  if(s2f < 0) return(-Inf)
  if(s2b<0) return(-Inf)
  if (bRange==1) {if(any(thdisc<0)) return(-Inf)}
  if (!is.null(priorUpBounds2b)) {
    upbound = min(priorUpBounds2b, s2upbound )
    if (s2b>upbound) return(-Inf)}
  if (s2f>s2upbound) return(-Inf)
    
  ## derive predictive distribution for XF paired with u
  XFU <- cbind(XF, matrix(rep(u, nrow(XF)), ncol=length(u), byrow=TRUE)) 
  p <- predict(GP, XFU, xprime=XFU)
  
  if (bRange==0) {Cdisc = Sigdisc * s2b} else {
    Cdisc = s2b * cov_gen(XF,theta=thdisc,type="Gaussian")
  }
  
  C <- s2f*diag(nrow(p$cov)) + (p$cov + t(p$cov))/2  +Cdisc  #discrepancy  + diag(p$nugs) # for variance of sto sim
  
  
 
  
  
  ## gaussian log density evaluation for yF under that predictive
  return(dmvnorm(yF, p$mean, C, log=TRUE) )
}




# function for calibration
#priorUpBounds2b for uniform prior on [0,priorUpBounds2b] for s2b variance discrepancy parameter 
# replicates
postcalibrationwithdiscrep = function(theta,XF,yF,GP,Sigdisc=NULL,priorUpBounds2b=NULL,logvar=FALSE,rep=1)
{
  
  XF = unique.matrix(XF)
  yFmat = matrix(yF,nrow=nrow(XF),ncol=rep,byrow=FALSE)
  yF = rowMeans(yFmat)
  
  # bound on variances 
  s2upbound = .1
  
  if (is.null(Sigdisc)) bRange = 1
  else bRange = 0
  
  ## input processing and checking
  if(length(theta) != (ncol(GP$X0) - ncol(XF))*(bRange+1) + 2) 
    stop("length(theta), ncol(XF), ncol(GP$X0) mismatch")
  u <- theta[1:ncol(XF)]
  
  
  
  if (logvar==TRUE)  
  {
    s2f <- exp(theta[ncol(XF)+1])
    s2b = exp(theta[ncol(XF)+2])
  }
  else {
    s2f <- theta[ncol(XF)+1]
    s2b = theta[ncol(XF)+2]
  }
  
  if (bRange==1) thdisc = theta[(ncol(XF)+3):length(theta)]
  
  
  
  
  ## prior checking  
  if(any(u < 0 | u > 1)) return (-Inf)
  if(s2f < 0) return(-Inf)
  if(s2b<0) return(-Inf)
  if (bRange==1) {if(any(thdisc<0)) return(-Inf)}
  if (!is.null(priorUpBounds2b)) {
    upbound = min(priorUpBounds2b, s2upbound )
    if (s2b>upbound) return(-Inf)}
  if (s2f>s2upbound) return(-Inf)
  
  ## derive predictive distribution for XF paired with u
  XFU <- cbind(XF, matrix(rep(u, nrow(XF)), ncol=length(u), byrow=TRUE)) 
  p <- predict(GP, XFU, xprime=XFU)
  
  if (bRange==0) {Cdisc = Sigdisc * s2b} else {
    Cdisc = s2b * cov_gen(XF,theta=thdisc,type="Gaussian")
  }
  
  C <- s2f/rep*diag(nrow(p$cov)) + (p$cov + t(p$cov))/2  +Cdisc  #discrepancy  + diag(p$nugs) # for variance of sto sim
  
  
  lik2sf = 0
  
  if (rep>1)
  {
    varobs = sum(apply(yFmat,1,var))
    liks2f = dchisq(varobs/s2f,(rep-1)*nrow(XF),log=TRUE) + log(s2f) 
  }
  
  
  
  ## gaussian log density evaluation for yF under that predictive
  return(dmvnorm(yF, p$mean, C, log=TRUE) + lik2sf)
}



postcalibrationwithoutdisc = function(theta,XF,yF,GP,logvar=FALSE)
{
  
  # bound on variances 
  s2upbound = .1
  
  
  u <- theta[1:ncol(XF)]

  if (logvar==TRUE)  
  {
    s2f <- exp(theta[ncol(XF)+1])
  }
  else {
    s2f <- theta[ncol(XF)+1]

  }
  
  
  ## prior checking  
  if(any(u < 0 | u > 1)) return (-Inf)
  if(s2f < 0) return(-Inf)
  if (s2f>s2upbound) return(-Inf)
  
  ## derive predictive distribution for XF paired with u
  XFU <- cbind(XF, matrix(rep(u, nrow(XF)), ncol=length(u), byrow=TRUE)) 
  p <- predict(GP, XFU, xprime=XFU)
  
 
  C <- s2f*diag(nrow(p$cov)) + (p$cov + t(p$cov))/2 #+ diag(p$nugs) for variance of stosim
  
  
  ## gaussian log density evaluation for yF under that predictive
  return(dmvnorm(yF, p$mean, C, log=TRUE)  )
}




# Sum of squares for L2 optimization
SumOfSquares = function(u,GP,XF,yF)
{
  XFU <- cbind(XF, matrix(rep(u, nrow(XF)), ncol=length(u), byrow=TRUE)) 
  p <- predict(GP, XFU)
  sum((yF - p$mean)^2)
}

# Functions for prediction from a fixed value for parameters
prednoncalfixed = function(k,GP,u,vareps,loc,Ym,Yv)
{
  nloc = nrow(loc)
  p = predict(GP, matrix(c(loc,rep(u,each=nloc)),nrow=nloc))
  Zprednoncal = (p$mean + rnorm(nloc,0,sqrt(p$sd2)))*sqrt(Yv) + Ym + rnorm(nloc,0,sqrt(vareps))
    #(p$mean + rnorm(nloc,0,sqrt(p$sd2+p$nugs)))*sqrt(Yv) + Ym + rnorm(nloc,0,sqrt(vareps))
  return(Zprednoncal)
}


# Function for prediction when theta is sampled uniformly in the domain (uniform prior distribution) 
prednoncal = function(k,GP,vareps,loc,Ym,Yv)
{
  nloc = nrow(loc)
  p = predict(GP, matrix(c(loc,rep(runif(2),each=nloc)),nrow=nloc))
  Zprednoncal = (p$mean + rnorm(nloc,0,sqrt(p$sd2)))*sqrt(Yv) + Ym + rnorm(nloc,0,sqrt(vareps))
    #(p$mean + rnorm(nloc,0,sqrt(p$sd2+p$nugs)))*sqrt(Yv) + Ym + rnorm(nloc,0,sqrt(vareps))
  return(Zprednoncal)
}




# Function for prediction incorporating discrepancy and posterior sample of theta
predcal = function(k,cal,GP,vareps,loc,YfN,Xfield,Ym,Yv,psi=NULL)
{
  # log scale
  cal[k,3] = exp(cal[k,3])
  cal[k,4] = exp(cal[k,4])
  if (is.null(psi)) {psi=cal[k,5:6]} # should be given in the post sample cal if not provided
  
  # prevent too small value of variance for numerical issues
  if (cal[k,3]<1e-8) {cal[k,3] = 1e-8}
  if (cal[k,4]<1e-8) {cal[k,4] = 1e-8}
  
  nloc = nrow(loc)
  u = c(cal[k,1],cal[k,2])
  testfield = rbind(cbind(loc,matrix(u,nrow(loc),2,byrow = T)),
                    cbind(Xfield[,1:2],matrix(u,nrow(Xfield),2,byrow=T)))
  pcal <- predict(GP, testfield ,xprime=testfield)
  #CGP = (pcal$cov + t(pcal$cov))/2 + diag(pcal$nugs)
  CGP = (pcal$cov + t(pcal$cov))/2 #+ diag(pcal$nugs)
  realCM = as.vector(pcal$mean + rmvnorm(1,rep(0,length(pcal$mean)),CGP))
  diff = YfN - realCM[-(1:nloc)]
  Cdisc = cal[k,4] * cov_gen(Xfield[,1:2],theta=psi,type="Gaussian") + cal[k,3]* diag(nrow(Xfield))
  Cdiscloc = cal[k,4] * cov_gen(X1=Xfield[,1:2],X2=loc,theta=psi,type="Gaussian")
  sigdiscloc =  cal[k,3]* diag(nrow(loc)) + cal[k,4] * cov_gen(loc,theta=psi,type="Gaussian") - t(Cdiscloc) %*% solve(Cdisc,Cdiscloc)
  
  # avoid non positive definite matrix
  sigdiscloceig = eigen(sigdiscloc)$values
  if (any(is.complex(sigdiscloceig))) {simdiscloc = 0}
  else {if (any(eigen(sigdiscloc)$values<=0)) simdiscloc = 0
  else simdiscloc = as.vector(rmvnorm(1,rep(0,nloc),sigdiscloc))}
  realdisc = as.vector(t(Cdiscloc) %*% solve(Cdisc,diff)) + simdiscloc
  return( (realCM[1:nloc] + realdisc)*sqrt(Yv) + Ym )
}





# Computing score from an empirical distribution
scoreEstDens = function(Ech,v)
{
  vscore=numeric(ncol(Ech))
  for (k in 1:ncol(Ech))
  {
    estd = density(Ech[,k],from=min(v),to=max(v))
    vscore[k] = log(max(1e-20,approx(estd$x,estd$y,xout=v[k])$y))
  }
  return(mean(vscore))
}