# function for initialization of the range parameters in the discrepancy
RangeEstim = function(XF,yF,GP,u,s2)
{
  XFU <- cbind(XF, matrix(rep(u, nrow(XF)), ncol=length(u), byrow=TRUE)) 
  p <- predict(GP, XFU, xprime=XFU)
  C <- s2*diag(nrow(p$cov)) + (p$cov + t(p$cov))/2 + diag(p$nugs)
  res = yF - p$mean
  GPdisc = mleHomGP(X=XF,Z=res,covtype = "Gaussian",known = list(g=s2,beta0=0))
  return(list(psi=GPdisc$theta,sigb=GPdisc$nu_hat))
}


# fonction for calibration
likcalibrationwithdisc = function(XF,yF,GP,theta,Sigdisc=NULL)
{
  
  
  if (is.null(Sigdisc)) bRange = 1
  else bRange = 0
  
  ## input processing and checking
  if(length(theta) != (ncol(GP$X0) - ncol(XF))*(bRange+1) + 2) 
    stop("length(theta), ncol(XF), ncol(GP$X0) mismatch")
  u <- theta[1:ncol(XF)]
  s2f <- theta[ncol(XF)+1]
  s2b = theta[ncol(XF)+2]
  if (bRange==1) thdisc = theta[(ncol(XF)+3):length(theta)]
  
  
  ## prior checking  
  if(any(u < 0 | u > 1)) return (-Inf)
  if(s2f < 0) return(-Inf)
  if(s2b<0) return(-Inf)
  if (bRange==1) {if(any(thdisc<0)) return(-Inf)}
  
  ## derive predictive distribution for XF paired with u
  XFU <- cbind(XF, matrix(rep(u, nrow(XF)), ncol=length(u), byrow=TRUE)) 
  p <- predict(GP, XFU, xprime=XFU)
  
  if (bRange==0) {Cdisc = Sigdisc * s2b} else {
    Cdisc = s2b * cov_gen(XF,theta=thdisc,type="Gaussian")
  }
  
  C <- s2f*diag(nrow(p$cov)) + (p$cov + t(p$cov))/2 + diag(p$nugs) +Cdisc  #discrepancy
  
  
  ## gaussian log density evaluation for yF under that predictive
  return(dmvnorm(yF, p$mean, C, log=TRUE) - log(s2f)-log(s2b))
}


# Sum of squares for L2 optimization
SumOfSquares = function(u,GP,XF,yF)
{
  XFU <- cbind(XF, matrix(rep(u, nrow(XF)), ncol=length(u), byrow=TRUE)) 
  p <- predict(GP, XFU)
  sum((yF - p$mean)^2)
}

# Functions for prediction from a fixed value for parameters
prednoncalfixed = function(k,GP,u,vareps,loc) #ATTENTION vareps ajoute par rapport au RMD
{
  p = predict(GP, matrix(c(loc,rep(u,each=nloc)),nrow=nloc))
  Zprednoncal = (p$mean + rnorm(nloc,0,sqrt(p$sd2+p$nugs)))*sqrt(Yv) + Ym + rnorm(nloc,0,sqrt(vareps)) # ajouter la variance d'obs ??
  return(Zprednoncal)
}


# Function for prediction when theta is sampled uniformly in the domain (uniform prior distribution) 
prednoncal = function(k,GP,vareps,loc)
{
  nloc = nrow(loc)
  p = predict(GP, matrix(c(loc,rep(runif(2),each=nloc)),nrow=nloc))
  Zprednoncal = (p$mean + rnorm(nloc,0,sqrt(p$sd2+p$nugs)))*sqrt(Yv) + Ym + rnorm(nloc,0,sqrt(vareps)) # ajouter la variance d'obs ??
  return(Zprednoncal)
}

# Function for prediction incorporating discrepancy and posterior sample of theta
predcal = function(k,cal,GP,vareps,loc,YfN,Xfield)
{
  nloc = nrow(loc)
  testfield = rbind(matrix(c(loc,rep(cal[k,1:2],each=nloc)),nrow=nloc),
                    matrix(c(Xfield[,1:2],rep(cal[k,1:2],each=nrow(Xfield))),nrow=nrow(Xfield)))
  pcal <- predict(GP, testfield ,xprime=testfield)
  CGP = (pcal$cov + t(pcal$cov))/2 + diag(pcal$nugs)
  realCM = as.vector(pcal$mean + rmvnorm(1,rep(0,length(pcal$mean)),CGP))
  diff = YfN - realCM[-(1:nloc)]
  Cdisc = cal[k,4] * cov_gen(Xfield[,1:2],theta=cal[k,5:6],type="Gaussian") + cal[k,3]* diag(nrow(Xfield))
  Cdiscloc = cal[k,4] * cov_gen(X1=Xfield[,1:2],X2=loc,theta=cal[k,5:6],type="Gaussian")
  sigdiscloc = cal[k,4] * cov_gen(loc,theta=cal[k,5:6],type="Gaussian") - t(Cdiscloc) %*% solve(Cdisc,Cdiscloc)
  realdisc = as.vector(t(Cdiscloc) %*% solve(Cdisc,diff)) + as.vector(rmvnorm(1,rep(0,nloc),sigdiscloc))
  return( (realCM[1:nloc] + realdisc)*sqrt(Yv) + Ym + rnorm(nloc,0,sqrt(Yv*cal[k,3])))
}

# Computing score from an empirical distribution
scoreEstDens = function(Ech,v)
{
  vscore=numeric(ncol(Ech))
  for (k in 1:ncol(Ech))
  {
    estd = density(Ech[,k])
    vscore[k] = log(approx(estd$x,estd$y,xout=v[k])$y)
  }
  return(mean(vscore))
}