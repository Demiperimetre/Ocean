rm(list=ls())
setwd("/home/pbarbillon/save/Ocean/")

## Load useful packages
library(ggplot2)
library(colorRamps)
library(gridExtra)
#library(DiceKriging)
library(DiceDesign)
library(hetGP)
library(lhs)
library(mvtnorm)
library(plyr)
library(reshape2)
library(parallel)

## Load the ocean simulator
source("SimulatorAndFunctions.R")


## number of emulation
Nrep = 100


## set seeds
seed=1234
set.seed(seed)



## GP settings
lower <- rep(0.01, 2)
upper <- rep(10, 2)
covtype <- "Matern5_2"
nc <- list(g_min=1e-6, g_bounds=c(1e-6, 1), lowerDelta=log(1e-6))
settings <- list(linkThetas="none", logN=TRUE, initStrategy="smoothed",
                 checkHom=TRUE, penalty=TRUE, trace=0, return.matrices=TRUE, return.hom=FALSE)
control <- list(tol_dist=1e-4, tol_diff=1e-4, multi.start=30)


## Load testdesign
test = read.csv("testdata2D.csv",sep=" ")
testdesign = as.matrix(test[,1:2])
Ztest.mean = test[,3]
Ztest.sd = test[,4]

## Do everything in parallel
RESrep = mclapply(1:Nrep,
function(k){
# space filling design
n=50#n = 80 #
X0 = randomLHS(n,2)
X0 = maximinSA_LHS(X0)$design
X <- rbind(X0, X0, X0, X0, X0, X0, X0, X0, X0, X0)
X <- rbind(X, X)
Z <- apply(X, 1, simulator)

# normalization
Zm <- mean(Z)
Zv <- var(Z)
Z <- (Z - Zm)/sqrt(Zv)


# GP fit
## Homoskedastic
Ghom <- mleHomGP(X, Z, lower=lower, upper=upper, covtype=covtype,
                known=list(beta0=0), maxit=1000)

## Heteroskedastic
Ghet <- mleHetGP(X, Z, lower=lower, upper=upper, covtype=covtype,
                noiseControl=nc, settings=settings, known=list(beta0=0), maxit=1000)



## -seqdesign
ninitseq = 5
n=nrow(X0)
Xseq <- X[1:(ninitseq*n),]
Y <- Z[1:(ninitseq*n)]
mod <- mleHetGP(Xseq, Y, lower=lower, upper=upper, covtype=covtype,
                noiseControl=nc, settings=settings, known=list(beta0=0), maxit=1000)

nadd = nrow(X)-nrow(Xseq)

h <- rep(NA, nadd)
## acquisitions
for(i in 1:nadd) {

  ## choose lookahead horizon and solve IMSPE
  h[i] <- horizon(mod)
  opt <- IMSPE_optim(mod, h[i], control=control)
  #cat("i=", i, ", h=", h[i], "\n", sep="")

  ## evaluate the simulator
  ynew <- (simulator(opt$par) -Zm)/sqrt(Zv)

  ## update the fit
  mod <- update(mod, Xnew=opt$par, Znew=ynew, ginit=mod$g*1.01)
  if(i %% 25 == 0){
    mod2 <- mleHetGP(list(X0=mod$X0, Z0=mod$Z0, mult=mod$mult),
                     Z=mod$Z, lower=lower, upper=upper, covtype=covtype,
                     noiseControl=nc, settings=settings, known=list(beta0=0),
                     maxit=1000)
    if(mod2$ll > mod$ll) mod <- mod2
  }
}
seqGhet=mod


## Predictions

#Hom GP
predhom = predict(x = testdesign, object = Ghom)
#Het GP
predhet = predict(x = testdesign, object = Ghet)
#Het GP seq
predHetseq = predict(x = testdesign, object = seqGhet)

# Normalize Ztest.mean
Ztest.mean.N = (Ztest.mean - Zm)/sqrt(Zv)
Ztest.sd.N = Ztest.sd/sqrt(Zv)

# MSE for mean
rmsehom = sqrt(mean(((predhom$mean-Ztest.mean.N))^2)*Zv)
rmsehet = sqrt(mean(((predhet$mean-Ztest.mean.N))^2)*Zv)
rmsehetseq = sqrt(mean(((predHetseq$mean-Ztest.mean.N))^2)*Zv)


#MSE for sd
rmsehomsd = sqrt(mean(((sqrt(predhom$nugs)-Ztest.sd.N))^2)*Zv)
rmsehetsd = sqrt(mean(((sqrt(predhet$nugs)-Ztest.sd.N))^2)*Zv)
rmsehetseqsd = sqrt(mean(((sqrt(predHetseq$nugs)-Ztest.sd.N))^2)*Zv)


# scores for the prediction of a single run of the simulator
Ztest = (apply(testdesign,1,simulator)-Zm)/sqrt(Zv) # for computing scores on a single realization of the simulator
schom = mean(-(Ztest-predhom$mean)^2/(predhom$sd2+predhom$nugs) -log(predhom$sd2+predhom$nugs) - log(Zv))
schet = mean(-(Ztest-predhet$mean)^2/(predhet$sd2+predhet$nugs) -log(predhet$sd2+predhet$nugs)- log(Zv))
schetseq = mean(-(Ztest-predHetseq$mean)^2/(predHetseq$sd2+predHetseq$nugs) -log(predHetseq$sd2+predHetseq$nugs)- log(Zv))

## export results
return(c(rmsehom=rmsehom,rmsehet=rmsehet,rmsehetseq=rmsehetseq,rmsehomsd=rmsehomsd,rmsehetsd=rmsehetsd,rmsehetseqsd=rmsehetseqsd,scorehom=schom,scorehet=schet,scorehetseq=schetseq))
},mc.cores=10)


RES = Reduce(rbind,RESrep)
colnames(RES) = c("rmsehom","rmsehet","rmsehetseq","rmsehomsd","rmsehetsd","rmsehetseqsd","scorehom","scorehet","scorehetseq")
rownames(RES) = 1:100

name="EmululationRep.csv"
write.table(RES,file=name)

