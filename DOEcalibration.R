rm(list=ls())

## Load useful packages
library(ggplot2)
library(colorRamps)
library(gridExtra)
library(DiceKriging)
library(DiceDesign)
library(hetGP)
library(MCMCpack)
library(lhs)
library(mvtnorm)




## Load the ocean simulator and the functions for DOE
source("SimulatorAndFunctions.R")
source("LHDcalibration.R")


## set seeds
seed=1234
set.seed(seed)

## Generate "field data" 
n = 150
X0 = lhs::randomLHS(n,2)
X0 = maximinSA_LHS(X0)$design
u1 = (700 -100)/900  
u2 = (200- 100)/900
Xfield = cbind(X0,u1,u2)
# one replicate hence two field data at the same location
Xfield = rbind(Xfield,Xfield)
### run the simulator
Yfield = sapply(1:n,function(i) simulator4d(Xfield[i,],NPATHS=1200))
### add a Gaussian withe noise
vareps = 2^2 
Yfieldnoise = Yfield + rnorm(length(Yfield),0,sd=sqrt(vareps)) 
### normalize the response data
Ym <- mean(Yfieldnoise)
Yv <- var(Yfieldnoise)
YfieldnoiseN <- (Yfieldnoise-Ym)/sqrt(Yv)


### Saving field data
write.table(cbind(Xfield,Yfieldnoise),file="FieldData.csv",col.names = FALSE,row.names = FALSE)


## Compute emulators for calibration purpose

### DOE with field input variables restricted to the field location data (in X0) and maximin LHD constraint for the theta (calibration parameters)
b_inf_theta=c(0,0)
b_sup_theta=c(1,1)
size_design = 500
out = GRID_design(Xfield[,1:2],b_inf_theta,b_sup_theta,size_design)
out_opt=maximinSA_LHS_grid(out,dim_xf=2,c=0.95,dim_theta=2,Imax=2000)

### 10 replicates for each of the 500 unique locations
nreplicStatic = 10
Xsim = matrix(rep(t(out_opt),nreplicStatic),ncol=4,byrow = T)
Ysim = apply(Xsim,1,simulator4d)
YsimN = (Ysim-Ym)/sqrt(Yv)

write.table(cbind(Xsim,Ysim),file="staticDOEcalibration.csv",col.names = FALSE,row.names = FALSE)




### HomGP and HetGP emulators
covtype <- "Matern5_2"
noiseControl <- list(g_min=1e-6, g_bounds=c(1e-6, 1), lowerDelta=log(1e-6))
settings <- list(linkThetas="none", initStrategy="smoothed", return.hom=TRUE)
lower <- c(0.01, 0.01, 0.001, 0.001) 
upper <- c(30, 30, 100, 100) 

het <- mleHetGP(Xsim, YsimN, lower=lower, upper=upper, covtype=covtype, noiseControl=noiseControl, 
                settings=settings, maxit=10000)
hom <- het$modHom


### Seq Het GP emulator with initial design from initial design with only 4 replicates
nreplicInit = 4
Xsim.init = Xsim[1:(nreplicInit*size_design),]
Ysim.init = YsimN[1:(nreplicInit*size_design)]

covtype <- "Matern5_2"
noiseControl <- list(g_min=1e-6, g_bounds=c(1e-6, 1), lowerDelta=log(1e-6))
settings <- list(linkThetas="none", initStrategy="smoothed", return.hom=FALSE)
lower <- c(0.01, 0.01, 0.001, 0.001) 
upper <- c(30, 30, 100, 100) 

mod <- mleHetGP(Xsim.init, Ysim.init, lower=lower, upper=upper, covtype=covtype, noiseControl=noiseControl, 
                settings=settings, maxit=10000)

Xseq = Xsim.init
Yseq = Ysim.init
nadd = nrow(Xsim)-nrow(Xseq)  
control <- list(tol_dist=1e-4, tol_diff=1e-4, multi.start=30)

h <- rep(NA, nadd)
## acquisitions
for(i in 1:nadd) { 
  
  ## choose lookahead horizon and solve IMSPE
  h[i] <- horizon(mod)
  opt <- IMSPE_optim(mod, h[i], control=control)
  cat("i=", i, ", h=", h[i], "\n", sep="")
  
  ## evaluate the simulator
  ynew <- (simulator4d(opt$par)-Ym) / sqrt(Yv)
  
  ## update the fit
  mod <- update(mod, Xnew=opt$par, Znew=ynew, ginit=mod$g*1.01)
  if(i %% 25 == 0){
    mod2 <- mleHetGP(list(X0=mod$X0, Z0=mod$Z0, mult=mod$mult),
                     Z=mod$Z, lower=lower, upper=upper, covtype=covtype, 
                     noiseControl=noiseControl, settings=settings, known=list(beta0=0), 
                     maxit=1000)
    if(mod2$ll > mod$ll) mod <- mod2  
  }
}
seqGhet = mod

### saving 

save.image(file="EmulatorCalib4Dall.Rdata")
save(hom,het,seqGhet,file="EmulatorCalib4D.Rdata")

Xseq = seqGhet$X0[rep(1:length(seqGhet$mult),times=seqGhet$mult),]
Yseq = seqGhet$Z*sqrt(Yv) + Ym

write.table(cbind(Xseq,Yseq),file="seqDOEcalibration.csv",col.names = FALSE,row.names = FALSE)

