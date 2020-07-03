Calibration of the Ocean Simulator
================
Pierre Barbillon
July 1, 2020

``` r
rm(list=ls())
library(ggplot2)
library(colorRamps)
library(gridExtra)
library(DiceKriging)
library(DiceDesign)
library(hetGP)
library(MCMCpack)
library(lhs)
library(mvtnorm)
library(parallel)
library(reshape2)
```

Loading simulator and GPEs (homGP, hetGP and seqhetGP which is hetGP
with a sequential design)

``` r
# functions
source("SimulatorAndFunctions.R")
source("FunctionCalibration.R")
```

Fixing seed

``` r
seed=1234
set.seed(seed)
```

## Retrieving field data

``` r
Field = read.table("FieldData.csv",sep=" ")
Xfield = as.matrix(Field[,1:4])
Yfieldnoise = Field[,5]
n = nrow(Xfield)
Ym <- mean(Yfieldnoise)
Yv <- var(Yfieldnoise)
YfieldnoiseN <- (Yfieldnoise-Ym)/sqrt(Yv)
vareps = 2^2 
```

Simulate discrepancy field

``` r
test = read.csv("testdata2D.csv",sep=" ")
testdesign = as.matrix(test[,1:2])
# simulate discrepancy
Xfielduniloc = unique.matrix(Xfield)
# simulate discrepancy and add it to field data
#disc = simdiscrepancy(rbind(Xfielduniloc[,1:2],testdesign))
#YfieldnoiseN <- (Yfieldnoise-Ym)/sqrt(Yv) + rep(disc[1:nrow(Xfielduniloc)],2)
# or reload it
disc = read.table("discrepancy.csv")[,1]
YfieldnoiseN = read.table("noisyfielddata.csv")[,1]
```

## Emulators

Either

``` r
#load("EmulatorCalib4D.Rdata")
```

Or

``` r
## From Static design
Static = read.table("staticDOEcalibration.csv")
Xsim = as.matrix(Static[,1:4])
Ysim = Static[,5]
YsimN = (Ysim-Ym)/sqrt(Yv)  
  
covtype <- "Matern5_2"
noiseControl <- list(g_min=1e-6, g_bounds=c(1e-6, 1), lowerDelta=log(1e-6))
settings <- list(linkThetas="none", initStrategy="smoothed", return.hom=TRUE)
lower <- c(0.01, 0.01, 0.001, 0.001) 
upper <- c(30, 30, 100, 100) 

het <- mleHetGP(Xsim, YsimN, lower=lower, upper=upper, covtype=covtype, noiseControl=noiseControl, 
                settings=settings, maxit=10000)
hom <- het$modHom

## Sequential design
Seq = read.table("seqDOEcalibration.csv")
Xseq = as.matrix(Seq[,1:4])
Yseq = Seq[,5]
YseqN = (Yseq-Ym)/sqrt(Yv)  
  
covtype <- "Matern5_2"
noiseControl <- list(g_min=1e-6, g_bounds=c(1e-6, 1), lowerDelta=log(1e-6))
settings <- list(linkThetas="none", initStrategy="smoothed", return.hom=TRUE)
lower <- c(0.01, 0.01, 0.001, 0.001) 
upper <- c(30, 30, 100, 100) 

seqhet <- mleHetGP(Xseq, YseqN, lower=lower, upper=upper, covtype=covtype, noiseControl=noiseControl, 
                settings=settings, maxit=10000)
```

## Calibration

Settings for MCMC

``` r
nMCMC = 1e5  #1e6
nburnin =  2e4 
# variance for random walks
Vrw <- diag(c(.05,.05,.1,.1))/2
Uprior =c(.5,.5)
Thin = 100
```

Calibration with homGP

``` r
homMLEdisc = RangeEstim(Xfield[,1:2],YfieldnoiseN,hom,Uprior,vareps,rep=2)
# Run MH algo for calibration is time consuming
# calHom <- MCMCmetrop1R(postcalibrationwithdiscrep, XF=Xfield[,1:2], yF=YfieldnoiseN, GP=hom,
#                        Sigdisc=homMLEdisc$Sigdisc,priorUpBounds2b=3/Yv,rep=2,logvar=T
#                              ,logfun=TRUE, theta.init=c(0.5,0.5,log(2/Yv),log(2/Yv)),
#                              burnin=nburnin, mcmc=nMCMC,V=Vrw,thin=Thin)
# write.table(calHom,file="calHom.csv")


# reloading saved calibration results
calHom = read.table("calHom.csv",h=T,sep=" ")
calHom = as.matrix(calHom)
dfcalhom = as.data.frame(calHom)
names(dfcalhom) = c("u1","u2","s2f","s2b")
dfcalhom[,1:2] = dfcalhom[,1:2] * 900 +100
dfcalhom[,3] = exp(dfcalhom[,3]) * Yv
dfcalhom[,4] = exp(dfcalhom[,4]) * Yv
hom1 = ggplot(dfcalhom,aes(x=u1,stat(density)))+ xlab(expression(K[x])) + geom_density() + theme_bw() + xlim(50,1000) + geom_vline(aes(xintercept=700),color="red")
hom2 = ggplot(dfcalhom,aes(x=u2,stat(density))) + xlab(expression(K[y])) + geom_density() + theme_bw() +xlim(50,1000)  + geom_vline(aes(xintercept=200),color="red")
hom3 = ggplot(dfcalhom,aes(x=s2f,stat(density)))+ xlab(expression(sigma[epsilon]^2))  + geom_density() + theme_bw() + xlim(0,7) + geom_vline(aes(xintercept=4),color="red")
hom4 = ggplot(dfcalhom,aes(x=s2b,stat(density)))+ xlab(expression(sigma[MD]^2)) + geom_density() + theme_bw() + xlim(0,3) +geom_vline(aes(xintercept=Yv*.02),color="red")
```

Calibration with hetGP

``` r
hetMLEdisc = RangeEstim(Xfield[,1:2],YfieldnoiseN,het,Uprior,vareps,rep=2)
# Run MH algo for calibration is time consuming
# calHet <- MCMCmetrop1R(postcalibrationwithdiscrep, XF=Xfield[,1:2], yF=YfieldnoiseN, GP=het,
#                              Sigdisc=hetMLEdisc$Sigdisc,priorUpBounds2b=3/Yv,rep=2,logvar=T
#                              ,logfun=TRUE, theta.init=c(0.5,0.5,log(2/Yv),log(2/Yv)),
#                              burnin=nburnin, mcmc=nMCMC,V=Vrw,thin=Thin)
# write.table(calHet,file="calHet.csv")

# reloading saved calibration results
calHet = read.table("calHet.csv",h=T,sep=" ")
calHet = as.matrix(calHet)
dfcalhet = as.data.frame(calHet)
names(dfcalhet) = c("u1","u2","s2f","s2b")
dfcalhet[,1:2] = dfcalhet[,1:2] * 900 +100
dfcalhet[,3] = exp(dfcalhet[,3]) * Yv
dfcalhet[,4] = exp(dfcalhet[,4]) * Yv
het1 = ggplot(dfcalhet,aes(x=u1,stat(density))) + xlab(expression(K[x]))+ geom_density() + theme_bw() + xlim(50,1000) + geom_vline(aes(xintercept=700),color="red")
het2 = ggplot(dfcalhet,aes(x=u2,stat(density)))+ xlab(expression(K[y])) + geom_density() + theme_bw() +xlim(50,1000)  + geom_vline(aes(xintercept=200),color="red")
het3 = ggplot(dfcalhet,aes(x=s2f,stat(density)))+ xlab(expression(sigma[epsilon]^2)) + geom_density() + theme_bw() + xlim(0,7) + geom_vline(aes(xintercept=4),color="red")
het4 = ggplot(dfcalhet,aes(x=s2b,stat(density))) + xlab(expression(sigma[MD]^2)) + geom_density() + theme_bw() + xlim(0,3) +geom_vline(aes(xintercept=Yv*.02),color="red")
```

Calibration with seqhetGP

``` r
seqhetMLEdisc = RangeEstim(Xfield[,1:2],YfieldnoiseN,seqhet,Uprior,vareps,rep=2)
# Run MH algo for calibration is time consuming
# calSeqHet <- MCMCmetrop1R(postcalibrationwithdiscrep, XF=Xfield[,1:2], yF=YfieldnoiseN, GP=seqhet,
#                                 Sigdisc=seqhetMLEdisc$Sigdisc,priorUpBounds2b=3/Yv,rep=2,logvar=T
#                              ,logfun=TRUE, theta.init=c(0.5,0.5,log(2/Yv),log(2/Yv)),
#                              burnin=nburnin, mcmc=nMCMC,V=Vrw,thin=Thin)
#write.table(calSeqHet,file="calSeqHet.csv")

# reloading saved calibration results
calSeqHet = read.table("calSeqHet.csv",h=T,sep=" ")
calSeqHet = as.matrix(calSeqHet)
dfcalseqhet = as.data.frame(calSeqHet)
names(dfcalseqhet) = c("u1","u2","s2f","s2b")
dfcalseqhet[,1:2] = dfcalseqhet[,1:2] * 900 +100
dfcalseqhet[,3] = exp(dfcalseqhet[,3]) * Yv
dfcalseqhet[,4] = exp(dfcalseqhet[,4]) * Yv
seqhet1 = ggplot(dfcalseqhet,aes(x=u1,stat(density))) + xlab(expression(K[x]))+ geom_density() + theme_bw() + xlim(50,1000) + geom_vline(aes(xintercept=700),color="red")
seqhet2 = ggplot(dfcalseqhet,aes(x=u2,stat(density)))+ xlab(expression(K[y])) + geom_density() + theme_bw() +xlim(50,1000)  + geom_vline(aes(xintercept=200),color="red")
seqhet3 = ggplot(dfcalseqhet,aes(x=s2f,stat(density))) + xlab(expression(sigma[epsilon]^2))+ geom_density() + theme_bw() + xlim(0,7) + geom_vline(aes(xintercept=4),color="red")
seqhet4 = ggplot(dfcalseqhet,aes(x=s2b,stat(density)))+ xlab(expression(sigma[MD]^2))  + geom_density() + theme_bw() + xlim(0,3) +geom_vline(aes(xintercept=Yv*.02),color="red")
```

``` r
grid.arrange(arrangeGrob(hom1,hom2,hom3,hom4,top="homGP",ncol=4),arrangeGrob(het1,het2,het3,het4,top="hetGP",ncol=4),
             arrangeGrob(seqhet1,seqhet2,seqhet3,seqhet4,top="seq hetGP",ncol=4),nrow=3,ncol=1)
```

    ## Warning: Removed 16 rows containing non-finite values (stat_density).

    ## Warning: Removed 44 rows containing non-finite values (stat_density).

    ## Warning: Removed 27 rows containing non-finite values (stat_density).

![](OceanCalibration_files/figure-gfm/plotcalib-1.png)<!-- -->

## Prediction

Loading test data

``` r
# test = read.csv("testdata2D.csv",sep=" ")
# testdesign = as.matrix(test[,1:2])
Ztest.mean = test[,3] + sqrt(Yv)*disc[-c(1:nrow(Xfielduniloc))]
Ztest.sim = Ztest.mean + rnorm(length(Ztest.mean),0,sqrt(vareps))
```

We run prediction from a uniform prior on the calibration parameter and
prediction from the posterior distribution of the calibration parameter
and the discrepancy.

Note that the seed is not fixed anymore since we use the parallel
computations.

``` r
npred = nMCMC/Thin
ZpredhomGPnoncal = mclapply(1:npred,prednoncal,GP=hom,vareps=vareps,loc=testdesign,Ym=Ym,Yv=Yv,mc.cores = 8)
ZpredhomGPcal = mclapply(1:npred,predcal,GP=hom,cal=calHom,vareps=vareps,
                         loc=testdesign,YfN=YfieldnoiseN,Xfield=Xfield,Ym=Ym,Yv=Yv,psi=homMLEdisc$psi
                         ,mc.cores = 8)
ZpredhomGPnoncal = Reduce(rbind,ZpredhomGPnoncal)
ZpredhomGPcal = Reduce(rbind,ZpredhomGPcal)
```

``` r
ZpredhetGPnoncal = mclapply(1:npred,prednoncal,GP=het,vareps=vareps,loc=testdesign,Ym=Ym,Yv=Yv,mc.cores = 8)
ZpredhetGPcal = mclapply(1:npred,predcal,GP=het,cal=calHet,
                         vareps=vareps,loc=testdesign,YfN=YfieldnoiseN,Xfield=Xfield,Ym=Ym,Yv=Yv,psi=hetMLEdisc$psi
                         ,mc.cores = 8)
ZpredhetGPnoncal = Reduce(rbind,ZpredhetGPnoncal)
ZpredhetGPcal = Reduce(rbind,ZpredhetGPcal)
```

``` r
ZpredseqhetGPnoncal = mclapply(1:npred,prednoncal,GP=seqhet,vareps=vareps,loc=testdesign,Ym=Ym,Yv=Yv,mc.cores = 8)
ZpredseqhetGPcal = mclapply(1:npred,predcal,GP=seqhet,cal=calSeqHet,
                            vareps=vareps,loc=testdesign,YfN=YfieldnoiseN,Xfield=Xfield,Ym=Ym,Yv=Yv,psi=seqhetMLEdisc$psi,
                            mc.cores = 8)
ZpredseqhetGPnoncal = Reduce(rbind,ZpredseqhetGPnoncal)
ZpredseqhetGPcal = Reduce(rbind,ZpredseqhetGPcal)
```

### Some plots

``` r
nlocplot=10
for (k in 1:nlocplot){
    dfpredhom = data.frame(cal=ZpredhomGPcal[,k],noncal = ZpredhomGPnoncal[,k] )
  df2predhom = melt(dfpredhom)
  ghom = ggplot(df2predhom,aes(x=value, fill=variable)) + geom_density(alpha=0.25)+ theme_bw()  + geom_vline(aes(xintercept=Ztest.mean[k]),color="red")+ggtitle("homGP")
  dfpredhet = data.frame(cal=ZpredhetGPcal[,k],noncal = ZpredhetGPnoncal[,k] )
  df2predhet = melt(dfpredhet)
  ghet = ggplot(df2predhet,aes(x=value, fill=variable)) + geom_density(alpha=0.25)+ theme_bw()  + geom_vline(aes(xintercept=Ztest.mean[k]),color="red")+ggtitle("hetGP")
  dfpredseqhet = data.frame(cal=ZpredseqhetGPcal[,k],noncal = ZpredseqhetGPnoncal[,k] )
  df2predsehet = melt(dfpredseqhet)
  ghetseq = ggplot(df2predsehet,aes(x=value, fill=variable)) + geom_density(alpha=0.25)+ theme_bw()  + geom_vline(aes(xintercept=Ztest.mean[k]),color="red")+ggtitle("seqhetGP")
  print(grid.arrange(ghom,ghet,ghetseq))
}
```

    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables

    ## TableGrob (3 x 1) "arrange": 3 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (2-2,1-1) arrange gtable[layout]
    ## 3 3 (3-3,1-1) arrange gtable[layout]

    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables

![](OceanCalibration_files/figure-gfm/plotspredic-1.png)<!-- -->

    ## TableGrob (3 x 1) "arrange": 3 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (2-2,1-1) arrange gtable[layout]
    ## 3 3 (3-3,1-1) arrange gtable[layout]

    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables

![](OceanCalibration_files/figure-gfm/plotspredic-2.png)<!-- -->

    ## TableGrob (3 x 1) "arrange": 3 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (2-2,1-1) arrange gtable[layout]
    ## 3 3 (3-3,1-1) arrange gtable[layout]

    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables

![](OceanCalibration_files/figure-gfm/plotspredic-3.png)<!-- -->

    ## TableGrob (3 x 1) "arrange": 3 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (2-2,1-1) arrange gtable[layout]
    ## 3 3 (3-3,1-1) arrange gtable[layout]

    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables

![](OceanCalibration_files/figure-gfm/plotspredic-4.png)<!-- -->

    ## TableGrob (3 x 1) "arrange": 3 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (2-2,1-1) arrange gtable[layout]
    ## 3 3 (3-3,1-1) arrange gtable[layout]

    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables

![](OceanCalibration_files/figure-gfm/plotspredic-5.png)<!-- -->

    ## TableGrob (3 x 1) "arrange": 3 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (2-2,1-1) arrange gtable[layout]
    ## 3 3 (3-3,1-1) arrange gtable[layout]

    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables

![](OceanCalibration_files/figure-gfm/plotspredic-6.png)<!-- -->

    ## TableGrob (3 x 1) "arrange": 3 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (2-2,1-1) arrange gtable[layout]
    ## 3 3 (3-3,1-1) arrange gtable[layout]

    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables

![](OceanCalibration_files/figure-gfm/plotspredic-7.png)<!-- -->

    ## TableGrob (3 x 1) "arrange": 3 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (2-2,1-1) arrange gtable[layout]
    ## 3 3 (3-3,1-1) arrange gtable[layout]

    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables

![](OceanCalibration_files/figure-gfm/plotspredic-8.png)<!-- -->

    ## TableGrob (3 x 1) "arrange": 3 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (2-2,1-1) arrange gtable[layout]
    ## 3 3 (3-3,1-1) arrange gtable[layout]

    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables
    ## No id variables; using all as measure variables

![](OceanCalibration_files/figure-gfm/plotspredic-9.png)<!-- -->![](OceanCalibration_files/figure-gfm/plotspredic-10.png)<!-- -->

    ## TableGrob (3 x 1) "arrange": 3 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (2-2,1-1) arrange gtable[layout]
    ## 3 3 (3-3,1-1) arrange gtable[layout]

### MSEs and scores

``` r
RMSEhomnoncal = sqrt(mean((rowMeans(ZpredhomGPnoncal)-Ztest.mean)^2))
RMSEhomcal = sqrt(mean((rowMeans(ZpredhomGPcal)-Ztest.mean)^2))
RMSEhetnoncal = sqrt(mean((rowMeans(ZpredhetGPnoncal)-Ztest.mean)^2))
RMSEhetcal = sqrt(mean((rowMeans(ZpredhetGPcal)-Ztest.mean)^2))
RMSEseqhetnoncal = sqrt(mean((rowMeans(ZpredseqhetGPnoncal)-Ztest.mean)^2))
RMSEseqhetcal = sqrt(mean((rowMeans(ZpredseqhetGPcal)-Ztest.mean)^2))
```

``` r
Scorehomnoncal = scoreEstDens(ZpredhomGPnoncal,Ztest.sim)
Scorehomcal = scoreEstDens(ZpredhomGPcal,Ztest.sim)
Scorehetnoncal = scoreEstDens(ZpredhetGPnoncal,Ztest.sim)
Scorehetcal = scoreEstDens(ZpredhetGPcal,Ztest.sim)
Scoreseqhetnoncal = scoreEstDens(ZpredseqhetGPnoncal,Ztest.sim)
Scoreseqhetcal = scoreEstDens(ZpredseqhetGPcal,Ztest.sim)
```

## Other pre-calibration

We try other three other calibrations from L2 minimization, and two
guesses (a good and a bad one) for the calibration parameter.

SS minimizer for homGP

``` r
L2minhom = optim(c(.5,.5),SumOfSquares,lower=0,upper=1,GP=hom,XF=Xfield[,1:2],yF=YfieldnoiseN)
```

    ## Warning in optim(c(0.5, 0.5), SumOfSquares, lower = 0, upper = 1, GP = hom, :
    ## bounds can only be used with method L-BFGS-B (or Brent)

``` r
L2minhom$par* 900 +100
```

    ## [1] 1000.0000  318.0109

``` r
L2minhom$value
```

    ## [1] 25.6542

``` r
SumOfSquares((c(700,200)-100)/900,hom,Xfield[,1:2],YfieldnoiseN)
```

    ## [1] 27.93046

SS minimizer for hetGP

``` r
L2minhet = optim(c(.5,.5),SumOfSquares,lower=0,upper=1,GP=het,XF=Xfield[,1:2],yF=YfieldnoiseN)
```

    ## Warning in optim(c(0.5, 0.5), SumOfSquares, lower = 0, upper = 1, GP = het, :
    ## bounds can only be used with method L-BFGS-B (or Brent)

``` r
L2minhet$par* 900 +100
```

    ## [1] 1000.00  337.22

``` r
L2minhet$value
```

    ## [1] 26.48761

``` r
SumOfSquares((c(700,200)-100)/900,het,Xfield[,1:2],YfieldnoiseN)
```

    ## [1] 28.55229

SS minimizer for seqhetGP

``` r
L2minseqhet = optim(c(.5,.5),SumOfSquares,lower=0,upper=1,GP=seqhet,XF=Xfield[,1:2],yF=YfieldnoiseN)
```

    ## Warning in optim(c(0.5, 0.5), SumOfSquares, lower = 0, upper = 1, GP = seqhet, :
    ## bounds can only be used with method L-BFGS-B (or Brent)

``` r
L2minseqhet$par* 900 +100
```

    ## [1] 605.2122 283.5692

``` r
L2minseqhet$value
```

    ## [1] 28.25017

``` r
SumOfSquares((c(700,200)-100)/900,seqhet,Xfield[,1:2],YfieldnoiseN)
```

    ## [1] 29.72664

We propose a “good guess” and a “bad guess”.

``` r
goodguess = (c(750, 215) -100)/900 
badguess = (c(200,600) -100)/900
```

Prediction for homGP

``` r
ZpredhomGPnoncalL2 = mclapply(1:npred,prednoncalfixed,GP=hom,u=L2minhom$par,vareps=vareps,Ym=Ym,Yv=Yv,loc=testdesign,mc.cores = 8)
ZpredhomGPnoncalL2 = Reduce(rbind,ZpredhomGPnoncalL2)
ZpredhomGPnoncalgoodguess = mclapply(1:npred,prednoncalfixed,GP=hom,u=goodguess,vareps=vareps,Ym=Ym,Yv=Yv,loc=testdesign,mc.cores = 8)
ZpredhomGPnoncalgoodguess = Reduce(rbind,ZpredhomGPnoncalgoodguess)
ZpredhomGPnoncalbadguess = mclapply(1:npred,prednoncalfixed,GP=hom,u=badguess,vareps=vareps,Ym=Ym,Yv=Yv,loc=testdesign,mc.cores = 8)
ZpredhomGPnoncalbadguess = Reduce(rbind,ZpredhomGPnoncalbadguess)
```

Prediction for hetGP

``` r
ZpredhetGPnoncalL2 = mclapply(1:npred,prednoncalfixed,GP=het,u=L2minhet$par,vareps=vareps,Ym=Ym,Yv=Yv,loc=testdesign,mc.cores = 8)
ZpredhetGPnoncalL2 = Reduce(rbind,ZpredhetGPnoncalL2)
ZpredhetGPnoncalgoodguess = mclapply(1:npred,prednoncalfixed,GP=het,u=goodguess,vareps=vareps,Ym=Ym,Yv=Yv,loc=testdesign,mc.cores = 8)
ZpredhetGPnoncalgoodguess = Reduce(rbind,ZpredhetGPnoncalgoodguess)
ZpredhetGPnoncalbadguess = mclapply(1:npred,prednoncalfixed,GP=het,u=badguess,vareps=vareps,Ym=Ym,Yv=Yv,loc=testdesign,mc.cores = 8)
ZpredhetGPnoncalbadguess = Reduce(rbind,ZpredhetGPnoncalbadguess)
```

Prediction for seqhetGP

``` r
ZpredseqhetGPnoncalL2 = mclapply(1:npred,prednoncalfixed,GP=seqhet,u=L2minseqhet$par,vareps=vareps,Ym=Ym,Yv=Yv,loc=testdesign,mc.cores = 8)
ZpredseqhetGPnoncalL2 = Reduce(rbind,ZpredseqhetGPnoncalL2)
ZpredseqhetGPnoncalgoodguess = mclapply(1:npred,prednoncalfixed,GP=seqhet,u=goodguess,vareps=vareps,Ym=Ym,Yv=Yv,loc=testdesign,mc.cores = 8)
ZpredseqhetGPnoncalgoodguess = Reduce(rbind,ZpredseqhetGPnoncalgoodguess)
ZpredseqhetGPnoncalbadguess = mclapply(1:npred,prednoncalfixed,GP=seqhet,u=badguess,vareps=vareps,Ym=Ym,Yv=Yv,loc=testdesign,mc.cores = 8)
ZpredseqhetGPnoncalbadguess = Reduce(rbind,ZpredseqhetGPnoncalbadguess)
```

``` r
RMSEhomnoncalL2 = sqrt(mean((rowMeans(ZpredhomGPnoncalL2)-Ztest.mean)^2))
RMSEhomnoncalgoodguess = sqrt(mean((rowMeans(ZpredhomGPnoncalgoodguess)-Ztest.mean)^2))
RMSEhomnoncalbadguess = sqrt(mean((rowMeans(ZpredhomGPnoncalbadguess)-Ztest.mean)^2))
RMSEhetnoncalL2 = sqrt(mean((rowMeans(ZpredhetGPnoncalL2)-Ztest.mean)^2))
RMSEhetnoncalgoodguess = sqrt(mean((rowMeans(ZpredhetGPnoncalgoodguess)-Ztest.mean)^2))
RMSEhetnoncalbadguess = sqrt(mean((rowMeans(ZpredhetGPnoncalbadguess)-Ztest.mean)^2))
RMSEseqhetnoncalL2 = sqrt(mean((rowMeans(ZpredseqhetGPnoncalL2)-Ztest.mean)^2))
RMSEseqhetnoncalgoodguess = sqrt(mean((rowMeans(ZpredseqhetGPnoncalgoodguess)-Ztest.mean)^2))
RMSEseqhetnoncalbadguess = sqrt(mean((rowMeans(ZpredseqhetGPnoncalbadguess)-Ztest.mean)^2))
```

1st row for homGP, 2nd row for hetGP, 3rd row for seqhetGP 1st column
for L2 estimate, 2nd column for bad guess, 3rd column for good guess,
4th for unif prior, 5th post calibration

``` r
rbind(c(RMSEhomnoncalL2,RMSEhomnoncalbadguess,RMSEhomnoncalgoodguess,RMSEhomnoncal,RMSEhomcal),
      c(RMSEhetnoncalL2,RMSEhetnoncalbadguess,RMSEhetnoncalgoodguess,RMSEhetnoncal,RMSEhetcal),
      c(RMSEseqhetnoncalL2,RMSEseqhetnoncalbadguess,RMSEseqhetnoncalgoodguess,RMSEseqhetnoncal,RMSEseqhetcal))
```

    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,] 9.244333 9.338416 9.246637 9.303135 9.244512
    ## [2,] 9.247127 9.352167 9.252128 9.325061 9.239035
    ## [3,] 9.239578 9.356559 9.251953 9.290245 9.245038

``` r
ScorehomnoncalL2 = scoreEstDens(ZpredhomGPnoncalL2,Ztest.sim)
Scorehomnoncalbadguess = scoreEstDens(ZpredhomGPnoncalbadguess,Ztest.sim)
Scorehomnoncalgoodguess = scoreEstDens(ZpredhomGPnoncalgoodguess,Ztest.sim)
ScorehetnoncalL2 = scoreEstDens(ZpredhetGPnoncalL2,Ztest.sim)
Scorehetnoncalbadguess = scoreEstDens(ZpredhetGPnoncalbadguess,Ztest.sim)
Scorehetnoncalgoodguess = scoreEstDens(ZpredhetGPnoncalgoodguess,Ztest.sim)
ScoreseqhetnoncalL2 = scoreEstDens(ZpredseqhetGPnoncalL2,Ztest.sim)
Scoreseqhetnoncalbadguess = scoreEstDens(ZpredseqhetGPnoncalbadguess,Ztest.sim)
Scoreseqhetnoncalgoodguess = scoreEstDens(ZpredseqhetGPnoncalgoodguess,Ztest.sim)
```

1st row for homGP, 2nd row for hetGP, 3rd row for seqhetGP 1st column
for L2 estimate, 2nd column for bad guess, 3rd column for good guess,
4th for unif prior, 5th post calibration

``` r
rbind(c(ScorehomnoncalL2,Scorehomnoncalbadguess,Scorehomnoncalgoodguess,Scorehomnoncal,Scorehomcal),
      c(ScorehetnoncalL2,Scorehetnoncalbadguess,Scorehetnoncalgoodguess,Scorehetnoncal,Scorehetcal),
      c(ScoreseqhetnoncalL2,Scoreseqhetnoncalbadguess,Scoreseqhetnoncalgoodguess,Scoreseqhetnoncal,Scoreseqhetcal))
```

    ##           [,1]      [,2]      [,3]      [,4]      [,5]
    ## [1,] -2.655429 -3.605982 -2.640154 -2.694576 -2.334756
    ## [2,] -2.646136 -3.664159 -2.704185 -2.696873 -2.351505
    ## [3,] -2.740549 -3.604130 -2.823942 -2.739051 -2.317552

<!--

```r
VarPredhomnoncal = mean(apply(ZpredhomGPnoncal,1,var))
VarPredhomcal = mean(apply(ZpredhomGPcal,1,var))
VarPredhetnoncal = mean(apply(ZpredhetGPnoncal,1,var))
VarPredhetcal = mean(apply(ZpredhetGPcal,1,var))
VarPredseqhetnoncal = mean(apply(ZpredseqhetGPnoncal,1,var))
VarPredseqhetcal = mean(apply(ZpredseqhetGPcal,1,var))
VarPredhomnoncalL2 = mean(apply(ZpredhomGPnoncalL2,1,var))
VarPredhomnoncalgoodguess = mean(apply(ZpredhomGPnoncalgoodguess,1,var))
VarPredhomnoncalbadguess = mean(apply(ZpredhomGPnoncalbadguess,1,var))
VarPredhetnoncalL2 = mean(apply(ZpredhetGPnoncalL2,1,var))
VarPredhetnoncalgoodguess = mean(apply(ZpredhetGPnoncalgoodguess,1,var))
VarPredhetnoncalbadguess = mean(apply(ZpredhetGPnoncalbadguess,1,var))
VarPredseqhetnoncalL2 = mean(apply(ZpredseqhetGPnoncalL2,1,var))
VarPredseqhetnoncalgoodguess = mean(apply(ZpredseqhetGPnoncalgoodguess,1,var))
VarPredseqhetnoncalbadguess = mean(apply(ZpredseqhetGPnoncalbadguess,1,var))
```


1st row for homGP, 2nd row for hetGP, 3rd row for seqhetGP
1st column for L2 estimate, 2nd column for bad guess, 3rd column for good guess, 4th for unif prior, 5th post calibration

```r
rbind(c(VarPredhomnoncalL2,VarPredhomnoncalbadguess,VarPredhomnoncalgoodguess,VarPredhomnoncal,VarPredhomcal),
      c(VarPredhetnoncalL2,VarPredhetnoncalbadguess,VarPredhetnoncalgoodguess,VarPredhetnoncal,VarPredhetcal),
      c(VarPredseqhetnoncalL2,VarPredseqhetnoncalbadguess,VarPredseqhetnoncalgoodguess,VarPredseqhetnoncal,VarPredseqhetcal))
```

```
##          [,1]     [,2]     [,3]     [,4]     [,5]
## [1,] 80.63582 75.35902 82.53045 76.73996 86.44405
## [2,] 77.27884 76.81310 80.46505 77.04646 86.14737
## [3,] 82.82664 79.22142 85.34760 80.69367 88.49830
```

-->
