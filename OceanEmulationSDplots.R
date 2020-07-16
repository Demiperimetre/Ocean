
library(ggplot2)
library(colorRamps)
library(gridExtra)
library(DiceDesign)
library(hetGP)
library(MCMCpack)
library(lhs)
library(mvtnorm)
library(reshape2)
library(plyr)
library(parallel)

source("SimulatorAndFunctions.R")



nseed = 40 

mclapply(1:nseed,function(i){
  
  
  
  
  seed= i  #1234 #234 #1 #12345
  set.seed(seed)
  
  
  ## ----space filling design------------------------------------------------------------------------------------------------
  n = 50
  X0 = randomLHS(n,2)
  X0 = maximinSA_LHS(X0)$design
  X <- rbind(X0, X0, X0, X0, X0, X0, X0, X0, X0, X0)
  X <- rbind(X, X)
  Z <- apply(X, 1, simulator)
  
  
  ## ----normalization-------------------------------------------------------------------------------------------------------
  Zm <- mean(Z)
  Zv <- var(Z)
  Z <- (Z - Zm)/sqrt(Zv)
  
  
  ## ----GP settings---------------------------------------------------------------------------------------------------------
  # settings for GPs
  lower <- rep(0.01, 2)
  upper <- rep(10, 2)
  covtype <- "Matern5_2"
  nc <- list(g_min=1e-6, g_bounds=c(1e-6, 1),  
             lowerDelta=log(1e-6))
  settings <- list(linkThetas="none", logN=TRUE, initStrategy="smoothed", 
                   checkHom=TRUE, penalty=TRUE, trace=0, return.matrices=TRUE, return.hom=FALSE)
  control <- list(tol_dist=1e-4, tol_diff=1e-4, multi.start=30)
  
  
  ## ----GP fit--------------------------------------------------------------------------------------------------------------
  ## Homoskedastic
  Ghom <- mleHomGP(X, Z, lower=lower, upper=upper, covtype=covtype,
                   known=list(beta0=0), maxit=1000) 
  
  ## Heteroskedastic 
  Ghet <- mleHetGP(X, Z, lower=lower, upper=upper, covtype=covtype, 
                   noiseControl=nc, settings=settings, known=list(beta0=0), maxit=1000)
  
  
  
  ## ----seqdesign,results="hide"--------------------------------------------------------------------------------------------
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
    cat("i=", i, ", h=", h[i], "\n", sep="")
    
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
  
  
  ## ----grid----------------------------------------------------------------------------------------------------------------
  de = as.matrix(expand.grid(seq(0,1,.01),seq(0,1,.01)))
  predHomde = predict(x = de, object = Ghom)
  gridHom = as.data.frame(cbind(de[,1:2],predHomde$mean,predHomde$sd2,predHomde$nugs,sqrt(predHomde$sd2 + predHomde$nugs)))
  names(gridHom) = c("long","lat","mean","varmean","nug","psd")
  predhetde = predict(x = de, object = Ghet)
  gridhet = as.data.frame(cbind(de[,1:2],predhetde$mean,predhetde$sd2,predhetde$nugs,sqrt(predhetde$sd2 + predhetde$nugs)))
  names(gridhet) = c("long","lat","mean","varmean","nug","psd")
  designseq = as.data.frame(cbind(mod$X0,mod$mult))
  names(designseq) = c("long","lat","rep")
  predhetseqde = predict(x = de, object = mod)
  gridhetseq = as.data.frame(cbind(de[,1:2],predhetseqde$mean,predhetseqde$sd2,predhetseqde$nugs, sqrt(predhetseqde$sd2 + predhetseqde$nugs)))
  names(gridhetseq) = c("long","lat","mean","varmean","nug","psd")
  
  # unnormalize the outputs
  gridHom$mean = gridHom$mean * sqrt(Zv) + Zm
  gridHom$varmean = gridHom$varmean * Zv
  gridHom$psd = gridHom$psd * sqrt(Zv)
  gridHom$nug = gridHom$nug * Zv
  gridhet$mean = gridhet$mean * sqrt(Zv) + Zm
  gridhet$varmean = gridhet$varmean * Zv
  gridhet$psd = gridhet$psd * sqrt(Zv)
  gridhet$nug = gridhet$nug * Zv
  gridhetseq$mean = gridhetseq$mean * sqrt(Zv) + Zm
  gridhetseq$varmean = gridhetseq$varmean * Zv
  gridhetseq$psd = gridhetseq$psd * sqrt(Zv)
  gridhetseq$nug = gridhetseq$nug * Zv
  
  
  ## ----plot homgp----------------------------------------------------------------------------------------------------------
  # ggplot(gridHom, aes(long, lat)) + geom_raster(aes(fill = mean), interpolate = TRUE) + scale_fill_gradientn(colours=matlab.like(10),limits=c(min(gridHom$mean,gridhet$mean),max(gridHom$mean,gridhet$mean)))
  # ggplot(gridHom, aes(long, lat)) + geom_raster(aes(fill = psd), interpolate = TRUE)+ scale_fill_gradientn(colours=matlab.like(10),limits=c(0,max(gridHom$psd,gridhet$psd,6.5)))
  
  
  ## ----plot hetgp----------------------------------------------------------------------------------------------------------
  # ggplot(gridhet, aes(long, lat)) + geom_raster(aes(fill = mean), interpolate = TRUE) +scale_fill_gradientn(colours=matlab.like(10),limits=c(min(gridHom$mean,gridhet$mean),max(gridHom$mean,gridhet$mean)))
  sdhet = ggplot(gridhet, aes(long, lat)) + geom_raster(aes(fill = psd), interpolate = TRUE)+ scale_fill_gradientn(colours=matlab.like(10),limits=c(0,max(gridHom$psd,gridhet$psd,6.5)))
  
  
  ## ----plot seqhetgp-------------------------------------------------------------------------------------------------------
  # ggplot(gridhetseq, aes(long, lat)) + geom_raster(aes(fill = mean), interpolate = TRUE) + scale_fill_gradientn(colours=matlab.like(10))
  # ggplot(gridhetseq, aes(long, lat)) + geom_raster(aes(fill = mean), interpolate = TRUE) + scale_fill_gradientn(colours=matlab.like(10))+geom_point(data=designseq,aes(x=long,y=lat,colour=rep))+scale_color_gradient(low="grey", high="black")
  # ggplot(gridhetseq, aes(long, lat)) + geom_raster(aes(fill = psd), interpolate = TRUE)+ scale_fill_gradientn(colours=matlab.like(10),limits=c(0,max(gridHom$psd,gridhet$psd,6.5)))+geom_point(data=designseq,aes(x=long,y=lat,colour=rep))+scale_color_gradient(low="grey", high="black")
  sdseqhet =  ggplot(gridhetseq, aes(long, lat)) + geom_raster(aes(fill = psd), interpolate = TRUE)+ scale_fill_gradientn(colours=matlab.like(10),limits=c(0,max(gridHom$psd,gridhet$psd,6.5)))+geom_text(data=designseq,aes(x=long,y=lat,label=rep))
  
  
  ## ----load data test2-----------------------------------------------------------------------------------------------------
  test = read.csv("testdata2D.csv",sep=" ")
  testdesign = as.matrix(test[,1:2])
  Ztest.mean = test[,3]
  #Hom GP
  predhom = predict(x = testdesign, object = Ghom)
  #Het GP
  predhet = predict(x = testdesign, object = Ghet)
  #Het GP seq
  predHetseq = predict(x = testdesign, object = seqGhet)
  
  # Normalize Ztest.mean
  Ztest.mean.N = (Ztest.mean - Zm)/sqrt(Zv)
  
  # MSE
  msehom = mean(((predhom$mean-Ztest.mean.N))^2)
  msehet = mean(((predhet$mean-Ztest.mean.N))^2)
  msehetseq = mean(((predHetseq$mean-Ztest.mean.N))^2)
  
  # scores for the mean prediction
  schom = mean(-(Ztest.mean.N-predhom$mean)^2/(predhom$sd2) -log(predhom$sd2))
  schet = mean(-(Ztest.mean.N-predhet$mean)^2/(predhet$sd2) -log(predhet$sd2))
  schetseq = mean(-(Ztest.mean.N-predHetseq$mean)^2/(predHetseq$sd2) -log(predHetseq$sd2))
  
  # scores for the prediction of a single run of the simulator
  Ztest = (apply(testdesign,1,simulator)-Zm)/sqrt(Zv) # for computing scores on a single realization of the simulator
  schom2 = mean(-(Ztest-predhom$mean)^2/(predhom$sd2+predhom$nugs) -log(predhom$sd2+predhom$nugs))
  schet2 = mean(-(Ztest-predhet$mean)^2/(predhet$sd2+predhet$nugs) -log(predhet$sd2+predhet$nugs))
  schetseq2 = mean(-(Ztest-predHetseq$mean)^2/(predHetseq$sd2+predHetseq$nugs) -log(predHetseq$sd2+predHetseq$nugs))
  
  
  ## ----results-------------------------------------------------------------------------------------------------------------
  c(msehom=msehom,msehet=msehet,msehetseq=msehetseq)
  c(scorehom=schom,scorehet=schet,scorehetseq=schetseq)
  c(scorehom=schom2,scorehet=schet2,scorehetseq=schetseq2)
  
  
  ## ----sdtest--------------------------------------------------------------------------------------------------------------
  #ggplot(test,aes(x=long,y=lat,col=m))+geom_point()+scale_color_gradientn(colours=matlab.like(10),limits=c(min(gridHom$mean,gridhet$mean),max(gridHom$mean,gridhet$mean)))+theme_bw()
  sdtrue = ggplot(test,aes(x=long,y=lat,col=sd))+geom_point()+scale_color_gradientn(colours=matlab.like(10),limits=c(0,max(gridHom$psd,gridhet$psd,6.5)))+theme_bw()
  
  
  

  name=paste0("Emul2Dseed",seed,".pdf")
  pdf(file=name)
  grid.arrange(sdhet,sdseqhet,sdtrue)
  dev.off()
},mc.cores = 10)
