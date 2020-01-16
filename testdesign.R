rm(list=ls())


# load the simulator
source("SimulatorAndFunctions.R")


# generate the test design
set.seed(1234)
testdesign = lhs::randomLHS(500,2)
testdesign = DiceDesign::maximinSA_LHS(testdesign)$design



vm = numeric(nrow(testdesign)) # mean of the simulator
vs = numeric(nrow(testdesign)) # sd between runs at the same location
vn = numeric(nrow(testdesign)) # number of runs for each location

# simulator is run for specified diffusion coefficients KX=700 and KY=200


# loop over location
for (i in 1:nrow(testdesign))
{
   res = as.vector(simulatorpar(testdesign[i,]))
   m = mean(res)
   s = sd(res)  
   n = length(res)
   #repeat until standard error is less than .1
   while( s/sqrt(n) >.1 )
   {
     res = c(res,as.vector(simulatorpar(testdesign[i,])))
     s = sd(res)
     n = length(res)
     m = mean(res)
   }
   vm[i] = m
   vs[i] = s
   vn[i] = length(res)
   print(i)
}


testdata = data.frame(long=testdesign[,1],lat=testdesign[,2],m=vm,sd=vs,num=vn)
write.table(testdata,file="testdata2D.csv")


