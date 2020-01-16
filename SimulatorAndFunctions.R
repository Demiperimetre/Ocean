require(parallel)


# useful settings for running the simulator
#load("fksettings.Rdata")

xy=c(-34.0, 11.0, -32.0, -3.0)
lat0=-17.5
domain=rep(0,4)
domain[1]=xy[1]*111e3*cos(lat0*pi/180.0)
domain[2]=xy[2]*111e3*cos(lat0*pi/180.0)
domain[3]=xy[3]*111e3
domain[4]=xy[4]*111e3
N = 37
M = 19
delx=(domain[2]-domain[1])/(N-1)
dely=(domain[4]-domain[3])/(M-1)

# velocity vectors
U = read.table('gridu.txt')
U =-U
V = read.table('gridv.txt')
V = -V


#this is needed only for the boundary values
OXY = as.matrix(read.table('oxygengrid.txt'))



# simulator with replications (NPATHS)
fk_simulator  = function(domain, sites, U, V, KX, KY, OXY, NPATHS){
  
  M = dim(OXY)[1]
  N = dim(OXY)[2]
  NOBS = dim(sites)[1]
  HH    = 1e7
  lam   = 1e-11
  
  OUT = matrix(0, nrow = NOBS, ncol = NPATHS)
  
  for(i in 1:NOBS){
    #print(i)
    for(kk in 1:NPATHS){
      
      
      xc = sites[i,1]
      yc = sites[i,2]
      
      tau=0;
      while( ((xc - domain[1])*(xc-domain[2])<0) & ( (yc-domain[3])*(yc-domain[4])<0)    ){
        
        #--------
        # find indices
        jx = ceiling( (xc - domain[1])/delx );
        ix = ceiling( (yc - domain[3])/dely );
        if (jx<=0){
          jx=1
        }
        else{
          if (jx>N){
            jx=N
          }
        }
        
        if (ix<=0){
          ix=1
        }
        else{
          if (ix>M) ix=M
        }
        #------------
        
        uc = U[ix,jx]
        vc = V[ix,jx]
        
        # Euler step
        xn = xc + HH * uc + sqrt(HH)*sqrt(2*KX)*rnorm(1);
        yn = yc + HH * vc + sqrt(HH)*sqrt(2*KY)*rnorm(1);    
        
        xc=xn;
        yc=yn;
        tau = tau + HH;
        
        
      }
      JJ = ceiling( (xc - domain[1])/delx );
      II = ceiling( (yc - domain[3])/dely );
      
      if (JJ<=0){
        JJ=1
      }
      else{
        if (JJ>N) JJ=N
      }
      
      if (II<=0){
        II=1
      }
      else{
        if (II>M) II=M
      }
      
      OUT[i,kk] = OXY[II,JJ]*exp(-lam*tau)
      
    }
  }
  
  return(OUT)
}


# simulator for emulation in 2D with diffusion parameters KX and KY fixed,
# to avoid non gaussianity of the stochastic simulator a default run is the average over 6 paths
# v is 2D normalized in [0,1]^2 longitutde and latitude
simulator = function(v,NPATHS=6) 
{
  #unnormalize
  long = v[1] * (domain[2]-domain[1]) + domain[1]
  lat = v[2] * (domain[4]-domain[3]) + domain[3]
  KX = 700
  KY = 200
  
  sim = fk_simulator(domain, matrix(c(long,lat),nrow=1), U, V, KX,KY, OXY, NPATHS)
  return(mean(sim))
}


# for calibration
# v has 4D:  longitude, latitude, KX and KY (KX and KY are diffusion coefficients)
# v is normalized in [0,1]^4
simulator4d = function(v,NPATHS=6) 
{
  #unnormalize
  long = v[1] * (domain[2]-domain[1]) + domain[1]
  lat = v[2] * (domain[4]-domain[3]) + domain[3]
  KX = v[3] * 900 + 100
  KY = v[4] * 900 + 100
  
  sim = fk_simulator(domain, matrix(c(long,lat),nrow=1), U, V, KX,KY, OXY, NPATHS)
  return(mean(sim))
}



# Run simulator function in parallel
library(parallel)
simulatorpar = function(v)
{
  Reduce(c,mclapply(1:100, function(i) {return(simulator(v))},mc.cores = 10))
}

