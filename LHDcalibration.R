require(lhs)


### DESIGN OPTIMIZED LHS RESTRICTED THE  FIELD VALUES###
maximinSA_LHS_grid=function (design,dim_xf,dim_theta, T0 = 10, c = 0.95, it = 2000, p = 50, profile = "GEOM", 
    Imax = 100) 
{
    crit <- NULL
    temp <- NULL
    proba <- NULL
    if (profile == "GEOM") {
        m <- design
        i <- 0
        T <- T0
        fi_p <- phiP(m, p)
        crit <- fi_p
        while (T > 0 & i < it) {
            G <- perturbationLHS(m,dim_xf,dim_theta)
            fi_p_ep <- phiP(G, p)
            diff <- min(exp((fi_p - fi_p_ep)/T), 1)
            if (diff == 1) {
                m <- G
                fi_p <- fi_p_ep
            }
            else {
                Bernoulli <- rbinom(1, 1, diff)
                if (Bernoulli == 1) {
                  m <- G
                  fi_p <- fi_p_ep
                }
            }
            i <- i + 1
            crit <- c(crit, fi_p)
            temp <- c(temp, T)
            proba <- c(proba, diff)
            T <- (c^i) * (T0)
      }
}
    return(m)
    
}


# For M an LHS, this function outputs another LHS as a small perturbation of M
perturbationLHS=function(M,dim_xf,dim_theta) 
{  
n=nrow(M)
deb=dim_xf+1
fin=dim_xf+dim_theta
G=M
u=trunc(runif(2,0,n))  ## choosing row
u=u+1
u1=trunc(runif(1,0,2)) ## choosing either x (field input variable) or theta (calibration parameter)
u1=u1+1                
if(u1==1)
 {x=G[u[1],1:dim_xf]
  G[u[1],1:dim_xf]=G[u[2],1:dim_xf]
  G[u[2],1:dim_xf]=x}
else {x=G[u[1],deb:fin]
      G[u[1],deb:fin]=G[u[2],deb:fin]
      G[u[2],deb:fin]=x}
return(G)
} 


GRID_design=function(Xf,b_inf_theta,b_sup_theta,size_design)
 #Xf=input field data (matrix or vector)
 #b_inf_theta=lower bound of the parameter to be calibrated
 #b_sup_theta=upper bound of the parameter to be calibrated
 #design_point=number of point in the design 
{ Xf=as.matrix(Xf)
  dim_theta=length(b_inf_theta)
  N=nrow(Xf) 
  number_bloc=size_design%/%N
  rest=size_design%%N
  Xff=Xf
  if(number_bloc>=2)
   {for(i in 1:(number_bloc-1))
      {Xff=rbind(Xff,Xf)
      }}
  if(rest>0)
   {if(ncol(Xf)>1)
     {Xf=Xf[1:rest,]
      Xff=rbind(Xff,Xf)}
    else {Xf=as.matrix(Xf[1:rest,])
          Xff=rbind(Xff,Xf)}
   }
  
  Theta=c()
  for(j in 1:dim_theta)
  {M=randomLHS(size_design,1)*(b_sup_theta[j]-b_inf_theta[j])+b_inf_theta[j]
  Theta=cbind(Theta,M)}
  out=cbind(Xff,Theta)
  return(out)

}
  