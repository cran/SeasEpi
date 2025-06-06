#' SeasEpi for a Simulation Study
#'
#' This function conducts a simulation study for spatiotemporal individual-level model of seasonal infectious disease transmission within the Susceptible-Exposed-Infectious-Recovered-Susceptible (SEIRS) framework, using a user-defined grid size. It applies a likelihood based Monte Carlo Expectation Conditional Maximization (MCECM) algorithm to estimate model parameters and compute the AIC.
#' @param GridDim1 First dimension of the grid
#' @param GridDim2 Second dimension of the grid
#' @param NPostPerGrid Number of postal codes per grid cell
#' @param MaxTimePand Last time point of the pandemic
#' @param tau0 Initial value for spatial precision
#' @param lambda0 Initial value for spatial dependence
#' @param alphaS0 Initial value for the susceptibility intercept
#' @param delta0 Initial value for the spatial decay parameter
#' @param zeta10 Initial value for the amplitude of the seasonal oscillation parameter (sin part)
#' @param zeta20 Initial value for the phase of the seasonal oscillation parameter (cos part)
#' @param alphaT0 Initial value for the infectivity intercept
#' @param PopMin Minimum population per postal code
#' @param PopMax Maximum population per postal code
#' @param InfFraction Fraction of each grid cell's population to be infected
#' @param InfPrd Infectious period that can be obtained either from the literature or by fitting an SEIRS model to the data
#' @param IncPrd Incubation period that can be obtained either from the literature or by fitting an SEIRS model to the data
#' @param NIterMC Number of MCMC iterations
#' @param NIterMCECM Number of MCECM iterations
#' @param T_cycle The duration of a complete seasonal cycle (e.g., 12 months for an annual cycle)
#' @return
#'   `alphaS` Estimate of alpha S
#'
#'   `BetaCovInf` Estimate of beta vector for the individual level infection covariate
#'
#'   `BetaCovSus` Estimate of beta vector for the areal susceptibility to first infection covariate
#'
#'   `alphaT` Estimate of alpha T
#'
#'   `delta` Estimate of delta
#'
#'   `zeta1` Estimate of zeta1
#'
#'   `zeta2` Estimate of zeta2
#'
#'   `tau1` Estimate of tau
#'
#'   `lambda1` Estimate of lambda
#'
#'   `AIC` AIC of the fitted GDILM SEIRS
#'
#' @export
#' @import MASS
#' @import mvtnorm
#' @import ngspatial
#' @import stats
#' @examples
#' \donttest{
#' SeasEpi_Sim_Par_Est(5,5,10,30,0.7, 0.7, -1, 0.1, 0,40, 50,0.6, 5, 5, 10, 3,0.2,0.2,5)
#' }
#'
SeasEpi_Sim_Par_Est=function(GridDim1,GridDim2,NPostPerGrid,MaxTimePand,tau0, lambda0, alphaS0, delta0, alphaT0,PopMin, PopMax,InfFraction, InfPrd, IncPrd, NIterMC, NIterMCECM,zeta10,zeta20,T_cycle){
 if(lambda0>1) stop("The spatial dependence parameter should be restricted to a range between 0 and 1.")
 if(lambda0==0) stop("Absence of spatial dependence: This model is designed for scenarios where spatial dependence is present.")
 if(delta0<0) stop("The spatial decay parameter must be greater than zero.")
 if(NIterMC<2) stop("The number of iterations must exceed 2.")
 if(InfPrd<0) stop("The infectious period must be greater than zero.")
 if(IncPrd<0) stop("The incubation period must be greater than zero.")
 if(InfFraction>1) stop("Fraction of each grid cell's population to be infected must be be restricted to a range between 0 and 1.")
 if(InfFraction<0) stop("Fraction of each grid cell's population to be infected must be be restricted to a range between 0 and 1.")
NTotalpost=NPostPerGrid*GridDim1^2
NTotalGrid=GridDim1^2
NAllPostPerGrid=rep(NPostPerGrid,GridDim1*GridDim1)
NTotalpost=sum(NAllPostPerGrid)
generate_grid_data <- function(NPostPerGrid, GridDim1, GridDim2) {
xxx <- c()
yyy <- c()
for (col in 1:GridDim1) {
x_col <- matrix(0, NPostPerGrid, GridDim2)
for (i in 1:GridDim2) {
x_col[, i] <- runif(NPostPerGrid, GridDim1 * (col - 1), GridDim1 * col)
}
xxx <- c(xxx, c(x_col))
}
for (row in 1:GridDim2) {
y_row <- matrix(0, NPostPerGrid, GridDim1)
for (i in 1:GridDim1) { y_row[, i] <- runif(NPostPerGrid, GridDim2 * (i - 1), GridDim2 * i) }
yyy <- c(yyy, c(y_row))
}
data <- data.frame(xxx = xxx, yyy = yyy)
return(data)
}
data=generate_grid_data(NPostPerGrid,GridDim1,GridDim2)
Lat=data[,1]
Long=data[,2]
################### Adjacency matrix #######################
I=diag(NTotalGrid)
A1=adjacency.matrix(GridDim1, GridDim1)
A2=colSums(A1)
D=-A1
diag(D)=A2
##################### Region labels ########################
NLableGrid1=list()
for(i in 1:NTotalGrid){NLableGrid1[[i]]=rep(i,NAllPostPerGrid[i])}
NLableGrid=unlist(NLableGrid1)
NewLabelGrid=matrix(0,NTotalpost,NTotalGrid)
for(GridIndic in 1:NTotalGrid){NewLabelGrid[,GridIndic]=rep(D[,GridIndic],NAllPostPerGrid)}
#################### Distance matrix #######################
Dist=matrix(0,NTotalpost,NTotalpost)
for(i in 1:NTotalpost){
for(j in 1:NTotalpost){
Dist[i,j]=sqrt((Lat[i]-Lat[j])^2+(Long[i]-Long[j])^2)
}
}
Dist=Dist*50
################## Population size #########################
################### Infected size ##########################
Pop=sample(PopMin:PopMax, NTotalpost, replace=TRUE)
NInf1=round(Pop*InfFraction)
NInf=ifelse(NInf1==0,1,NInf1)
################### Spatial rendom effects #################
mu=rep(0,NTotalGrid)
Sigma0=solve(tau0^2*(lambda0*D+(1-lambda0)*I))
phi=mvrnorm(1, mu, Sigma0, tol = 1e-6)
######################### Covaraites #######################
CovSus=as.matrix(cbind(rnorm(NTotalGrid,1,1),runif(NTotalGrid,1,2)))
CovInf=as.matrix(cbind(rnorm(NTotalpost,0,1),runif(NTotalpost,0,1)))
DimCovInf=dim(CovInf)[2]
DimCovSus=dim(CovSus)[2]
BetaCovInf0=rep(1,DimCovInf)
BetaCovSus0=rep(2,DimCovSus)
######################## Infected MaxTimePand #####################
D1=c()
for(GridIndic in 1:NTotalGrid){D1[GridIndic]=sample(NAllPostPerGrid[GridIndic],1,replace=F)}
ExpoTime1=list()
for(i in 1:NTotalGrid){
      ExpoTime1[[i]]=rep(0,NAllPostPerGrid[i])
      ExpoTime1[[i]][D1[i]]=1
    }
    ExpoTime=unlist(ExpoTime1)
    InfPeriod=rep(InfPrd,NTotalpost) ############# ifected period ###############
    IncPeriod=rep(IncPrd,NTotalpost)    ############ Incubation  period ############
    IncPeriod1=IncPeriod[1]
    InfTime=ifelse(ExpoTime>0,ExpoTime+IncPeriod1,ExpoTime)

    for(t in 1:MaxTimePand){
      for(i in 1:NTotalpost){
        for(GridIndic in 1:NTotalGrid){
          if (NLableGrid[i]==GridIndic){
            if(InfTime[i]== 0){
              dx1=rep(0,NTotalpost)
              for(j in 1:NTotalpost){
                if (NewLabelGrid[j,GridIndic]!=0){
                  if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
                    dx1[j]=exp(alphaT0+CovInf[j,]%*%BetaCovInf0-delta0*Dist[i,j]+zeta10*sin(2*pi*t/T_cycle)+zeta20*cos(2*pi*t/T_cycle))
                  }
                }
              }
              dx1=replace(dx1,is.infinite(dx1),0)
              dx=sum(dx1)
              P=1-exp(-exp(alphaS0+CovSus[GridIndic,]%*%BetaCovSus0+phi[GridIndic])*dx)
              u=runif(1,0,1)
              if(P>u){
                InfTime[i]=t+1
              }
            }
          }
        }
      }
    }
    InfTime
    ExpoTime=InfTime-IncPeriod1


    ##############################################################################
    SumH=function(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT){
      SumH1=rep(0,NTotalpost)
      for(j in 1:NTotalpost){
        if (NewLabelGrid[j,GridIndic]!=0){
          if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
            SumH1[j]=exp(alphaT+CovInf[j,]%*%BetaCovInf-delta*Dist[i,j]+zeta1*sin(2*pi*t/T_cycle)+zeta2*cos(2*pi*t/T_cycle))
          }
        }
      }
      SumH1=replace(SumH1,is.infinite(SumH1),0)
      SumH2=sum(SumH1)
      SumH3=exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus)*SumH2
      return(SumH3)
    }

##############################################################################
    SumW=function(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT){
      SumW1=array(0,c(DimCovInf,1,NTotalpost))
      for(j in 1:NTotalpost){
        if (NewLabelGrid[j,GridIndic]!=0){
          if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
            SumW1[,,j]=CovInf[j,]*as.numeric(exp(alphaT+CovInf[j,]%*%BetaCovInf)-delta*Dist[i,j]+zeta1*sin(2*pi*t/T_cycle)+zeta2*cos(2*pi*t/T_cycle))
          }
        }
      }
      SumW1[!is.finite(SumW1)]=NA
      SumW2=apply(SumW1,c(1,2),sum,na.rm=T)
      SumW3=SumW2*as.numeric(exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus))
      return(SumW3)
    }
##############################################################################
    SumNewW1=function(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT){
      SumW11=array(0,c(DimCovInf,DimCovInf,NTotalpost))
      for(j in 1:NTotalpost){
        if (NewLabelGrid[j,GridIndic]!=0){
          if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
            SumW11[,,j]=CovInf[j,]%*%t(CovInf[j,])*as.numeric(exp(alphaT+CovInf[j,]%*%BetaCovInf)-delta*Dist[i,j]+zeta1*sin(2*pi*t/T_cycle)+zeta2*cos(2*pi*t/T_cycle))
          }
        }
      }
      SumW11[!is.finite(SumW11)]=NA
      SumW21=apply(SumW11,c(1,2),sum,na.rm=T)
      SumW31=SumW21*as.numeric(exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus))
      return(SumW31)
    }
##############################################################################
    SumHnew1=function(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT){
      SumH4=rep(0,NTotalpost)
      for(j in 1:NTotalpost){
        if (NewLabelGrid[j,GridIndic]!=0){
          if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
            SumH4[j]=exp(alphaT+CovInf[j,]%*%BetaCovInf)*Dist[i,j]^(-delta)*(log(Dist[i,j]))
          }
        }
      }
      SumH4=replace(SumH4,is.infinite(SumH4),0)
      SumH5=sum(SumH4)
      SumH6=exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus)*SumH5
      return(SumH6)
    }
##############################################################################
    SumHnew2=function(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT){
      SumH7=rep(0,NTotalpost)
      for(j in 1:NTotalpost){
        if (NewLabelGrid[j,GridIndic]!=0){
          if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
            SumH7[j]=exp(alphaT+CovInf[j,]%*%BetaCovInf-delta*Dist[i,j]+zeta1*sin(2*pi*t/T_cycle)+zeta2*cos(2*pi*t/T_cycle))
          }
        }
      }
      SumH7=replace(SumH7,is.infinite(SumH7),0)
      SumH8=sum(SumH7)
      SumH9=exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus)*SumH8
      return(SumH9)
    }


    ##############################################################################
    SumHnew1z1=function(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT){
      SumH4=rep(0,NTotalpost)
      for(j in 1:NTotalpost){
        if (NewLabelGrid[j,GridIndic]!=0){
          if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
            SumH4[j]=sin(2*pi*t/T_cycle)*exp(alphaT+CovInf[j,]%*%BetaCovInf-delta*Dist[i,j]+zeta1*sin(2*pi*t/T_cycle)+zeta2*cos(2*pi*t/T_cycle))
          }
        }
      }
      SumH4=replace(SumH4,is.infinite(SumH4),0)
      SumH5=sum(SumH4)
      SumH6=exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus)*SumH5
      return(SumH6)
    }
    ##############################################################################
    SumHnew2z1=function(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT){
      SumH7=rep(0,NTotalpost)
      for(j in 1:NTotalpost){
        if (NewLabelGrid[j,GridIndic]!=0){
          if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
            SumH7[j]=(sin(2*pi*t/T_cycle))^2*exp(alphaT+CovInf[j,]%*%BetaCovInf-delta*Dist[i,j]+zeta1*sin(2*pi*t/T_cycle)+zeta2*cos(2*pi*t/T_cycle))
          }
        }
      }
      SumH7=replace(SumH7,is.infinite(SumH7),0)
      SumH8=sum(SumH7)
      SumH9=exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus)*SumH8
      return(SumH9)
    }


    ##############################################################################
    SumHnew1z2=function(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT){
      SumH4=rep(0,NTotalpost)
      for(j in 1:NTotalpost){
        if (NewLabelGrid[j,GridIndic]!=0){
          if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
            SumH4[j]=cos(2*pi*t/T_cycle)*exp(alphaT+CovInf[j,]%*%BetaCovInf-delta*Dist[i,j]+zeta1*sin(2*pi*t/T_cycle)+zeta2*cos(2*pi*t/T_cycle))
          }
        }
      }
      SumH4=replace(SumH4,is.infinite(SumH4),0)
      SumH5=sum(SumH4)
      SumH6=exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus)*SumH5
      return(SumH6)
    }
    ##############################################################################
    SumHnew2z2=function(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT){
      SumH7=rep(0,NTotalpost)
      for(j in 1:NTotalpost){
        if (NewLabelGrid[j,GridIndic]!=0){
          if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
            SumH7[j]=(cos(2*pi*t/T_cycle))^2*exp(alphaT+CovInf[j,]%*%BetaCovInf-delta*Dist[i,j]+zeta1*sin(2*pi*t/T_cycle)+zeta2*cos(2*pi*t/T_cycle))
          }
        }
      }
      SumH7=replace(SumH7,is.infinite(SumH7),0)
      SumH8=sum(SumH7)
      SumH9=exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus)*SumH8
      return(SumH9)
    }

    ################################ Expectations ##########################
    Fy1=function(phi,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,tau1,BetaCovInf,BetaCovSus,alphaT){
      fy1=array(0,c(NTotalpost,MaxTimePand,NTotalGrid))
      for(i in 1:NTotalpost){
        for(t in 1:MaxTimePand){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                dx1=rep(0,NTotalpost)
                for(j in 1:NTotalpost){
                  if (NewLabelGrid[j,GridIndic]!=0){
                    if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
                      dx1[j]=exp(alphaT+CovInf[j,]%*%BetaCovInf-delta*Dist[i,j]+zeta1*sin(2*pi*t/T_cycle)+zeta2*cos(2*pi*t/T_cycle))
                    }
                  }
                }
                dx1=replace(dx1,is.infinite(dx1),0)
                dx=sum(dx1)
                prob1=1-exp(-exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus+phi[GridIndic])*dx)
                fy1[i,t,GridIndic]=(1-prob1)
              }

              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                dx4=rep(0,NTotalpost)
                for(j in 1:NTotalpost){
                  if (NewLabelGrid[j,GridIndic]!=0){
                    if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
                      dx4[j]=exp(alphaT+CovInf[j,]%*%BetaCovInf-delta*Dist[i,j]+zeta1*sin(2*pi*t/T_cycle)+zeta2*cos(2*pi*t/T_cycle))
                    }
                  }
                }
                dx4=replace(dx4,is.infinite(dx4),0)
                dx5=sum(dx4)
                prob2=1-exp(-exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus+phi[GridIndic])*dx5)
                fy1[i,t,GridIndic]=(prob2)
              }
            }
          }
        }
      }
      PP=c()
      H=list()
      fy2=c()
      for(GridIndic in 1:NTotalGrid){
        for(i in 1:NTotalpost){
          PP[i]=round(prod(fy1[i,,GridIndic][fy1[i,,GridIndic]>0]),10)
        }
        H[[GridIndic]]=which(PP!=1)
        fy2[GridIndic]=prod(PP[H[[GridIndic]]][PP[H[[GridIndic]]]>0])
      }
      return(fy2)
    }
##############################################################################
    alphaS=alphaS0
    delta=delta0
    zeta1=zeta10
    zeta2=zeta20
    tau1=tau0
    lambda1=lambda0
    BetaCovInf=BetaCovInf0
    BetaCovSus=BetaCovSus0
    alphaT=alphaT0
    ##############################################################################
    estfun=function(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,tau1,BetaCovInf,BetaCovSus,alphaT){
      Pos=matrix(0,NIterMC+1,NTotalGrid)
      mu=rep(0,NTotalGrid)
      Sigma1=solve(tau1^2*(lambda1*D+(1-lambda1)*I))
      phi0=mvrnorm(1, mu, Sigma1, tol = 1e-6)

      Pos[1,]=phi0
      Uni=c()
      MPHphi=c()
      for(L in 2:NIterMC){
        phi=mvrnorm(1, mu, Sigma1, tol = 1e-6)
        UEST1=Fy1(phi,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,tau1,BetaCovInf,BetaCovSus,alphaT)
        UEST2=Fy1(Pos[L-1,],alphaS,delta,zeta1,zeta2,T_cycle,lambda1,tau1,BetaCovInf,BetaCovSus,alphaT)
        Uni[L]=runif(1,0,1)
        MPHphi[L]=min(1,prod(UEST1/UEST2))
        if(Uni[L]<MPHphi[L]){
          Pos[L,]=phi
        }
        if(Uni[L]>=MPHphi[L]){
          Pos[L,]=Pos[L-1,]
        }
      }

 mean1=function(Pos,GridIndic){
        d1=c()
        for(L in 1:NIterMC){
          d1[L]=exp(Pos[L,GridIndic])
        }
        SMean1=mean(d1)
        return(SMean1)
      }
mean2=function(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L){
        dx1=rep(0,NTotalpost)
        for(j in 1:NTotalpost){
          if (NewLabelGrid[j,GridIndic]!=0){
            if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
              dx1[j]=exp(alphaT+CovInf[j,]%*%BetaCovInf-delta*Dist[i,j]+zeta1*sin(2*pi*t/T_cycle)+zeta2*cos(2*pi*t/T_cycle))
            }
          }
        }
        dx1=replace(dx1,is.infinite(dx1),0)
        dx=sum(dx1)
        prob1=1-exp(-exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus+Pos[L,GridIndic])*dx)
        if(prob1==0){SMean2=0}
        if(prob1!=0){
          SMean2=(1-prob1)/prob1*exp(Pos[L,GridIndic])
        }
        return(SMean2)
      }
mean3=function(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L){
      dx1=rep(0,NTotalpost)
        for(j in 1:NTotalpost){
          if (NewLabelGrid[j,GridIndic]!=0){
            if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
              dx1[j]=exp(alphaT+CovInf[j,]%*%BetaCovInf-delta*Dist[i,j]+zeta1*sin(2*pi*t/T_cycle)+zeta2*cos(2*pi*t/T_cycle))
            }
          }
        }
        dx1=replace(dx1,is.infinite(dx1),0)
        dx=sum(dx1)
        prob1=1-exp(-exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus+Pos[L,GridIndic])*dx)
        if(prob1==0){SMean3=0}
        if(prob1!=0){
          SMean3=(1-prob1)/prob1^2*exp(2*Pos[L,GridIndic])
        }
        return(SMean3)
      }

      ######################################## alpha #######################################
      A1=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        A2=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if (NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                A2[i]=-as.numeric(SumH(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT)*mean1(Pos,GridIndic))
              }
            }
          }
        }
        A1[t]=sum(A2)
      }
      SusA1=sum(A1)
      ########################################
      A3=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        A4=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                SA4=c()
                for(L in 1:NIterMC){
                  SA4[L]=mean2(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }

                A4[i]=as.numeric(SumH(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT)*mean(SA4))
              }
            }
          }
        }
        A3[t]=sum(A4)
      }
      InfA3=sum(A3,na.rm=T)
      EndA3=SusA1+InfA3
      ######################### Second Der...########################
    A5=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        A6=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if (NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                SA6=c()
                for(L in 1:NIterMC){
                  SA6[L]=mean3(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                A6[i]=-1*as.numeric((SumH(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT))^2*mean(SA6))
              }
            }
          }
        }
        A5[t]=sum(A6)
      }
      InfA5=sum(A5,na.rm=T)
      EndA5=EndA3+InfA5
      EstAlphaS=alphaS-EndA3/EndA5
      #################################BetaCovSus##################################
      Y1=array(0,c(DimCovSus,1,MaxTimePand))
      for(t in 1:MaxTimePand){
        Y2=array(0,c(DimCovSus,1,NTotalpost))
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if (NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                Y2[,,i]=-CovSus[GridIndic,]*as.numeric(SumH(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT)*mean1(Pos,GridIndic))
              }
            }
          }
        }
        Y1[,,t]=apply(Y2,c(1,2),sum)
      }
      SusY1=apply(Y1,c(1,2),sum)
      ##############################################################################
      Y3=array(0,c(DimCovSus,1,MaxTimePand))
      for(t in 1:MaxTimePand){
        Y4=array(0,c(DimCovSus,1,NTotalpost))
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                YA4=c()
                for(L in 1:NIterMC){
                  YA4[L]=mean2(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                Y4[,,i]=CovSus[GridIndic,]*as.numeric(SumH(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT)*mean(YA4))
              }
            }
          }
        }
        Y3[,,t]=apply(Y4,c(1,2),sum,na.rm=T)
      }
      InfY3=apply(Y3,c(1,2),sum)
      EndY3=SusY1+InfY3
      ######################### Second Der...########################
 YY1=array(0,c(DimCovSus,DimCovSus,MaxTimePand))
      for(t in 1:MaxTimePand){
        YY2=array(0,c(DimCovSus,DimCovSus,NTotalpost))
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if (NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                YY2[,,i]=-CovSus[GridIndic,]%*%t(CovSus[GridIndic,])*as.numeric(SumH(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT)*mean1(Pos,GridIndic))
              }
            }
          }
        }
        YY1[,,t]=apply(YY2,c(1,2),sum)
      }
      SusYY1=apply(YY1,c(1,2),sum)
      ##############################################################################
      YY3=array(0,c(DimCovSus,DimCovSus,MaxTimePand))
      for(t in 1:MaxTimePand){
        YY4=array(0,c(DimCovSus,DimCovSus,NTotalpost))
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                YA45=c()
                for(L in 1:NIterMC){
                  YA45[L]=mean2(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                YY4[,,i]=CovSus[GridIndic,]%*%t(CovSus[GridIndic,])*as.numeric(SumH(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT)*mean(YA45))
              }
            }
          }
        }
        YY3[,,t]=apply(YY4,c(1,2),sum)
      }
      InfYY3=apply(YY3,c(1,2),sum,na.rm=T)

      Y5=array(0,c(DimCovSus,DimCovSus,MaxTimePand))
      for(t in 1:MaxTimePand){
        Y6=array(0,c(DimCovSus,DimCovSus,NTotalpost))
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if (NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                YA65=c()
                for(L in 1:NIterMC){
                  YA65[L]=mean3(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                Y6[,,i]=-1*CovSus[GridIndic,]%*%t(CovSus[GridIndic,])*as.numeric((SumH(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT))^2*mean(YA65))
              }
            }
          }
        }
        Y5[,,t]=apply(Y6,c(1,2),sum,na.rm=T)
      }
      InfY5=apply(Y5,c(1,2),sum)
      EndY5=SusYY1+InfYY3+InfY5

      epsilon <- 1e-3
      if (det(EndY5) < epsilon) {
        EndY5 <- EndY5 + diag(epsilon, nrow(EndY5))
      }

      EstBetaCovSus=BetaCovSus-solve(EndY5)%*%EndY3



      ######################################## alphaT #######################################

      TA1=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        TA2=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if (NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                TA2[i]=-as.numeric(SumH(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,EstBetaCovSus,alphaT)*mean1(Pos,GridIndic))
              }
            }
          }
        }
        TA1[t]=sum(TA2)
      }
      TSusA1=sum(TA1)
      ##############################################################################
      TA3=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        TA4=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                IA4=c()
                for(L in 1:NIterMC){
                  IA4[L]=mean2(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                TA4[i]=as.numeric(SumH(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,EstBetaCovSus,alphaT)*mean(IA4))
              }
            }
          }
        }
        TA3[t]=sum(TA4)
      }
      TInfA3=sum(TA3,na.rm=T)
      TEndA3=TSusA1+TInfA3
      ######################### Second Der...########################

      TA5=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        TA6=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if (NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                IA5=c()
                for(L in 1:NIterMC){
                  IA5[L]=mean3(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                TA6[i]=-1*as.numeric((SumH(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,EstBetaCovSus,alphaT))^2*mean(IA5))
              }
            }
          }
        }
        TA5[t]=sum(TA6)
      }
      TInfA5=sum(TA5,na.rm=T)
      TEndA5=TEndA3+TInfA5
      EstAlphaT=alphaT-TEndA3/TEndA5

      #################################BetaCovInf##################################

      B1=array(0,c(DimCovInf,1,MaxTimePand))
      for(t in 1:MaxTimePand){
        B2=array(0,c(DimCovInf,1,NTotalpost))
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if (NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                B2[,,i]=-SumW(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,EstBetaCovSus,EstAlphaT)*mean1(Pos,GridIndic)
              }
            }
          }
        }
        B1[,,t]=apply(B2,c(1,2),sum,na.rm=T)
      }
      SusB1=apply(B1,c(1,2),sum,na.rm=T)

      ##############################################################################

      B3=array(0,c(DimCovInf,1,MaxTimePand))
      for(t in 1:MaxTimePand){
        B4=array(0,c(DimCovInf,1,NTotalpost))
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                IAW4=c()
                for(L in 1:NIterMC){
                  IAW4[L]=mean2(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                B4[,,i]=SumW(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,EstBetaCovSus,EstAlphaT)*as.numeric(mean(IAW4))
              }
            }
          }
        }
        B3[,,t]=apply(B4,c(1,2),sum,na.rm=T)
      }
      InfB3=apply(B3,c(1,2),sum,na.rm=T)
      EndB3=SusB1+InfB3
      ######################### Second Der...########################

      BB1=array(0,c(DimCovInf,DimCovInf,MaxTimePand))
      for(t in 1:MaxTimePand){
        BB2=array(0,c(DimCovInf,DimCovInf,NTotalpost))
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if (NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                BB2[,,i]=-SumNewW1(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,EstBetaCovSus,EstAlphaT)*as.numeric(mean1(Pos,GridIndic))
              }
            }
          }
        }
        BB1[,,t]=apply(BB2,c(1,2),sum,na.rm=T)
      }
      SusBB1=apply(BB1,c(1,2),sum,na.rm=T)

      ##############################################################################

      BB3=array(0,c(DimCovInf,DimCovInf,MaxTimePand))
      for(t in 1:MaxTimePand){
        BB4=array(0,c(DimCovInf,DimCovInf,NTotalpost))
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){

                IAW5=c()
                for(L in 1:NIterMC){
                  IAW5[L]=mean2(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                BB4[,,i]=SumNewW1(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,EstBetaCovSus,EstAlphaT)*as.numeric(mean(IAW5))
              }
            }
          }
        }
        BB3[,,t]=apply(BB4,c(1,2),sum,na.rm=T)
      }
      InfBB3=apply(BB3,c(1,2),sum,na.rm=T)
      EndBB3=SusBB1+InfBB3

      B5=array(0,c(DimCovInf,DimCovInf,MaxTimePand))
      for(t in 1:MaxTimePand){
        B6=array(0,c(DimCovInf,DimCovInf,NTotalpost))
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if (NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                IAW6=c()
                for(L in 1:NIterMC){
                  IAW6[L]=mean3(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                B6[,,i]=-1*SumW(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,EstBetaCovSus,EstAlphaT)%*%t(SumW(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,EstBetaCovSus,EstAlphaT))*as.numeric(mean(IAW6))
              }
            }
          }
        }
        B5[,,t]=apply(B6,c(1,2),sum,na.rm=T)
      }
      InfB5=apply(B5,c(1,2),sum,na.rm=T)
      EndB5=EndBB3+InfB5
      EstBetaCovInf=BetaCovInf-solve(EndB5)%*%EndB3


      ################################# Delta #############################

      S8=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        S7=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                S7[i]=as.numeric(SumHnew1(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT)*mean1(Pos,GridIndic))
              }
            }
          }
        }
        S8[t]=sum(S7,na.rm=T)
      }
      SusNew8=sum(S8,na.rm=T)

      ##############################################################################

      I10=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        I9=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                ID1=c()
                for(L in 1:NIterMC){
                  ID1[L]=mean2(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                I9[i]=-as.numeric(SumHnew1(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT)*mean(ID1))
              }
            }
          }
        }
        I10[t]=sum(I9)
      }
      InfeNew10=sum(I10,na.rm=T)
      EndNew10=SusNew8+InfeNew10
      ##############################################################################

      S10=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        S9=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                S9[i]=-as.numeric(SumHnew2(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT)*mean1(Pos,GridIndic))
              }
            }
          }
        }
        S10[t]=sum(S9)
      }
      SusNew10=sum(S10)
      ##############################################################################

      I12=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        I11=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                ID2=c()
                for(L in 1:NIterMC){
                  ID2[L]=mean2(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                ID3=c()
                for(L in 1:NIterMC){
                  ID3[L]=mean3(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                I11[i]=as.numeric(SumHnew2(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT)*mean(ID2))-1*as.numeric((SumHnew1(NLableGrid,Dist,EstAlphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT))^2*mean(ID3))
              }
            }
          }
        }
        I12[t]=sum(I11)
      }
      InfeNew12=sum(I12,na.rm=T)
      EndNew12=SusNew10+InfeNew12
      Estdelta=delta-EndNew10/EndNew12





      ################################# zeta1 #############################

      S8z1=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        S7z1=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                S7z1[i]=as.numeric(SumHnew1z1(NLableGrid,Dist,EstAlphaS,Estdelta,zeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT)*mean1(Pos,GridIndic))
              }
            }
          }
        }
        S8z1[t]=sum(S7z1,na.rm=T)
      }
      SusNew8z1=sum(S8z1,na.rm=T)

      ##############################################################################

      I10z1=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        I9z1=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                ID1z1=c()
                for(L in 1:NIterMC){
                  ID1z1[L]=mean2(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                I9z1[i]=-as.numeric(SumHnew1z1(NLableGrid,Dist,EstAlphaS,Estdelta,zeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT)*mean(ID1z1))
              }
            }
          }
        }
        I10z1[t]=sum(I9z1)
      }
      InfeNew10z1=sum(I10z1,na.rm=T)
      EndNew10z1=SusNew8z1+InfeNew10z1
      ##############################################################################

      S10z1=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        S9z1=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                S9z1[i]=-as.numeric(SumHnew2z1(NLableGrid,Dist,EstAlphaS,Estdelta,zeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT)*mean1(Pos,GridIndic))
              }
            }
          }
        }
        S10z1[t]=sum(S9z1)
      }
      SusNew10z1=sum(S10z1)
      ##############################################################################

      I12z1=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        I11z1=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                ID2z1=c()
                for(L in 1:NIterMC){
                  ID2z1[L]=mean2(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                ID3z1=c()
                for(L in 1:NIterMC){
                  ID3z1[L]=mean3(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                I11z1[i]=as.numeric(SumHnew2z1(NLableGrid,Dist,EstAlphaS,Estdelta,zeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT)*mean(ID2z1))-1*as.numeric((SumHnew1z1(NLableGrid,Dist,EstAlphaS,Estdelta,zeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT))^2*mean(ID3z1))
              }
            }
          }
        }
        I12z1[t]=sum(I11z1)
      }
      InfeNew12z1=sum(I12z1,na.rm=T)
      EndNew12z1=SusNew10z1+InfeNew12z1
      Estzeta1=zeta1-EndNew10z1/EndNew12z1

      ################################# zeta2 #############################

      S8z2=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        S7z2=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                S7z2[i]=as.numeric(SumHnew1z2(NLableGrid,Dist,EstAlphaS,Estdelta,Estzeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT)*mean1(Pos,GridIndic))
              }
            }
          }
        }
        S8z2[t]=sum(S7z2,na.rm=T)
      }
      SusNew8z2=sum(S8z2,na.rm=T)

      ##############################################################################

      I10z2=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        I9z2=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                ID1z2=c()
                for(L in 1:NIterMC){
                  ID1z2[L]=mean2(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                I9z2[i]=-as.numeric(SumHnew1z2(NLableGrid,Dist,EstAlphaS,Estdelta,Estzeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT)*mean(ID1z2))
              }
            }
          }
        }
        I10z2[t]=sum(I9z2)
      }
      InfeNew10z2=sum(I10z2,na.rm=T)
      EndNew10z2=SusNew8z2+InfeNew10z2
      ##############################################################################

      S10z2=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        S9z2=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                S9z2[i]=-as.numeric(SumHnew2z2(NLableGrid,Dist,EstAlphaS,Estdelta,Estzeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT)*mean1(Pos,GridIndic))
              }
            }
          }
        }
        S10z2[t]=sum(S9z2)
      }
      SusNew10z2=sum(S10z2)
      ##############################################################################

      I12z2=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        I11z2=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if(NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                ID2z2=c()
                for(L in 1:NIterMC){
                  ID2z2[L]=mean2(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                ID3z2=c()
                for(L in 1:NIterMC){
                  ID3z2[L]=mean3(NLableGrid,Pos,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                I11z2[i]=as.numeric(SumHnew2z2(NLableGrid,Dist,EstAlphaS,Estdelta,Estzeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT)*mean(ID2z2))-1*as.numeric((SumHnew1z2(NLableGrid,Dist,EstAlphaS,Estdelta,Estzeta1,zeta2,T_cycle,i,GridIndic,t,EstBetaCovInf,EstBetaCovSus,EstAlphaT))^2*mean(ID3z2))
              }
            }
          }
        }
        I12z2[t]=sum(I11z2)
      }
      InfeNew12z2=sum(I12z2,na.rm=T)
      EndNew12z2=SusNew10z2+InfeNew12z2
      Estzeta2=zeta2-EndNew10z2/EndNew12z2

      ##################################### Sigma_U #####################

      Loglik1=function(par){
        lambda1=par[[1]]
        tau1=par[[2]]
        Log1=c()
        for(L in 1:NIterMC){
          Log1[L]=dmvnorm(Pos[L,], rep(0,NTotalGrid), solve(tau1^2*(as.numeric(lambda1)*D+(1-as.numeric(lambda1))*I)), log = t)
        }
        LLL=-mean(Log1)
      }

      init=c(lambda1,tau1)
      m1=optim(init,Loglik1)
      EstU1=m1$par

      EstGammau=EstU1[1]
      HatSigmmaU=EstU1[2]

      ##############################################

      result=list(Pos=Pos,BetaCovInf=EstBetaCovInf,BetaCovSus=EstBetaCovSus,Uhat=EstU1,alphaS=EstAlphaS,alphaT=EstAlphaT,delta=Estdelta,zeta1=Estzeta1,zeta2=Estzeta2,tau1=HatSigmmaU,lambda1=EstGammau)
      result
    }




    OUT2=list()
    LA=numeric()
    Loglik=function(NLableGrid,Pos1,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,tau1,BetaCovInf,BetaCovSus,alphaT){

      mean1L=function(Pos1,GridIndic){
        d1=c()
        for(L in 1:NIterMC){
          d1[L]=exp(Pos1[L,GridIndic])
        }
        SMean1L=mean(d1)
        return(SMean1L)
      }

      mean2L=function(NLableGrid,Pos1,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L){
        dx1=rep(0,NTotalpost)
        for(j in 1:NTotalpost){
          if (NewLabelGrid[j,GridIndic]!=0){
            if(InfTime[j]<=t & (InfTime[j]+InfPeriod[j])>=t & InfTime[j]!=0){
              dx1[j]=exp(alphaT+CovInf[j,]%*%BetaCovInf-delta*Dist[i,j]+zeta1*sin(2*pi*t/T_cycle)+zeta2*cos(2*pi*t/T_cycle))
            }
          }
        }
        dx1=replace(dx1,is.infinite(dx1),0)
        dx=sum(dx1)

        prob1=1-exp(-exp(alphaS+CovSus[GridIndic,]%*%BetaCovSus+Pos1[L,GridIndic])*dx)
        if(prob1==0){SMean2L=0}
        if(prob1!=0){
          SMean2L=log(prob1)
        }
        return(SMean2L)
      }

      A1L=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        A2L=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if (NLableGrid[i]==GridIndic){
              if(ExpoTime[i]>t|ExpoTime[i]==0){
                A2L[i]=-as.numeric(SumH(NLableGrid,Dist,alphaS,delta,zeta1,zeta2,T_cycle,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT)*mean1L(Pos1,GridIndic))
              }
            }
          }
        }
        A1L[t]=sum(A2L)
      }
      SusA1L=sum(A1L)

      A3L=rep(0,MaxTimePand)
      for(t in 1:MaxTimePand){
        A4L=rep(0,NTotalpost)
        for(i in 1:NTotalpost){
          for(GridIndic in 1:NTotalGrid){
            if (NLableGrid[i]==GridIndic){
              if(ExpoTime[i]<=t & (ExpoTime[i]+IncPeriod[i])>t & ExpoTime[i]!=0){
                MM1=c()
                for(L in 1:NIterMC){
                  MM1[L]=mean2L(NLableGrid,Pos1,Dist,alphaS,delta,zeta1,zeta2,T_cycle,lambda1,i,GridIndic,t,BetaCovInf,BetaCovSus,alphaT,L)
                }
                A4L[i]=as.numeric(mean(MM1))
              }
            }
          }
        }
        A3L[t]=sum(A4L)
      }
      INFA1=sum(A3L)

      Log1=c()
      for(L in 1:NIterMC){
        Log1[L]=dmvnorm(Pos1[L,], rep(0,NTotalGrid), solve(tau1^2*(as.numeric(lambda1)*D+(1-as.numeric(lambda1))*I)), log = T)
      }
      LLL=mean(Log1)

      LogLIK1=LLL+INFA1+SusA1L
      return(LogLIK1)
    }


    ##############################################################################
    ##############################################################################
    ##############################################################################
    ##############################################################################


    est0=estfun(NLableGrid,Dist,alphaS0,delta0,zeta10,zeta20,T_cycle,lambda0,tau0,BetaCovInf0,BetaCovSus0,alphaT0)
    alphaS=est0$alphaS
    alphaS
    delta=est0$delta
    delta
    zeta1=est0$zeta1
    zeta1
    zeta2=est0$zeta2
    zeta2
    lambda1=est0$lambda1
    lambda1
    tau1=est0$tau1
    tau1
    BetaCovInf=est0$BetaCovInf
    BetaCovInf
    BetaCovSus=est0$BetaCovSus
    BetaCovSus
    alphaT=est0$alphaT
    alphaT
    Uhat=est0$Uhat
    Uhat
    Pos1=est0$Pos
    Pos1


    ########################################################################
    ########################################################################
    ########################################################################
    E1 <- list()
    AIC <- numeric()
    mes <- numeric()

    tolerance <- 0.1


    for (ss in 1:NIterMCECM) {

      est <- estfun(NLableGrid, Dist, alphaS, delta, zeta1,zeta2,T_cycle,lambda1, tau1, BetaCovInf, BetaCovSus, alphaT)

      alphaS <- est$alphaS
      BetaCovInf <- est$BetaCovInf
      BetaCovSus <- est$BetaCovSus
      delta <- est$delta
      zeta1<- est$zeta1
      zeta2<- est$zeta2
      lambda1 <- est$lambda1
      tau1 <- est$tau1
      alphaT <- est$alphaT
      Uhat <- est$Uhat
      Pos <- est$Pos

      LA[ss] <- Loglik(NLableGrid, Pos, Dist, alphaS, delta, zeta1,zeta2,T_cycle,lambda1, tau1, BetaCovInf, BetaCovSus,alphaT)

      AIC[ss] <- -2 * LA[ss] + 11

      out1 <- list(
        alphaS = alphaS,
        BetaCovInf = BetaCovInf,
        BetaCovSus = BetaCovSus,
        alphaT = alphaT,
        delta = delta,
        zeta1 = zeta1,
        zeta2 = zeta2,
        tau1 = tau1,
        lambda1 = lambda1,
        AIC = AIC[ss]
      )

      E1[[ss]] <- c(alphaS, BetaCovInf, BetaCovSus,alphaT, delta, zeta1,zeta2,tau1, lambda1)

      if (ss > 1) {
        mes[ss] <- sqrt(sum((E1[[ss]] - E1[[ss-1]])^2))
      } else {
        mes[ss] <- Inf
      }

      if (mes[ss] < tolerance) {
        message("MCECM Converged at iteration ", ss, " with parameter change: ", mes[ss])
        break
      }

      if (ss %% 10 == 0) {
        message("Iteration ", ss, ": Parameter change = ", mes[ss])
      }
    }

    out1

  }
