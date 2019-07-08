library(MASS)
library(tidyverse) # watch out MASS::select

pk1coma = function(CL, V, Ka, BioA=1, DosingHistory, Time, PropE=0, AddE=0, LLoQ=0, Jitter=0)
{
  nObs = length(Time)
  Conc = rep(0, nObs)
  ke = CL/V
  
  if (Jitter > 0) Time = round(jitter(Time), Jitter)
  Time[Time < 0] = 0
  
  TERM1 = BioA*Ka/(Ka - ke)/V
  DH = DosingHistory[DosingHistory[,2] > 0,,drop=FALSE]
  nAmt = nrow(DH)
  for (i in 1:nAmt) {
    TERM2 = DH[i, 2]*TERM1
    dTime = Time - DH[i, 1]
    dTime2 = dTime[dTime >= 0]
    Conc = Conc + c(rep(0, length(dTime) - length(dTime2)), TERM2*(exp(-ke*dTime2) - exp(-Ka*dTime2)))
  }
  
  Err1 = rnorm(nObs, mean=0, sd=PropE)
  Err2 = rnorm(nObs, mean=0, sd=AddE)
  
  Conc = Conc + Conc*Err1 + Err2
  Conc[Conc < LLoQ] = 0
  return(cbind(Time, Conc))
}


Var <- list()
CL = 20
V = 100
Ka = 2
mu <- c(CL, V, Ka)

Time = 0:48
nObs = length(Time)

DH1 = matrix(c(0, 100000), nrow=1, ncol=2, byrow=TRUE)
# DH2 = matrix(c(0, 100000, 12, 100000, 24, 100000), nrow=3, ncol=2, byrow=TRUE)

nSubj = 36
nTrt = 1

FullCov <- matrix(c(0.1,0,0,0,0,0,0,0,0), nrow = 3, byrow=TRUE)

rpk = MASS::mvrnorm(nSubj*nTrt, rep(0, 3), FullCov) ; rpk
iPK = matrix(rep(mu, nSubj*nTrt), nrow=nSubj*nTrt, byrow=TRUE) * exp(rpk) ; iPK

BA=1

individual_Conc <- list()
for (i in 1:nSubj) {
  cSUBJ = i
  cCL = iPK[i, 1]
  cV  = iPK[i, 2]
  cKa = iPK[i, 3]
  BA = 1
  iConc = as.data.frame(pk1coma(cCL, cV, cKa, BioA=BA, DH1, Time, PropE=0.1, AddE=1, LLoQ=0, Jitter=2)) %>% 
    mutate(SUBJ = i)
  individual_Conc[[i]] <- iConc
}
Conc <- bind_rows(individual_Conc)
