library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(MASS)
library(tibble)
library(NonCompart)
# library(ncar)
library(markdown)

# source("function.R")

# Define server logic required for pk1c ----

shinyServer(function(input, output) {
  
  # concTable ----
  
  output$concTable <- renderTable({
    Vars <- Init(input$nSubj, input$CL, input$V, input$Ka, input$DH1, input$FullCov, input$Time, input$PropE, input$AddE, input$Jitter)
    return(Vars$Conc)
  })
 
  # concTimePlot ----
  
  output$concTimePlot <- renderPlot({
    Vars <- Init(input$nSubj, input$CL, input$V, input$Ka, input$DH1, input$FullCov, input$Time,input$PropE, input$AddE, input$Jitter)
    p <- ggplot(Vars$Conc, aes(x=TIME, y=CONC, group=SUBJ, color=SUBJ)) + 
      geom_line() +
      geom_point(size = 3)
    
    if (input$concLog == FALSE) print(p) else print(p + scale_y_log10())
  })
  
  output$concTimeFacet <- renderPlot({
    Vars <- Init(input$nSubj, input$CL, input$V, input$Ka, input$DH1, input$FullCov, input$Time,input$PropE, input$AddE, input$Jitter)
    p <- ggplot(Vars$Conc, aes(x=TIME, y=CONC, group=SUBJ, color=SUBJ)) + 
      geom_line() +
      geom_point(size = 3) +
      facet_wrap(~ SUBJ, ncol = 4)
    
    if (input$concLog == FALSE) print(p) else print(p + scale_y_log10())
  })
  
  output$ncarTable <- renderTable({
    Vars <- Init(input$nSubj, input$CL, input$V, input$Ka, input$DH1, input$FullCov, input$Time, input$PropE, input$AddE, input$Jitter)
    
    steady_state_conc <- Vars$Conc %>% 
      filter(TIME >= 60) %>% 
      mutate(TIME = TIME - 60) %>% 
      as.data.frame() 
    
    resNCA = NonCompart::tblNCA(steady_state_conc, key="SUBJ", colTime="TIME", colConc="CONC") %>% 
      dplyr::select(SUBJ, CMAX, TMAX, AUCLST, LAMZHL)
    return(resNCA)
  })
})

#Init_oral <- function(nSubj, CLo, Vo, Ka, DosingHistory, FullCov, Time, PropE=0, AddE=0, Jitter=0)

#Vars <- Init(input$nSubj, input$CL, input$V, input$Ka, input$DH1, input$FullCov, input$Time,input$PropE, input$AddE, input$Jitter)

Init <- function(nSubj, CL, V, Ka, DH1, FullCov, Time, PropE, AddE, Jitter)
{
  Var <- list()
  
  mu <- c(CL, V, Ka)
  
  Time <-eval(parse(text = paste0("c(", Time, ")"))) 
  
  Time = seq(0, 100)
  nObs <- length(Time)
  
  #DH1 <- eval(parse(text = paste0("c(", DH1, ")")))
  #DH1 <- matrix(DH1, ncol=2, byrow=TRUE) # c(0, 100000)
  
  dosing_amt <- 6000
  
  DH1 <-  matrix(c(0, dosing_amt, 12, dosing_amt, 24, dosing_amt, 36, dosing_amt, 
                   48, dosing_amt, 60, dosing_amt), nrow=6, ncol=2, byrow=TRUE)
  
  #DH1 = matrix(c(0, 100000), nrow=1, ncol=2, byrow=TRUE)
  # DH2 = matrix(c(0, 100000, 12, 100000, 24, 100000), nrow=3, ncol=2, byrow=TRUE)
  
  FullCov <- eval(parse(text = paste0("c(", FullCov, ")")))
  FullCov <- matrix(FullCov, nrow=3)  # CL, V, Ka
  #FullCov <- matrix(c(0,0,0,0,,0,0,0,1), nrow = 3, byrow=TRUE)
  
  nTrt = 1
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
    iConc = as.data.frame(pk1coma(cCL, cV, cKa, BioA=BA, DH1, Time, PropE=PropE, AddE=AddE, LLoQ=0, Jitter=Jitter)) %>% 
      mutate(SUBJ = i)
    individual_Conc[[i]] <- iConc
  }
  
  Conc <- bind_rows(individual_Conc) %>% 
    mutate(SUBJ = as.factor(SUBJ)) %>% 
    dplyr::select(CONC = 2, TIME = 1, SUBJ = 3)
  
  Var$DH1 <- DH1
  Var$Time <- Time
  Var$Conc <- Conc
  
  return(Var)
}
#nSubj <- 36
# Init(10, 2, 10, 1, DH1, FullCov, 0:24, 0, 0, 0)
#pk1coma(2, 10, 1, BioA=1, DH1, 0:24, PropE=0, AddE=0, LLoQ=0, Jitter=0)
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



