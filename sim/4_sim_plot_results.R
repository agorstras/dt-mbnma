################################################################################
# Description: Display simulation results
# - Figures: S3-S5 (supplement)
# - Tables: 2-4 (main manuscript)
################################################################################

rm(list = ls())
library(ggplot2)
library(tidyverse)

# set location
loc = "sim"

################################################################################
getParm = function(run, model, scenario, path){
  
  lst = lapply(1:length(run), function(j){
    
    temp = readRDS( paste0(path, run[j]) )
    # ITCs
    resList2 = lapply(1:100, function(i){
      tmp = as.data.frame(temp[[i]][[1]]) 
      tmp$rep = i
      tmp$parameter = rownames(tmp)
      tmp
    })  
    parms = do.call("rbind", resList2)
    parms$Model = model[j]
    parms$Scenario = scenario[j]
    parms
  })
  
  plotDat = do.call("rbind", lst)
  return(plotDat)
}  

################################################################################
getDat = function(run, model, scenario, path){
  lst = lapply(1:length(run), function(j){
    
    temp = readRDS( paste0(path, run[j]) )
    # ITCs
    resList2 = lapply(1:100, function(i){
      if(j == 1) tmp = as.data.frame(temp[[i]][[2]]) else tmp = as.data.frame(temp[[i]]) 
      tmp$rep = i
      tmp
    })  
    releff = do.call("rbind", resList2)
    names(releff) = c("parameter", "Estimate", "Est.Error", "Q2.5", "Q97.5", "rep")
    releff = subset(releff, parameter %in% c("Estimate_AvsC", "Estimate_BvsC") )
    releff$true = rep( c( -10*1/(exp(-2) + 1) - -20*2/( exp(-0.75) + 2), -15*3/( exp(-1.5) + 3) - -20*2/( exp(-0.75) + 2) ), times=100)
    releff$Model = model[j]
    releff
  })
  
  plotDat = do.call("rbind", lst)
  plotDat$Scenario = scenario
  return(plotDat)
}  


################################################################################
# Estimation model: univariate likelihood

# set path to simulation results
path  = paste0(loc, "/sim_results/univariate_ll/")

# FIGURE S3 (univariate)

run   = c("dosetimeNMAList.rds", "dosetimeNMAList2.rds", "dosetimeNMAList3.rds") 
scenario = c("(1)", "(2)", "(3)")
model = c("Time-course: Exponential \n Dose-response: Emax", 
          "Time-course: Exponential \n Dose-response: Emax",
          "Time-course: Exponential \n Dose-response: Emax")

plotDat = getParm(run = run, model = model, scenario = scenario, path = path)
plotDat = subset(plotDat, parameter %in% c("emaxA_Intercept", "emaxB_Intercept", "emaxC_Intercept", 
                                           "ed50A_Intercept", "ed50B_Intercept", "ed50C_Intercept",
                                           "sA_Intercept"   , "sB_Intercept"   , "sC_Intercept"))
plotDat$true = rep(c(-10,-15,-20,-2,-1.5,-0.75,-0.2,-0.4,-0.8), times=100)
labbrm = as_labeller(c(emaxA_Intercept="EmaxD^(A)",emaxB_Intercept="EmaxD^(B)",emaxC_Intercept="EmaxD^(C)", 
                     ed50A_Intercept="ED50^(A)", ed50B_Intercept="ED50^(B)", ed50C_Intercept="ED50^(C)",
                     sA_Intercept   ="b^(A)",sB_Intercept       ="b^(B)",sC_Intercept       ="b^(C)"), default = label_parsed)
  
ggplot(plotDat, aes(x=Scenario, y=Estimate)) + facet_wrap(~ parameter, labeller = labbrm, scales = "free_y") + geom_boxplot() +
     geom_hline(data = plotDat, aes(yintercept = true), linetype=2, colour="darkgrey") + theme_bw() + ylab("Parameter Estimates")


# TABLE 2 in main manuscript

# Scenario 1  
run   = c("dosetimeNMAList.rds", "timeNMAList.rds", "NMAList.rds") 
model = c("Time-course Dose-response", 
          "Time-course NMA",
          "Standard NMA")
plotDat1 = getDat(run = run, model = model, scenario = "(1)", path = path)

# Scenario 2
run   = c("dosetimeNMAList2.rds", "timeNMAList.rds", "NMAList.rds") 
plotDat2 = getDat(run = run, model = model, scenario = "(2)", path = path)

# Scenario 3
run      = c("dosetimeNMAList3.rds") 
plotDat3 = getDat(run = run, model = "Time-course Dose-response", scenario = "(3)", path = path)

plotDat = rbind(plotDat1, plotDat2)

evalres = rbind(plotDat, plotDat3) %>% 
  group_by(Scenario, parameter, Model) %>% 
  summarise(
    mean = mean(Estimate),
    true = unique(true),
    bias = mean( (Estimate - true ) ),
    rmse = sqrt( mean( ( Estimate - true )**2) ),
    mad  = mean( abs(Estimate - true ) )
  )

out = as.data.frame(evalres)
out$true = format(round(out$true, 2), nsmall = 2)
out$mean = format(round(out$mean, 2), nsmall = 2)
out$bias = format(round(out$bias, 2), nsmall = 2)
out$rmse = format(round(out$rmse, 2), nsmall = 2)
out$mad  = format(round(out$mad, 2), nsmall = 2)

out


################################################################################
# Estimation model: multi-variate likelihood - AR1 within-arm correlations

# set path
path  = paste0(loc, "/sim_results/AR1/")

# FIGURE S4 (AR-1)

run   = c("dosetimeNMAList.rds", "dosetimeNMAList2.rds", "dosetimeNMAList3.rds") 
scenario = c("(1)", "(2)", "(3)")
model = c("Time-course: Exponential \n Dose-response: Emax", 
          "Time-course: Exponential \n Dose-response: Emax",
          "Time-course: Exponential \n Dose-response: Emax")

plotDat = getParm(run = run, model = model, scenario = scenario, path = path)
plotDat = subset(plotDat, parameter %in% c("emaxA_Intercept", "emaxB_Intercept", "emaxC_Intercept", 
                                           "ed50A_Intercept", "ed50B_Intercept", "ed50C_Intercept",
                                           "sA_Intercept"   , "sB_Intercept"   , "sC_Intercept"))
plotDat$true = rep(c(-10,-15,-20,-2,-1.5,-0.75,-0.2,-0.4,-0.8), times=100)

ggplot(plotDat, aes(x=Scenario, y=Estimate)) + facet_wrap(~ parameter, labeller = labbrm, scales = "free_y") + geom_boxplot() +
  geom_hline(data = plotDat, aes(yintercept = true), linetype=2, colour="darkgrey") + theme_bw() + ylab("Parameter Estimates")

# Table 3 in main manscript

run   = c("dosetimeNMAList.rds", "timeNMAList.rds", "NMAList.rds") 
model = c("Time-course Dose-response", 
          "Time-course NMA",
          "Standard NMA")
plotDat1 = getDat(run = run, model = model, scenario = "(1)", path = path)

# Scenario 2
run   = c("dosetimeNMAList2.rds", "timeNMAList.rds", "NMAList.rds") 
plotDat2 = getDat(run = run, model = model, scenario = "(2)", path = path)

# Scenario 3
run      = c("dosetimeNMAList3.rds") 
plotDat3 = getDat(run = run, model = "Time-course Dose-response", scenario = "(3)", path = path)

plotDat = rbind(plotDat1, plotDat2)

evalres = rbind(plotDat, plotDat3) %>% 
  group_by(Scenario, parameter, Model) %>% 
  summarise(
    mean = mean(Estimate),
    true = unique(true),
    bias = mean( (Estimate - true ) ),
    rmse = sqrt( mean( ( Estimate - true )**2) ),
    mad  = mean( abs(Estimate - true ) )
  )

out = as.data.frame(evalres)
out$true = format(round(out$true, 2), nsmall = 2)
out$mean = format(round(out$mean, 2), nsmall = 2)
out$bias = format(round(out$bias, 2), nsmall = 2)
out$rmse = format(round(out$rmse, 2), nsmall = 2)
out$mad  = format(round(out$mad, 2), nsmall = 2)

out

################################################################################
# Estimation model: multivariate likelihood - use within-arm correlations

# set path 
path  = paste0(loc, "/sim_results/obs/")

# FIGURE S5 (observed correlations)

run   = c("dosetimeNMAList.rds", "dosetimeNMAList2.rds", "dosetimeNMAList3.rds") 
scenario = c("(1)", "(2)", "(3)")
model = c("Time-course: Exponential \n Dose-response: Emax", 
          "Time-course: Exponential \n Dose-response: Emax",
          "Time-course: Exponential \n Dose-response: Emax")

plotDat = getParm(run = run, model = model, scenario = scenario, path = path)

plotDat = subset(plotDat, parameter %in% c("b_emaxA[1]", "b_emaxB[1]", "b_emaxC[1]", 
                                           "b_ed50A[1]", "b_ed50B[1]", "b_ed50C[1]", 
                                           "b_sA[1]",    "b_sB[1]",    "b_sC[1]"))

plotDat$true = rep(c(-10,-15,-20,-2,-1.5,-0.75,-0.2,-0.4,-0.8), times=100)

labrstan = as_labeller(c("b_emaxA[1]"="EmaxD^(A)", "b_emaxB[1]"="EmaxD^(B)", "b_emaxC[1]" ="EmaxD^(C)", 
                     "b_ed50A[1]"="ED50^(A)", "b_ed50B[1]" ="ED50^(B)",  "b_ed50C[1]" ="ED50^(C)",
                     "b_sA[1]"   ="b^(A)",    "b_sB[1]"    ="b^(B)",     "b_sC[1]"   ="b^(C)" ), default = label_parsed)

ggplot(plotDat, aes(x=Scenario, y=mean)) + facet_wrap(~ parameter, labeller = labrstan, scales = "free_y") + geom_boxplot() +
  geom_hline(data = plotDat, aes(yintercept = true), linetype=2, colour="darkgrey") + theme_bw() + ylab("Parameter Estimates")

#TABLE 4  

run   = c("dosetimeNMAList.rds", "timeNMAList.rds", "NMAList.rds") 
model = c("Time-course Dose-response", 
          "Time-course NMA",
          "Standard NMA")
plotDat1 = getDat(run = run, model = model, scenario = "(1)", path = path)

# Scenario 2
run   = c("dosetimeNMAList2.rds", "timeNMAList.rds", "NMAList.rds") 
plotDat2 = getDat(run = run, model = model, scenario = "(2)", path = path)

# Scenario 3
run      = c("dosetimeNMAList3.rds") 
plotDat3 = getDat(run = run, model = "Time-course Dose-response", scenario = "(3)", path = path)

plotDat = rbind(plotDat1, plotDat2)

evalres = rbind(plotDat, plotDat3) %>% 
  group_by(Scenario, parameter, Model) %>% 
  summarise(
    mean = mean(Estimate),
    true = unique(true),
    bias = mean( (Estimate - true ) ),
    rmse = sqrt( mean( ( Estimate - true )**2) ),
    mad  = mean( abs(Estimate - true ) )
  )

out = as.data.frame(evalres)
out$true = format(round(out$true, 2), nsmall = 2)
out$mean = format(round(out$mean, 2), nsmall = 2)
out$bias = format(round(out$bias, 2), nsmall = 2)
out$rmse = format(round(out$rmse, 2), nsmall = 2)
out$mad  = format(round(out$mad, 2), nsmall = 2)

out