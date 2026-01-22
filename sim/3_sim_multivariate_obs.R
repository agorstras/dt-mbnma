################################################################################
# Description: Estimation of Standard NMA, time course, and dose-response time-course
#   - loops over 100 replicas, saved in ./sim_data
#   - dose-response time-course assumes multivariate likelihood, with observed correlations
#   - under 3 different data scenarios
#   - estimates are saved in ./sim_results/obs/
#
################################################################################

rm(list = ls())
options(mc.cores = parallel::detectCores())

library(tidyverse)
library(MASS)
library(brms)
library(Matrix)
library(haven)
library(tidybayes)
library(future.apply)
plan(multisession)

set.seed(1234)

# set location
scer = F
loc  = ifelse(scer, "~/hdrive/mbnma", "H:/mmbna")
file.stan = paste0(loc,"/sim_stan/")

export = T # save estimation datasets?

# set seed
set.seed(1234)

# utility functions

# correlation matrix for correlated random effects
mat.fun = function(i){
  tmp = matrix(NA, i, i)
  diag(tmp) = 1
  tmp[lower.tri(tmp)] = 1/2
  tmp[upper.tri(tmp)] = 1/2
  tmp
  
}

# load data
dataLst = readRDS(paste0(loc,"/sim_data/dataLst.rds"))
varLst = readRDS(paste0(loc,"/sim_data/varLst.rds"))

################################################################################
# loop through 100 simulated datasets - Use standard NMA 

NMAList = future_lapply(future.seed = T, X = seq_along(dataLst), FUN = function(i){
  
  dat = dataLst[[i]]
  dat = subset(dat, week == 52 & study %in% c(2,3,4,6,7,9,10))
  dat$treatment = paste(dat$agent, dat$dose, "mg") 
  dat$treatment[dat$agent == "soc"] = "soc"
  
  # Setup covariance matrix for study effects
  dat$stcom = rep(1:length(unique(paste0(dat$study,dat$arm))), times=table(paste0(dat$study,dat$arm)))
  # Check treatment allocation
  dat2 = dat[!duplicated(paste0(dat$study, dat$arm), fromLast = T), ]
  # Setup covariance matrix, handling trials with multiple arms
  lst = split(dat2, dat2$study)
  lst = lapply(1:length(lst), function(j){
    tmp = lst[[j]]
    tmp = tmp[tmp$arm != 1, ]
    as.matrix(bdiag(matrix(1E-20,1,1), mat.fun(length(tmp$study)))) 
  } )
  st_mat = as.matrix(bdiag(lst))
  rownames(st_mat) <- colnames(st_mat) <- 1:nrow(dat2) 
  dat$study = factor(dat$study)
  
  # Setup design matrix for treatment effects
  d = unique(dat$treatment[dat$treatment!="soc"]); dname = paste0("d",1:length(d))
  for(i in 1:length(d)) dat[,dname[i]] = ifelse(dat$treatment == d[i], 1, 0) 
  
  # Fit model:
  f =  as.formula(paste("y | se(se) ~ 0 + study + ", paste(dname, collapse= "+"), "+ (1 | gr(stcom, cov = st_mat))" ))
  prior1 <- prior(normal(0,100), class = b) + prior(normal(0,5), lb=0, class = sd)
  arm_brms <- brm(f, 
                  data = dat,
                  data2 = list(st_mat=st_mat),
                  prior = prior1,
                  control = list(adapt_delta = 0.999, max_treedepth = 25),
                  iter = 2000, cores = 2, chains = 2)
  
  out <- spread_draws(arm_brms, b_d2, b_d3, b_d6, ndraws = NULL) %>% 
    mutate(Estimate_AvsC = b_d2 - b_d6, 
           Estimate_BvsC = b_d3 - b_d6 ) %>% 
    dplyr::select( .chain, .iteration, .draw,  Estimate_AvsC, Estimate_BvsC ) %>%
    pivot_longer(cols=c("Estimate_AvsC", "Estimate_BvsC"),  names_to='parameter', values_to='values') %>% 
    group_by(parameter) %>% 
    summarise(Estimate = mean(values), 
              Est.Error= sd(values),       
              Q2.5 = quantile(values, probs = 0.025),
              Q97.5= quantile(values, probs = 0.975) )  
  return(out)  
  
})

if(export) saveRDS(NMAList, paste0(loc,"/sim_results/obs/NMAList.rds"))

################################################################################
# time course NMA
#  Include data phase 3 for all agents  

timeNMAList = future_lapply(future.seed = T, X = seq_along(dataLst), FUN = function(i){
  datList = readRDS(paste0(loc,"/sim_data/tempListTime.RDS"))
  dat = dataLst[[i]]
  dat = subset(dat, study %in% c(2,3,4,6,7,9,10))
  Vlst= varLst[[i]]
  Vlst=Vlst[c(2,3,4,6,7,9,10)]
  V   = as.matrix(bdiag( Vlst )) 
  datList$Mfcor <- V
  datList$Y <- dat$y
  datList$se <- dat$se
  nonlin <- stan(data=datList, file = paste0(file.stan, "exp_Time_NMA.stan"), model_name = "exp.NMA.time", 
                 iter = 2000, cores = 2, chains = 2, control = list(max_treedepth = 20, adapt_delta = 0.999), seed=1234 )
  
  out <- spread_draws(nonlin,  b_emaxT[i],  ndraws = NULL) %>%
    filter( i %in% c(9,10,13)) %>%
    pivot_wider(names_from='i', 
                names_glue = "{i}_{.value}",
                values_from='b_emaxT') %>%
    rename( b_emaxA = '9_b_emaxT', b_emaxB='10_b_emaxT', b_emaxC='13_b_emaxT') %>%
    mutate(Estimate_AvsC = b_emaxA - b_emaxC, 
           Estimate_BvsC = b_emaxB - b_emaxC ) %>% 
    dplyr::select( .chain, .iteration, .draw,  Estimate_AvsC, Estimate_BvsC ) %>%
    pivot_longer(cols=c("Estimate_AvsC", "Estimate_BvsC"),  names_to='parameter', values_to='values') %>% 
    group_by(parameter) %>% 
    summarise(Estimate = mean(values), 
              Est.Error= sd(values),       
              Q2.5 = quantile(values, probs = 0.025),
              Q97.5= quantile(values, probs = 0.975) )  
  return(out)
  
})

if(export) saveRDS(timeNMAList, paste0(loc, "/sim_results/obs/timeNMAList.rds"))

################################################################################
# loop through the simulated datasets - Use Dose-Response, Time-Course NMA 
# scenario (1)
#

# Read data with priors and estimate
 dosetimeNMAList = future_lapply(future.seed = T, X = seq_along(dataLst), FUN = function(i){
   datList <- readRDS(paste0(loc,"/sim_data/tempList.RDS")) # helper to ensure correct data format, values overwritten below
   dat = dataLst[[i]]
   Vlst= varLst[[i]]
   V   = as.matrix(bdiag( Vlst )) 
   datList$Mfcor <- V
   datList$Y <- dat$y
   datList$se <- dat$se
   nonlin2 <- stan(data=datList, file = paste0(file.stan, "exp_NMA_emax_log.stan"), model_name = "exp.NMA.emax.log", 
                   iter = 2000, cores = 2, chains = 2, control = list(max_treedepth = 20, adapt_delta = 0.999), seed=1234 )
   
   est = summary(nonlin2, pars = c( "b_emaxA[1]","b_emaxB[1]","b_emaxC[1]","b_ed50A[1]","b_ed50B[1]","b_ed50C[1]","b_sA[1]","b_sB[1]","b_sC[1]","sd_1[2]","sd_1[1]","Cor_1[2,1]") )$summary

   out <- spread_draws(nonlin2, b_emaxA[1], b_emaxB[1], b_emaxC[1], b_ed50A[1], b_ed50B[1], b_ed50C[1], ndraws = NULL) %>% 
     rename( b_emaxA_Intercept=b_emaxA, b_emaxB_Intercept=b_emaxB, b_emaxC_Intercept=b_emaxC, b_ed50A_Intercept=b_ed50A, b_ed50B_Intercept=b_ed50B, b_ed50C_Intercept=b_ed50C ) 
   for(i in names(out)[4:9])  out[,i] = as.numeric(unlist(out[,i]))
   out <- out %>%
     mutate(Estimate_A    = b_emaxA_Intercept*1/(exp(b_ed50A_Intercept) + 1), 
            Estimate_B    = b_emaxB_Intercept*3/( exp(b_ed50B_Intercept) + 3),
            Estimate_C    = b_emaxC_Intercept*2/( exp(b_ed50C_Intercept) + 2),
            Estimate_AvsC = b_emaxA_Intercept*1/(exp(b_ed50A_Intercept) + 1) -  b_emaxC_Intercept*2/( exp(b_ed50C_Intercept) + 2), 
            Estimate_BvsC = b_emaxB_Intercept*3/( exp(b_ed50B_Intercept) + 3) - b_emaxC_Intercept*2/( exp(b_ed50C_Intercept) + 2) ) %>% 
     dplyr::select( .chain, .iteration, .draw,  Estimate_A, Estimate_B, Estimate_C, Estimate_AvsC, Estimate_BvsC ) %>%
     pivot_longer(cols=c("Estimate_A", "Estimate_B", "Estimate_C", "Estimate_AvsC", "Estimate_BvsC"),  names_to='parameter', values_to='values') %>% 
     group_by(parameter) %>% 
     summarise(Estimate = mean(values), 
               Est.Error= sd(values),       
               Q2.5 = quantile(values, probs = 0.025),
               Q97.5= quantile(values, probs = 0.975) )
   
   return( list(Est=est, releff=out) )
 
 })

if( export ) saveRDS(dosetimeNMAList,paste0(loc,"/sim_results/obs/dosetimeNMAList.rds"))

################################################################################
# loop through 100 simulated datasets - Use Dose-Response, Time-Course NMA 
# scenario (2): Exclude data after 28 week for phase 2 for agent C
#

dosetimeNMAList2 = future_lapply(future.seed = T, X = seq_along(dataLst), FUN = function(i){
   datList = readRDS(paste0(loc,"/sim_data/tempList2.RDS"))
   dat = subset(dataLst[[i]], !(study == 8 & week > 28) )
   Vlst= varLst[[i]]
   V8  = Vlst[[8]][c(1:7, 11:17, 21:27, 31:37, 41:47),c(1:7, 11:17, 21:27, 31:37, 41:47)]
   Vlst[[8]] = V8
   V   = as.matrix(bdiag( Vlst )) 
   datList$Mfcor <- V
   datList$Y <- dat$y
   datList$se <- dat$se
   nonlin2 <- stan(data=datList, file = paste0(file.stan, "exp_NMA_emax_log2.stan"), model_name = "exp.NMA.emax.log", 
                   iter = 2000, cores = 2, chains = 2, control = list(max_treedepth = 20, adapt_delta = 0.999), seed=1234 )
   
   est = summary(nonlin2, pars = c( "b_emaxA[1]","b_emaxB[1]","b_emaxC[1]","b_ed50A[1]","b_ed50B[1]","b_ed50C[1]","b_sA[1]","b_sB[1]","b_sC[1]","sd_1[2]","sd_1[1]","Cor_1[2,1]") )$summary
   
   out <- spread_draws(nonlin2, b_emaxA[1], b_emaxB[1], b_emaxC[1], b_ed50A[1], b_ed50B[1], b_ed50C[1], ndraws = NULL) %>% 
     rename( b_emaxA_Intercept=b_emaxA, b_emaxB_Intercept=b_emaxB, b_emaxC_Intercept=b_emaxC, b_ed50A_Intercept=b_ed50A, b_ed50B_Intercept=b_ed50B, b_ed50C_Intercept=b_ed50C ) 
   for(i in names(out)[4:9])  out[,i] = as.numeric(unlist(out[,i]))
   out <- out %>%
     mutate(Estimate_A    = b_emaxA_Intercept*1/(exp(b_ed50A_Intercept) + 1), 
            Estimate_B    = b_emaxB_Intercept*3/( exp(b_ed50B_Intercept) + 3),
            Estimate_C    = b_emaxC_Intercept*2/( exp(b_ed50C_Intercept) + 2),
            Estimate_AvsC = b_emaxA_Intercept*1/(exp(b_ed50A_Intercept) + 1) -  b_emaxC_Intercept*2/( exp(b_ed50C_Intercept) + 2), 
            Estimate_BvsC = b_emaxB_Intercept*3/( exp(b_ed50B_Intercept) + 3) - b_emaxC_Intercept*2/( exp(b_ed50C_Intercept) + 2) ) %>% 
     dplyr::select( .chain, .iteration, .draw,  Estimate_A, Estimate_B, Estimate_C, Estimate_AvsC, Estimate_BvsC ) %>%
     pivot_longer(cols=c("Estimate_A", "Estimate_B", "Estimate_C", "Estimate_AvsC", "Estimate_BvsC"),  names_to='parameter', values_to='values') %>% 
     group_by(parameter) %>% 
     summarise(Estimate = mean(values), 
               Est.Error= sd(values),       
               Q2.5 = quantile(values, probs = 0.025),
               Q97.5= quantile(values, probs = 0.975) )
   
   return( list(Est=est, releff=out) )
   
 })

if( export ) saveRDS(dosetimeNMAList2,paste0(loc,"/sim_results/obs/dosetimeNMAList2.rds"))


################################################################################
# loop through 100 simulated datasets - Use Dose-Response, Time-Course NMA 
# scenario (3): Exclude phase 3 for agent C, and data after 28 week for phase 2
#

dosetimeNMAList3 = future_lapply(future.seed = T, X = seq_along(dataLst), FUN = function(i){
  datList = readRDS(paste0(loc,"/sim_data/tempList3.RDS"))
  dat = subset(dataLst[[i]], !(study == 8 & week > 28 | study %in% c(9,10)) )
  Vlst= varLst[[i]]
  V8  = Vlst[[8]][c(1:7, 11:17, 21:27, 31:37, 41:47),c(1:7, 11:17, 21:27, 31:37, 41:47)]
  Vlst[[8]] = V8
  Vlst=Vlst[1:8]
  V   = as.matrix(bdiag( Vlst )) 
  datList$Mfcor <- V
  datList$Y <- dat$y
  datList$se <- dat$se
  nonlin2 <- stan(data=datList, file = paste0(file.stan, "exp_NMA_emax_log3.stan"), model_name = "exp.NMA.emax.log", 
                  iter = 2000, cores = 2, chains = 2, control = list(max_treedepth = 20, adapt_delta = 0.999), seed=1234 )
  
  est = summary(nonlin2, pars = c( "b_emaxA[1]","b_emaxB[1]","b_emaxC[1]","b_ed50A[1]","b_ed50B[1]","b_ed50C[1]","b_sA[1]","b_sB[1]","b_sC[1]","sd_1[2]","sd_1[1]","Cor_1[2,1]") )$summary
  
  out <- spread_draws(nonlin2, b_emaxA[1], b_emaxB[1], b_emaxC[1], b_ed50A[1], b_ed50B[1], b_ed50C[1], ndraws = NULL) %>% 
    rename( b_emaxA_Intercept=b_emaxA, b_emaxB_Intercept=b_emaxB, b_emaxC_Intercept=b_emaxC, b_ed50A_Intercept=b_ed50A, b_ed50B_Intercept=b_ed50B, b_ed50C_Intercept=b_ed50C ) 
  for(i in names(out)[4:9])  out[,i] = as.numeric(unlist(out[,i]))
  out <- out %>%
    mutate(Estimate_A    = b_emaxA_Intercept*1/(exp(b_ed50A_Intercept) + 1), 
           Estimate_B    = b_emaxB_Intercept*3/( exp(b_ed50B_Intercept) + 3),
           Estimate_C    = b_emaxC_Intercept*2/( exp(b_ed50C_Intercept) + 2),
           Estimate_AvsC = b_emaxA_Intercept*1/(exp(b_ed50A_Intercept) + 1) -  b_emaxC_Intercept*2/( exp(b_ed50C_Intercept) + 2), 
           Estimate_BvsC = b_emaxB_Intercept*3/( exp(b_ed50B_Intercept) + 3) - b_emaxC_Intercept*2/( exp(b_ed50C_Intercept) + 2) ) %>% 
    dplyr::select( .chain, .iteration, .draw,  Estimate_A, Estimate_B, Estimate_C, Estimate_AvsC, Estimate_BvsC ) %>%
    pivot_longer(cols=c("Estimate_A", "Estimate_B", "Estimate_C", "Estimate_AvsC", "Estimate_BvsC"),  names_to='parameter', values_to='values') %>% 
    group_by(parameter) %>% 
    summarise(Estimate = mean(values), 
              Est.Error= sd(values),       
              Q2.5 = quantile(values, probs = 0.025),
              Q97.5= quantile(values, probs = 0.975) )
  
  return( list(Est=est, releff=out) )
  
})

if( export ) saveRDS(dosetimeNMAList2,paste0(loc,"/sim_results/obs/dosetimeNMAList3.rds"))