################################################################################
# Description: Estimation of Standard NMA, time course, and dose-response time-course
#   - loops over 100 replicas, saved in ./sim_data
#   - dose-response time-course assumes multivariate likelihood, with AR-1 correlations
#   - under 3 different data scenarios
#   - estimates are saved in ./sim_results/AR1/
#
################################################################################
#
rm(list = ls()) # clear memory
options(mc.cores = parallel::detectCores()) # setup for parallel 

# load required packages
library(tidyverse)
library(MASS)
library(brms)
library(Matrix)
library(haven)
library(tidybayes)
library(future.apply)

plan(multisession)

# set location
scer = F
loc  = ifelse(scer, "~/hdrive/mbnma", "H:/mmbna")

# set seed for reproducibility
set.seed(1234)
export = T # set to T if datasets are to be saved

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

if(export) saveRDS(NMAList, paste0(loc,"/sim_results/AR1/NMAList.rds"))

################################################################################
# loop through 100 simulated datasets - Use Time-Course NMA 

# Specify prior for random effects:  
exp.prior =
  c(
    prior(normal(-2.5, 10), coef="study2", nlpar = "emaxT"),
    prior(normal(-2.5, 10), coef="study3", nlpar = "emaxT"),
    prior(normal(-2.5, 10), coef="study4", nlpar = "emaxT"),
    prior(normal(-2.5, 10), coef="study6", nlpar = "emaxT"),
    prior(normal(-2.5, 10), coef="study7", nlpar = "emaxT"),
    prior(normal(-2.5, 10), coef="study9", nlpar = "emaxT"),
    prior(normal(-2.5, 10), coef="study10",nlpar = "emaxT"),
    prior(normal(0, 100),  coef="d1", nlpar = "emaxT"),
    prior(normal(0, 100),  coef="d2", nlpar = "emaxT"),
    prior(normal(0, 100),  coef="d3",nlpar = "emaxT"),
    prior(normal(0, 100),  coef="d4", nlpar = "emaxT"),
    prior(normal(0, 100),  coef="d5", nlpar = "emaxT"),
    prior(normal(0, 100),  coef="d6",  nlpar = "emaxT"),
    # Priors on kT
    prior(normal(-2.5, 10), coef="study2", nlpar = "kT"),
    prior(normal(-2.5, 10), coef="study3", nlpar = "kT"),
    prior(normal(-2.5, 10), coef="study4", nlpar = "kT"),
    prior(normal(-2.5, 10), coef="study6", nlpar = "kT"),
    prior(normal(-2.5, 10), coef="study7", nlpar = "kT"),
    prior(normal(-2.5, 10), coef="study9", nlpar = "kT"),
    prior(normal(-2.5, 10), coef="study10",nlpar = "kT"),
    prior(normal(0, 100),  coef="d1", nlpar = "kT"),
    prior(normal(0, 100),  coef="d2", nlpar = "kT"),
    prior(normal(0, 100),  coef="d3",nlpar = "kT"),
    prior(normal(0, 100),  coef="d4", nlpar = "kT"),
    prior(normal(0, 100),  coef="d5", nlpar = "kT"),
    prior(normal(0, 100),  coef="d6",  nlpar = "kT"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "emaxT"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "kT")
  )


timeNMAList = future_lapply(future.seed = T, X = seq_along(dataLst), FUN = function(i){
  
  dat = dataLst[[i]]
  dat = subset(dat, study %in% c(2,3,4,6,7,9,10))
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
  
  # Setup design matrix for treatment effects
  d = unique(dat$treatment[dat$treatment!="soc"]); dname = paste0("d",1:length(d))
  for(i in 1:length(d)) dat[,dname[i]] = ifelse(dat$treatment == d[i], 1, 0) 
  
  dat$study = factor(dat$study)
  
  p1 = as.formula(paste("emaxT ~ 0 + study + ", paste(dname, collapse= "+"), "+ (1 | study | gr(stcom, cov = st_mat))" ))
  p2 = as.formula(paste("kT ~ 0 + study + ", paste(dname, collapse= "+"), "+ (1 | study | gr(stcom, cov = st_mat))" ))
  
  exp.NMA.emax <- brm(bf(y | se(se) ~ emaxT*( 1 - exp( -exp(kT)*week ) ),
                         p1, p2,
                         nl=TRUE),
                      data = dat,  
                      data2 = list(st_mat=st_mat),
                      prior = exp.prior,
                      control = list(adapt_delta = 0.999, max_treedepth = 25),
                      iter = 2000, chains = 2,
                      backend = "cmdstanr",
                      threads = threading(2) )
  
  out <- spread_draws(exp.NMA.emax, b_emaxT_d2, b_emaxT_d3, b_emaxT_d6, ndraws = NULL) %>% 
    mutate(Estimate_AvsC = b_emaxT_d2 - b_emaxT_d6, 
           Estimate_BvsC = b_emaxT_d3 - b_emaxT_d6 ) %>% 
    dplyr::select( .chain, .iteration, .draw,  Estimate_AvsC, Estimate_BvsC ) %>%
    pivot_longer(cols=c("Estimate_AvsC", "Estimate_BvsC"),  names_to='parameter', values_to='values') %>% 
    group_by(parameter) %>% 
    summarise(Estimate = mean(values), 
              Est.Error= sd(values),       
              Q2.5 = quantile(values, probs = 0.025),
              Q97.5= quantile(values, probs = 0.975) )  
  
  return(out)  
  
} )

if(export) saveRDS(timeNMAList, paste0(loc, "/sim_results/AR1/timeNMAList.rds"))


################################################################################
# loop through the simulated datasets - Use Dose-Response, Time-Course NMA 
# scenario (1)
#

# Specify prior   
exp.prior =
  c(
    prior(normal(-2.5, 10), coef="study1", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study2", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study3", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study4", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study5", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study6", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study7", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study8", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study9", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study10",nlpar = "MUemaxT"),
    prior(normal(-10, 10),  nlpar = "emaxA"),
    prior(normal(-15, 10),  nlpar = "emaxB"),
    prior(normal(-20, 10),  nlpar = "emaxC"),
    prior(normal(-2.0, 10),  nlpar = "ed50A"),
    prior(normal(-1.5, 10),  nlpar = "ed50B"),
    prior(normal(-0.75,10),  nlpar = "ed50C"),
    # Priors on kT
    prior(normal(-2.5, 10), coef="study1", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study2", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study3", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study4", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study5", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study6", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study7", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study8", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study9", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study10",nlpar = "MUkT"),
    prior(normal(0, 10),  nlpar="sA"),
    prior(normal(0, 10),  nlpar="sB"),
    prior(normal(0, 10),  nlpar="sC"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "MUemaxT"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "MUkT")
  )

# Analysis model = True model
dosetimeNMAList = future_lapply(future.seed = T, X = seq_along(dataLst), FUN = function(i){
  dat = dataLst[[i]]
  # Setup covariance matrix for study effects                       
  dat$stcom = rep(1:length(unique(paste0(dat$study,dat$arm))), times=table(factor(paste0(dat$study,dat$arm), levels=unique(paste0(dat$study,dat$arm)))) )
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
  
  # Estimate parameters:
  exp.NMA.emax.log <- brm(bf(y | se(se, sigma=TRUE) ~ (MUemaxT + (emaxA*dA + emaxB*dB + emaxC*dC)*dose/(exp( ed50A*dA + ed50B*dB + ed50C*dC ) + dose))*( 1 - exp( -exp( MUkT + sA*dA*log(dose + 1) + sB*dB*log(dose + 1) + sC*dC*log(dose + 1) )*week ) ),
                             MUemaxT ~ 0 + study + (1 | study | gr(stcom, cov = st_mat)), 
                             emaxA + emaxB + emaxC + ed50A + ed50B + ed50C ~ 1,
                             MUkT ~ 0 + study + (1 | study | gr(stcom, cov = st_mat)), 
                             sA + sB + sC ~ 1,
                             nl=TRUE, 
                             autocor = cor_ar(~ week | stcom, cov = TRUE) ),
                          data = dat,
                          data2 = list(st_mat=st_mat),
                          prior = exp.prior,
                          control = list(adapt_delta = 0.999, max_treedepth = 20), 
                          iter = 200, chains = 2,
                          backend = "cmdstanr",
                          seed = 1234)


  out <- spread_draws(exp.NMA.emax.log, b_emaxA_Intercept, b_emaxB_Intercept, b_emaxC_Intercept, b_ed50A_Intercept, b_ed50B_Intercept, b_ed50C_Intercept, ndraws = NULL) %>%   
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
  
  
  return( list(Est=as.data.frame( cbind(fixef(exp.NMA.emax.log), rhat=rhat(exp.NMA.emax.log)[1:29]) ), releff=out) )
  
} )

if(export) saveRDS(dosetimeNMAList, paste0(loc, "/sim_results/AR1/dosetimeNMAList.rds"))

################################################################################
# loop through 100 simulated datasets - Use Dose-Response, Time-Course NMA 
# scenario (2): Exclude data after 28 week for phase 2 for agent C
#

# Analysis model = True model
dosetimeNMAList2 = future_lapply(future.seed = T, X = seq_along(dataLst), FUN = function(i){
  dat = subset(dataLst[[i]], !(study == 8 & week > 28) )
  # Setup covariance matrix for study effects
  dat$stcom = rep(1:length(unique(paste0(dat$study,dat$arm))), times=table(factor(paste0(dat$study,dat$arm), levels=unique(paste0(dat$study,dat$arm)))) )
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
  
  # Estimate parameters:
  exp.NMA.emax.log <- brm(bf(y | se(se, sigma=TRUE) ~ (MUemaxT + (emaxA*dA + emaxB*dB + emaxC*dC)*dose/(exp( ed50A*dA + ed50B*dB + ed50C*dC ) + dose))*( 1 - exp( -exp( MUkT + sA*dA*log(dose + 1) + sB*dB*log(dose + 1) + sC*dC*log(dose + 1) )*week ) ),
                             MUemaxT ~ 0 + study + (1 | study | gr(stcom, cov = st_mat)), 
                             emaxA + emaxB + emaxC + ed50A + ed50B + ed50C ~ 1,
                             MUkT ~ 0 + study + (1 | study | gr(stcom, cov = st_mat)), 
                             sA + sB + sC ~ 1,
                             nl=TRUE, 
                             autocor = cor_ar(~ week | stcom, cov = TRUE) ),
                          data = dat,
                          data2 = list(st_mat=st_mat),
                          prior = exp.prior,
                          control = list(adapt_delta = 0.999, max_treedepth = 20), 
                          iter = 2000, chains = 2,
                          backend = "cmdstanr",
                          seed = 1234)
  
  
  out <- spread_draws(exp.NMA.emax.log, b_emaxA_Intercept, b_emaxB_Intercept, b_emaxC_Intercept, b_ed50A_Intercept, b_ed50B_Intercept, b_ed50C_Intercept, ndraws = NULL) %>%   
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
  
  
  return( list(Est=as.data.frame( cbind(fixef(exp.NMA.emax.log), rhat=rhat(exp.NMA.emax.log)[1:29]) ), releff=out) )
  
} )

if(export) saveRDS(dosetimeNMAList2, paste0(loc, "/sim_results/AR1/dosetimeNMAList2.rds"))

################################################################################
# loop through 100 simulated datasets - Use Dose-Response, Time-Course NMA 
# scenario (3): Exclude phase 3 for agent C, and data after 28 week for phase 2
#

# Specify prior   
exp.prior =
  c(
    prior(normal(-2.5, 10), coef="study1", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study2", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study3", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study4", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study5", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study6", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study7", nlpar = "MUemaxT"),
    prior(normal(-2.5, 10), coef="study8", nlpar = "MUemaxT"),
#    prior(normal(-2.5, 10), coef="study9", nlpar = "MUemaxT"),
#    prior(normal(-2.5, 10), coef="study10",nlpar = "MUemaxT"),
    prior(normal(-10, 10),  nlpar = "emaxA"),
    prior(normal(-15, 10),  nlpar = "emaxB"),
    prior(normal(-20, 10),  nlpar = "emaxC"),
    prior(normal(-2.0, 10),  nlpar = "ed50A"),
    prior(normal(-1.5, 10),  nlpar = "ed50B"),
    prior(normal(-0.75,10),  nlpar = "ed50C"),
    # Priors on kT
    prior(normal(-2.5, 10), coef="study1", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study2", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study3", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study4", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study5", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study6", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study7", nlpar = "MUkT"),
    prior(normal(-2.5, 10), coef="study8", nlpar = "MUkT"),
#    prior(normal(-2.5, 10), coef="study9", nlpar = "MUkT"),
#    prior(normal(-2.5, 10), coef="study10",nlpar = "MUkT"),
    prior(normal(0, 10),  nlpar="sA"),
    prior(normal(0, 10),  nlpar="sB"),
    prior(normal(0, 10),  nlpar="sC"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "MUemaxT"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "MUkT")
  )


dosetimeNMAList3 = future_lapply(future.seed = T, X = seq_along(dataLst), FUN = function(i){
  dat = subset(dataLst[[i]], !(study == 8 & week > 28 | study %in% c(9,10)) )
  # Setup covariance matrix for study effects
  dat$stcom = rep(1:length(unique(paste0(dat$study,dat$arm))), times=table(factor(paste0(dat$study,dat$arm), levels=unique(paste0(dat$study,dat$arm)))) )
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
  
  # Estimate parameters:
  exp.NMA.emax.log <- brm(bf(y | se(se, sigma=TRUE) ~ (MUemaxT + (emaxA*dA + emaxB*dB + emaxC*dC)*dose/(exp( ed50A*dA + ed50B*dB + ed50C*dC ) + dose))*( 1 - exp( -exp( MUkT + sA*dA*log(dose + 1) + sB*dB*log(dose + 1) + sC*dC*log(dose + 1) )*week ) ),
                             MUemaxT ~ 0 + study + (1 | study | gr(stcom, cov = st_mat)), 
                             emaxA + emaxB + emaxC + ed50A + ed50B + ed50C ~ 1,
                             MUkT ~ 0 + study + (1 | study | gr(stcom, cov = st_mat)), 
                             sA + sB + sC ~ 1,
                             nl=TRUE, 
                             autocor = cor_ar(~ week | stcom, cov = TRUE) ),
                          data = dat,
                          data2 = list(st_mat=st_mat),
                          prior = exp.prior,
                          control = list(adapt_delta = 0.999, max_treedepth = 20), 
                          iter = 2000, chains = 2,
                          backend = "cmdstanr",
                          seed = 1234)
  
  
  out <- spread_draws(exp.NMA.emax.log, b_emaxA_Intercept, b_emaxB_Intercept, b_emaxC_Intercept, b_ed50A_Intercept, b_ed50B_Intercept, b_ed50C_Intercept, ndraws = NULL) %>%   
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
  
  
  return( list(Est=as.data.frame( cbind(fixef(exp.NMA.emax.log), rhat=rhat(exp.NMA.emax.log)[1:25]) ), releff=out) )
  
} )

if(export) saveRDS(dosetimeNMAList3, paste0(loc, "/sim_results/AR1/dosetimeNMAList3.rds"))