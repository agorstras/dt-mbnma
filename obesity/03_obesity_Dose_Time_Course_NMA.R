################################################################################
# Description: Dose response, time course Network meta-analysis, on obesity data set, independent and AR-1  
# - estimation results saved in: /obesity_results/DoseTime/
#
################################################################################
rm(list = ls())

library(tidyverse)
library(gridExtra)
library(stringr)
library(brms)
library(cmdstanr)
library(Matrix)
library(tidybayes)    
library(ggplot2)
options(mc.cores = parallel::detectCores())

export = F

#set location
loc = "obesity"

set.seed(1234)

# Setup covariance matrix for random study effects
mat.fun = function(i){
  tmp = matrix(NA, i, i)
  diag(tmp) = 1
  tmp[lower.tri(tmp)] = 1/2
  tmp[upper.tri(tmp)] = 1/2
  tmp
  
}

# load data
dat = readRDS(paste0(loc, "/obesity_data/datDoseTime.rds"))

dat2 = dat[!duplicated(paste0(dat$author, dat$trt)), ]

# Setup design matrix for compound effects
d = unique(dat$compound[dat$trt!="Placebo"]); dname = paste0("d",1:length(d))
for(i in 1:length(d)) dat[,dname[i]] = ifelse(dat$compound == d[i], 1, 0) 

# Setup covariance matrix, handling trials with multiple arms
lst = split(dat2, dat2$author)
lst = lapply(1:length(lst), function(j){
  tmp = lst[[j]]
  tmp = tmp[tmp$trt != "Placebo", ]
  as.matrix(bdiag(matrix(1E-20,1,1), mat.fun(length(tmp$author)))) 
} )

st_mat = as.matrix(bdiag(lst))
rownames(st_mat) <- colnames(st_mat) <- 1:nrow(dat2) 

#################################################################################
# Time-course model: Exponential 
# Dose-response: Emax parameter ~ log-linear model and Kt ~ log-linear model

# Derive dose effects, i.e. regression variables
dat$dLira = dat$d1*log(dat$mdose + 1)  
dat$dSema = dat$d2*log(dat$mdose + 1)

exp.prior =
  c(
    prior(normal(-2.5, 100), coef="authorAstrup",     nlpar = "emaxT"),
    prior(normal(-2.5, 100), coef="authorGarvey",     nlpar = "emaxT"),
    prior(normal(-2.5, 100), coef="authorKadowaki",   nlpar = "emaxT"),
    prior(normal(-2.5, 100), coef="authorONeil",      nlpar = "emaxT"),
    prior(normal(-2.5, 100), coef="authorPiMSunyer",  nlpar = "emaxT"),
    prior(normal(-2.5, 100), coef="authorRubino",     nlpar = "emaxT"),
    prior(normal(-2.5, 100), coef="authorWilding",    nlpar = "emaxT"),
    prior(normal(0, 100),  coef="dLira", nlpar = "emaxT"),
    prior(normal(0, 100),  coef="dSema", nlpar = "emaxT"),
    # Priors on kT
    prior(normal(-2.5, 100), coef="authorAstrup",     nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorGarvey",     nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorKadowaki",   nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorONeil",      nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorPiMSunyer",  nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorRubino",     nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorWilding",    nlpar = "kT"),
    prior(normal(0, 100),  coef="dLira", nlpar = "kT"),
    prior(normal(0, 100),  coef="dSema", nlpar = "kT"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "emaxT"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "kT")
  )

# Estimate parameters:
exp.NMA.log.1 <- brm(bf(mean | se(se) ~ emaxT*( 1 - exp( -exp(kT)*week ) ),
                          emaxT ~ 0 + author + dLira + dSema + (1 | author | gr(stcom, cov = st_mat)), 
                          kT    ~ 0 + author + dLira + dSema + (1 | author | gr(stcom, cov = st_mat)),
                          nl = TRUE), 
                          data = dat,
                          data2 = list(st_mat=st_mat),
                          prior = exp.prior,
                          control = list(adapt_delta = 0.999, max_treedepth = 20), 
                          iter = 200, cores = 4,
                          backend = "cmdstanr",
                          threads = threading(6),
                          save_pars = save_pars(all = TRUE), seed = 1234)

exp.NMA.log.1

if(export) saveRDS(exp.NMA.log.1, paste0(loc, "/obesity_results/DoseTime/exp.NMA.log.1.rds"))

# Estimate parameters, AR-1 correlations

exp.NMA.log.1.ar <- brm(bf(mean | se(se, sigma=TRUE) ~ emaxT*( 1 - exp( -exp(kT)*week ) ),
                           emaxT ~ 0 + author + dLira + dSema + (1 | author | gr(stcom, cov = st_mat)), 
                           kT    ~ 0 + author + dLira + dSema + (1 | author | gr(stcom, cov = st_mat)),
                           nl = TRUE, autocor = cor_ar(~ week | stcom, cov = TRUE) ),
                     data = dat,
                     data2 = list(st_mat=st_mat),
                     prior = exp.prior,
                     control = list(adapt_delta = 0.999, max_treedepth = 20), 
                     iter = 5000, cores = 4,
                     seed = 1234)

exp.NMA.log.1.ar 

if(export) saveRDS(exp.NMA.log.1.ar, paste0(loc, "/obesity_results/DoseTime/exp.NMA.log.1.ar.rds") )

################################################################################
# Time-course model: Exponential 
# Dose-response: Emax parameter ~ Emax model and Kt ~ emax model

exp.prior =
  c(
    prior(normal(-2.5, 100), coef="authorAstrup",     nlpar = "MUemaxT"),
    prior(normal(-2.5, 100), coef="authorGarvey",     nlpar = "MUemaxT"),
    prior(normal(-2.5, 100), coef="authorKadowaki",   nlpar = "MUemaxT"),
    prior(normal(-2.5, 100), coef="authorONeil",      nlpar = "MUemaxT"),
    prior(normal(-2.5, 100), coef="authorPiMSunyer",  nlpar = "MUemaxT"),
    prior(normal(-2.5, 100), coef="authorRubino",     nlpar = "MUemaxT"),
    prior(normal(-2.5, 100), coef="authorWilding",    nlpar = "MUemaxT"),
    prior(normal(-15, 100),  nlpar = "emaxLira"),
    prior(normal(-25, 100),  nlpar = "emaxSema"),
    prior(normal(1.5, 100),  nlpar = "ed50Lira"),
    prior(normal(0.75, 100),  nlpar = "ed50Sema"),
    # Priors on kT
    prior(normal(-2.5, 100), coef="authorAstrup",     nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorGarvey",     nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorKadowaki",   nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorONeil",      nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorPiMSunyer",  nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorRubino",     nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorWilding",    nlpar = "kT"),
    prior(normal(-1, 100),  nlpar = "kemaxLira"),
    prior(normal(-2, 100),  nlpar = "kemaxSema"),
    prior(normal(1.5, 100),  nlpar = "ked50Lira"),
    prior(normal(0.75,100),  nlpar = "ked50Sema"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "MUemaxT"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "kT")
  )

exp.NMA.emax.1 <- brm(bf(mean | se(se) ~ ( MUemaxT + (emaxLira*d1 + emaxSema*d2 )*mdose/(exp( ed50Lira*d1 + ed50Sema*d2 ) + mdose) )*( 1 - exp( -exp(   MUkT + (kemaxLira*d1 + kemaxSema*d2 )*mdose/(exp( ked50Lira*d1 + ked50Sema*d2) + mdose) + deltakT   )*week ) ),
                         MUemaxT ~ 0 + author + (1 | author | gr(stcom, cov = st_mat)), 
                         emaxLira + emaxSema + ed50Lira + ed50Sema ~ 1,
                         MUkT ~ 0 + author + (1 | author | gr(stcom, cov = st_mat)),
                         kemaxLira + kemaxSema + ked50Lira + ked50Sema ~ 1,
                         deltakT ~ 0 , 
                         nl=TRUE),
                      data = dat,  
                      data2 = list(st_mat=st_mat),
                      control = list(adapt_delta = 0.999, max_treedepth = 20),
                      iter = 5000, cores = 4, chains = 4,
                      backend = "cmdstanr",
                      threads = threading(6),
                      save_pars = save_pars(all = TRUE), seed = 1234)

if(export) saveRDS(exp.NMA.emax.1, paste0(loc, "/obesity_results/DoseTime/exp.NMA.emax.1.rds") )


# Dose-response: Emax parameter ~ Emax model and Kt ~ log-linear model

exp.prior =
  c(
    prior(normal(-2.5, 100), coef="authorAstrup",     nlpar = "MUemaxT"),
    prior(normal(-2.5, 100), coef="authorGarvey",     nlpar = "MUemaxT"),
    prior(normal(-2.5, 100), coef="authorKadowaki",   nlpar = "MUemaxT"),
    prior(normal(-2.5, 100), coef="authorONeil",      nlpar = "MUemaxT"),
    prior(normal(-2.5, 100), coef="authorPiMSunyer",  nlpar = "MUemaxT"),
    prior(normal(-2.5, 100), coef="authorRubino",     nlpar = "MUemaxT"),
    prior(normal(-2.5, 100), coef="authorWilding",    nlpar = "MUemaxT"),
    prior(normal(-15, 100),  nlpar = "emaxLira"),
    prior(normal(-25, 100),  nlpar = "emaxSema"),
    prior(normal(1.5, 100),  nlpar = "ed50Lira"),
    prior(normal(0.75, 100),  nlpar = "ed50Sema"),
    # Priors on kT
    prior(normal(-2.5, 100), coef="authorAstrup",     nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorGarvey",     nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorKadowaki",   nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorONeil",      nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorPiMSunyer",  nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorRubino",     nlpar = "kT"),
    prior(normal(-2.5, 100), coef="authorWilding",    nlpar = "kT"),
    prior(normal(0, 100),  coef="dLira", nlpar = "kT"),
    prior(normal(0, 100),  coef="dSema", nlpar = "kT"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "MUemaxT"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "kT")
  )

exp.NMA.emax.log <- brm(bf(mean | se(se) ~ (MUemaxT + (emaxLira*d1 + emaxSema*d2 )*mdose/(exp( ed50Lira*d1 + ed50Sema*d2  ) + mdose))*( 1 - exp( -exp(kT)*week ) ),
                             MUemaxT ~ 0 + author + (1 | author | gr(stcom, cov = st_mat)), 
                             emaxLira + emaxSema +  ed50Lira + ed50Sema  ~ 1,
                             kT ~ 0 + author + dLira + dSema + (1 | author | gr(stcom, cov = st_mat)), 
                             nl=TRUE),
                          data = dat,  
                          data2 = list(st_mat=st_mat),
                          prior = exp.prior,
                          control = list(adapt_delta = 0.999, max_treedepth = 20),
                          iter = 5000, cores = 4, chains = 4,
                          backend = "cmdstanr",
                          threads = threading(6),
                          save_pars = save_pars(all = TRUE), seed = 1234)

exp.NMA.emax.log

if(export) saveRDS(exp.NMA.emax.log, paste0(loc, "/obesity_results/DoseTime/exp.NMA.emax.log.rds") )

# AR1 structure fitted
exp.NMA.emax.log.ar <- brm(bf(mean | se(se, sigma=TRUE) ~ (MUemaxT + (emaxLira*d1 + emaxSema*d2 )*mdose/(exp( ed50Lira*d1 + ed50Sema*d2  ) + mdose))*( 1 - exp( -exp(kT)*week ) ),
                           MUemaxT ~ 0 + author + (1 | author | gr(stcom, cov = st_mat)), 
                           emaxLira + emaxSema +  ed50Lira + ed50Sema  ~ 1,
                           kT ~ 0 + author + dLira + dSema + (1 | author | gr(stcom, cov = st_mat)), 
                           nl=TRUE, autocor = cor_ar(~ week | stcom, p=1, cov = TRUE)),
                        data = dat,  
                        data2 = list(st_mat=st_mat),
                        prior = exp.prior,
                        control = list(adapt_delta = 0.999, max_treedepth = 20),
                        iter = 2500, cores = 4, chains = 4,
                        save_pars = save_pars(all = TRUE), seed = 12345)
exp.NMA.emax.log.ar
if(export) saveRDS(exp.NMA.emax.log.ar, paste0(loc, "/obesity_results/DoseTime/exp.NMA.emax.log.ar.rds") )


#################################################################################
# Time-course model: Exponential 
# Dose-response model: Linear

# Derive dose effects, i.e. regression variables
dat$dLira = dat$d1*dat$mdose    
dat$dSema = dat$d2*dat$mdose  

exp.prior =
  c(
    prior(normal(-2.5, 10), coef="authorAstrup",     nlpar = "emaxT"),
    prior(normal(-2.5, 10), coef="authorGarvey",     nlpar = "emaxT"),
    prior(normal(-2.5, 10), coef="authorKadowaki",   nlpar = "emaxT"),
    prior(normal(-2.5, 10), coef="authorONeil",      nlpar = "emaxT"),
    prior(normal(-2.5, 10), coef="authorPiMSunyer",  nlpar = "emaxT"),
    prior(normal(-2.5, 10), coef="authorRubino",     nlpar = "emaxT"),
    prior(normal(-2.5, 10), coef="authorWilding",    nlpar = "emaxT"),
    prior(normal(0, 100),  coef="dLira", nlpar = "emaxT"),
    prior(normal(0, 100),  coef="dSema", nlpar = "emaxT"),
    # Priors on kT
    prior(normal(-2.5, 10), coef="authorAstrup",     nlpar = "kT"),
    prior(normal(-2.5, 10), coef="authorGarvey",     nlpar = "kT"),
    prior(normal(-2.5, 10), coef="authorKadowaki",   nlpar = "kT"),
    prior(normal(-2.5, 10), coef="authorONeil",      nlpar = "kT"),
    prior(normal(-2.5, 10), coef="authorPiMSunyer",  nlpar = "kT"),
    prior(normal(-2.5, 10), coef="authorRubino",     nlpar = "kT"),
    prior(normal(-2.5, 10), coef="authorWilding",    nlpar = "kT"),
    prior(normal(0, 100),  coef="dLira", nlpar = "kT"),
    prior(normal(0, 100),  coef="dSema", nlpar = "kT"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "emaxT"),
    prior(normal(0,5), lb=0, class = "sd", nlpar = "kT")
  )

exp.NMA.lin.1 <- brm(bf(mean | se(se) ~ emaxT*( 1 - exp( -exp(kT)*week ) ),
                           emaxT ~ 0 + author + dLira + dSema + (1 | author | gr(stcom, cov = st_mat)), 
                           kT    ~ 0 + author + dLira + dSema + (1 | author | gr(stcom, cov = st_mat)),
                           nl = TRUE),
                        data = dat,
                        data2 = list(st_mat=st_mat),
                        prior = exp.prior,
                        control = list(adapt_delta = 0.9999, max_treedepth = 20),
                        iter = 5000, cores = 4, chains = 4,
                        save_pars = save_pars(all = TRUE), seed = 1234)

exp.NMA.lin.1

if(export) saveRDS(exp.NMA.lin.1, paste0(loc, "/obesity_results/DoseTime/exp.NMA.lin.1.rds") )

# Add AR-1 process
exp.NMA.lin.1.ar <- brm(bf(mean | se(se, sigma=TRUE) ~ emaxT*( 1 - exp( -exp(kT)*week ) ),
                        emaxT ~ 0 + author + dLira + dSema + (1 | author | gr(stcom, cov = st_mat)), 
                        kT    ~ 0 + author + dLira + dSema + (1 | author | gr(stcom, cov = st_mat)),
                        nl = TRUE, autocor = cor_ar(~ week | stcom, p=1, cov = TRUE) ),
                     data = dat,
                     data2 = list(st_mat=st_mat),
                     prior = exp.prior,
                     control = list(adapt_delta = 0.9999, max_treedepth = 20),
                     iter = 5000, cores = 4, chains = 4,
                     save_pars = save_pars(all = TRUE), seed = 1234)
exp.NMA.lin.1.ar
if(export) saveRDS(exp.NMA.lin.1.ar, paste0(loc, "/obesity_results/DoseTime/exp.NMA.lin.1.ar.rds") )
