################################################################################
# Description: Dose-response, time-course Network meta-analysis, on obesity data set, using observed correlations   
# - results saved in /obesity_results/DoseTimeRCOV/
#
################################################################################
rm(list = ls())

library(tidyverse)
library(gridExtra)
library(stringr)
library(brms)
library(Matrix)
library(tidybayes)    
library(ggplot2)
library(rstan)
options(mc.cores = parallel::detectCores())

export = F 

#set location
scer = T 
loc  = ifelse(scer, "~/hdrive/mbnma", "H:/mmbna")
stan.loc = paste0(loc, "/obesity_stan/")


# utility functions
# Setup covariance matrix for study effects
mat.fun = function(i){
  tmp = matrix(NA, i, i)
  diag(tmp) = 1
  tmp[lower.tri(tmp)] = 1/2
  tmp[upper.tri(tmp)] = 1/2
  tmp
  
}

# load data
dat = readRDS(paste0(loc, "/obesity_data/datDoseTimeRCOV.rds"))
v_mat = readRDS(paste0(loc,"/obesity_data/v_mat.rds")) # variance matrices


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

# Derive dose effects, i.e. regression variables
dat$dLira = dat$d1*log(dat$mdose + 1)  
dat$dSema = dat$d2*log(dat$mdose + 1)

################################################################################
# Estimate parameters
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

################################################################################
# Estimate parameters - univariate likelihood
# in order to pass the data structure

exp.NMA.log.1 <- brm(bf(mean | se(se) ~ emaxT*( 1 - exp( -exp(kT)*week ) ),
                        emaxT ~ 0 + author + dLira + dSema + (1 | author | gr(stcom, cov = st_mat)), 
                        kT    ~ 0 + author + dLira + dSema + (1 | author | gr(stcom, cov = st_mat)),
                        nl = TRUE),
                     data = dat,
                     data2 = list(st_mat=st_mat),
                     prior = exp.prior,
                     control = list(adapt_delta = 0.999, max_treedepth = 20), 
                     iter = 2000, cores = 2, chains = 2, 
                     save_pars = save_pars(all = TRUE), seed = 1234)

################################################################################
# Estimate parameters - multivariate likelihood:
# Dose-response: Emax parameter ~ log-linear and Kt ~ log-linear model
# 

datlist  = standata(exp.NMA.log.1) 
datlist$Mfcor = v_mat


exp.NMA.log.2 <- stan(data=datlist, file = paste0(stan.loc, "exp_NMA_log.stan"), model_name = "exp.NMA.log.log", 
                      iter = 2000, cores = 2, chains = 2, refresh=100, 
                      control = list(max_treedepth = 20, adapt_delta = 0.999), seed=1234 )

summary(exp.NMA.log.2, pars=c("b_emaxT", "b_kT", "sd_1", "Cor_1[2,1]"), probs=c(0.025, 0.975))$summary
names(exp.NMA.log.2)
saveRDS(exp.NMA.log.2, paste0(loc,"/obesity_results/DoseTimeRCOV/exp.NMA.log.log.rds"))

################################################################################
# Estimate parameters - univariate likelihood:
# Dose-response: Emax parameter ~ Emax model and Kt ~ log-linear model
# to get the right data structure
#

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

exp.NMA.emax.1 <- brm(bf(mean | se(se) ~ (MUemaxT + (emaxLira*d1 + emaxSema*d2 )*mdose/(exp( ed50Lira*d1 + ed50Sema*d2  ) + mdose))*( 1 - exp( -exp(kT)*week ) ),
                           MUemaxT ~ 0 + author + (1 | author | gr(stcom, cov = st_mat)), 
                           emaxLira + emaxSema +  ed50Lira + ed50Sema  ~ 1,
                           kT ~ 0 + author + dLira + dSema + (1 | author | gr(stcom, cov = st_mat)), 
                           nl=TRUE),
                        data = dat,  
                        data2 = list(st_mat=st_mat),
                        prior = exp.prior,
                        control = list(adapt_delta = 0.999, max_treedepth = 20),
                        iter = 200, cores = 2, chains = 2,
                        save_pars = save_pars(all = TRUE), seed = 1234)

################################################################################
# Estimate parameters - multivariate likelihood:

datlist  = standata(exp.NMA.emax.1) 
datlist$Mfcor = v_mat

exp.NMA.log.2.1 <- stan(data=datlist, file = paste0(stan.loc, "exp_NMA_emax.stan"), model_name = "exp.NMA.emax.log", 
                        iter = 5000, cores = 1, chains = 2, refresh=25, 
                        control = list(max_treedepth = 20, adapt_delta = 0.999), seed=5678 )

saveRDS(exp.NMA.log.2.1, paste0(loc, "/obesity_results/DoseTimeRCOV/exp.NMA.emax.log.rds"))







