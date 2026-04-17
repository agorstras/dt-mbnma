################################################################################
# Description: Time course network meta-analysis, on obesity data set   
# - Figure 2 (main manuscript)
# - results saved to /obesity_results/Time
#
#
################################################################################

rm(list = ls())

library(stringr)
library(brms)
library(cmdstanr)
library(tidybayes)    
library(gridExtra)
library(Matrix)
library(ggplot2)
library(dplyr)

options(mc.cores = parallel::detectCores())

#set location
loc = "obesity"

#save output?
export = F

# utility functions
#
# correlation matrix for correlated random effects
mat.fun = function(i){
  tmp = matrix(NA, i, i)
  diag(tmp) = 1
  tmp[lower.tri(tmp)] = 1/2
  tmp[upper.tri(tmp)] = 1/2
  tmp
  
}

# independent random effects (for UME)

mat.fun.indep = function(i){
  tmp = matrix(NA, i, i)
  diag(tmp) = 1
  tmp[lower.tri(tmp)] = 0
  tmp[upper.tri(tmp)] = 0
  tmp
}


# MAP treatments

Map.trt = function(object, pars, dname, d){
  out = fixef(object)
  for(i in 1:length(pars)){
    rownames(out)[rownames(out) %in% paste0(pars[i],"_",dname)] = paste0(pars[i],"[",d,"]")  
  }
  print(out, digits=2)
}

# load time course data
dat = readRDS(paste0(loc, "/obesity_data/datTime.rds"))


################################################################################
# Plot of time course data
# Figure 2 (main manuscript)
#

 ggplot(dat, aes(x=week, y=mean, group=trt, colour = Treatment)) + geom_point() + 
       geom_line() + facet_wrap(~Study) +
       scale_color_manual(values=c("darkorange1","darkorange2","darkorange3","darkorange4","#999999",
                                   "blue","blue1","blue2","blue3", "green3", "blue4", "green4",
                                   "blue1", "blue4")) + 
      theme_bw() + theme(legend.position = "bottom") + scale_y_continuous(limits = c(-20,0)) +
      xlab("Time from randomisation (weeks)") + ylab("Body weight (%) change from baseline") 


# Check treatment allocation
dat2 = dat[!duplicated(paste0(dat$author, dat$trt)), ]

################################################################################

# Setup covariance matrix, handling trials with multiple arms
lst = split(dat2, dat2$author)
lst = lapply(1:length(lst), function(j){
  tmp = lst[[j]]
  tmp = tmp[tmp$trt != "Placebo", ]
  as.matrix(bdiag(matrix(1E-20,1,1), mat.fun(length(tmp$author)))) 
} )
st_mat = as.matrix(bdiag(lst))
rownames(st_mat) <- colnames(st_mat) <- 1:nrow(dat2) 


################################################################################
# Model estimation: Exponential model
################################################################################

# Setup design matrix for treatment effects
d = unique(dat$trt[dat$trt!="Placebo"]); dname = paste0("d",1:length(d))
for(i in 1:length(d)) dat[,dname[i]] = ifelse(dat$trt == d[i], 1, 0) 


# Exponential model - fixed effect:
exp.prior = c(
  prior(normal(-3, 100), coef="authorAstrup",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorGarvey",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorKadowaki",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorONeil",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorPiMSunyer",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorRubino",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorWilding",   nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d1", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d2", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d3",nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d4", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d5", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d6",  nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d7", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d8", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d9", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d10", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d11", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d12", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d13",nlpar = "emaxT"),
# Priors on kT
  prior(normal(-2.5, 100), coef="authorAstrup",   nlpar = "kT"),
  prior(normal(-2.5, 100), coef="authorGarvey",   nlpar = "kT"),
  prior(normal(-2.5, 100), coef="authorKadowaki",   nlpar = "kT"),
  prior(normal(-2.5, 100), coef="authorONeil",   nlpar = "kT"),
  prior(normal(-2.5, 100), coef="authorPiMSunyer",   nlpar = "kT"),
  prior(normal(-2.5, 100), coef="authorRubino",   nlpar = "kT"),
  prior(normal(-2.5, 100), coef="authorWilding",   nlpar = "kT"),
  prior(normal(0, 100),  coef="d1", nlpar = "kT"),
  prior(normal(0, 100),  coef="d2", nlpar = "kT"),
  prior(normal(0, 100),  coef="d3",nlpar = "kT"),
  prior(normal(0, 100),  coef="d4", nlpar = "kT"),
  prior(normal(0, 100),  coef="d5", nlpar = "kT"),
  prior(normal(0, 100),  coef="d6",  nlpar = "kT"),
  prior(normal(0, 100),  coef="d7", nlpar = "kT"),
  prior(normal(0, 100),  coef="d8", nlpar = "kT"),
  prior(normal(0, 100),  coef="d9", nlpar = "kT"),
  prior(normal(0, 100),  coef="d10", nlpar = "kT"),
  prior(normal(0, 100),  coef="d11", nlpar = "kT"),
  prior(normal(0, 100),  coef="d12", nlpar = "kT"),
  prior(normal(0, 100),  coef="d13",nlpar = "kT")
)

# Estimate parameters - fixed model:
p1 = as.formula(paste("emaxT ~ -1 + author + ", paste(dname, collapse= "+")))
p2 = as.formula(paste("kT ~ -1 + author + ", paste(dname, collapse= "+")))
exp.NMA.fixed <- brm(bf(mean | se(se) ~ emaxT*( 1 - exp( -exp(kT)*week ) ),
                            p1, p2,
                            nl = TRUE),
                         data = dat,
                         prior = exp.prior,
                     control = list(adapt_delta = 0.999, max_treedepth = 20),
                     iter = 5000, cores = 4, chains = 4,
                     backend = "cmdstanr",
                     threads = threading(6),
                     save_pars = save_pars(all = TRUE), seed = 1234)


exp.NMA.fixed

Map.trt(exp.NMA.fixed, pars = c("emaxT", "kT"), dname, d )

if(export) saveRDS(exp.NMA.fixed, paste0(loc, "/obesity_results/Time/exp.NMA.fixed.rds"))

####################################
# Exponential model - random effect:

exp.prior = c(
  exp.prior,
  prior(normal(0,5), lb=0, class = "sd", nlpar = "emaxT"), 
  prior(normal(0,5), lb=0, class = "sd", nlpar = "kT")
)

# Estimate parameters:
p1 = as.formula(paste("emaxT ~ -1 + author + ", paste(dname, collapse= "+"), "+ (1 | author | gr(stcom, cov = st_mat))" ))
p2 = as.formula(paste("kT ~ -1 + author + ", paste(dname, collapse= "+"), "+ (1 | author | gr(stcom, cov = st_mat))" ))

exp.NMA.ran <- brm(bf(mean | se(se) ~ emaxT*( 1 - exp( -exp(kT)*week ) ),
                      p1, p2, 
                      nl = TRUE),
                   data = dat,
                   data2 = list(st_mat=st_mat),
                   prior = exp.prior,
                   control = list(adapt_delta = 0.999, max_treedepth = 20),
                   iter = 5000, cores = 4, chains = 4,
                   backend = "cmdstanr",
                   threads = threading(6),
                   save_pars = save_pars(all = TRUE), seed = 1234)

exp.NMA.ran

Map.trt(exp.NMA.ran, pars = c("emaxT", "kT"), dname, d )

if(export) saveRDS(exp.NMA.ran, paste0(loc, "/obesity_results/Time/exp.NMA.ran.rds"))

####################################
# Exponential model - random effect with AR1 process:

exp.NMA.ran.ar <- brm(bf(mean | se(se, sigma=TRUE) ~ emaxT*( 1 - exp( -exp(kT)*week ) ),
                      p1, p2, 
                      nl = TRUE, autocor = cor_ar(~ week | stcom, p=1, cov = TRUE)),
                   data = dat,
                   data2 = list(st_mat=st_mat),
                   prior = exp.prior,
                   control = list(adapt_delta = 0.999, max_treedepth = 20),
                   iter = 5000, cores = 4, chains = 4,
                   save_pars = save_pars(all = TRUE), seed = 12345)

exp.NMA.ran.ar

Map.trt(exp.NMA.ran, pars = c("emaxT", "kT"), dname, d )

if(export) saveRDS(exp.NMA.ran.ar, paste0(loc, "/obesity_results/Time/exp.NMA.ran.ar.rds"))

##########################################
# ------------------- Unrelated mean model for model checking:
# When the network includes multiarm trials,
# the UME-Dias model estimates the same vector effects with
# the NMA model in multi-arm trial. However, contrary
# to the NMA model, the UME model treats the random
# study effects as separate univariate normal distributions 
# ------------------

# Setup covariance matrix, handling trials with multiple arms
lst = split(dat2, dat2$author)
lst = lapply(1:length(lst), function(j){
  tmp = lst[[j]]
  tmp = tmp[tmp$trt != "Placebo", ]
  as.matrix(bdiag(matrix(1E-20,1,1), mat.fun.indep(length(tmp$author)))) 
} )
st_mat = as.matrix(bdiag(lst))
rownames(st_mat) <- colnames(st_mat) <- 1:nrow(dat2) 

exp.NMA.ume <- brm(bf(mean | se(se) ~ emaxT*( 1 - exp( -exp(kT)*week ) ),
                      p1, p2, 
                      nl = TRUE),
                   data = dat,
                   data2 = list(st_mat=st_mat),
                   prior = exp.prior,
                   control = list(adapt_delta = 0.9999, max_treedepth = 20), 
                   iter = 5000, cores = 4, chains=4)

exp.NMA.ume

if(export) saveRDS(exp.NMA.ume, paste0(loc, "/obesity_results/Time/expModUme.rds"))

# AR1 process added
exp.NMA.ume.ar <- brm(bf(mean | se(se, sigma=TRUE) ~ emaxT*( 1 - exp( -exp(kT)*week ) ),
                      p1, p2, 
                      nl = TRUE, autocor = cor_ar(~ week | stcom, p=1, cov = TRUE) ),
                   data = dat,
                   data2 = list(st_mat=st_mat),
                   prior = exp.prior,
                   control = list(adapt_delta = 0.9999, max_treedepth = 20), 
                   iter = 5000, cores = 4, chains=4)

exp.NMA.ume.ar

if(export) saveRDS(exp.NMA.ume.ar, paste0(loc, "/obesity_results/Time/exp.NMA.ume.ar.rds"))


################################################################################
# Model estimation: Emax model
################################################################################

# Setup design matrix for treatment effects
d = unique(dat$trt[dat$trt!="Placebo"]); dname = paste0("d",1:length(d))
for(i in 1:length(d)) dat[,dname[i]] = ifelse(dat$trt == d[i], 1, 0) 

# Emax model - fixed effect:
emax.prior = c(
  prior(normal(-3, 100), coef="authorAstrup",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorGarvey",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorKadowaki",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorONeil",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorPiMSunyer",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorRubino",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorWilding",   nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d1", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d2", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d3",nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d4", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d5", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d6",  nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d7", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d8", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d9", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d10", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d11", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d12", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d13",nlpar = "emaxT"),
  # Priors on eT50
  prior(normal(2.5, 100), coef="authorAstrup",    nlpar = "eT50"),
  prior(normal(2.5, 100), coef="authorGarvey",    nlpar = "eT50"),
  prior(normal(2.5, 100), coef="authorKadowaki",  nlpar = "eT50"),
  prior(normal(2.5, 100), coef="authorONeil",     nlpar = "eT50"),
  prior(normal(2.5, 100), coef="authorPiMSunyer", nlpar = "eT50"),
  prior(normal(2.5, 100), coef="authorRubino",    nlpar = "eT50"),
  prior(normal(2.5, 100), coef="authorWilding",   nlpar = "eT50"),
  prior(normal(0, 100),  coef="d1", nlpar = "eT50"),
  prior(normal(0, 100),  coef="d2", nlpar = "eT50"),
  prior(normal(0, 100),  coef="d3",nlpar = "eT50"),
  prior(normal(0, 100),  coef="d4", nlpar = "eT50"),
  prior(normal(0, 100),  coef="d5", nlpar = "eT50"),
  prior(normal(0, 100),  coef="d6",  nlpar = "eT50"),
  prior(normal(0, 100),  coef="d7", nlpar = "eT50"),
  prior(normal(0, 100),  coef="d8", nlpar = "eT50"),
  prior(normal(0, 100),  coef="d9", nlpar = "eT50"),
  prior(normal(0, 100),  coef="d10", nlpar = "eT50"),
  prior(normal(0, 100),  coef="d11", nlpar = "eT50"),
  prior(normal(0, 100),  coef="d12", nlpar = "eT50"),
  prior(normal(0, 100),  coef="d13",nlpar = "eT50")
)

# Estimate parameters - fixed model:
p1 = as.formula(paste("emaxT ~ 0 + author + ", paste(dname, collapse= "+")))
p2 = as.formula(paste("eT50 ~  0 + author + ", paste(dname, collapse= "+")))

Emax.NMA.fixed <- brm(bf(mean | se(se) ~ emaxT*( week/( exp(eT50) + week) ),
                        p1, p2,
                        nl = TRUE),
                      data = dat,
                      prior = emax.prior,
                      control = list(adapt_delta = 0.999, max_treedepth = 20),
                      iter = 5000, cores = 4, chains = 4,
                      backend = "cmdstanr",
                      threads = threading(6),
                      save_pars = save_pars(all = TRUE), seed = 1234)

Emax.NMA.fixed

Map.trt(Emax.NMA.fixed, pars = c("emaxT", "eT50"), dname, d )

if(export) saveRDS(Emax.NMA.fixed, paste0(loc, "/obesity_results/Time/emax.NMA.fixed.rds"))

####################################
# Emax model - random effect:
emax.prior = c(
  emax.prior,
  prior(normal(0,5), lb=0, class = "sd", nlpar = "emaxT"), 
  prior(normal(0,5), lb=0, class = "sd", nlpar = "eT50")
)
#
# Setup covariance matrix, handling trials with multiple arms
lst = split(dat2, dat2$author)
lst = lapply(1:length(lst), function(j){
  tmp = lst[[j]]
  tmp = tmp[tmp$trt != "Placebo", ]
  as.matrix(bdiag(matrix(1E-20,1,1), mat.fun(length(tmp$author)))) 
} )
st_mat = as.matrix(bdiag(lst))
rownames(st_mat) <- colnames(st_mat) <- 1:nrow(dat2) 

# Estimate parameters:
p1 = as.formula(paste("emaxT ~ -1 + author + ", paste(dname, collapse= "+"), "+ (1 | author | gr(stcom, cov = st_mat))" ))
p2 = as.formula(paste("eT50  ~ -1 + author + ", paste(dname, collapse= "+"), "+ (1 | author | gr(stcom, cov = st_mat))" ))

emax.NMA.ran <- brm(bf(mean | se(se) ~ emaxT*( week/( exp(eT50) + week) ),
                      p1, p2, 
                      nl = TRUE),
                   data = dat,
                   data2 = list(st_mat=st_mat),
                   prior = emax.prior,
                   control = list(adapt_delta = 0.999, max_treedepth = 20),
                   iter = 5000, cores = 4, chains = 4,
                   backend = "cmdstanr",
                   threads = threading(6),
                   save_pars = save_pars(all = TRUE), seed = 1234)

out = summary(emax.NMA.ran)

out

Map.trt(emax.NMA.ran, pars = c("emaxT", "eT50"), dname, d )

if(export) saveRDS(emax.NMA.ran, paste0(loc, "/obesity_results/Time/emax.NMA.ran.rds") )

# ar 1 process

emax.NMA.ran.ar <- brm(bf(mean | se(se, sigma = TRUE) ~ emaxT*( week/( exp(eT50) + week) ),
                       p1, p2,
                       nl = TRUE, autocor = cor_ar(~ week | stcom, p=1, cov = TRUE) ),
                    data = dat,
                    data2 = list(st_mat=st_mat), # same issue here
                    prior = emax.prior,
                    control = list(adapt_delta = 0.9999, max_treedepth = 20),
                    iter = 5000, cores = 4, chains=4)

out = emax.NMA.ran.ar

out

Map.trt(out, pars = c("emaxT", "eT50"), dname, d )

if(export) saveRDS(emax.NMA.ran.ar, paste0(loc, "/obesity_results/Time/emax.NMA.ran.ar.rds") )

################################################################################
# Model estimation: Bateman model
################################################################################

# Setup design matrix for treatment effects
d = unique(dat$trt[dat$trt!="Placebo"]); dname = paste0("d",1:length(d))
for(i in 1:length(d)) dat[,dname[i]] = ifelse(dat$trt == d[i], 1, 0) 

# Fixed effect model specification - Bateman model:
bateman.prior = c(
  prior(normal(-3, 100), coef="authorAstrup",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorGarvey",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorKadowaki",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorONeil",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorPiMSunyer",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorRubino",   nlpar = "emaxT"),
  prior(normal(-3, 100), coef="authorWilding",   nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d1", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d2", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d3",nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d4", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d5", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d6",  nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d7", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d8", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d9", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d10", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d11", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d12", nlpar = "emaxT"),
  prior(normal(0, 100),  coef="d13",nlpar = "emaxT"),
  # Priors on kT
  prior(normal(-2.5, 100), coef="authorAstrup",   nlpar = "kT"),
  prior(normal(-2.5, 100), coef="authorGarvey",   nlpar = "kT"),
  prior(normal(-2.5, 100), coef="authorKadowaki",   nlpar = "kT"),
  prior(normal(-2.5, 100), coef="authorONeil",   nlpar = "kT"),
  prior(normal(-2.5, 100), coef="authorPiMSunyer",   nlpar = "kT"),
  prior(normal(-2.5, 100), coef="authorRubino",   nlpar = "kT"),
  prior(normal(-2.5, 100), coef="authorWilding",   nlpar = "kT"),
  prior(normal(0, 100),  coef="d1", nlpar = "kT"),
  prior(normal(0, 100),  coef="d2", nlpar = "kT"),
  prior(normal(0, 100),  coef="d3",nlpar = "kT"),
  prior(normal(0, 100),  coef="d4", nlpar = "kT"),
  prior(normal(0, 100),  coef="d5", nlpar = "kT"),
  prior(normal(0, 100),  coef="d6",  nlpar = "kT"),
  prior(normal(0, 100),  coef="d7", nlpar = "kT"),
  prior(normal(0, 100),  coef="d8", nlpar = "kT"),
  prior(normal(0, 100),  coef="d9", nlpar = "kT"),
  prior(normal(0, 100),  coef="d10", nlpar = "kT"),
  prior(normal(0, 100),  coef="d11", nlpar = "kT"),
  prior(normal(0, 100),  coef="d12", nlpar = "kT"),
  prior(normal(0, 100),  coef="d13",nlpar = "kT"),

  # Priors on k
  prior(normal(-5, 100), coef="authorAstrup",   nlpar = "k"),
  prior(normal(-5, 100), coef="authorGarvey",   nlpar = "k"),
  prior(normal(-5, 100), coef="authorKadowaki",   nlpar = "k"),
  prior(normal(-5, 100), coef="authorONeil",   nlpar = "k"),
  prior(normal(-5, 100), coef="authorPiMSunyer",   nlpar = "k"),
  prior(normal(-5, 100), coef="authorRubino",   nlpar = "k"),
  prior(normal(-5, 100), coef="authorWilding",   nlpar = "k"),
  prior(normal(0, 100),  coef="d1", nlpar = "k"),
  prior(normal(0, 100),  coef="d2", nlpar = "k"),
  prior(normal(0, 100),  coef="d3",nlpar = "k"),
  prior(normal(0, 100),  coef="d4", nlpar = "k"),
  prior(normal(0, 100),  coef="d5", nlpar = "k"),
  prior(normal(0, 100),  coef="d6",  nlpar = "k"),
  prior(normal(0, 100),  coef="d7", nlpar = "k"),
  prior(normal(0, 100),  coef="d8", nlpar = "k"),
  prior(normal(0, 100),  coef="d9", nlpar = "k"),
  prior(normal(0, 100),  coef="d10", nlpar = "k"),
  prior(normal(0, 100),  coef="d11", nlpar = "k"),
  prior(normal(0, 100),  coef="d12", nlpar = "k"),
  prior(normal(0, 100),  coef="d13",nlpar = "k")
 )

# Estimate parameters:
p1 = as.formula(paste("emaxT ~ 0 + author + ", paste(dname, collapse= "+")))
p2 = as.formula(paste("kT ~ 0 + author + ", paste(dname, collapse= "+")))
p3 = as.formula(paste("k ~ 0 + author + ", paste(dname, collapse= "+")))

bateman.NMA.fixed <- brm(bf(mean | se(se) ~ emaxT*( exp(kT)/( exp(kT) - exp(k) ) )*( exp(-exp(k)*week) - exp(-exp(kT)*week) ),
                            p1, p2, p3,
                            nl = TRUE),
                         data = dat,
                         prior = bateman.prior,
                         control = list(adapt_delta = 0.999, max_treedepth = 20),
                         iter = 15000, cores = 4, chains = 4,
                         backend = "cmdstanr",
                         threads = threading(6),
                         save_pars = save_pars(all = TRUE), seed = 1234)

bateman.NMA.fixed
Map.trt(bateman.NMA.fixed, pars = c("emaxT", "kT", "k"), dname, d )
if(export) saveRDS(bateman.NMA.fixed, paste0(loc, "/obesity_results/Time/bateman.NMA.fixed.rds"))