# Time course models:
#
# - Table 5 main manuscript
#
rm(list = ls()) # clear memory
options(mc.cores = parallel::detectCores()) # setup for parallel

library(gridExtra)
library(ggplot2)

loc = "obesity"

set.seed(1234)
# load time course data
dat = readRDS(paste0(loc, "/obesity_data/datTime.rds"))

# utility functions
################################################################################
Map.trt = function(object, pars, dname, d){
  out = fixef(object)
  for(i in 1:length(pars)){
    rownames(out)[rownames(out) %in% paste0(pars[i],"_",dname)] = paste0(pars[i],"[",d,"]")  
  }
  print(out, digits=2)
}

# ---- Set study and plot
new.dat = function(data, pub, end){
  new = subset(data, author == pub)
  new = split(new, new$trt)
  new = lapply(seq_along(new), function(i){
    tmp=new[[i]]
    tmp=tmp[,c("author","trt",dname,"stcom","se")]
    tmp=data.frame(week = 0:end, tmp[1,])
    tmp
  })
  new = do.call("rbind", new)
  return(new)  
}


# ---- Deviance plots
dev.fun = function(object, dat, model){
  res = residual_draws(object, newdata = dat, seed=1234) %>%
    mutate( deviance = .residual**2/se**2 ) %>%
    group_by(author, trt, week) %>%
    summarise(res = median(deviance), q25 = quantile(deviance, probs=0.025), q975= quantile(deviance, probs=0.975) )
  res$model = model
  return(res)
}

################################################################################
# Goodness of statistics

# Load fixed effect models:
exp.fix = readRDS(paste0(loc,"/obesity_results/Time/exp.NMA.fixed.rds"))
emax.fix = readRDS(paste0(loc,"/obesity_results/Time/emax.NMA.fixed.rds"))
bateman.fix = readRDS(paste0(loc,"/obesity_results/Time/bateman.NMA.fixed.rds"))  
loo.exp.fix     = loo(exp.fix, moment_match = TRUE) 
loo.emax.fix    = loo(emax.fix, moment_match = TRUE) 
loo.bateman.fix = loo(bateman.fix, moment_match = TRUE)

# Load random effect models:
exp.ran     = readRDS(paste0(loc,"/obesity_results/Time/exp.NMA.ran.rds"))
emax.ran    = readRDS(paste0(loc,"/obesity_results/Time/emax.NMA.ran.rds"))
bateman.ran = readRDS(paste0(loc,"/obesity_results/Time/batemanModRan.rds"))  
loo.exp.ran    = loo(exp.ran, moment_match = TRUE) 
loo.emax.ran    = loo(emax.ran, moment_match = TRUE) 
loo.bateman.ran = loo(bateman.ran, moment_match = TRUE)

# Load random effect with AR(1)
exp.ar      = readRDS(paste0(loc,"/obesity_results/Time/exp.NMA.ran.ar.rds"))
emax.ar     = readRDS(paste0(loc,"/obesity_results/Time/emax.NMA.ran.ar.rds"))
loo.exp.ar  = loo(exp.ar) 
loo.emax.ar = loo(emax.ar) 

out = loo_compare(loo.exp.fix, loo.emax.fix, loo.bateman.fix, loo.exp.ran, loo.emax.ran, loo.exp.ar, loo.emax.ar)
out = as.data.frame(out)
# row.names(out):       "exp.ar"       "emax.ar"      "emax.ran" "emax.fix"     "exp.ran"       "exp.fix"  "bateman.fix" 
out$Function   = c("Exponential" ,        "Emax",         "Emax",   "Emax",  "Exponential", "Exponential", "Bateman")
out$Effect     = c("Exchangeable","Exchangeable", "Exchangeable", "Common", "Exchangeable",     "Common" ,  "Common")
out$Within.arm = c("AR(1)"       ,"AR(1)"       ,         "None",   "None",         "None",       "None" ,  "None"  )
out$name = row.names(out)
print(out[, c("Function", "Effect", "Within.arm", "elpd_loo", "se_elpd_loo", "elpd_diff", "se_diff")], row.names = F, digits=3)

for(i in c("elpd_loo", "se_elpd_loo", "elpd_diff", "se_diff") ) out[,i] = format(round(out[,i], 2), nsmall = 2)

out

################################################################################
# Final model

final = exp.ar

# Setup design matrix for treatment effects
d = unique(dat$trt[dat$trt!="Placebo"]); dname = paste0("d",1:length(d))
for(i in 1:length(d)) dat[,dname[i]] = ifelse(dat$trt == d[i], 1, 0) 

Map.trt(final, pars = c("emaxT", "kT"), dname, d )
