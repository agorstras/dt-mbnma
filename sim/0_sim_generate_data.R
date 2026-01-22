################################################################################
# Description: Generate simulation data
#   - generates 100 replicates of data set (10 trials)
#   - saved in ./sim_data/
#   - plots of mean trajectories of trials, and dose response (Figures S2A, S2B)
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

# create directories for storing /results
#
dir.create(paste0(loc, "/sim_results/AR1"),showWarnings = FALSE,recursive = TRUE)
dir.create(paste0(loc, "/sim_results/obs"),showWarnings = FALSE,recursive = TRUE)
dir.create(paste0(loc, "/sim_results/univariate_ll"),showWarnings = FALSE,recursive = TRUE)


# set seed for reproducibility
set.seed(1234)
export = T # set to T if datasets are to be saved

# utility functions
emax = function(dose, ed50, emax) dose*emax/(dose + ed50)

#########################
# Simulation design
design = readxl::read_xlsx(path = paste0(loc,"/sim_design/Simulation_Design.xlsx"))
design$dose = as.numeric(design$dose)

# observed within-arm correlations 
cov    = read_sas( paste0(loc,"/sim_design/covparms_pchg_bw.sas7bdat") )
visit  = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40, 44, 48, 52)

covMat = matrix(NA, 18, 18)
covMat[upper.tri(covMat, diag = T)] = cov$Estimate
covMat[lower.tri(covMat, diag = F)] = t(covMat)[lower.tri(covMat, diag = F)]
rownames(covMat) = colnames(covMat) = visit
covSim = covMat[rownames(covMat) %in% c(2, 4, 8, 12, 16, 20, 28, 36, 44, 52), colnames(covMat) %in% c(2, 4, 8, 12, 16, 20, 28, 36, 44, 52)]


# Simulation function mean dose-response

mu = function(dA, dB, dC, dose, time, arm, etaE, etaK){
  
  if( arm == 1 ){
    emaxT = -2.5  
    kT    = -2.0     
  } else{
    emaxT = -2.5 + (-10*dA + -15*dB + -20*dC)*dose/( exp( -2*dA + -1.5*dB + -0.75*dC ) + dose ) + etaE
    kT    = -2.0 + -0.2*dA*log(dose + 1) + -0.4*dB*log(dose + 1) + -0.8*dC*log(dose + 1) + etaK
  }
  
  return( emaxT * (1 - exp(-exp(kT) * time)) ) 
  
}

# correlation matrix for correlated random effects

mat.fun = function(i){
  tmp = matrix(NA, i, i)
  diag(tmp) = 1
  tmp[lower.tri(tmp)] = 1/2
  tmp[upper.tri(tmp)] = 1/2
  tmp
  
}

## Main Simulation Step
#
designLst = split(design, design$study)
dataLst   = vector("list", length=100) # number of data sets
varLst    = vector("list", length=100) 

for(i in 1:length(dataLst)){
  
  tmplst = lapply(1:length(designLst), function(j){
    
    tmp = designLst[[j]]
    
    # Study effects
    tau = matrix( c(0.5**2, 0.25*0.5*0.2, 0.25*0.5*0.2, 0.2**2), 2, 2)
    n.cntr = nrow(tmp)-1
    omega = kronecker(tau, mat.fun(n.cntr))
    eta  = mvrnorm(n = 1, mu = rep(0,2*n.cntr), omega)
    etaE = c(NA, eta[1:n.cntr])
    etaK = c(NA, eta[(n.cntr+1):(2*n.cntr)])
    
    # Generate study aggregate level data
    datlst = lapply(1:nrow(tmp), function(i){
      ipd = mvrnorm(n = tmp$N[i], mu = mu(dA = tmp$dA[i], dB = tmp$dB[i], dC = tmp$dC[i], dose=tmp$dose[i], time = c(2, 4, 8, 12, 16, 20, 28, 36, 44, 52), arm = tmp$arm[i], etaE = etaE[i], etaK[i]), Sigma = covSim) 
      R   = cor(ipd)
      M   = data.frame(study=tmp$study[i], arm=tmp$arm[i], week=c(2, 4, 8, 12, 16, 20, 28, 36, 44, 52), y=apply(ipd, 2, mean), se=apply(ipd, 2, function(x) sd(x)/sqrt(length(x))) )  
      return( list( M=M, V = diag(x = M$se, ncol=length(M$se), nrow=length(M$se)) %*% R %*% t(diag(x = M$se, ncol=length(M$se), nrow=length(M$se))) ) )  
    })
    Mlst = lapply(1:nrow(tmp), function(i){ datlst[[i]]$M })
    M    = do.call("rbind", Mlst)
    Vlst = lapply(1:nrow(tmp), function(i){ datlst[[i]]$V })
    V    = as.matrix(bdiag( Vlst ))
    return(list(M=M,V=V))
  })
  
  Mlst = lapply(1:length(designLst), function(j){ tmplst[[j]]$M  })
  dat = do.call("rbind", Mlst )   
  dat = merge(dat, design, by = c("study", "arm"), all.x = T, sort=F)
  Vlst= lapply(1:length(designLst), function(j){ tmplst[[j]]$V  })
  
  
  dataLst[[i]] = dat
  varLst[[i]]  = Vlst
}

if( export ) saveRDS(dataLst,paste0(loc,"/sim_data/dataLst.rds"))
if( export ) saveRDS(varLst,paste0(loc,"/sim_data/varLst.rds"))

## Visualise

# Treatment network (Figure S1)
design$treatment = paste(design$agent, design$dose, " mg") 
design$treatment[design$agent == "soc"] = "soc"
design$estimate = design$se = 1


net <-  set_agd_arm(design, 
                    study = study,
                    trt = treatment,
                    y = estimate, 
                    se = se,
                    sample_size = N,
                    trt_ref = "soc")

plot(net, weight_edges = TRUE, weight_nodes = TRUE)

# mean trajectories by study (Figure S2A)
lst = lapply(1:nrow(design), function(i){
  tmp = design[i,]
  data.frame(study=tmp$study, arm=tmp$arm, agent=tmp$agent, dose=tmp$dose, time=0:52,
             bw = mu(dA=tmp$dA, dB=tmp$dB, dC=tmp$dC, dose=tmp$dose, time=0:52, arm=tmp$arm, etaE = 0, etaK = 0) )
})

plt = do.call("rbind", lst)

ggplot(plt, aes(x=time, y=bw, group=arm, colour=agent)) + 
  geom_line(linewidth=1) + 
  facet_wrap(~study, nrow = 3, ncol = 4) + 
  theme_bw() + 
  xlab("Time from randomisation") + ylab("Absolute response")


# Allocation of doses and dose-response profiles (figure S2B) 
dose=seq(0,4,by=0.01)
plot(x=c(0,4), y=c(-20,0), ylab="Relative effect (%-point)", xlab="Dose (mg)", type="n")
lines(x=dose, y=emax(dose, ed50 = exp(-2), emax = -10), lwd=2); points(x=design$dose[design$agent=="A"], y= emax(design$dose[design$agent=="A"], ed50 = exp(-2), emax = -10), pch=16 )
lines(x=dose, y=emax(dose, ed50 = exp(-1.5), emax = -15), lwd=2, col=2); points(x=design$dose[design$agent=="B"], y=emax(design$dose[design$agent=="B"], ed50 = exp(-1.5), emax = -15), col=2, pch=16 )
lines(x=dose, y=emax(dose, ed50 = exp(-0.75), emax = -20), lwd=2, col=4); points(x=design$dose[design$agent=="C"], y=emax(design$dose[design$agent=="C"], ed50 = exp(-0.75), emax = -20), col=4, pch=16 )
legend("topright", paste("Agent", LETTERS[1:3]), col=c(1,2,4), pch = 16 )

