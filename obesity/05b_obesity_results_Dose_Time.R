######################################
# Figure 3: main manuscript
# Figures S7-S12: supplement
# Table 6: main manuscript
# Figure S6: supplement
# Table 7: main manuscript
#

library(gridExtra)
library(ggplot2)
library(multinma)

rm(list = ls()) # clear memory
options(mc.cores = parallel::detectCores()) # setup for parallel

loc = "obesity"

set.seed(1234)
# utility functions
################################################################################
Map.trt = function(object, pars, dname, d){
  out = fixef(object)
  for(i in 1:length(pars)){
    rownames(out)[rownames(out) %in% paste0(pars[i],"_",dname)] = paste0(pars[i],"[",d,"]")  
  }
  print(out, digits=2)
}

new.dat = function(data, pub, end, predictors){
  new = subset(data, author == pub)
  new = split(new, new$trt)
  new = lapply(seq_along(new), function(i){
    tmp=new[[i]]
    tmp=tmp[1,c("author","trt",predictors,"stcom","se")]
    tmp=merge(data.frame(week = 0:end), tmp, all.x=F)
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

mat.fun = function(i){
  tmp = matrix(NA, i, i)
  diag(tmp) = 1
  tmp[lower.tri(tmp)] = 1/2
  tmp[upper.tri(tmp)] = 1/2
  tmp
  
}

# load time course data
dat = readRDS(paste0(loc, "/obesity_data/datDoseTime.rds"))
#
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

dat$dLira = dat$d1*log(dat$mdose + 1)  
dat$dSema = dat$d2*log(dat$mdose + 1)
#
# set Final model
final = readRDS( paste0(loc, "/obesity_results/DoseTime/exp.NMA.emax.log.ar.rds") )

################################################################################                
# Change plot settings to align with data plot 

studyEff = NULL # Predicted outcomes for trial arms, including random effects in structural parameters of time-course models
# studyEff = NA # Predicted outcomes for trial arms, excluding random effects in structural parameters of time-course models

preds = c("d1","d2","dLira","dSema","mdose", "Study","Treatment") 
new.data = rbind(new.dat(dat, "O'Neil", 52, predictors=preds),
                 new.dat(dat, "Wilding", 68, predictors=preds),
                 new.dat(dat, "Kadowaki", 68, predictors=preds), 
                 new.dat(dat, "Rubino", 68, predictors=preds),
                 new.dat(dat, "Garvey", 104, predictors=preds),
                 new.dat(dat, "Astrup", 20, predictors=preds),
                 new.dat(dat, "Pi-Sunyer", 56, predictors=preds) )

pred <- final %>% 
  epred_draws( newdata = new.data, re_formula = studyEff )

ggplot(pred, aes(x = week, y = .epred, group = Treatment, colour=Treatment)) + facet_wrap(~Study) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Greys") +
  labs(x = "Time from randomisation (weeks)", y = "Body weight (%) change from baseline", fill = "Credible interval") +
  geom_pointrange(data=dat, aes(x=week, y=mean, ymin = mean - se, ymax=mean + se), size=0.25) +
  theme_bw() +  theme(legend.position = "bottom") + 
  scale_color_manual(values=c("darkorange1","darkorange2","darkorange3","darkorange4","#999999",
                              "blue","blue1","blue2","blue3", "green3", "blue4", "green4",
                              "blue1", "blue4"))

################################################################################
# Semaglutide

studyEff = NULL # Predicted outcomes for trial arms, including random effects in structural parameters of time-course models
# studyEff = NA   # Predicted outcomes for trial arms, excluding random effects in structural parameters of time-course models

x = "O'Neil"
pred <- final %>% 
  epred_draws( newdata = new.dat(dat, x, 52, predictors=c("d1","d2","dLira","dSema","mdose")), re_formula = studyEff )

sema = ggplot(pred, aes(x = week, y = .epred, group = trt, colour=trt)) + # facet_wrap(~trt, scales = "free_y") +
  stat_lineribbon() + ggtitle( paste("Author:", x) ) + 
  scale_fill_brewer(palette = "Greys") +
  labs(x = "Time from randomisation (weeks)", y = "Body weight (%) change from baseline", fill = "Credible interval") +
  geom_pointrange(data=subset(dat, author == x), aes(x=week, y=mean, ymin = mean - se, ymax=mean + se), size=0.5) +
  theme_bw() +  theme(legend.position = "bottom") +
  scale_color_manual(values=c("darkorange4","#999999","blue","blue1","blue2","blue3", "green3", "blue4", "green4"))
print(sema)


x = "Wilding"
pred <- final %>% 
  epred_draws( newdata = new.dat(dat, x, 68, predictors=c("d1","d2","dLira","dSema","mdose")), re_formula = studyEff )

step1 = ggplot(pred, aes(x = week, y = .epred, group = trt, colour=trt)) + #facet_wrap(~trt, scales = "free_y") +
  stat_lineribbon() + ggtitle( paste("Author:", x) ) +
  scale_fill_brewer(palette = "Greys") +
  labs(x = "Time from randomisation (weeks)", y = "Body weight (%) change from baseline", fill = "Credible interval") +
  geom_pointrange(data=subset(dat, author == x), aes(x=week, y=mean, ymin = mean - se, ymax=mean + se), size=0.5) +
  theme_bw() +  theme(legend.position = "bottom") + scale_color_manual(values=c("#999999", "blue4"))
step1

x = "Kadowaki"
pred <- final %>% 
  epred_draws( newdata = new.dat(dat, x, 68, predictors=c("d1","d2","dLira","dSema","mdose")), re_formula = studyEff )

step6 = ggplot(pred, aes(x = week, y = .epred, group = trt, colour=trt)) + #facet_wrap(~trt, scales = "free_y") +
  stat_lineribbon() + ggtitle( paste("Author:", x) ) +
  scale_fill_brewer(palette = "Greys") +
  labs(x = "Time from randomisation (weeks)", y = "Body weight (%) change from baseline", fill = "Credible interval") +
  geom_pointrange(data=subset(dat, author == x), aes(x=week, y=mean, ymin = mean - se, ymax=mean + se), size=0.5) +
  theme_bw() +  theme(legend.position = "bottom") + scale_color_manual(values=c("#999999", "blue1", "blue4"))
step6

x = "Rubino"
pred <- final %>% 
  epred_draws( newdata = new.dat(dat, x, 68, predictors=c("d1","d2","dLira","dSema","mdose")), re_formula = studyEff )

step8 = ggplot(pred, aes(x = week, y = .epred, group = trt, colour=trt)) + # facet_wrap(~trt, scales = "free_y") +
  stat_lineribbon() + ggtitle( paste("Author:", x) ) + 
  scale_fill_brewer(palette = "Greys") +
  labs(x = "Time from randomisation (weeks)", y = "Body weight (%) change from baseline", fill = "Credible interval") +
  geom_pointrange(data=subset(dat, author == x), aes(x=week, y=mean, ymin = mean - se, ymax=mean + se), size=0.5) +
  theme_bw() +  theme(legend.position = "bottom") + scale_color_manual(values=c("darkorange4", "#999999", "blue4"))
step8

x = "Garvey"
pred <- final %>% 
  epred_draws( newdata = new.dat(dat, x, 104, predictors=c("d1","d2","dLira","dSema","mdose")), re_formula = studyEff )

step5 = ggplot(pred, aes(x = week, y = .epred, group = trt, colour=trt)) + #facet_wrap(~trt, scales = "free_y") +
  stat_lineribbon() + ggtitle( paste("Author:", x) ) + 
  scale_fill_brewer(palette = "Greys") +
  labs(x = "Time from randomisation (weeks)", y = "Body weight (%) change from baseline", fill = "Credible interval") +
  geom_pointrange(data=subset(dat, author == x), aes(x=week, y=mean, ymin = mean - se, ymax=mean + se), size=0.5) +
  theme_bw() +  theme(legend.position = "bottom") + scale_color_manual(values=c("#999999", "blue4"))
step5

#print( grid.arrange(step1, step5, step6, step8, nrow = 2, ncol=2) )

################################################################################
# Liraglutide

x = "Astrup"
pred <- final %>% 
  epred_draws( newdata = new.dat(dat, x, 20, predictors=c("d1","d2","dLira","dSema","mdose")), re_formula = studyEff )

lira = ggplot(pred, aes(x = week, y = .epred, group = trt, colour=trt)) + # facet_wrap(~trt, scales = "free_y") +
  stat_lineribbon() + ggtitle( paste("Author:", x) ) +
  scale_fill_brewer(palette = "Greys") +
  labs(x = "Time from randomisation (weeks)", y = "Body weight (%) change from baseline", fill = "Credible interval") +
  geom_pointrange(data=subset(dat, author == x), aes(x=week, y=mean, ymin = mean - se, ymax=mean + se), size=0.5) +
  theme_bw() +  theme(legend.position = "bottom") + scale_color_manual(values=c("darkorange1","darkorange2","darkorange3","darkorange4","#999999"))
lira

x = "Pi-Sunyer"
pred <- final %>% 
  epred_draws( newdata = new.dat(dat, x, 56, predictors=c("d1","d2","dLira","dSema","mdose")), re_formula = studyEff )

scale = ggplot(pred, aes(x = week, y = .epred, group = trt, colour=trt)) + # facet_wrap(~trt, scales = "free_y") +
  stat_lineribbon() + ggtitle( paste("Author:", x) ) +
  scale_fill_brewer(palette = "Greys") +
  labs(x = "Time from randomisation (weeks)", y = "Body weight (%) change from baseline", fill = "Credible interval") +
  geom_pointrange(data=subset(dat, author == x), aes(x=week, y=mean, ymin = mean - se, ymax=mean + se), size=0.5) +
  theme_bw() +  theme(legend.position = "bottom") + scale_color_manual(values=c("darkorange4", "#999999"))
scale

print( grid.arrange(lira, scale, nrow = 2, ncol=1) )

################################################################################
# Compare models

# Table 6
# Load random effect with AR(1)
NMA.ar      = readRDS(paste0(loc,"/obesity_results/Time/exp.NMA.ran.ar.rds"))
loo.NMA.ar  = loo(NMA.ar) 

UME.ar      = readRDS(paste0(loc,"/obesity_results/Time/exp.NMA.ume.ar.rds"))
loo.ume.ar  = loo(UME.ar) 

Emax.log.ar  = readRDS(paste0(loc,"/obesity_results/DoseTime/exp.NMA.emax.log.ar.rds"))
loo.emax.log.ar = loo(Emax.log.ar) 

Log.log.ar  = readRDS(paste0(loc,"/obesity_results/DoseTime/exp.NMA.log.1.ar.rds"))
loo.log.log.ar = loo(Log.log.ar) 

Lin.lin.ar  = readRDS(paste0(loc,"/obesity_results/DoseTime/exp.NMA.lin.1.ar.rds"))
loo.Lin.lin.ar = loo(Lin.lin.ar) 

# Rank models
out = loo_compare(loo.NMA.ar, loo.ume.ar, loo.emax.log.ar, loo.log.log.ar, loo.Lin.lin.ar)
out = as.data.frame(out)
# row.names(out):       
out$EmaxT= c("Linear","Log-linear", "Treatment level", "UME","Emax")
out$kT   = c("Linear","Log-linear", "Treatment level", "UME","Log-linear")
out$name = row.names(out)
print(out[, c("EmaxT", "kT", "elpd_loo", "se_elpd_loo", "elpd_diff", "se_diff")], row.names = F, digits=3)

out$Between.EmaxT[1] = paste0(format(round(VarCorr(Lin.lin.ar)$stcom$sd[1,1], 2), nsmall = 2), " [", format(round(VarCorr(Lin.lin.ar)$stcom$sd[1,3], 2), nsmall = 2), " ; ", format(round(VarCorr(Lin.lin.ar)$stcom$sd[1,4], 2), nsmall = 2), "]")   
out$Between.kT[1]    = paste0(format(round(VarCorr(Lin.lin.ar)$stcom$sd[2,1], 2), nsmall = 2), " [", format(round(VarCorr(Lin.lin.ar)$stcom$sd[2,3], 2), nsmall = 2), " ; ", format(round(VarCorr(Lin.lin.ar)$stcom$sd[2,4], 2), nsmall = 2), "]")   
out$Corr.EmaxT.kT[1] = paste0(format(round(VarCorr(Lin.lin.ar)$stcom$cor[[2]], 2), nsmall = 2), " [", format(round(VarCorr(Lin.lin.ar)$stcom$cor[[6]], 2), nsmall = 2), " ; ", format(round(VarCorr(Lin.lin.ar)$stcom$cor[[8]], 2), nsmall = 2), "]")   
out$ar.corr[1]       = paste0(format(round(posterior_summary(Lin.lin.ar, variable =  "ar[1]")[1], 2), nsmall = 2), " [", format(round(posterior_summary(Lin.lin.ar, variable =  "ar[1]")[3], 2), nsmall = 2), " ; ", format(round(posterior_summary(Lin.lin.ar, variable =  "ar[1]")[4], 2), nsmall = 2), "]" )

out$Between.EmaxT[2] = paste0(format(round(VarCorr(Log.log.ar)$stcom$sd[1,1], 2), nsmall = 2), " [", format(round(VarCorr(Log.log.ar)$stcom$sd[1,3], 2), nsmall = 2), " ; ", format(round(VarCorr(Log.log.ar)$stcom$sd[1,4], 2), nsmall = 2), "]")   
out$Between.kT[2]    = paste0(format(round(VarCorr(Log.log.ar)$stcom$sd[2,1], 2), nsmall = 2), " [", format(round(VarCorr(Log.log.ar)$stcom$sd[2,3], 2), nsmall = 2), " ; ", format(round(VarCorr(Log.log.ar)$stcom$sd[2,4], 2), nsmall = 2), "]")   
out$Corr.EmaxT.kT[2] = paste0(format(round(VarCorr(Log.log.ar)$stcom$cor[[2]], 2), nsmall = 2), " [", format(round(VarCorr(Log.log.ar)$stcom$cor[[6]], 2), nsmall = 2), " ; ", format(round(VarCorr(Log.log.ar)$stcom$cor[[8]], 2), nsmall = 2), "]")   
out$ar.corr[2]       = paste0(format(round(posterior_summary(Log.log.ar, variable =  "ar[1]")[1], 2), nsmall = 2), " [", format(round(posterior_summary(Log.log.ar, variable =  "ar[1]")[3], 2), nsmall = 2), " ; ", format(round(posterior_summary(Log.log.ar, variable =  "ar[1]")[4], 2), nsmall = 2), "]" )

out$Between.EmaxT[3] = paste0(format(round(VarCorr(NMA.ar)$stcom$sd[1,1], 2), nsmall = 2), " [", format(round(VarCorr(NMA.ar)$stcom$sd[1,3], 2), nsmall = 2), " ; ", format(round(VarCorr(NMA.ar)$stcom$sd[1,4], 2), nsmall = 2), "]")   
out$Between.kT[3]    = paste0(format(round(VarCorr(NMA.ar)$stcom$sd[2,1], 2), nsmall = 2), " [", format(round(VarCorr(NMA.ar)$stcom$sd[2,3], 2), nsmall = 2), " ; ", format(round(VarCorr(NMA.ar)$stcom$sd[2,4], 2), nsmall = 2), "]")   
out$Corr.EmaxT.kT[3] = paste0(format(round(VarCorr(NMA.ar)$stcom$cor[[2]], 2), nsmall = 2), " [", format(round(VarCorr(NMA.ar)$stcom$cor[[6]], 2), nsmall = 2), " ; ", format(round(VarCorr(NMA.ar)$stcom$cor[[8]], 2), nsmall = 2), "]")   
out$ar.corr[3]       = paste0(format(round(posterior_summary(NMA.ar, variable =  "ar[1]")[1], 2), nsmall = 2), " [", format(round(posterior_summary(NMA.ar, variable =  "ar[1]")[3], 2), nsmall = 2), " ; ", format(round(posterior_summary(NMA.ar, variable =  "ar[1]")[4], 2), nsmall = 2), "]" )

out$Between.EmaxT[4] = paste0(format(round(VarCorr(UME.ar)$stcom$sd[1,1], 2), nsmall = 2), " [", format(round(VarCorr(UME.ar)$stcom$sd[1,3], 2), nsmall = 2), " ; ", format(round(VarCorr(UME.ar)$stcom$sd[1,4], 2), nsmall = 2), "]")   
out$Between.kT[4]    = paste0(format(round(VarCorr(UME.ar)$stcom$sd[2,1], 2), nsmall = 2), " [", format(round(VarCorr(UME.ar)$stcom$sd[2,3], 2), nsmall = 2), " ; ", format(round(VarCorr(UME.ar)$stcom$sd[2,4], 2), nsmall = 2), "]")   
out$Corr.EmaxT.kT[4] = paste0(format(round(VarCorr(UME.ar)$stcom$cor[[2]], 2), nsmall = 2), " [", format(round(VarCorr(UME.ar)$stcom$cor[[6]], 2), nsmall = 2), " ; ", format(round(VarCorr(UME.ar)$stcom$cor[[8]], 2), nsmall = 2), "]")   
out$ar.corr[4]       = paste0(format(round(posterior_summary(UME.ar, variable =  "ar[1]")[1], 2), nsmall = 2), " [", format(round(posterior_summary(UME.ar, variable =  "ar[1]")[3], 2), nsmall = 2), " ; ", format(round(posterior_summary(UME.ar, variable =  "ar[1]")[4], 2), nsmall = 2), "]" )

out$Between.EmaxT[5] = paste0(format(round(VarCorr(Emax.log.ar)$stcom$sd[1,1], 2), nsmall = 2), " [", format(round(VarCorr(Emax.log.ar)$stcom$sd[1,3], 2), nsmall = 2), " ; ", format(round(VarCorr(Emax.log.ar)$stcom$sd[1,4], 2), nsmall = 2), "]")   
out$Between.kT[5]    = paste0(format(round(VarCorr(Emax.log.ar)$stcom$sd[2,1], 2), nsmall = 2), " [", format(round(VarCorr(Emax.log.ar)$stcom$sd[2,3], 2), nsmall = 2), " ; ", format(round(VarCorr(Emax.log.ar)$stcom$sd[2,4], 2), nsmall = 2), "]")   
out$Corr.EmaxT.kT[5] = paste0(format(round(VarCorr(Emax.log.ar)$stcom$cor[[2]], 2), nsmall = 2), " [", format(round(VarCorr(Emax.log.ar)$stcom$cor[[6]], 2), nsmall = 2), " ; ", format(round(VarCorr(Emax.log.ar)$stcom$cor[[8]], 2), nsmall = 2), "]")   
out$ar.corr[5]       = paste0(format(round(posterior_summary(Emax.log.ar, variable =  "ar[1]")[1], 2), nsmall = 2), " [", format(round(posterior_summary(Emax.log.ar, variable =  "ar[1]")[3], 2), nsmall = 2), " ; ", format(round(posterior_summary(Emax.log.ar, variable =  "ar[1]")[4], 2), nsmall = 2), "]" )

print(out[, c("EmaxT", "kT", "elpd_loo", "se_elpd_loo", "elpd_diff", "se_diff", "Between.EmaxT", "Between.kT", "Corr.EmaxT.kT", "ar.corr")], row.names = F, digits=3)

################################################################################        
#
# Figure S6
#
# Time course: Exponential | Dose-Response: Treatment-level NMA
# Univariate models
NMA = readRDS(paste0(loc,"/obesity_results/Time/exp.NMA.ran.rds"))

# Setup design matrix for treatment effects - all relative to placebo
d = unique(dat$trt[dat$trt!="Placebo"]); dname = paste0("d",1:length(d))
for(i in 1:length(d)) dat[,dname[i]] = ifelse(dat$trt == d[i], 1, 0) 

exp.NMA.res = residual_draws(NMA, newdata = dat) %>%
  group_by(author, trt, week) %>%
  summarise(res = mean(.residual**2/se**2), 
            q25 = quantile(.residual**2/se**2, probs=0.025),
            q975= quantile(.residual**2/se**2, probs=0.975) )
exp.NMA.res$model = "Dose-Response: Treatment-level"

# Time-course model: Exponential & Dose-response model: UME
UME  = readRDS(paste0(loc,"/obesity_results/Time/expModUme.rds"))
exp.NMA.ume.res = residual_draws(UME, newdata = dat) %>%
  group_by(author, trt, week) %>%
  summarise(res = mean(.residual**2/se**2), 
            q25 = quantile(.residual**2/se**2, probs=0.025),
            q975= quantile(.residual**2/se**2, probs=0.975) )
exp.NMA.ume.res$model = "Dose-Response: UME"

# Time-course model: Exponential & Dose-response model: Linear
Lin.lin  = readRDS(paste0(loc,"/obesity_results/DoseTime/exp.NMA.lin.1.rds"))

# Setup design matrix for compound effects
d = unique(dat$compound[dat$trt!="Placebo"]); dname = paste0("d",1:length(d))
for(i in 1:length(d)) dat[,dname[i]] = ifelse(dat$compound == d[i], 1, 0) 

# Derive dose effects, i.e. regression variables
dat$dLira = dat$d1*dat$mdose    
dat$dSema = dat$d2*dat$mdose  

exp.NMA.lin.res = residual_draws(Lin.lin, newdata = dat) %>%
  group_by(author, trt, week) %>%
  summarise(res = mean(.residual**2/se**2), 
            q25 = quantile(.residual**2/se**2, probs=0.025),
            q975= quantile(.residual**2/se**2, probs=0.975) )
exp.NMA.lin.res$model = "Dose-Response: Linear"

# Dose-response: log-linear model

# Derive dose effects, i.e. regression variables
dat$dLira = dat$d1*log(dat$mdose + 1)  
dat$dSema = dat$d2*log(dat$mdose + 1)

Log.log = readRDS(paste0(loc,"/obesity_results/DoseTime/exp.NMA.log.1.rds"))

exp.NMA.log.res = residual_draws(Log.log, newdata = dat) %>%
  group_by(author, trt, week) %>%
  summarise(res = mean(.residual**2/se**2), 
            q25 = quantile(.residual**2/se**2, probs=0.025),
            q975= quantile(.residual**2/se**2, probs=0.975) )
exp.NMA.log.res$model = "Dose-Response: Log-linear"

# Dose-response: Emax model
Emax.log = readRDS(paste0(loc,"/obesity_results/DoseTime/exp.NMA.emax.log.rds"))

exp.NMA.emax.res = residual_draws(Emax.log, newdata = dat) %>%
  group_by(author, trt, week) %>%
  summarise(res = mean(.residual**2/se**2), 
            q25 = quantile(.residual**2/se**2, probs=0.025),
            q975= quantile(.residual**2/se**2, probs=0.975) )
exp.NMA.emax.res$model = "Dose-Response: Emax"

################################################################################
res = rbind(exp.NMA.res,  exp.NMA.ume.res, exp.NMA.log.res, exp.NMA.emax.res)
# Residual Deviance vs time
p = ggplot(res, aes(x=week, y=res, ymin=q25, ymax=q975) ) + facet_wrap( ~ model, nrow = 2, ncol = 2) + 
  geom_point() +
  ylab("Deviance") + xlab("Time (week)") + theme_bw() +
  scale_y_continuous(limits = c(0,50))

p

# Relative to UME model (Figure S6) 
res = rbind(exp.NMA.res, exp.NMA.log.res, exp.NMA.emax.res)
ref = exp.NMA.ume.res; names(ref)[4:7] = paste0(names(ref)[4:7], ".ref") 
res = merge(ref, res, all.y = T)
p = ggplot(res, aes(y=res.ref, x=res)) + facet_wrap( ~ model, nrow = 2, ncol = 2) + 
  geom_point(size=3) + geom_abline(intercept = 0, slope = 1, linewidth=1, colour="grey") + 
  geom_pointrange(aes(ymin = q25, ymax = q975)) +
  geom_pointrange(aes(xmin = q25.ref, xmax=q975.ref)) + 
  ylab("Deviance: inconsistency model") + xlab("Deviance: consistency model") + theme_bw()

p

################################################################################
# Compute ITCs from Time course, and Dose-response Time-course MBNMAs
################################################################################
# Table 7


# Log-linear 
out <- spread_draws(Log.log.ar, b_emaxT_dLira, b_emaxT_dSema, ndraws = NULL) %>%   
  mutate( sema_17vslira30 = b_emaxT_dSema*log(1.7 + 1) - b_emaxT_dLira*log(3 + 1),
          sema_24vslira30 = b_emaxT_dSema*log(2.4 + 1) - b_emaxT_dLira*log(3 + 1),
          sema_17vssema24 = b_emaxT_dSema*log(1.7 + 1) - b_emaxT_dSema*log(2.4 + 1)
  ) %>% 
  dplyr::select( .chain, .iteration, .draw, sema_17vslira30, sema_24vslira30, sema_17vssema24 ) %>%
  pivot_longer(cols=c("sema_17vslira30", "sema_24vslira30", "sema_17vssema24"),  names_to='parameter', values_to='values') %>% 
  group_by(parameter) %>% 
  summarise(Estimate = mean(values), 
            Est.Error= sd(values),       
            Q2.5 = quantile(values, probs = 0.025),
            Q97.5= quantile(values, probs = 0.975) )  
out
out$Method = "Dose-response Time-course NMA: Log-linear"
outres = out

# Emax
out <- spread_draws(Emax.log.ar, b_emaxLira_Intercept, b_ed50Lira_Intercept, b_emaxSema_Intercept, b_ed50Sema_Intercept, ndraws = NULL) %>%   
  mutate( sema_17vslira30 = b_emaxSema_Intercept*1.7 / (exp(b_ed50Sema_Intercept) + 1.7) - b_emaxLira_Intercept*3/(exp(b_ed50Lira_Intercept) + 3),
          sema_24vslira30 = b_emaxSema_Intercept*2.4 / (exp(b_ed50Sema_Intercept) + 2.4) - b_emaxLira_Intercept*3/(exp(b_ed50Lira_Intercept) + 3),
          sema_17vssema24 = b_emaxSema_Intercept*1.7 / (exp(b_ed50Sema_Intercept) + 1.7) - b_emaxSema_Intercept*2.4 / (exp(b_ed50Sema_Intercept) + 2.4)
  ) %>% 
  dplyr::select( .chain, .iteration, .draw, sema_17vslira30, sema_24vslira30, sema_17vssema24 ) %>%
  pivot_longer(cols=c("sema_17vslira30", "sema_24vslira30", "sema_17vssema24"),  names_to='parameter', values_to='values') %>% 
  group_by(parameter) %>% 
  summarise(Estimate = mean(values), 
            Est.Error= sd(values),       
            Q2.5 = quantile(values, probs = 0.025),
            Q97.5= quantile(values, probs = 0.975) )  
out
out$Method = "Dose-response Time-course NMA: Emax"
outres = rbind(outres, out)

out <- spread_draws(NMA.ar, b_emaxT_d4, b_emaxT_d5, b_emaxT_d6, ndraws = NULL) %>%   
  mutate( sema_17vslira30 = b_emaxT_d6 - b_emaxT_d4,
          sema_24vslira30 = b_emaxT_d5 - b_emaxT_d4,
          sema_17vssema24 = b_emaxT_d6 - b_emaxT_d5 ) %>% 
  dplyr::select( .chain, .iteration, .draw, sema_17vslira30, sema_24vslira30, sema_17vssema24 ) %>%
  pivot_longer(cols=c("sema_17vslira30", "sema_24vslira30", "sema_17vssema24"),  names_to='parameter', values_to='values') %>% 
  group_by(parameter) %>% 
  summarise(Estimate = mean(values), 
            Est.Error= sd(values),       
            Q2.5 = quantile(values, probs = 0.025),
            Q97.5= quantile(values, probs = 0.975) )  

out
out$Method = "Time-course NMA"
outres = rbind(outres, out)

################################################################################
# Compute ITCs from Dose-response Time-course MBNMAs RCOV
################################################################################
Log.log.mv  = readRDS(paste0(loc,file = "/obesity_results/DoseTimeRCOV/exp.NMA.log.log.rds"))
summary(Log.log.mv, pars=c("b_emaxT", "b_kT", "sd_1", "Cor_1[2,1]"), probs=c(0.025, 0.975))$summary

# Log-linear 
out <- spread_draws(Log.log.mv, b_emaxT[i], ndraws = NULL) %>%
  filter( i %in% c(8,9)) %>%
  pivot_wider(names_from='i', 
              names_glue = "{i}_{.value}",
              values_from='b_emaxT') %>%      
  rename( b_emaxT_dLira="8_b_emaxT", b_emaxT_dSema = "9_b_emaxT") %>%
  mutate( sema_17vslira30 = b_emaxT_dSema*log(1.7 + 1) - b_emaxT_dLira*log(3 + 1),
          sema_24vslira30 = b_emaxT_dSema*log(2.4 + 1) - b_emaxT_dLira*log(3 + 1),
          sema_17vssema24 = b_emaxT_dSema*log(1.7 + 1) - b_emaxT_dSema*log(2.4 + 1)
  ) %>% 
  dplyr::select( .chain, .iteration, .draw, sema_17vslira30, sema_24vslira30, sema_17vssema24 ) %>%
  pivot_longer(cols=c("sema_17vslira30", "sema_24vslira30", "sema_17vssema24"),  names_to='parameter', values_to='values') %>% 
  group_by(parameter) %>% 
  summarise(Estimate = mean(values), 
            Est.Error= sd(values),       
            Q2.5 = quantile(values, probs = 0.025),
            Q97.5= quantile(values, probs = 0.975) )  
out
out$Method = "Dose-response Time-course NMA: Log-linear, observed corr"
outres = rbind(outres, out)
outres
# 
# Standard NMA results
#
nma.p = readRDS(paste0(loc, "/obesity_results/NMA/NMAp.rds"))

x = as.data.frame(relative_effects(nma.p, trt_ref = "lira 3.0 mg", probs = c(0.025, 0.975)))
y = as.data.frame(relative_effects(nma.p, trt_ref = "sema 2.4 mg", probs = c(0.025, 0.975)))
out.NMA = rbind(x[2:3,1:7], y[3,1:7])
out.NMA$Method = "Standard NMA (P)" 
out.NMA$ITC    = paste(out.NMA$.trtb, " vs ", out.NMA$.trta)
out.NMA        = out.NMA[,4:9]
out.NMA$Estimate = out.NMA$Estimate
out.NMA$Q2.5  = out.NMA$Q2.5
out.NMA$Q97.5 = out.NMA$Q97.5
#

nma.e = readRDS(paste0(loc, "/obesity_results/NMA/NMAe.rds"))

x = as.data.frame(relative_effects(nma.e, trt_ref = "lira 3.0 mg", probs = c(0.025, 0.975)))
y = as.data.frame(relative_effects(nma.e, trt_ref = "sema 2.4 mg", probs = c(0.025, 0.975)))
out = rbind(x[2:3,1:7], y[3,1:7])
out$Method = "Standard NMA (E)" 
out$ITC    = paste(out$.trtb, " vs ", out$.trta)
out        = out[,4:9]
out$Estimate = out$Estimate
out$Q2.5  = out$Q2.5
out$Q97.5 = out$Q97.5
#
out.NMA = rbind(out.NMA, out)
#
out.NMA
