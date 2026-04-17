################################################################################
# Description: Network meta-analysis, on obesity data set   
# - Figure 1A, 1B (main manuscript)
# - saves results of analyses to ./obesity_results
################################################################################
#
rm(list = ls())

library(multinma)
library(stringr)
library(ggplot2)
library(brms)
library(Matrix)
library(tidybayes)

#set location
loc = "obesity"

set.seed(1234)
# load data
datNMAe = readRDS(paste0(loc, "/obesity_data/datNMAeotobs.rds"))
datNMAp = readRDS(paste0(loc, "/obesity_data/datNMAprimary.rds"))

################################################################################
# Use multinma to visualise full network

arm_net <- set_agd_arm(datNMAp,
                       study = author,
                       trt = treatment,
                       y = y,
                       se = y.se,
                       sample_size = randomised,
                       trt_ref = "placebo")

# Figure 1A
plot(arm_net, weight_edges = TRUE, weight_nodes = TRUE)


################################################################################
# Use multinma to fit model on phase 3 data only

datNMAp = subset(datNMAp, !(author %in% c("Astrup", "O'Neil")))

arm_net <- set_agd_arm(datNMAp,
                       study = author,
                       trt = treatment,
                       y = y,
                       se = y.se,
                       sample_size = randomised,
                       trt_ref = "placebo")

# Figure 1B
plot(arm_net, weight_edges = TRUE, weight_nodes = TRUE)

# fit random effects NMA; note some divergent transitions 
arm_fit_RE <- nma(arm_net,
                  trt_effects = "random",
                  prior_intercept = normal(scale = 100),
                  prior_trt = normal(scale = 100),
                  prior_het = half_normal(scale = 5),
                  control=list(adapt_delta = 0.9999, max_treedepth = 25), iter = 5000, cores = 4, seed=1234)

summary(arm_fit_RE)

saveRDS(arm_fit_RE, paste0(loc, "/obesity_results/NMA/NMAp.rds"))

# end of trial, observed data
arm_net <- set_agd_arm(datNMAe,
                       study = author,
                       trt = trt,
                       y = mean,
                       se = mean.se,
                       sample_size = N,
                       trt_ref = "placebo")

arm_fit_RE <- nma(arm_net,
                  trt_effects = "random",
                  prior_intercept = normal(scale = 100),
                  prior_trt = normal(scale = 100),
                  prior_het = half_normal(scale = 5),
                  control=list(adapt_delta = 0.9999, max_treedepth = 25), iter = 5000, cores = 4, seed=1234)

summary(arm_fit_RE)

saveRDS(arm_fit_RE, paste0(loc, "/obesity_results/NMA/NMAe.rds"))