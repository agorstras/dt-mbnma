# Model-Based Network Meta-Analysis (MB-NMA) Repository

This repository contains code and data for two separate lines of work on model-based network meta-analysis:

1. **Simulation Study** (`_sim_`): Simulation-based evaluation of MB-NMA methods
2. **Obesity Case Study** (`_obesity_`): Application of MB-NMA to weight management treatment data

## Repository Structure

The repository is organized with a modular approach where each task follows a numbered workflow and has dedicated subdirectories for data, Stan models, and results. 

### Simulation Study

**Main Scripts** (prefix: `0_sim_` through `4_sim_`)
- `0_sim_generate_data.R` - Generate 100 replicates of simulation data (10 trials each) with dose-response trajectories
- `1_sim_univariate_ll.R` - Univariate likelihood analysis
- `2_sim_multivariate_AR1.R` - Multivariate analysis with AR1 correlation structure
- `3_sim_multivariate_obs.R` - Multivariate analysis with observed correlation structure
- `4_sim_plot_results.R` - Visualization of simulation results

**Subdirectories:**
- `sim_data/` - Generated simulation datasets (`.rds` files)
  - `dataLst.rds`, `varLst.rds` - Main simulation data and variance lists
  - Various temporary lists for different analyses
- `sim_design/` - Simulation design specifications
- `sim_stan/` - Stan model files for simulation analyses
  - `exp_NMA_emax_log.stan`, `exp_NMA_emax_log2.stan`, `exp_NMA_emax_log3.stan`
  - `exp_Time_NMA.stan`
- `sim_results/` - Analysis outputs and results

### Obesity Case Study

**Main Scripts** (prefix: `01_obesity_` through `05b_obesity_`)
- `01_obesity_NMA.R` - Network meta-analysis (Figures 1A, 1B)
- `02_obesity_Time_course_NMA.R` - Time course network meta-analysis
- `03_obesity_Dose_Time_Course_NMA.R` - Combined dose-time course NMA
- `04_obesity_Dose_Time_Course_NMA_RCOV.R` - Dose-time course NMA with random covariance structure
- `05a_obesity_results_Time.R` - Time-specific results
- `05b_obesity_results_Dose_Time.R` - Dose-time specific results

**Subdirectories:**
- `obesity_data/` - Obesity treatment datasets (`.rds` files)
  - `datNMAeotobs.rds` - End of trial observed data
  - `datNMAprimary.rds` - Primary analysis data
  - `datDoseTime.rds`, `datDoseTimeRCOV.rds` - Dose-time data
  - `datTime.rds` - Time course data
  - `v_mat.rds` - Variance-covariance matrix
- `obesity_stan/` - Stan model files for obesity analyses
  - `exp_NMA_emax.stan`
  - `exp_NMA_log.stan`
- `obesity_results/` - Analysis outputs and results

## Workflow

### Simulation Study Workflow
1. Generate simulated trial data with known dose-response relationships
2. Fit univariate models to establish baseline
3. Fit multivariate models with different correlation structures (AR1 vs observed)
4. Compare and visualize model performance

### Obesity Case Study Workflow
1. Conduct standard network meta-analysis on phase 3 data
2. Extend to time course analysis
3. Incorporate dose-response relationships
4. Refine with random covariance structures
5. Generate comprehensive results and figures

## Dependencies

Key R packages required:
- `multinma` - Network meta-analysis
- `brms` - Bayesian regression models with Stan
- `tidyverse` - Data manipulation and visualization
- `tidybayes` - Bayesian analysis utilities
- `MASS` - Multivariate normal distributions
- `Matrix` - Sparse matrix operations
- `rstan` / `cmdstanr` - Stan interface

## Key Features

- **Model-based approach**: Incorporates dose-response and time course relationships
- **Multivariate modeling**: Accounts for within-study correlations across time points
- **Bayesian framework**: Uses Stan for flexible hierarchical modeling
- **Simulation validation**: Extensive simulation study to evaluate method performance
- **Real-world application**: Applied to obesity treatment trial network

## Notes

- File paths in scripts should be updated to reflect your directory structure
- Initial runs supporting the publication are already saved in the /results folders, such that it is possible to execute script 4 (sim), and 5a/5b (obesity) without rerunning analyses 
- Parallel processing enabled where appropriate (`mc.cores`, `future.apply`)
- Reproducible results in simulation scripts

## Contact

Please contact the authors: Anders Strathe (aqss@novonordisk.com), Martin Bøg (axbq@novonordisk.com), or Anders Gorst-Rasmussen (agtr@novonordisk.com) in case of questions or comments.

---

**Last Updated**: January 2026
