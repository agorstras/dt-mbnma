# Model-Based Network Meta-Analysis (MB-NMA) Repository

This repository contains code and data for two separate lines of work on model-based network meta-analysis with joint dose-response and time-course modeling:

1. **Simulation Study** (`_sim_`): Simulation-based evaluation of MB-NMA methods
2. **Obesity Case Study** (`_obesity_`): Application of MB-NMA to weight management treatment data

## Repository Structure

The repository is organized with a modular approach where each task follows a numbered workflow and has dedicated subdirectories for data, Stan models, and results. 

### Simulation Study

**Main Scripts** (prefix: `0_sim_` through `4_sim_`)
- `0_sim_generate_data.R` - Generate 100 replicates of simulation data (10 trials each) with dose-response trajectories. Also generates Figures S1 and S2 (treatment network and mean trajectories/dose-response profiles)
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

### Obesity Case Study

**Main Scripts** (prefix: `00_obesity_` through `05b_obesity_`)
- `00_obesity_setup.R` - Setup directory structure 
- `01_obesity_NMA.R` - Network meta-analysis 
- `02_obesity_Time_course_NMA.R` - Time course network meta-analysis
- `03_obesity_Dose_Time_Course_NMA.R` - Combined dose-time course NMA
- `04_obesity_Dose_Time_Course_NMA_RCOV.R` - Dose-time course NMA with covariance structure based on observed correlations
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
4. Extend to cover observed correlations
5. Generate results and figures

## Dependencies

Key R packages required:
- `brms` - Bayesian regression models with Stan
- `tidyverse` - Data manipulation and visualization
- `tidybayes` - Bayesian analysis utilities
- `MASS` - Multivariate normal distributions
- `Matrix` - Sparse matrix operations
- `multinma` - Network meta-analysis
- `rstan` / `cmdstanr` - Stan interface

### Installation Notes

The `cmdstanr` package requires special installation. After installing the package, you need to install the CmdStan toolchain:

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
library(cmdstanr)
check_cmdstan_toolchain()
install_cmdstan(cores = 4)
cmdstan_path()
cmdstan_version()
```

## Key Features

- **Model-based approach**: Incorporates dose-response and time course relationships
- **Multivariate modeling**: Accounts for within-study correlations across time points
- **Bayesian framework**: Uses Stan for flexible hierarchical modeling
- **Simulation validation**: Extensive simulation study to evaluate method performance
- **Real-world application**: Applied to weight management treatments trial network

## Notes

- All scripts that perform random sampling (model estimation, MCMC) initialize a fixed seed (1234) for reproducibility. This ensures deterministic results across runs.
- File paths in scripts use relative paths (`sim` and `obesity` directories) for cross-machine reproducibility
- Parallel processing enabled where appropriate (`mc.cores`, `future.apply`)

## Contact

Please contact the authors: Anders Strathe (aqss@novonordisk.com), Martin Bøg (axbq@novonordisk.com), or Anders Gorst-Rasmussen (agtr@novonordisk.com) in case of questions or comments.

---

**Last Updated**: 27 March 2026
