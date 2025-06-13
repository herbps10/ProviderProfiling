# Replication Materials for Doubly Robust Nonparametric Efficient Estimation for Healthcare Provider Evaluation
*Authors: Herbert P. Susmann, Yiting Li, Mara A. McAdams-DeMarco, Iván Díaz, and Wenbo Wu*

This repository includes the R code used for the two simulation studies in the paper _Doubly Robust Nonparametric Efficient Estimation for Healthcare Provider Evaluation_. 

The simulations depend on the [`TargetedRisk`](https://github.com/herbps10/TargetedRisk). This package can be installed directly from Github:
```{r}
remotes::install_github("herbps10/TargetedRisk")
```

Code for the two simulation studies included in the paper are divided into the folders `simulation_study_1` and `simulation_study_2`, respectively. Each folder contains the following files:
- `simulate_data.R`: contains the function `simulate_data` for generating a dataset from the simulation data-generating process.
- `simulation_study.R`: runs a single replicate of the simulation study (see below for further details).
- `collect_results.R`: collates cached results from the simulation study into a single results file.
- `analyze_results.R`: calculates error and coverage metrics and generates results tables.

The additional file `R/glm.R` provides a helper function for fitting the comparison `GLM` based method for direct and indirect standardization.

The simulation study code is written for execution on a parallel computing cluster using [Slurm](https://en.wikipedia.org/wiki/Slurm_Workload_Manager). Slurm [job arrays](https://slurm.schedmd.com/job_array.html) are used to generate jobs, with each job in the array executing one replicate of the simulation study. Slurm provides an environment variable to running jobs, `SLURM_ARRAY_TASK_ID`, specifying the job array index.

The simulation code caches intermediate results to a temporary directory specified in the environment variable `SIMULATION_CACHE_PATH`. Collated results are saved to the path in the environment variable `SIMULATION_RESULT_PATH`.

If you would like to run the simulation study outside of a Slurm environment, you first need to set several environment variables. For example, to run one replicate of the first simulation study in `R`:
```{r}
# Run a single replicate 
Sys.setenv(SLURM_ARRAY_TASK_ID = 1)

# Set cache path
Sys.setenv(SIMULATION_CACHE_PATH = "/path/to/cache")

# Set results path
Sys.setenv(SIMULATION_RESULTS_PATH = "/path/to/results")

source("R/simulation_study_1/simulation_study.R")
```

This will execute one replicate of the first simulation study, and save the results to the cache. You can then collect all cached results into a single results file by running the `R/collect_results.R` script(make sure the same environment variables as before are set):
```{r}
source("R/simulation_study_1/collect_results.R")
```

Then, to analyze the results and generate results tables, run
```{r}
source("R/simulation_study_1/analyze_simulation_study.R")
```
