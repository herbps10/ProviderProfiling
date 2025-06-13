#
# Simulation Study 2
#
library(tidyverse)
library(TargetedRisk)
library(mlr3extralearners)

source("simulation_study_2/simulate.R")
source("R/glm.R")

wrapper <- function(index, seed, N, m, parameter, P, outcome_type, data) {
  # Check cache
  cache <- glue::glue("{cache_path}/{N}/{parameter}_{seed}.rds")
  if(file.exists(cache)) {
    return(read_rds(cache))
  }

  
  # Calculate true parameter values
  largedat <- simulate_data(N = 1e6, m = m, seed = seed, P = P, outcome_type = outcome_type)
  truth <- largedat |> group_by(A) |> 
    summarize(psi1 = mean(gamma1), 
              psi2 = mean(gamma2),
        n = n()) |>
    ungroup() |>
    mutate(a = as.character(A)) |>
    mutate(smr = psi1 / psi2,
           er = psi1 - psi2)

  truth_direct <- largedat |>
    select(starts_with("Qbar")) |>
    summarize_all(mean) |>
    pivot_longer(everything()) |>
    mutate(a = str_replace(name, "Qbar", "")) |>
    select(-name) |>
    rename(direct = value)

  rm(largedat)

  # Simulate incorrect nuisance parameter estimates
  if(parameter == "direct") {
    Qtilde <- matrix(0.5, ncol = m, nrow = N)
    colnames(Qtilde) <- 1:m
  }
  else {
    Qtilde <- matrix(0.5, ncol = 1, nrow = N)
  }
  g <- matrix(runif(m * N, 0.1, 0.9), ncol = m, nrow = N)
  g <- g / matrix(rowSums(g), ncol = m, nrow = N, byrow = FALSE)
  colnames(g) <- 1:m

# SuperLearner algorithms
  correct_outcome <- c("mean", "glm", "glmnet", "ranger", "knn")
  correct_trt <- c("mean", "glm", "glmnet", "ranger")

  learners_outcome <- correct_outcome
  learners_trt <- correct_trt


  # Set up estimators
  if(parameter == "direct") {
    f_tmle <- direct_tmle
    f_onestep <- direct_onestep
  }
  else {
    f_tmle <- indirect_tmle
    f_onestep <- indirect_onestep
  }

  baseline <- paste0("W", 1:P)
  folds <- 5

  start_time <- Sys.time()

  # Scenario 1
  fit1_tmle <- f_tmle(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    outcome_type = outcome_type,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    folds = folds,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  ) 

  fit1_onestep <- f_onestep(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    outcome_type = outcome_type,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    folds = folds,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )

  # Scenario 2
  fit2_tmle <- f_tmle(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    outcome_type = outcome_type,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    g = g,
    Qtilde = fit1_tmle$Qtilde$predicted_outcomes,
    folds = folds,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )

  fit2_onestep <- f_onestep(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    outcome_type = outcome_type,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    g = g,
    Qtilde = fit1_tmle$Qtilde$predicted_outcomes,
    folds = folds,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )

  # Scenario 3
  fit3_tmle <- f_tmle(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    outcome_type = outcome_type,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    g = fit1_tmle$g$treatment_probs,
    Qtilde = Qtilde,
    folds = folds,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )

  fit3_onestep <- f_onestep(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    outcome_type = outcome_type,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    g = fit1_tmle$g$treatment_probs,
    Qtilde = Qtilde,
    folds = folds,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )

  # Scenario 4
  fit4_tmle <- f_tmle(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    outcome_type = outcome_type,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    g = g,
    Qtilde = Qtilde,
    folds = folds,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )

  fit4_onestep <- f_onestep(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    outcome_type = outcome_type,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    g = g,
    Qtilde = Qtilde,
    folds = folds,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )

  # Fit GLM and WeightIt based estimators
  tidy_fit_weightit <- NULL
  if(parameter == "direct") {
    fit_glm <- direct_glm(data, "Y", "A", baseline) |>
      rename(trt = a, estimate = glm_direct) |>
      mutate(estimator = "GLM", trt = as.character(trt), parameter = "direct")

    fit_weightit <- direct_weightit(
      data,
      outcome = "Y",
      trt = "A",
      baseline = baseline,
      method = "ebal"
    )

    tidy_fit_weightit <- tidy(fit_weightit)
  }
  else {
    fit_glm <- indirect_glm(data, "Y", "A", baseline) |>
      rename(trt = a) |>
      mutate(estimator = "GLM", trt = as.character(trt)) |>
      pivot_longer(c(glm_psi1, glm_psi2, glm_ER, glm_SMR)) |>
      mutate(parameter = str_replace(name, "glm_", ""), estimate = value) |>
      select(-name, -value)
  }

  # Collate results
  res <- tibble(
    index = rep(index, 10),
    N = rep(N, 10),
    seed = rep(seed, 10), 
    scenario = c(1:4, 1:4, 1, 1),
    outcome_type = rep(outcome_type, 10),
    parameter = rep(parameter, 10),
    fit = list(tidy(fit1_tmle), tidy(fit2_tmle), tidy(fit3_tmle), tidy(fit4_tmle), tidy(fit1_onestep), tidy(fit2_onestep), tidy(fit3_onestep), tidy(fit4_onestep), fit_glm, tidy_fit_weightit), 
    truth = list(truth),
    truth_direct = list(truth_direct),
    m = rep(m, 10)
  )

  # Write to cache
  write_rds(res, cache)

  return(res)
}

cache_path <- Sys.getenv("SIMULATION_CACHE_PATH")
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

if(cache_path == "") stop("Please set SIMULATION_CACHE_PATH environment variable.")
if(task_id == "") stop("Task id not set. Please set SLURM_ARRAY_TASK_ID, or run simulations through a Slurm job array, which will set this environment variable for you.")

# Set up simulations
N_simulations <- 1
simulations <- expand_grid(
  index = as.numeric(task_id),
  m = 75,
  #N = c(5e3, 1e4, 2e4),
  N = 5e3,
  outcome_type = c("binomial", "continuous"),
  parameter = c("direct", "smr")
) 

# Make sure cache folders exist
for(N in unique(simulations$N)) {
  dir <- glue::glue("{cache_path}/{N}")
  if(!dir.exists(dir)) {
    dir.create(dir)
  }
}

# Run simulations
simulations <- simulations |>
  mutate(
    seed = index * 1e5 + N
  ) |>
  mutate(
    data = pmap(list(seed, N, m, outcome_type), simulate_data, P = 10),
    fits = pmap(list(index, seed, N, m, parameter, outcome_type, data), wrapper, P = 10)
  )
