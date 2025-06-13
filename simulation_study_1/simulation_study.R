#
# Simulation Study 1
#

library(tidyverse)
library(TargetedRisk)
library(mlr3extralearners)

source("simulation_study_1/simulate.R")
source("R/glm.R")

wrapper <- function(index, seed, N, m, parameter, P, outcome_type, data) {
  # Check cache
  cache <- glue::glue("{cache_path}/{N}/{parameter}_{outcome_type}_{seed}.rds")
  if(file.exists(cache)) {
    return(read_rds(cache))
  }

  # Calculate true parameter values
  largedat <- simulate_data(N = 5e6, m = m, seed = seed, P = P, outcome_type = outcome_type)
  truth <- largedat %>% group_by(A) %>% 
    summarize(psi1 = mean(gamma1), 
              psi2 = mean(gamma2),
        n = n()) %>%
    ungroup() %>%
    mutate(a = as.character(A)) %>%
    mutate(smr = psi1 / psi2,
           er = psi1 - psi2)

  truth_direct <- data %>%
    distinct(A, hospital_effect1) %>%
    arrange(A) %>%
    mutate(
      direct = case_when(
        outcome_type == "continuous" & hospital_effect1 == 1 ~ 0.87,
        outcome_type == "continuous" & hospital_effect1 == 0 ~ 0.43,
        outcome_type == "binomial" & hospital_effect1 == 1 ~ 0.4 * plogis(0.3) + 0.3 * plogis(1) + 0.3 * plogis(1.5),
        outcome_type == "binomial" & hospital_effect1 == 0 ~ 0.4 * plogis(0.7) + 0.3 * plogis(0.5) + 0.3 * plogis(0)
      )
    )

  rm(largedat)

  # Set up algorithms
  if(parameter == "direct") {
    f_tmle <- direct_tmle
    f_onestep <- direct_onestep
  }
  else {
    f_tmle <- indirect_tmle
    f_onestep <- indirect_onestep
  }

  correct_outcome <- list(list("ranger", num.trees = 100, max.depth = 2), list("ranger", num.trees = 50, max.depth = 2))
  correct_trt <- list(list("ranger", num.trees = 100, max.depth = 2), list("ranger", num.trees = 50, max.depth = 2))

  learners_outcome <- correct_outcome
  learners_trt <- correct_trt

  baseline <- c("W1")
  folds <- 5

  start_time <- Sys.time()

  # TMLE
  fit_tmle <- f_tmle(
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

  fit_onestep <- f_onestep(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    outcome_type = outcome_type,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    Qtilde = fit_tmle$Qtilde$predicted_outcomes,
    g = fit_tmle$g$treatment_probs,
    folds = folds,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )

  fit_weightit <- NULL
  if(parameter == "direct") {
    fit_weightit <- direct_weightit(
      data,
      outcome = "Y",
      trt = "A",
      baseline = baseline,
      method = "ebal"
    )
  }

  # One-step

  # Fit GLM model
  if(parameter == "direct") {
    fit_glm <- direct_glm(data, "Y", "A", baseline, outcome_type)
  }
  else {
    fit_glm <- indirect_glm(data, "Y", "A", baseline, outcome_type)
  }

  tidy_weightit <- NULL
  if(parameter == "direct") {
    tidy_weightit <- tidy(fit_weightit)
  }

  res <- tibble(
    index = index,
    N = N,
    seed = seed,
    parameter = parameter,
    tmle = list(tidy(fit_tmle)),
    onestep = list(tidy(fit_onestep)),
    weightit = list(tidy_weightit),
    glm = list(fit_glm),
    truth = list(truth),
    truth_direct = list(truth_direct),
    m = m,
    outcome_type = outcome_type
  )

  write_rds(res, cache)

  return(res)
}

cache_path <- Sys.getenv("SIMULATION_CACHE_PATH")
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

if(cache_path == "") stop("Please set SIMULATION_CACHE_PATH environment variable.")
if(task_id == "") stop("Task id not set. Please set SLURM_ARRAY_TASK_ID, or run simulations through a Slurm job array, which will set this environment variable for you.")

N_simulations <- 1
simulations <- expand_grid(
  index = as.numeric(task_id),
  m = 10,
  N = c(1e3, 2.5e3, 5e3),
  outcome_type = c("binomial", "continuous"),
  parameter = c("direct", "smr")
) 

for(N in unique(simulations$N)) {
  dir <- glue::glue("{cache_path}/{N}")
  if(!dir.exists(dir)) {
    dir.create(dir)
  }
}

simulations <- simulations %>%
  mutate(
    seed = index * 1e5 + N
  ) %>%
  mutate(
    data = pmap(list(seed, N, m, outcome_type), simulate_data, P = 1),
    fits = pmap(list(index, seed, N, m, parameter, outcome_type, data), wrapper, P = 1)
  )
