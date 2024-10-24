dyn.load("/gpfs/share/apps/gcc/11.2.0/lib64/libstdc++.so.6")
library(torch)

library(tidyverse)
library(TargetedSMR)
library(SuperRiesz)
library(mlr3extralearners)

source("R/simulate.R")
source("R/glm.R")

wrapper <- function(index, seed, N, m, estimator, trt_method, parameter, P, data) {
  cache <- glue::glue("/gpfs/scratch/susmah01/ProviderProfiling/cache/{N}/{parameter}_{trt_method}_{estimator}_{seed}.rds")

  if(file.exists(cache)) {
    return(read_rds(cache))
  }

  #correct_outcome <- list(list("knn", k = 10), list("knn", k = 25), list("knn", k = 50), list("knn", k = 100))
  #correct_outcome <- list(list("lightgbm", num_iterations = 5), list("lightgbm", num_iterations = 10), list("lightgbm", num_iterations = 25), list("lightgbm", num_iterations = 50), list("lightgbm", num_iterations = 100), list("lightgbm", num_iterations = 200))
  correct_outcome <- list(list("ranger", num.trees = 100, max.depth = 2), list("ranger", num.trees = 50, max.depth = 2))

  if(trt_method == "SuperRiesz") {
    correct_trt <- list(list("nn",
      layers = 2,
      hidden = 20,
      learning_rate = 0.01,
      constrain_positive = TRUE,
      epochs = 250,
      dropout = 0.05,
      verbose = TRUE
    ))
  }
  else if(trt_method == "torch") {
    correct_trt = "torch"
  }
  else {
    #correct_trt <- list("randomforest", list("knn", k = 5), list("knn", k = 7), list("knn", k = 9))
    #correct_trt <- list(list("knn", k = 10), list("knn", k = 25), list("knn", k = 50), list("knn", k = 100))
    #correct_trt <- list(list("lightgbm", num_iterations = 5), list("lightgbm", num_iterations = 10), list("lightgbm", num_iterations = 25), list("lightgbm", num_iterations = 50), list("lightgbm", num_iterations = 100), list("lightgbm", num_iterations = 200))
    correct_trt <- list(list("ranger", num.trees = 100, max.depth = 2), list("ranger", num.trees = 50, max.depth = 2))
  }

  learners_outcome <- correct_outcome
  learners_trt <- correct_trt

  file <- glue::glue("/gpfs/scratch/susmah01/ProviderProfiling/logs/log-{Sys.getpid()}.txt")
  write(glue::glue("{Sys.time()}: starting {estimator} trt_method={trt_method} parameter={parameter} index={index} seed={seed} N={N}"), file = file, append = TRUE)
  start <- Sys.time()

  if(parameter == "direct") {
    f <- direct_tmle
  }
  else {
    f <- indirect_tmle
  }


  largedat <- simulate_data(N = 5e6, m = m, seed = seed, P = P)
  truth <- largedat %>% group_by(A) %>% 
    summarize(psi1 = mean(gamma1), 
              psi2 = mean(gamma2),
        n = n()) %>%
    ungroup() %>%
    mutate(a = as.character(A)) %>%
    mutate(smr = psi1 / psi2,
           er = psi1 - psi2)

  #truth_direct <- largedat %>%
  #  select(starts_with("Qbar")) %>%
  #  summarize_all(mean) %>%
  #  pivot_longer(everything()) %>%
  #  mutate(a = str_replace(name, "Qbar", "")) %>%
  #  select(-name) %>%
  #  rename(direct = value)

  truth_direct <- data %>%
    distinct(A, hospital_effect1) %>%
    arrange(A) %>%
    mutate(direct = ifelse(hospital_effect1 == 1, 0.87, 0.43))

  rm(largedat)

  #baseline <- c("W1", "W2", "W3", "W4", "W5")
  baseline <- c("W1")

  start_time <- Sys.time()

  g <- data[, paste0("g", 1:m)]
  colnames(g) <- 1:m

  folds <- 5

  fit1 <- f(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    trt_method = trt_method,
    outcome_type = "continuous",
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    #g = g,
    folds = folds,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  ) 
  elapsed <- Sys.time() - start_time

  #Qtilde <- matrix(runif(N), ncol = 1, nrow = N)
  if(parameter == "direct") {
    Qtilde <- matrix(0.5, ncol = m, nrow = N)
    colnames(Qtilde) <- 1:m
  }
  else {
    Qtilde <- matrix(0.5, ncol = 1, nrow = N)
  }
  g <- matrix(runif(m * N), ncol = m, nrow = N)
  g <- g / matrix(rowSums(g), ncol = m, nrow = N, byrow = FALSE)
  colnames(g) <- 1:m

  fit2 <- f(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    trt_method = trt_method,
    outcome_type = "continuous",
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    g = g,
    Qtilde = fit1$Qtilde$predicted_outcomes,
    folds = folds,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )

  fit3 <- f(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    trt_method = trt_method,
    outcome_type = "continuous",
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    g = fit1$g$treatment_probs,
    Qtilde = Qtilde,
    folds = folds,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )

  fit4 <- f(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    trt_method = trt_method,
    outcome_type = "continuous",
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    g = g,
    Qtilde = Qtilde,
    folds = folds,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )


  if(parameter == "direct") {
    fit_glm <- direct_glm(data, "Y", "A", baseline)
  }
  else {
    fit_glm <- indirect_glm(data, "Y", "A", baseline)
  }

  write(glue::glue("{Sys.time()}: ending {estimator} trt_method={trt_method} parameter={parameter} index={index} seed={seed} N={N} duration={hms::as.hms(Sys.time() - start)}"), file = file, append = TRUE)

  res <- tibble(
    index = rep(index, 4),
    N = rep(N, 4),
    estimator = rep(estimator, 4),
    seed = rep(seed, 4), 
    scenario = 1:4,
    trt_method = rep(trt_method, 4),
    parameter = rep(parameter, 4),
    fit = list(fit1, fit2, fit3, fit4), 
    fit_glm = list(fit_glm),
    truth = list(truth),
    truth_direct = list(truth_direct),
    m = rep(m, 4),
    duration = rep(elapsed, 4)
  )

  write_rds(res, cache)

  return(res)
}

N_simulations <- 1
simulations <- expand_grid(
  index = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")),
  m = 10,
  N = c(1e3, 2.5e3, 5e3),
  trt_method = "default",
  parameter = c("direct", "smr")
) 

for(N in unique(simulations$N)) {
  dir <- glue::glue("/gpfs/scratch/susmah01/ProviderProfiling/cache/{N}")
  if(!dir.exists(dir)) {
    dir.create(dir)
  }
}

simulations <- simulations %>%
  mutate(
    seed = index * 1e5 + N
  ) %>%
  mutate(
    data = pmap(list(seed, N, m), simulate_data, P = 1),
    fit_tmle = pmap(list(index, seed, N, m, trt_method, parameter, data), wrapper, P = 1, estimator = "TMLE")
  )
