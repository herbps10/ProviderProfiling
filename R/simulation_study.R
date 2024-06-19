dyn.load("/gpfs/share/apps/gcc/11.2.0/lib64/libstdc++.so.6")
library(torch)

library(tidyverse)
library(TargetedSMR)
library(SuperRiesz)

source("R/simulate.R")

wrapper <- function(index, seed, N, m, estimator, trt_method, parameter, data) {
  cache <- glue::glue("/gpfs/scratch/susmah01/ProviderProfiling/cache/{N}/{parameter}_{trt_method}_{estimator}_{seed}.rds")

  if(file.exists(cache)) {
    return(read_rds(cache))
  }

  m <- length(unique(data$A))

  correct_outcome <- c("mean", "glm", "glmnet", "ranger")

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
    correct_trt <- c("mean", "glm", "nnet", "glmnet", "ranger")
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

  baseline <- c("W1", "W2", "W3", "W4", "W5")

  start_time <- Sys.time()

  fit1 <- f(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    trt_method = trt_method,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    folds = 5,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  ) 
  elapsed <- Sys.time() - start_time

  #Qtilde <- matrix(runif(N), ncol = 1, nrow = N)
  Qtilde <- matrix(0.5, ncol = 1, nrow = N)
  g <- matrix(runif(m * N), ncol = m, nrow = N)
  g <- g / matrix(rowSums(g), ncol = m, nrow = N, byrow = FALSE)
  colnames(g) <- 1:m

  fit2 <- f(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    trt_method = trt_method,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    g = g,
    Qtilde = fit1$Qtilde$predicted_outcomes,
    folds = 5,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )

  fit3 <- f(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    trt_method = trt_method,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    g = NULL,
    Qtilde = Qtilde,
    folds = 5,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )

  fit4 <- f(
    data,
    outcome = "Y",
    trt = "A",
    baseline = baseline,
    trt_method = trt_method,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    g = g,
    Qtilde = Qtilde,
    folds = 5,
    control = standardization_control(.learners_trt_folds = 5, .learners_outcome_folds = 5)
  )


  largedat <- simulate_data(1e6, m = m, seed = seed)
  truth <- largedat %>% group_by(A) %>% 
    summarize(psi1 = mean(Y), 
              psi2 = mean(YA),
        n = n()) %>%
    ungroup() %>%
    mutate(a = as.character(A)) %>%
    mutate(smr = psi1 / psi2,
           er = psi1 - psi2)

  truth_direct <- largedat %>%
    select(starts_with("Y")) %>%
    select(-Y, -YA) %>%
    summarize_all(mean) %>%
    pivot_longer(everything()) %>%
    mutate(a = str_replace(name, "Y", "")) %>%
    select(-name) %>%
    rename(direct = value)


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
  m = 25,
  N = c(2500, 5e3, 1e4),
  trt_method = "torch",
  #parameter = c("direct", "smr")
  parameter = "smr"
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
    data = pmap(list(seed, N, m), simulate_data),
    fit_tmle = pmap(list(index, seed, N, m, trt_method, parameter, data), wrapper, estimator = "TMLE")
  )
