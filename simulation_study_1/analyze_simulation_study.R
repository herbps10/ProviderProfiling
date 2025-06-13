library(tidyverse)
library(TargetedRisk)
library(broom)

results_path <- Sys.getenv("SIMULATION_RESULTS_PATH")
if(results_path == "") stop("Please set SIMULATION_RESULTS_PATH environment variable.")

simulations <- read_rds(glue::glue("{results_path}/simulation_results.rds"))

truth_indirect <- simulations |>
  filter(parameter == "smr") |>
  select(index, N, truth, outcome_type) |> 
  unnest(c(truth)) |> 
  select(index, N, A, outcome_type, psi1, psi2, smr, er) |>
  mutate(A = as.character(A)) |>
  distinct()

truth_direct <- simulations |>
  select(index, N, truth_direct, outcome_type) |> 
  unnest(c(truth_direct)) |> 
  select(index, N, A, direct, outcome_type) |>
  mutate(A = as.character(A)) |>
  distinct()

glm_estimates <- simulations |>
  select(index, N, outcome_type, glm) |> 
  unnest(glm) |> 
  mutate(a = as.character(a)) |> 
  pivot_longer(c(glm_direct, glm_psi1, glm_psi2, glm_SMR, glm_ER)) |> 
  #pivot_longer(c(glm_psi1, glm_psi2, glm_SMR, glm_ER)) |> 
  mutate(name = str_replace(name, "glm_", "")) |>
  filter(!is.na(value)) |>
  rename(estimate = value, trt = a, parameter = name) |>
  mutate(estimator = "GLM")

tmle <- simulations |>
  ungroup() |>
  select(index, N, tmle, outcome_type) |>
  unnest(cols = c("tmle")) |>
  mutate(estimator = "TMLE")

onestep <- simulations |>
  ungroup() |>
  select(index, N, onestep, outcome_type) |>
  unnest(cols = c("onestep")) |>
  mutate(estimator = "Onestep")

weightit <- simulations |>
  ungroup() |>
  select(index, N, weightit, outcome_type) |>
  unnest(cols = c("weightit")) |>
  mutate(estimator = "WeightIt")

simulation_results <-  bind_rows(tmle, onestep, weightit, glm_estimates) |>
  left_join(truth_indirect, by = c("index", "N", "outcome_type", trt = "A")) |>
  left_join(truth_direct, by = c("index", "N", "outcome_type", trt = "A")) |>
  mutate(truth = case_when(
    parameter == "psi1" ~ psi1,
    parameter == "psi2" ~ psi2,
    parameter == "SMR" ~ smr,
    parameter == "ER" ~ er,
    parameter == "direct" ~ direct
  )) |>
  mutate(
    covered = round(conf.low, 2) <= round(truth, 2) & round(conf.high, 2) >= round(truth, 2), error = truth - estimate,
    covered = conf.low < truth & conf.high > truth, 
    error = truth - estimate
  )

table <- simulation_results |>
  group_by(N, parameter, outcome_type, estimator) |>
  summarize(
    coverage = mean(covered), 
    n = n(), 
    mae = mean(abs(error)), 
    me = mean(error)
  ) |>
  arrange(N)

table_formatted <- table |> 
  select(estimator, N, outcome_type, parameter, mae, me, coverage) |> 
  mutate_at(vars(c(mae, me)), scales::number_format(accuracy = 0.01, scale = 100)) |> 
  mutate_at(vars(coverage), scales::percent_format(accuracy = 0.1)) |> 
  pivot_wider(names_from = c("estimator", "parameter"), values_from = c("me", "mae", "coverage")) |>
  ungroup() |>
  select(outcome_type, N, ends_with("psi1"), ends_with("psi2"), ends_with("ER"), ends_with("SMR"), ends_with("direct")) 

table_formatted_1 <- table_formatted |> 
  ungroup() |>
  select(N, ends_with("psi1")) |>
  knitr::kable(format = "latex")

table_formatted_direct <- table_formatted |>
  ungroup() |>
  select(outcome_type, N, ends_with("direct")) |>
  select(outcome_type, N, contains("me"), contains("mae"), contains("coverage")) |>
  select(-coverage_GLM_direct, -coverage_WeightIt_direct)
  #knitr::kable(format = "latex")

table_formatted_psi2 <- table_formatted |>
  ungroup() |>
  select(outcome_type, N, ends_with("psi2")) |>
  select(outcome_type, N, contains("me"), contains("mae"), contains("coverage")) |>
  select(-coverage_GLM_psi2)

table_formatted_smr <- table_formatted |>
  ungroup() |>
  select(outcome_type, N, ends_with("SMR")) |>
  select(outcome_type, N, contains("me"), contains("mae"), contains("coverage")) |>
  select(-coverage_GLM_SMR)

table_formatted_er <- table_formatted |>
  ungroup() |>
  select(outcome_type, N, ends_with("ER")) |>
  select(outcome_type, N, contains("me"), contains("mae"), contains("coverage")) |>
  select(-coverage_GLM_ER)
