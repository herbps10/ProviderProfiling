library(tidyverse)
library(TargetedSMR)
library(broom)

results_path <- Sys.getenv("SIMULATION_RESULTS_PATH")
if(results_path == "") stop("Please set SIMULATION_RESULTS_PATH environment variable.")

simulations <- read_rds(glue::glue("{results_path}/simulation_results.rds"))

truth_indirect <- simulations |>
  filter(parameter == "smr") |>
  select(index, N, scenario, outcome_type, truth) |> 
  unnest(c(truth)) |> 
  select(index, outcome_type, N, A, scenario, psi1, psi2, smr, er) |>
  mutate(A = as.character(A)) |>
  distinct()

truth_direct <- simulations |>
  select(index, N, scenario, outcome_type, truth_direct) |> 
  unnest(c(truth_direct)) |> 
  rename(A = a) |>
  select(index, outcome_type, N, A, scenario, direct) |>
  mutate(A = as.character(A)) |>
  distinct()

simulation_results <- simulations |>
  ungroup() |>
  select(index, N, outcome_type, fit, scenario) |>
  unnest(cols = c("fit")) |>
  left_join(truth_indirect, by = c("index", "N", "scenario", trt = "A")) |>
  left_join(truth_direct, by = c("index", "N", "scenario", trt = "A")) |>
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
  group_by(outcome_type, N, parameter, scenario, estimator) |>
  summarize(
    coverage = mean(covered), 
    n = n(), 
    mae = mean(abs(error)), 
    me = mean(error),
  ) |>
  arrange(scenario, N)

table_formatted <- table |> 
  select(outcome_type, N, parameter, estimator, scenario, mae, me, coverage) |> 
  mutate_at(vars(c(mae, me)), scales::number_format(accuracy = 0.01, scale = 100)) |> 
  mutate_at(vars(coverage), scales::percent_format(accuracy = 0.1)) |> 
  pivot_wider(names_from = c("parameter"), values_from = c("me", "mae", "coverage")) |>
  select(estimator, scenario, N, ends_with("psi1"), ends_with("psi2"), ends_with("ER"), ends_with("SMR"), ends_with("direct")) 

table_formatted_direct <- table_formatted |>
  ungroup() |>
  select(estimator, scenario, N, me_direct, mae_direct, coverage_direct) |>
  arrange(estimator, scenario, N) |>
  knitr::kable(format = "latex")

table_formatted_psi2 <- table_formatted |>
  ungroup() |>
  select(estimator, scenario, N, me_psi2, mae_psi2, coverage_psi2) |>
  filter(!is.na(me_psi2)) |>
  arrange(estimator, scenario, N) |>
  knitr::kable(format = "latex")

table_formatted_smr <- table_formatted |>
  ungroup() |>
  select(estimator, scenario, N, me_SMR, mae_SMR, coverage_SMR) |>
  filter(!is.na(me_SMR)) |>
  arrange(estimator, scenario, N) |>
  knitr::kable(format = "latex")

table_formatted_er <- table_formatted |>
  ungroup() |>
  select(estimator, scenario, N, me_ER, mae_ER, coverage_ER) |>
  filter(!is.na(me_ER)) |>
  arrange(estimator, scenario, N) |>
  knitr::kable(format = "latex")
