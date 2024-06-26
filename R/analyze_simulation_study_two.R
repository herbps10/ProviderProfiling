dyn.load("/gpfs/share/apps/gcc/11.2.0/lib64/libstdc++.so.6")
library(torch)
library(tidyverse)
library(TargetedSMR)
library(broom)

source("R/simulate.R")

simulations <- read_rds("results/simulation_results_two.rds")

truth_indirect <- simulations %>%
  filter(parameter == "smr") %>%
  select(index, N, scenario, truth) %>% 
  unnest(c(truth)) %>% 
  select(index, N, A, scenario, psi1, psi2, smr, er) %>%
  mutate(A = as.character(A)) %>%
  distinct()

truth_direct <- simulations %>%
  select(index, N, scenario, truth_direct) %>% 
  unnest(c(truth_direct)) %>% 
  rename(A = a) %>%
  select(index, N, A, scenario, direct) %>%
  mutate(A = as.character(A)) %>%
  distinct()

glm_estimates <- simulations %>%
  select(index, N, estimator, scenario, trt_method, fit_glm) %>% 
  unnest(fit_glm) %>% 
  mutate(a = as.character(a)) %>% 
  pivot_longer(c(glm_direct, glm_psi1, glm_psi2, glm_SMR, glm_ER)) %>% 
  #pivot_longer(c(glm_psi1, glm_psi2, glm_SMR, glm_ER)) %>% 
  mutate(name = str_replace(name, "glm_", "")) %>%
  filter(!is.na(value)) %>%
  rename(glm = value)

simulation_results <- simulations %>%
  ungroup() %>%
  mutate(tmle = map(fit, tidy)) %>%
  select(index, N, tmle, scenario, trt_method) %>%
  unnest(cols = c("tmle")) %>%
  left_join(truth_indirect, by = c("index", "N", "scenario", trt = "A")) %>%
  left_join(truth_direct, by = c("index", "N", "scenario", trt = "A")) %>%
  left_join(glm_estimates, by = c("index", "N", "estimator", "parameter" = "name", "scenario", "trt_method", trt = "a")) %>%
  mutate(truth = case_when(
    parameter == "psi1" ~ psi1,
    parameter == "psi2" ~ psi2,
    parameter == "SMR" ~ smr,
    parameter == "ER" ~ er,
    parameter == "direct" ~ direct
  )) %>%
  mutate(
    covered = round(conf.low, 2) <= round(truth, 2) & round(conf.high, 2) >= round(truth, 2), error = truth - estimate,
    covered = conf.low < truth & conf.high > truth, error = truth - estimate,
    glm_error = truth - glm
  )

table <- simulation_results %>%
  group_by(N, parameter, scenario, trt_method) %>%
  summarize(
    coverage = mean(covered), 
    n = n(), 
    mae = mean(abs(error)), 
    me = mean(error), 
    #ci_width = mean(conf.high - conf.low)
    glm_me = mean(glm_error),
    glm_mae = mean(abs(glm_error))
  ) %>%
  arrange(scenario, N)

table_formatted <- table %>% 
  filter(scenario == 1) %>%
  select(N, trt_method, parameter, mae, me, coverage, glm_me, glm_mae) %>% 
  mutate_at(vars(c(mae, me, glm_me, glm_mae)), scales::number_format(accuracy = 0.001)) %>% 
  mutate_at(vars(coverage), scales::percent_format(accuracy = 0.1)) %>% 
  #mutate_at(vars(ci_width), scales::number_format(accuracy = 0.01)) %>% 
  #pivot_wider(names_from = c("parameter"), values_from = c("me", "mae", "coverage", "ci_width")) %>%
  pivot_wider(names_from = c("parameter"), values_from = c("me", "mae", "coverage", "glm_me", "glm_mae")) %>%
  select(scenario, N, ends_with("psi1"), ends_with("psi2"), ends_with("ER"), ends_with("SMR"), ends_with("direct")) 

table_formatted_1 <- table_formatted %>% 
  ungroup() %>%
  select(scenario, N, ends_with("psi1")) %>%
  knitr::kable(format = "latex")

table_formatted_direct <- table_formatted %>%
  ungroup() %>%
  select(N, me_direct, glm_me_direct, mae_direct, glm_mae_direct, coverage_direct) %>%
  #pivot_longer(c(me_direct, glm_me_direct, mae_direct, glm_mae_direct, coverage_direct)) %>%
  #mutate(estimator = ifelse(str_detect(name, "glm"), "GLM", "TMLE"),
  #       metric = str_replace_all(name, "_direct", ""),
  #       metric = str_replace_all(metric, "glm_", "")) %>%
  #pivot_wider(names_from = "metric", values_from = "value") %>%
  #select(scenario, N, estimator, me, mae, coverage) %>%
  knitr::kable(format = "latex")

table_formatted_psi2 <- table_formatted %>%
  ungroup() %>%
  select(N, me_psi2, glm_me_psi2, mae_psi2, glm_mae_psi2, coverage_psi2) %>%
  knitr::kable(format = "latex")

table_formatted_smr <- table_formatted %>%
  ungroup() %>%
  select(N, me_SMR, glm_me_SMR, mae_SMR, glm_mae_SMR, coverage_SMR) %>%
  knitr::kable(format = "latex")

table_formatted_er <- table_formatted %>%
  ungroup() %>%
  select(N, me_ER, glm_me_ER, mae_ER, glm_mae_ER, coverage_ER) %>%
  knitr::kable(format = "latex")
