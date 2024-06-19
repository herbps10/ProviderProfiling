dyn.load("/gpfs/share/apps/gcc/11.2.0/lib64/libstdc++.so.6")
library(torch)
library(tidyverse)
library(TargetedSMR)
library(broom)

source("R/simulate.R")

simulations <- read_rds("results/simulation_results.rds")

truth_indirect <- simulations %>%
  select(index, N, scenario, truth) %>% 
  unnest(c(truth)) %>% 
  select(index, N, a, scenario, psi1, psi2, smr, er)

truth_direct <- simulations %>%
  select(index, N, scenario, truth_direct) %>% 
  unnest(c(truth_direct)) %>% 
  select(index, N, a, scenario, direct)

simulation_results <- simulations %>%
  ungroup() %>%
  mutate(tmle = map(fit, tidy)) %>%
  select(index, N, tmle, scenario, trt_method) %>%
  unnest(cols = c("tmle")) %>%
  left_join(truth_indirect, by = c("index", "N", "scenario", trt = "a")) %>%
  left_join(truth_direct, by = c("index", "N", "scenario", trt = "a")) %>%
  mutate(truth = case_when(
    parameter == "psi1" ~ psi1,
    parameter == "psi2" ~ psi2,
    parameter == "SMR" ~ smr,
    parameter == "ER" ~ er,
    parameter == "direct" ~ direct
  )) %>%
  mutate(covered = conf.low < truth & conf.high > truth, error = truth - estimate)

table <- simulation_results %>%
  group_by(N, parameter, scenario, trt_method) %>%
  summarize(
    coverage = mean(covered), 
    n = n(), 
    mae = mean(abs(error)), 
    me = mean(error), 
    #ci_width = mean(conf.high - conf.low)
  ) %>%
  arrange(scenario, N)

table_formatted <- table %>% 
  select(N, scenario, trt_method, parameter, mae, me, coverage) %>% 
  mutate_at(vars(c(mae, me)), scales::number_format(accuracy = 0.001)) %>% 
  mutate_at(vars(coverage), scales::percent_format(accuracy = 0.1)) %>% 
  #mutate_at(vars(ci_width), scales::number_format(accuracy = 0.01)) %>% 
  #pivot_wider(names_from = c("parameter"), values_from = c("me", "mae", "coverage", "ci_width")) %>%
  pivot_wider(names_from = c("parameter"), values_from = c("me", "mae", "coverage")) %>%
  select(scenario, N, ends_with("psi1"), ends_with("psi2"), ends_with("ER"), ends_with("SMR"), ends_with("direct")) 

table_formatted_1 <- table_formatted %>% 
  ungroup() %>%
  select(scenario, N, ends_with("psi1"), ends_with("psi2")) %>%
  knitr::kable(format = "latex")

table_formatted_2 <- table_formatted %>%
  ungroup() %>%
  select(scenario, N, ends_with("direct"), ends_with("ER"), ends_with("SMR")) %>%
  knitr::kable(format = "latex")
