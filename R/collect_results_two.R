library(tidyverse)
files <- Sys.glob("/gpfs/scratch/susmah01/ProviderProfiling/cache/*/sim2_*.rds")

results <- map_df(files, read_rds) %>% dplyr::bind_rows()

write_rds(results, file = "/gpfs/home/susmah01/ProviderProfiling/results/simulation_results_two.rds")
