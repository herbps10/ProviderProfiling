library(tidyverse)
files <- Sys.glob("/gpfs/scratch/susmah01/ProviderProfiling/cache/*/*.rds")
files <- files[!str_detect(files, "sim2")]

results <- map_df(files, read_rds) %>% dplyr::bind_rows()

write_rds(results, file = "/gpfs/home/susmah01/ProviderProfiling/results/simulation_results.rds")
