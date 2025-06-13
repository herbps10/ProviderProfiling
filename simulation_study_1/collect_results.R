#
# Combine cached results into single results file
#
library(tidyverse)

cache_path <- Sys.getenv("SIMULATION_CACHE_PATH")
results_path <- Sys.getenv("SIMULATION_RESULTS_PATH")

if(cache_path == "") stop("Please set SIMULATION_CACHE_PATH environment variable.")
if(results_path == "") stop("Please set SIMULATION_RESULTS_PATH environment variable.")

files <- Sys.glob(glue::glue("{cache_path}/*/*.rds"))

results <- map_df(files, read_rds) %>% dplyr::bind_rows()

write_rds(results, file = glue::glue("{results_path}/simulation_results.rds")
