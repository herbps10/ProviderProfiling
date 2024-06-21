direct_glm <- function(data, outcome, trt, baseline) {
  trt_levels <- sort(unique(data[[trt]]))
  data[[trt]] <- factor(data[[trt]], levels = trt_levels)
  f <- as.formula(paste0(outcome, "~ ", trt, " + ", paste0(baseline, collapse = " + "))) 
  fit <- glm(f, data = data, family = binomial)
  tibble(
    a = trt_levels
  ) %>%
    mutate(direct = map_dbl(a, \(trt_level) {
      newdata <- data
      newdata[[trt]] <- factor(trt_level, levels = trt_levels)
      mean(predict(fit, newdata, type = "response"), na.rm = TRUE)
    }))
}