direct_glm <- function(data, outcome, trt, baseline, outcome_type = "continuous") {
  trt_levels <- sort(unique(data[[trt]]))
  data[[trt]] <- factor(data[[trt]], levels = trt_levels)
  f <- as.formula(paste0(outcome, "~ ", trt, " + ", paste0(baseline, collapse = " + "), " + ", paste0(trt, "*", baseline, collapse = " +"))) 
  if(outcome_type == "binomial") {
    fit <- glm(f, data = data, family = binomial(link = "logit"))
  }
  else {
    fit <- glm(f, data = data)
  }
  
  tibble(
    a = trt_levels
  ) %>%
    mutate(glm_direct = map_dbl(a, \(trt_level) {
      newdata <- data
      newdata[[trt]] <- factor(trt_level, levels = trt_levels)
      mean(predict(fit, newdata, type = "response"), na.rm = TRUE)
    }))
}

indirect_glm <- function(data, outcome, trt, baseline, outcome_type = "continuous") {
  trt_levels <- sort(unique(data[[trt]]))
  f <- as.formula(paste0(outcome, "~ ", paste0(baseline, collapse = " + ")))

  if(outcome_type == "binomial") {
    fit <- glm(f, data = data, family = binomial(link = "logit"))
  }
  else {
    fit <- glm(f, data = data)
  }

  tibble(
    a = trt_levels
  ) %>%
    mutate(
      glm_psi1 = map_dbl(a, \(trt_level) mean(data[data[[trt]] == trt_level, outcome])),
      glm_psi2 = map_dbl(a, \(trt_level) {
        mean(predict(fit, newdata = data[data[[trt]] == trt_level, ], type = "response"), na.rm = TRUE)
      }),
      glm_ER = glm_psi1 - glm_psi2,
      glm_SMR = glm_psi1 / glm_psi2
    )
}
