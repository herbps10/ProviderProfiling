#' Simulate data for simulation study 1.
#'
#' @param seed random number seed
#' @param N number of observations
#' @param m number of providers
#' @param outcome_type set to "continuous" to generate a continuous outcome variable, or "binomial" to generate a binary outcome
#' @param P number of baseline covariates
#'
#' @return data frame of simulated data. 
#' 
simulate_data <- function(seed = 1, N = 1e3, m = 5, P = 5, outcome_type = "continuous") {
  set.seed(seed)
  hospital_effect1 <- rbinom(m, 1, 0.5)
  W <- matrix(runif(N * P, 0, 1), ncol = P, nrow = N)
  colnames(W) <- paste0("W", 1:P)
  
  g <- plogis(
    2 * matrix(ifelse(W[, 1] > 0.7, 1, -1), ncol = m, nrow = N, byrow = FALSE) * matrix(hospital_effect1 * 2 - 1, ncol = m, nrow = N, byrow = TRUE)
  )
  g <- g / rowSums(g)
  colnames(g) <- paste0("g", 1:m)
  
  A <- map_int(1:N, \(i) sample(1:m, 1, prob = g[i, ]))
  
  trt_indicator <- matrix(0, ncol = m, nrow = N)
  for(a in 1:m) {
    trt_indicator[, a] <- A == a
  }
  
  Qbar <- matrix(0, ncol = m, nrow = N)
  for(i in 1:m) {
    Qbar[, i] <- case_when(
      W[, 1] <= 0.4 & hospital_effect1[i] == 1 ~ 0.3,
      W[, 1]  > 0.7 & hospital_effect1[i] == 1 ~ 1,
      W[, 1]  > 0.4 & hospital_effect1[i] == 1 ~ 1.5,
      W[, 1] <= 0.4 & hospital_effect1[i] == 0 ~ 0.7,
      W[, 1]  > 0.7 & hospital_effect1[i] == 0 ~ 0.5,
      W[, 1]  > 0.4 & hospital_effect1[i] == 0 ~ 0
    )
  }
  colnames(Qbar) <- paste0("Qbar", 1:m)

  if(outcome_type == "binomial") {
    Qbar <- plogis(Qbar)
  }
  
  Qtilde <- rowSums(Qbar * g)

  if(outcome_type == "continuous") {
    Ya <- matrix(rnorm(N * m, as.vector(Qbar), 0.1), ncol = m, nrow = N, byrow = FALSE)
  }
  else {
    Ya <- matrix(rbinom(N * m, 1, as.vector(Qbar)), ncol = m, nrow = N, byrow = FALSE)
  }
  colnames(Ya) <- paste0("Y", 1:m)
  
  Y  <- rowSums(Ya * trt_indicator)
  YA <- rowSums(Ya * g)
  
  gamma1 <- rowSums(Qbar * trt_indicator)
  gamma2 <- rowSums(Qbar * g)
  
  cbind(
    as.data.frame(W),
    as.data.frame(g),
    as.data.frame(Qbar),
    as.data.frame(Ya)
  ) %>%
    mutate(Y = Y, A = A, YA = YA, Qtilde = Qtilde, hospital_effect1 = hospital_effect1[A], gamma1 = gamma1, gamma2 = gamma2)
}
