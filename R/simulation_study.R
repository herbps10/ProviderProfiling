simulate_data <- function(seed = 5462, N = 1e3, m = 5) {
  set.seed(5)
  hospital_effects <- rbinom(m, 1, 0.1)
  set.seed(seed)
  W <- matrix(runif(N * 5, 0, 1), ncol = 5, nrow = N)
  A <- map2_int(W[, 1], W[, 2], \(w1, w2) {
    sample(1:m, 1, prob = rep(1, m) + ifelse(rep(w1, m) > 0.5, 10 * hospital_effects, 0))
  })
  Ya <- pmap(
    list(W[, 1], W[, 2], W[, 3], W[, 4], W[, 5]), 
    \(w1, w2, w3, w4, w5) {
      tibble(
        a = 1:m, 
        #mu = -2 + -2 * w1 + sqrt(w2) - w3 + 5 * hospital_effects[a], 
        mu = -2 - 10 * w1 + w5 + 10 * hospital_effects[a],
        Ya = rbinom(m, size = 1, prob = plogis(mu))
      ) 
    }
  )
  
  tibble(
    index = 1:N,
    W1 = W[, 1],
    W2 = W[, 2],
    W3 = W[, 3],
    W4 = W[, 4],
    W5 = W[, 5],
    A = A,
    Ya = Ya,
    hospital_effect = hospital_effects[A]
  ) %>%
    mutate(Y = map2_dbl(Ya, A, \(ya, a) ya$Ya[a]),
           YA = pmap_dbl(list(Ya, W1, W2), \(ya, w1, w2) weighted.mean(ya$Ya, rep(1, m) + ifelse(rep(w1, m) > 0.5, 10 * hospital_effects, 0))))
}

largedat <- simulate_data(1e4, m = 5, seed = 5)
truth <- largedat %>% group_by(A, hospital_effect) %>% 
  summarize(psi1 = mean(Y), 
            psi2 = mean(YA)) %>%
  ungroup() %>%
  mutate(a = as.character(A)) %>%
  mutate(smr = psi1 / psi2)

truth

hist(truth$psi2)

dat <- simulate_data(5e3, m = 5, seed = 5) %>%
  select(-Ya, -YA)

dat %>% group_by(A) %>% summarize(Y = mean(Y))

fit1 <- smr_tmle(dat, "A", "Y", c("W1", "W2", "W3", "W4", "W5"), learners_trt = c("mean"), learners_outcome = c("mean"), folds = 1)
fit2 <- smr_tmle(dat, "A", "Y", c("W1", "W2", "W3", "W4", "W5"), learners_trt = c("glm", "nnet"), learners_outcome = c("glm", "nnet"), folds = 1)

plot(fit1$g$treatment_probs[, 1], fit2$g$treatment_probs[, 1])

plot(truth$psi2, fit1$estimates[, 2])
plot(truth$psi2, fit2$estimates[, 2], col = "red")

plot(truth$smr, fit1$estimates[, 3])
plot(fit1$estimates[, 2], fit2$estimates[, 2])

plot(fit1$estimates[, 3], fit2$estimates[, 3], xlim = c(0, 2), ylim = c(0, 2))
truth$smr

table(dat$Y)
table(dat$A)

incorrect <- c("mean")
correct <- c("mean", "glm", "nnet")

N_simulations <- 200
simulations <- expand_grid(
  index = 1:N_simulations, 
  #N = c(5e2, 1e3, 2500, 5e3),
  N = c(1e3, 2500, 5e3),
  learners_trt = list(incorrect, correct),
  learners_outcome = list(incorrect, correct)
) %>%
  mutate(
    seed = 1:n(),
  ) %>%
  .[1236,] %>%
  mutate(
    data = map2(seed, N, simulate_data, m = 5),
    fit_tmle = pmap(list(data, learners_trt, learners_outcome), \(data, learners_trt, learners_outcome) {
      smr_tmle(
        data,
        outcome = "Y",
        trt = "A",
        baseline = c("W1", "W2", "W3", "W4", "W5"),
        learners_trt = learners_trt,
        learners_outcome = learners_outcome,
        control = smr_control(.learners_trt_folds = 3, .learners_outcome_folds = 3)
      )
    })
  )

simulation_results <- simulations %>%
  mutate(tmle = map(fit_tmle, tidy), learners_trt = map_chr(learners_trt, paste, collapse = ", "), learners_outcome = map_chr(learners_outcome, paste, collapse = ", ")) %>%
  select(index, N, tmle, starts_with("learners")) %>%
  unnest(cols = c("tmle")) %>%
  left_join(select(truth, a, smr), by = c(trt = "a")) %>%
  mutate(covered = conf.low < smr & conf.high > smr, error = smr - estimate)

simulation_results %>%
  group_by(N, learners_trt, learners_outcome) %>%
  summarize(coverage = mean(covered), mae = mean(abs(error)), me = mean(error)) %>%
  arrange(learners_trt, learners_outcome, N)
