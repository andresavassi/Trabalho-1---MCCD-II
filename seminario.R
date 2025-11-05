# Simulated Annealing para estimação de parâmetros (R)
# Rodar tudo de uma vez. Usa apenas pacotes base + graphics.

set.seed(20252)

# ---------------------------
# Funções: negative log-likelihoods
# ---------------------------

nll_weibull <- function(params, data) {
  # params: vector in original scale: c(shape, scale)
  k <- params[1]
  lambda <- params[2]
  if (k <= 0 || lambda <= 0) {
    return(1e10)
  }
  -sum(dweibull(data, shape = k, scale = lambda, log = TRUE))
}

nll_gamma <- function(params, data) {
  # params: c(shape, scale) (R's rgamma uses shape, scale)
  a <- params[1]
  b <- params[2]
  if (a <= 0 || b <= 0) {
    return(1e10)
  }
  -sum(dgamma(data, shape = a, scale = b, log = TRUE))
}

nll_lognorm <- function(params, data) {
  # params: c(meanlog, sdlog) (sdlog > 0)
  mu <- params[1]
  s <- params[2]
  if (s <= 0) {
    return(1e10)
  }
  -sum(dlnorm(data, meanlog = mu, sdlog = s, log = TRUE))
}

# ---------------------------
# Função genérica de SA
# ---------------------------
simulated_annealing <- function(
  nll_fn, # function(params, data) -> numeric (to minimize)
  init_params, # initial values in original scale
  data,
  transform = identity, # function to map original -> working (e.g., log)
  inv_transform = identity, # inverse transform working -> original
  proposal_sd = rep(0.1, length(init_params)), # sd of Gaussian proposals in working space
  n_iter = 5000,
  T0 = 1.0,
  alpha = 0.995,
  verbose = TRUE
) {
  # Working parameterization
  curr_work <- transform(init_params)
  curr_params <- inv_transform(curr_work)
  curr_cost <- nll_fn(curr_params, data)
  best_work <- curr_work
  best_params <- curr_params
  best_cost <- curr_cost

  cost_trace <- numeric(n_iter)
  T <- T0

  for (k in 1:n_iter) {
    # propose in working space
    prop_work <- curr_work +
      rnorm(length(curr_work), mean = 0, sd = proposal_sd)
    prop_params <- inv_transform(prop_work)
    prop_cost <- nll_fn(prop_params, data)

    delta <- prop_cost - curr_cost
    if (delta <= 0 || runif(1) < exp(-delta / T)) {
      # accept
      curr_work <- prop_work
      curr_params <- prop_params
      curr_cost <- prop_cost
      if (curr_cost < best_cost) {
        best_cost <- curr_cost
        best_work <- curr_work
        best_params <- curr_params
      }
    }
    cost_trace[k] <- curr_cost
    # cool down
    T <- T * alpha
    # optional: adjust proposal_sd adaptively (not implemented)
  }

  list(
    best_params = best_params,
    best_cost = best_cost,
    final_params = curr_params,
    final_cost = curr_cost,
    cost_trace = cost_trace
  )
}

# ---------------------------
# Simulação de dados (três distribuições)
# ---------------------------

n <- 500

# 1) Weibull: shape k, scale lambda
true_weib_k <- 1.8
true_weib_lambda <- 2.5
data_weib <- rweibull(n, shape = true_weib_k, scale = true_weib_lambda)

# 2) Gamma: shape a, scale b
true_gamma_a <- 2.2
true_gamma_b <- 1.3
data_gamma <- rgamma(n, shape = true_gamma_a, scale = true_gamma_b)

# 3) Log-normal: meanlog mu, sdlog s
true_logn_mu <- 0.5
true_logn_s <- 0.7
data_logn <- rlnorm(n, meanlog = true_logn_mu, sdlog = true_logn_s)

# ---------------------------
# Estimativas via SA para cada distribuição
# ---------------------------

# ---- Weibull ----
# optimize shape>0, scale>0 -> use log-transform
transform_weib <- function(p) log(p) # original -> working
inv_transform_weib <- function(w) exp(w) # working -> original

init_weib <- c(shape = 1.0, scale = 1.0) # starting guess
res_weib <- simulated_annealing(
  nll_fn = nll_weibull,
  init_params = init_weib,
  data = data_weib,
  transform = transform_weib,
  inv_transform = inv_transform_weib,
  proposal_sd = c(0.08, 0.08),
  n_iter = 8000,
  T0 = 1.0,
  alpha = 0.998,
  verbose = FALSE
)

# ---- Gamma ----
transform_gamma <- function(p) log(p)
inv_transform_gamma <- function(w) exp(w)
init_gamma <- c(shape = 1.0, scale = 1.0)
res_gamma <- simulated_annealing(
  nll_fn = nll_gamma,
  init_params = init_gamma,
  data = data_gamma,
  transform = transform_gamma,
  inv_transform = inv_transform_gamma,
  proposal_sd = c(0.08, 0.08),
  n_iter = 8000,
  T0 = 1.0,
  alpha = 0.998,
  verbose = FALSE
)

# ---- Log-normal ----
# mu in R, s > 0. We'll transform only sd via log, keep mu unconstrained (identity)
transform_logn <- function(p) c(p[1], log(p[2]))
inv_transform_logn <- function(w) c(w[1], exp(w[2]))
init_logn <- c(mu = 0.0, sd = 1.0)
res_logn <- simulated_annealing(
  nll_fn = nll_lognorm,
  init_params = init_logn,
  data = data_logn,
  transform = transform_logn,
  inv_transform = inv_transform_logn,
  proposal_sd = c(0.02, 0.06),
  n_iter = 6000,
  T0 = 1.0,
  alpha = 0.997,
  verbose = FALSE
)

# ---------------------------
# Resultados: comparar estimados vs verdadeiros
# ---------------------------

cat("=== Weibull (true)     ===\n")
cat("true shape =", true_weib_k, " true scale =", true_weib_lambda, "\n")
cat("SA estimate (best)     :", round(res_weib$best_params, 4), "\n\n")

cat("=== Gamma (true)       ===\n")
cat("true shape =", true_gamma_a, " true scale =", true_gamma_b, "\n")
cat("SA estimate (best)     :", round(res_gamma$best_params, 4), "\n\n")

cat("=== Log-Normal (true)  ===\n")
cat("true mu =", true_logn_mu, " true sd =", true_logn_s, "\n")
cat("SA estimate (best)     :", round(res_logn$best_params, 4), "\n\n")

# ---------------------------
# Plots de convergência
# ---------------------------
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
plot(
  res_weib$cost_trace,
  type = "l",
  main = "Convergência (Weibull)",
  ylab = "NLL",
  xlab = "Iter"
)
abline(h = res_weib$best_cost, col = "red")
plot(
  res_gamma$cost_trace,
  type = "l",
  main = "Convergência (Gamma)",
  ylab = "NLL",
  xlab = "Iter"
)
abline(h = res_gamma$best_cost, col = "red")
plot(
  res_logn$cost_trace,
  type = "l",
  main = "Convergência (LogNormal)",
  ylab = "NLL",
  xlab = "Iter"
)
abline(h = res_logn$best_cost, col = "red")

# ---------------------------
# Densidade comparativa (visual)
# ---------------------------
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
# Weibull density comparison
hist(
  data_weib,
  breaks = 40,
  freq = FALSE,
  main = "Weibull: dados vs fit",
  xlab = ""
)
xx <- seq(0, max(data_weib) * 1.1, length.out = 200)
lines(xx, dweibull(xx, shape = true_weib_k, scale = true_weib_lambda), lwd = 2)
lines(
  xx,
  dweibull(
    xx,
    shape = res_weib$best_params[1],
    scale = res_weib$best_params[2]
  ),
  lwd = 2,
  lty = 2
)

# Gamma
hist(
  data_gamma,
  breaks = 40,
  freq = FALSE,
  main = "Gamma: dados vs fit",
  xlab = ""
)
xx <- seq(0, max(data_gamma) * 1.1, length.out = 200)
lines(xx, dgamma(xx, shape = true_gamma_a, scale = true_gamma_b), lwd = 2)
lines(
  xx,
  dgamma(
    xx,
    shape = res_gamma$best_params[1],
    scale = res_gamma$best_params[2]
  ),
  lwd = 2,
  lty = 2
)

# Log-normal
hist(
  data_logn,
  breaks = 40,
  freq = FALSE,
  main = "LogNormal: dados vs fit",
  xlab = ""
)
xx <- seq(0, max(data_logn) * 1.1, length.out = 200)
lines(xx, dlnorm(xx, meanlog = true_logn_mu, sdlog = true_logn_s), lwd = 2)
lines(
  xx,
  dlnorm(
    xx,
    meanlog = res_logn$best_params[1],
    sdlog = res_logn$best_params[2]
  ),
  lwd = 2,
  lty = 2
)

par(mfrow = c(1, 1))
