# sa_mistura_normais.R
# Estimação de parâmetros da mistura de duas Normais via Simulated Annealing
# p * N(mu1, sigma1^2) + (1-p) * N(mu2, sigma2^2)

set.seed(2025)
library(ggplot2)
library(dplyr)
library(gridExtra)

# ---------------------------
# Função de log-verossimilhança
# ---------------------------
loglik_mistura <- function(params, x) {
  p <- params[1]
  mu1 <- params[2]
  s1 <- params[3]
  mu2 <- params[4]
  s2 <- params[5]
  if (p <= 0 || p >= 1 || s1 <= 0 || s2 <= 0) {
    return(-Inf)
  }
  dens <- p * dnorm(x, mu1, s1) + (1 - p) * dnorm(x, mu2, s2)
  if (any(dens <= 0)) {
    return(-Inf)
  }
  sum(log(dens))
}

# ---------------------------
# Função para gerar vizinho
# ---------------------------
propose_neighbor <- function(curr, sds = c(0.05, 0.2, 0.1, 0.2, 0.1)) {
  p <- curr[1] + rnorm(1, 0, sds[1])
  mu1 <- curr[2] + rnorm(1, 0, sds[2])
  s1 <- curr[3] + rnorm(1, 0, sds[3])
  mu2 <- curr[4] + rnorm(1, 0, sds[4])
  s2 <- curr[5] + rnorm(1, 0, sds[5])
  p <- min(max(p, 1e-3), 1 - 1e-3) # garantir 0 < p < 1
  s1 <- abs(s1) + 1e-6
  s2 <- abs(s2) + 1e-6
  c(p, mu1, s1, mu2, s2)
}

# ---------------------------
# Algoritmo Simulated Annealing
# ---------------------------
sa_mistura <- function(
  x,
  init = NULL,
  T0 = 10,
  Tf = 0.001,
  cooling_rate = 0.98,
  L = 10,
  sds = c(0.05, 0.2, 0.1, 0.2, 0.1)
) {
  if (is.null(init)) {
    init <- c(0.5, mean(x) - sd(x) / 2, sd(x) / 2, mean(x) + sd(x) / 2, sd(x))
  }
  curr <- init
  curr_ll <- loglik_mistura(curr, x)
  best <- curr
  best_ll <- curr_ll
  T <- T0
  eval <- 1
  history <- data.frame(eval = eval, ll = curr_ll)

  while (T > Tf) {
    for (i in 1:L) {
      prop <- propose_neighbor(curr, sds)
      prop_ll <- loglik_mistura(prop, x)
      if (prop_ll > curr_ll) {
        curr <- prop
        curr_ll <- prop_ll
      } else {
        if (runif(1) < exp((prop_ll - curr_ll) / T)) {
          curr <- prop
          curr_ll <- prop_ll
        }
      }
      if (curr_ll > best_ll) {
        best <- curr
        best_ll <- curr_ll
      }
      eval <- eval + 1
      history <- rbind(history, data.frame(eval = eval, ll = curr_ll))
    }
    T <- cooling_rate * T
  }
  list(best = best, best_ll = best_ll, history = history)
}

# ---------------------------
# Simulação de dados
# ---------------------------
gerar_mistura <- function(n, p, mu1, s1, mu2, s2) {
  z <- rbinom(n, 1, p)
  rnorm(n, mean = ifelse(z == 1, mu1, mu2), sd = ifelse(z == 1, s1, s2))
}

# ---------------------------
# Execução principal
# ---------------------------
true_params <- c(p = 0.4, mu1 = 0, s1 = 1, mu2 = 4, s2 = 1.2)
n <- 1000
x <- gerar_mistura(
  n,
  true_params[1],
  true_params[2],
  true_params[3],
  true_params[4],
  true_params[5]
)

# Rodar SA
t0 <- Sys.time()
res <- sa_mistura(x, T0 = 15, Tf = 0.0001, cooling_rate = 0.98, L = 15)
t1 <- Sys.time()
runtime <- as.numeric(difftime(t1, t0, units = "secs"))

# ---------------------------
# Resultados
# ---------------------------
cat("\nParâmetros verdadeiros:\n")
print(true_params)
cat("\nEstimativas via SA:\n")
print(res$best)
cat(sprintf(
  "\nLog-verossimilhança (verdadeiros): %.4f",
  loglik_mistura(true_params, x)
))
cat(sprintf("\nLog-verossimilhança (estimados): %.4f", res$best_ll))
cat(sprintf("\nTempo: %.2f s\n", runtime))

# ---------------------------
# Gráfico da evolução da log-verossimilhança
# ---------------------------
p <- ggplot(res$history, aes(x = eval, y = ll)) +
  geom_line(color = "blue") +
  geom_hline(
    yintercept = loglik_mistura(true_params, x),
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    title = "Evolução da log-verossimilhança (SA - Mistura de Normais)",
    x = "Iteração",
    y = "Log-verossimilhança"
  ) +
  theme_minimal(base_size = 13)
print(p)
