# sa_gumbel.R
# Estimação dos parâmetros (mu, beta) da distribuição de Gumbel
# via Simulated Annealing (maximização da log-verossimilhança)

set.seed(2025)
library(ggplot2)

# -----------------------------
# Função log-verossimilhança
# -----------------------------
loglik_gumbel <- function(params, x) {
  mu <- params[1]
  beta <- params[2]
  if (beta <= 0) {
    return(-Inf)
  }
  n <- length(x)
  ll <- -n * log(beta) - sum((x - mu) / beta) - sum(exp(-(x - mu) / beta))
  return(ll)
}

# -----------------------------
# Gerador de vizinho
# -----------------------------
propose_neighbor <- function(curr, sds = c(0.2, 0.1)) {
  mu1 <- curr[1] + rnorm(1, 0, sds[1])
  beta1 <- curr[2] + rnorm(1, 0, sds[2])
  if (beta1 <= 0) {
    beta1 <- abs(beta1) + 1e-6
  }
  c(mu1, beta1)
}

# -----------------------------
# Algoritmo SA
# -----------------------------
sa_gumbel <- function(
  x,
  init = NULL,
  T0 = 10,
  Tf = 0.001,
  cooling_rate = 0.98,
  L = 15,
  sds = c(0.2, 0.1)
) {
  if (is.null(init)) {
    init <- c(mean(x), sd(x)) # chute inicial
  }

  curr <- init
  curr_ll <- loglik_gumbel(curr, x)
  best <- curr
  best_ll <- curr_ll
  T <- T0
  eval <- 1
  history <- data.frame(eval = eval, ll = curr_ll)

  while (T > Tf) {
    for (i in 1:L) {
      prop <- propose_neighbor(curr, sds)
      prop_ll <- loglik_gumbel(prop, x)
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

# -----------------------------
# Geração de dados Gumbel
# -----------------------------
rgumbel <- function(n, mu, beta) {
  u <- runif(n)
  mu - beta * log(-log(u))
}

# -----------------------------
# Execução do experimento
# -----------------------------
true_params <- c(mu = 2, beta = 1.5)
n <- 1000
x <- rgumbel(n, true_params[1], true_params[2])

# Rodar SA
t0 <- Sys.time()
res <- sa_gumbel(
  x,
  T0 = 10,
  Tf = 1e-4,
  cooling_rate = 0.99,
  L = 30,
  sds = c(0.1, 0.05)
)
t1 <- Sys.time()
runtime <- as.numeric(difftime(t1, t0, units = "secs"))

# -----------------------------
# Resultados
# -----------------------------
cat("\nParâmetros verdadeiros:\n")
print(true_params)

cat("\nEstimativas via SA:\n")
print(res$best)

cat(sprintf(
  "\nLog-verossimilhança (verdadeiros): %.4f",
  loglik_gumbel(true_params, x)
))
cat(sprintf("\nLog-verossimilhança (estimados): %.4f", res$best_ll))
cat(sprintf("\nTempo: %.2f s\n", runtime))

# -----------------------------
# Gráfico da evolução da log-verossimilhança
# -----------------------------
p <- ggplot(res$history, aes(x = eval, y = ll)) +
  geom_line(color = "blue") +
  geom_hline(
    yintercept = loglik_gumbel(true_params, x),
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    title = "Evolução da log-verossimilhança (SA - Gumbel)",
    x = "Iterações",
    y = "Log-verossimilhança"
  ) +
  theme_minimal(base_size = 13)
print(p)
