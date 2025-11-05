# =====================================================
# Simulated Annealing para Ajuste de Modelo Logístico
# =====================================================

set.seed(42)

# ----- 1. Simulação de dados -----
n <- 100
x <- seq(0, 10, length.out = n)
A_true <- 10
k_true <- 1.5
x0_true <- 5
y_true <- A_true / (1 + exp(-k_true * (x - x0_true)))
y_obs <- y_true + rnorm(n, sd = 0.4)

# ----- 2. Função logística -----
mu <- function(x, params) {
  A <- params[1]
  k <- params[2]
  x0 <- params[3]
  A / (1 + exp(-k * (x - x0)))
}

# ----- 3. Função de perda (Soma dos Quadrados) -----
ssq <- function(params) {
  sum((y_obs - mu(x, params))^2)
}

# ----- 4. Implementação do Simulated Annealing -----
simulated_annealing <- function(fn, init, T0 = 1, alpha = 0.95, iter = 5000) {
  current <- init
  best <- init
  current_val <- fn(current)
  best_val <- current_val
  T <- T0

  hist_params <- matrix(NA, nrow = iter, ncol = 3)
  hist_loss <- numeric(iter)
  hist_T <- numeric(iter)

  for (i in 1:iter) {
    # Vizinhança (perturbação aleatória)
    proposal <- current + rnorm(length(init), 0, 0.1)
    proposal_val <- fn(proposal)

    # Critério de aceitação
    delta <- proposal_val - current_val
    if (delta < 0 || runif(1) < exp(-delta / T)) {
      current <- proposal
      current_val <- proposal_val
    }

    # Atualiza o melhor encontrado
    if (current_val < best_val) {
      best <- current
      best_val <- current_val
    }

    # Armazenar histórico
    hist_params[i, ] <- current
    hist_loss[i] <- current_val
    hist_T[i] <- T

    # Resfriamento
    T <- alpha * T
  }

  return(list(
    best_params = best,
    best_val = best_val,
    hist_params = hist_params,
    hist_loss = hist_loss,
    hist_T = hist_T
  ))
}

# ----- 5. Rodar o SA -----
init_guess <- c(5, 0.5, 2) # chute inicial (A, k, x0)
result_SA <- simulated_annealing(
  ssq,
  init_guess,
  T0 = 1,
  alpha = 0.99,
  iter = 8000
)
result_SA$best_params

# ----- 6. Comparar com nls() -----
model_nls <- nls(
  y_obs ~ A / (1 + exp(-k * (x - x0))),
  start = list(A = 5, k = 0.5, x0 = 2)
)
coef(model_nls)

# ----- 7. Gráficos -----

par(mfrow = c(2, 2))

# (1) Dados + Ajustes
plot(
  x,
  y_obs,
  pch = 16,
  col = "gray",
  main = "Ajuste Logístico via SA e nls",
  xlab = "x",
  ylab = "y"
)
lines(x, mu(x, result_SA$best_params), col = "blue", lwd = 2)
lines(x, predict(model_nls), col = "red", lwd = 2, lty = 2)
legend(
  "bottomright",
  legend = c("Dados", "SA", "nls"),
  col = c("gray", "blue", "red"),
  pch = c(16, NA, NA),
  lty = c(NA, 1, 2),
  lwd = c(NA, 2, 2)
)

# (2) Evolução da Função de Perda
plot(
  result_SA$hist_loss,
  type = "l",
  lwd = 1.5,
  col = "darkgreen",
  main = "Evolução da Função de Perda (SSQ)",
  xlab = "Iterações",
  ylab = "Soma dos Quadrados"
)
abline(h = result_SA$best_val, col = "red", lty = 2)

# (3) Evolução dos Parâmetros
matplot(
  result_SA$hist_params,
  type = "l",
  lwd = 1.2,
  main = "Trajetória dos Parâmetros (A, k, x0)",
  xlab = "Iterações",
  ylab = "Valor",
  col = c("blue", "purple", "orange")
)
legend(
  "topright",
  legend = c("A", "k", "x0"),
  col = c("blue", "purple", "orange"),
  lty = 1
)

# (4) Resfriamento da Temperatura
plot(
  result_SA$hist_T,
  type = "l",
  lwd = 1.5,
  col = "red",
  main = "Esquema de Resfriamento",
  xlab = "Iterações",
  ylab = "Temperatura (T)"
)

par(mfrow = c(1, 1))

# --- Superfície de perda para A e k fixando x0 ---
A_seq <- seq(8, 12, length.out = 50)
k_seq <- seq(1, 2, length.out = 50)
loss_surface <- outer(A_seq, k_seq, Vectorize(function(A, k) ssq(c(A, k, 5))))

persp(
  A_seq,
  k_seq,
  loss_surface,
  xlab = "A",
  ylab = "k",
  zlab = "Soma dos Quadrados",
  main = "Superfície da Função de Perda (x0 fixo = 5)",
  theta = 45,
  phi = 30,
  col = "lightblue"
)
