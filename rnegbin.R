rnegbinom <- function(m, n, mu, theta) {
  f_negbinom <- function(x, mu, theta) {
    # Cálculo da densidade
    log_f <- lgamma(x + theta) -
      lgamma(theta) -
      lgamma(x + 1) +
      theta * log(theta / (mu + theta)) +
      x * log(mu / (mu + theta))
    return(exp(log_f))
  }

  g_pois <- function(k, lambda) {
    return(exp(k * log(lambda) - lambda - lgamma(k + 1)))
  }

  source("rpoisson.R")

  y <- rpois(n = m, lambda = 0.9 * mu)

  pesos <- f_negbinom(y, mu = mu, theta = theta) / g_pois(y, lambda = 0.9 * mu)
  pesos <- pesos / sum(pesos)

  amostra <- sample(y, n, F, pesos)

  return(amostra)
}



# --- Parâmetros e amostra ---
set.seed(123)
mu <- 5
theta <- 2
m <- 200000
n <- 10000
amostra <- rnegbinom(m = m, n = n, mu = mu, theta = theta)


#dir.create("figuras", showWarnings = FALSE)
#png("figuras/graf_bn.png", width = 800, height = 600, res = 120)
# --- Histograma das frequências relativas ---
hist(amostra,
     breaks = seq(-0.5, max(amostra) + 0.5, 1),
     freq = FALSE,
     col = "darkblue",
     border = "white",
     main = sprintf("Distribuição Negativa Binomial (μ=%.1f, θ=%.1f): Teórico vs Estimado (n=%d)", mu, theta, n),
     xlab = "x",
     ylab = "Probabilidade / Frequência relativa"
)

# --- Probabilidades teóricas ---
x_vals <- 0:max(amostra)
prob_teo <- dnbinom(x_vals, size = theta, mu = mu)

# --- Sobreposição das probabilidades teóricas ---
points(x_vals, prob_teo, col = "red", pch = 19)
lines(x_vals, prob_teo, col = "red", lwd = 2)

# --- Legenda ---
legend("topright",
       legend = c("Frequência empírica", "Probabilidade teórica"),
       col = c("darkblue", "red"),
       pch = c(15, 19),
       pt.cex = c(1.5, 1.2),
       bty = "n")
#dev.off()

