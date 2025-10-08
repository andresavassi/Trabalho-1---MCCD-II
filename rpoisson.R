poisson_ar_1 <- function(lambda) {
  f_pois <- function(k, lambda) {
    return(exp(k * log(lambda) - lambda - lgamma(k + 1)))
  }
  fmax <- f_pois(floor(lambda), lambda)
  kmax <- qpois(0.9999, lambda = 5)
  repeat {
    y <- sample(0:kmax, 1)
    fy <- f_pois(y, lambda)
    gy <- 1 / (kmax + 1)
    u <- runif(1)
    if (u < fy / (fmax)) {
      return(y)
    }
  }
}
poisson_ar <- function(lambda, n = 1) {
  replicate(n, poisson_ar_1(lambda))
}

#Conferindo resultados
#lambda <- 5
#amostra <- poisson_ar(n = 10000, lambda = lambda)

#tab <- rbind(
#  round(prop.table(table(amostra)), 3),
#  round(dpois(0:max(amostra), lambda), 3)
#)
#colnames(tab) <- 0:max(amostra)
#rownames(tab) <- c("freq", "prob")

#tab <- tab[, 1:10] %>%
#  kbl(
#    align = c(rep("c", 5)),
#    caption = "Comparação de frequências relativas com probabilidades teóricas"
#  )
#tab %>%
#  kable_styling(
#    latex_options = "hold_position",
#    bootstrap_options = c("bordered"),
#    font_size = 8
#  )



# --- Função auxiliar: densidade Poisson ---
f_pois <- function(k, lambda) {
  exp(k * log(lambda) - lambda - lgamma(k + 1))
}

# --- Geração por aceitação e rejeição ---
poisson_ar_1 <- function(lambda) {
  fmax <- f_pois(floor(lambda), lambda)  # valor máximo da f(k)
  kmax <- qpois(0.9999, lambda = lambda) # valor máximo de amostragem
  repeat {
    y <- sample(0:kmax, 1)
    fy <- f_pois(y, lambda)
    u <- runif(1)
    if (u < fy / fmax) {
      return(y)
    }
  }
}

poisson_ar <- function(lambda, n = 1) {
  replicate(n, poisson_ar_1(lambda))
}

# --- Parâmetros e amostra ---
set.seed(20252)
lambda <- 5
n <- 10000
amostra <- poisson_ar(lambda = lambda, n = n)
dir.create("figuras", showWarnings = FALSE)
png("figuras/graf_poisson.png", width = 800, height = 600, res = 120)
# --- Histograma das frequências relativas ---
hist(amostra,
     breaks = seq(-0.5, max(amostra) + 0.5, 1),
     freq = FALSE,
     col = "darkblue",
     border = "white",
     main = sprintf("Distribuição Poisson(λ = %.1f): Teórico vs Estimado (n = %d)", lambda, n),
     xlab = "x",
     ylab = "Probabilidade / Frequência relativa"
)

# --- Probabilidades teóricas ---
x_vals <- 0:max(amostra)
prob_teo <- dpois(x_vals, lambda)

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

dev.off()
