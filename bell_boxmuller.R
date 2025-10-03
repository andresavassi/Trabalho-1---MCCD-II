# Número de amostras (pares de Normais)
n <- 10000

u1 <- runif(n)
u2 <- runif(n)

z1 <- sqrt(-2 * log(u1)) * cos(2 * pi * u2)
z2 <- sqrt(-2 * log(u1)) * sin(2 * pi * u2)

# Junta as duas amostras
x_bm <- c(z1, z2)

# Histograma
hist(
  x_bm,
  breaks = 50,
  probability = TRUE,
  main = "Normal(0,1) via Box-Muller",
  col = "lightgreen",
  border = "white"
)

# Curva teórica
curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
