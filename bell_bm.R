set.seed(20252)
# Número de amostras (pares de Normais)
n <- 10000

u1 <- runif(n)
u2 <- runif(n)

z1 <- sqrt(-2 * log(u1)) * cos(2 * pi * u2)
z2 <- sqrt(-2 * log(u1)) * sin(2 * pi * u2)

# Junta as duas amostras
x_bm <- c(z1, z2)


dir.create("figuras", showWarnings = FALSE)
png("figuras/graf_bell_bm.png", width = 800, height = 600, res = 120)
# Histograma
hist(
  x_bm,
  breaks = 50,
  probability = TRUE,
  main = "Normal(0,1) via Box-Muller",
  col = "darkblue",
  border = "white"
)

# Curva teórica
curve(dnorm(x), add = TRUE, col = "red", lwd = 2)

dev.off()
