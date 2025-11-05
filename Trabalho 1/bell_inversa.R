set.seed(20252)
# Número de amostras
n <- 10000

# Gera Uniformes(0,1)
u <- runif(n)

# Transformação inversa via função quantil
x_inv <- qnorm(u, mean = 0, sd = 1)

dir.create("figuras", showWarnings = FALSE)
png("figuras/graf_bell_inv.png", width = 800, height = 600, res = 120)
# Histograma
hist(
  x_inv,
  breaks = 50,
  probability = TRUE,
  main = "Normal(0,1) via Transformação Inversa (qnorm)",
  col = "darkblue",
  border = "white"
)

# Curva teórica
curve(dnorm(x), add = TRUE, col = "red", lwd = 2)


dev.off()
