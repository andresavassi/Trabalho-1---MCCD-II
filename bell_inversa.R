# Número de amostras
n <- 10000

# Gera Uniformes(0,1)
u <- runif(n)

# Transformação inversa via função quantil
x_inv <- qnorm(u, mean = 0, sd = 1)

# Histograma
hist(
  x_inv,
  breaks = 50,
  probability = TRUE,
  main = "Normal(0,1) via Transformação Inversa (qnorm)",
  col = "lightblue",
  border = "white"
)

# Curva teórica
curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
