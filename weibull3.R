# Instalar se necessário
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

# --- Função densidade da Weibull 3-parâmetros ---
dweibull3 <- function(x, shape, scale, location = 0) {
  x_adj <- x - location
  ifelse(
    x_adj > 0,
    (shape / scale) * (x_adj / scale)^(shape - 1) * exp(-(x_adj / scale)^shape),
    0
  )
}

# --- Parâmetros ---
alpha <- 1.2 # shape
beta <- 2.0 # scale
gamma <- 1.5 # location

# --- Grade de x ---
x <- seq(gamma - 0.5, gamma + 10, length.out = 1000)
dens <- dweibull3(x, shape = alpha, scale = beta, location = gamma)
df <- data.frame(x, dens)

# --- Gráfico ---
ggplot(df, aes(x, y = dens)) +
  geom_line(size = 1.2, color = "blue") +
  labs(
    title = "Função Densidade - Weibull Tri-Parâmetrica",
    subtitle = paste0(
      "shape = ",
      alpha,
      ", scale = ",
      beta,
      ", location = ",
      gamma
    ),
    x = "x",
    y = "f(x)"
  ) +
  theme_minimal(base_size = 14)
