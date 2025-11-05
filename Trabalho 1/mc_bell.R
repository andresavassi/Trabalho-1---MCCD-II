MC <- 1000
nsizes <- c(50, 100, 500, 1000)


library(dplyr)
library(kableExtra)
library(tibble)
library(foreach)
library(doParallel)
library(ggplot2)

source("rbell.R")

# ------------------- Gera teste
set.seed(9999)
df_teste <- data.frame(
  x1 = sample(c(0, 1), 10, replace = T, prob = c(0.7, 0.3)),
  x2 = rnorm(10),
  x3 = rnorm(10),
  x4 = sample(c(0, 1), 10, replace = T, prob = c(0.2, 0.8))
)
# -----------------------

run_MC <- function(r, teste, nsize) {
  set.seed(r)
  n <- nsize
  # Parametros reais
  betas <- c(1.2, 0.25, -0.08, 0.15, -0.12)
  X_teste <- model.matrix(~., teste)
  eta_teste_real <- X_teste %*% betas
  fitted_teste_real <- exp(eta_teste_real)

  # Gerando Covariáveis
  df <- data.frame(
    x1 = sample(c(0, 1), n, replace = T, prob = c(0.7, 0.3)),
    x2 = rnorm(n),
    x3 = rnorm(n),
    x4 = sample(c(0, 1), n, replace = T, prob = c(0.2, 0.8))
  )
  X <- model.matrix(~., df)

  # Calculo dos preditores lineares e médias
  eta <- X %*% betas
  mu <- exp(eta)

  # Gerando y
  Y <- rbell(1, lamW::lambertW0(mu))

  # Ajusta modelo com dados simulados
  df <- cbind(df, Y)
  fit <- glm(Y ~ ., family = "poisson", data = df)

  # Recuperando estimativa dos parametros e intervalos de confiança
  betas_est <- coef(fit)
  sigma <- vcov(fit)

  vars <- diag(X_teste %*% sigma %*% t(X_teste)) # Metodo delta

  d <- stats::qnorm(0.975) * sqrt(vars)

  eta_teste <- X_teste %*% betas_est
  fitted_teste <- exp(eta_teste)
  upper_teste <- exp(eta_teste + d)
  lower_teste <- exp(eta_teste - d)

  par <- paste0("y", 1:nrow(X_teste))

  tb <- tibble::tibble(
    par = par,
    real = fitted_teste_real,
    estimates = fitted_teste,
    se = vars,
    rb = 100 * (estimates - real) / abs(real),
    li = lower_teste,
    ls = upper_teste,
    cp = real > li & real < ls,
    nsize = nsize,
    rep = r
  ) |>
    dplyr::relocate(nsize, .before = "par")

  return(tb)
}


# Criar todas as combinações de MC e nsizes
params <- expand.grid(
  i = 1:MC,
  n = nsizes
)

# Paralelismo
cl <- makeCluster(detectCores())
registerDoParallel(cl)

t1 <- Sys.time()
tb <- foreach(
  row = 1:nrow(params),
  .combine = "rbind",
  .packages = c("dplyr", "tibble")
) %dopar%
  {
    run_MC(
      r = params$i[row],
      teste = df_teste,
      nsize = params$n[row]
    )
  }
t2 <- Sys.time()

stopCluster(cl)

# Tempo de execução
t2 - t1


names(tb) <- c(
  "nsize",
  "par",
  "real",
  "estimate",
  "se",
  "RB",
  "lwr",
  "upr",
  "cp",
  "rep"
)


tb <- tb |>
  mutate(
    nsize = as.factor(nsize)
  )


tbl <- tb |>
  group_by(nsize, par) |>
  summarise(
    sde = sd(estimate),
    across(
      c("real", "estimate", "se", "RB", "lwr", "upr", "cp"),
      ~ mean(.x, na.rm = TRUE)
    )
  ) |>
  relocate(sde, .after = se) |>
  arrange(nsize)


# Tabela para n = 100

tabela_kbl <- tbl[11:20, ] %>%
  kable(
    format = "html",
    caption = NULL,
    col.names = c(
      "Tamanho da Amostra",
      "Parâmetro",
      "Valor Real",
      "Estimativa Média",
      "Erro Padrão",
      "Desvio Padrão",
      "Viés Relativo",
      "Limite Inferior",
      "Limite Superior",
      "Cobertura"
    ),
    digits = 3,
    align = "c"
  )


# Gráficos

# ----------------------
# Viés Relativo

p1 <- ggplot(df, aes(x = par, y = RB, fill = nsize)) +
  geom_boxplot(alpha = 0.7, outlier.color = "red", outlier.shape = 21) +
  geom_hline(yintercept = 0, col = "blue", linetype = "dashed", linewidth = 1) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Distribuição do Viés Relativo (RB) por Parâmetro",
    subtitle = "Para Modelo log-linear Poisson (Variável resposta gerada da Bel): Comparação entre tamanhos de amostra",
    x = "Parâmetro",
    y = "Viés Relativo (RB)",
    fill = "Tamanho da Amostra"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray40"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

# -----------------------

# Viés relativo dividido em tamanhos amostrais
p2 <- ggplot(df, aes(x = par, y = RB, fill = par)) +
  geom_boxplot(alpha = 0.8, outlier.color = "red", outlier.shape = 21) +
  geom_hline(yintercept = 0, col = "blue", linetype = "dashed", linewidth = 1) +
  scale_fill_brewer(palette = "Pastel1") +
  facet_wrap(~nsize, ncol = 2) +
  labs(
    title = "Viés Relativo (RB) por Parâmetro e Tamanho da Amostra",
    subtitle = "Para Modelo log-linear Poisson (Variável resposta gerada da Bell): Boxplots separados por tamanho de amostra",
    x = "Parâmetro",
    y = "Viés Relativo (RB)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray40"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#2E86C1", color = NA),
    strip.text = element_text(color = "white", face = "bold")
  )

# --------------------------

# Salva os gráficos em .png

ggsave("grafico_bell1.png", plot = p1, width = 8, height = 6, dpi = 300)
ggsave("grafico_bell2.png", plot = p2, width = 8, height = 6, dpi = 300)

# --------------------------
