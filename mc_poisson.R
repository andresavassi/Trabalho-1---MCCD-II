MC <- 1000
nsize <- 100

library(foreach)
library(dplyr)
library(tibble)
library(doParallel)


set.seed(9999)
df_teste <- data.frame(
  x1 = sample(c(0, 1), 10, replace = T, prob = c(0.7, 0.3)),
  x2 = rnorm(10),
  x3 = rnorm(10),
  x4 = sample(c(0, 1), 10, replace = T, prob = c(0.2, 0.8))
)

X_teste <- model.matrix(~., df_teste)

results <- data.frame()

run_MC <- function(r, nsize) {
  set.seed(r)
  n <- nsize
  betas <- c(1.2, 0.25, -0.08, 0.15, -0.12)

  df <- data.frame(
    x1 = sample(c(0, 1), n, replace = T, prob = c(0.7, 0.3)),
    x2 = rnorm(n),
    x3 = rnorm(n),
    x4 = sample(c(0, 1), n, replace = T, prob = c(0.2, 0.8))
  )
  X <- model.matrix(~., df)
  eta <- X %*% betas
  mu <- exp(eta)

  Y <- rpois(n = n, lambda = mu)

  df <- cbind(df, Y)

  fit <- glm(Y ~ ., family = "poisson", data = df)

  betas_est <- coef(fit)
  sigma <- vcov(fit)

  # Metodo delta
  vars <- diag(X_teste %*% sigma %*% t(X_teste))

  d <- stats::qnorm(0.975) * sqrt(vars)

  eta_teste <- X_teste %*% betas_est
  fitted_teste <- exp(eta_teste)
  upper_teste <- exp(eta_teste + d)
  lower_teste <- exp(eta_teste - d)

  eta_teste_real <- X_teste %*% betas
  fitted_teste_real <- exp(eta_teste_real)

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


cl <- makeCluster(detectCores()) # parallel
registerDoParallel(cl)

t1 <- Sys.time()
tb <- foreach(i = 1:MC, .combine = "rbind") %dopar%
  {
    run_MC(i, nsize = n)
  }
t2 <- Sys.time()
t2 - t1

tb
