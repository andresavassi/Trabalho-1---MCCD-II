# sa_weibull.R
# Simulated Annealing para estimar os parâmetros (b, g, c) da Weibull 3-param
# Implementação segundo Abbasi et al. (2006).
# Gera tabelas e gráficos para os três exemplos do artigo.

set.seed(20252) # reprodutibilidade

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# --- Funções úteis ---------------------------------------------------------

# log-verossimilhança (soma) para Weibull 3-param: b (shape), g (scale), c (location)
# x: vetor de observações. retorna -Inf se qualquer xi <= c or b<=0 or g<=0
loglik_weibull3 <- function(params, x) {
  b <- params[1]
  g <- params[2]
  c <- params[3]
  if (b <= 0 || g <= 0 || any(x <= c)) {
    return(-Inf)
  }
  z <- (x - c) / g
  ll <- length(x) * (log(b) - log(g)) + (b - 1) * sum(log(z)) - sum(z^b)
  return(ll)
}

# Gera vizinho (proposta) garantindo g>0, b>0, c < min(x) (senão rejeita)
propose_neighbor <- function(curr, x, sds = c(0.1, 0.1, 0.1)) {
  # curr: c(b,g,c)
  attempt <- 0
  repeat {
    attempt <- attempt + 1
    b1 <- curr[1] + rnorm(1, mean = 0, sd = sds[1])
    g1 <- curr[2] + rnorm(1, mean = 0, sd = sds[2])
    c1 <- curr[3] + rnorm(1, mean = 0, sd = sds[3])
    # enforce positivity of b,g by reflecting small negatives
    if (b1 <= 0) {
      b1 <- abs(b1) + 1e-6
    }
    if (g1 <= 0) {
      g1 <- abs(g1) + 1e-6
    }
    # require c1 < min(x) - tiny
    if (c1 < min(x) - 1e-8) {
      return(c(b1, g1, c1))
    }
    # else try again; after many attempts increase sd for c
    if (attempt > 50) {
      sds[3] <- sds[3] * 1.5
    }
    if (attempt > 1000) {
      # extremely unlikely; fallback: set c1 = min(x) - epsilon
      return(c(b1, g1, min(x) - 1e-6))
    }
  }
}

# SA algorithm (maximização da log-verossimilhança)
sa_weibull3 <- function(
  x,
  init = NULL,
  T0 = 100,
  Tf = 0.001,
  cooling_rate = 0.99,
  L = 5,
  sds = c(0.1, 0.1, 0.1),
  verbose = FALSE
) {
  n <- length(x)
  # initial solution: if not provided, use method of moments-ish guesses
  if (is.null(init)) {
    # use rweibull fit approximations: estimate c0 slightly below min(x)
    c0 <- min(x) - 0.1 * sd(x)
    if (c0 >= min(x)) {
      c0 <- min(x) - 1e-3
    }
    # fit two-parameter Weibull to x - c0 via linearization (quick guess)
    y <- x - c0
    y[y <= 0] <- min(y[y > 0]) * 0.1
    # estimate using log-log regression
    est <- try(
      {
        ln_x <- log(y)
        ln_ln <- log(-log(1 - (1:n) / (n + 1)))
        fit <- lm(ln_x ~ ln_ln)
        b0 <- 1 / coef(fit)[2]
        g0 <- exp(coef(fit)[1])
        c(b0, g0, c0)
      },
      silent = TRUE
    )
    if (inherits(est, "try-error") || any(!is.finite(est))) {
      est <- c(1.5, sd(x), min(x) - 0.1)
    }
    curr <- est
  } else {
    curr <- init
  }

  curr_ll <- loglik_weibull3(curr, x)
  best <- curr
  best_ll <- curr_ll

  T <- T0
  history <- data.frame(eval = 1, ll = curr_ll)
  eval_count <- 1

  while (T > Tf) {
    ninner <- 1
    while (ninner <= L) {
      prop <- propose_neighbor(curr, x, sds = sds)
      prop_ll <- loglik_weibull3(prop, x)
      # If better, accept; else accept with probability exp((prop_ll - curr_ll)/T)
      if (prop_ll > curr_ll) {
        curr <- prop
        curr_ll <- prop_ll
      } else {
        delta <- prop_ll - curr_ll
        if (runif(1) < exp(delta / T)) {
          curr <- prop
          curr_ll <- prop_ll
        }
      }
      if (curr_ll > best_ll) {
        best <- curr
        best_ll <- curr_ll
      }
      eval_count <- eval_count + 1
      history <- rbind(history, data.frame(eval = eval_count, ll = curr_ll))
      ninner <- ninner + 1
    }
    T <- cooling_rate * T
  }
  return(list(
    best = best,
    best_ll = best_ll,
    history = history,
    evals = eval_count
  ))
}

# --- Experimentos como no artigo ------------------------------------------

run_experiment <- function(
  true_params,
  sample_sizes = c(2500, 1000, 500, 100),
  T0 = 100,
  Tf = 0.001,
  C = 0.99,
  L = 5,
  sds = c(0.1, 0.1, 0.1),
  seed_base = 20252,
  plot_prefix = "Fig"
) {
  # true_params: c(b, g, c)
  results_list <- list()
  plots <- list()
  for (n in sample_sizes) {
    set.seed(seed_base + n) # differ seeds by sample size
    b_true <- true_params[1]
    g_true <- true_params[2]
    c_true <- true_params[3]
    # generate three-parameter Weibull sample: x = c + g * rweibull(n, shape = b, scale=1)
    x <- c_true + g_true * rweibull(n, shape = b_true, scale = 1)
    # initialize starting point randomly away from truth
    init <- c(
      abs(true_params[1] * runif(1, 0.6, 1.4)),
      abs(true_params[2] * runif(1, 0.6, 1.4)),
      min(x) - runif(1, 0.05, 0.2) * sd(x)
    )
    t0 <- Sys.time()
    sa_res <- sa_weibull3(
      x,
      init = init,
      T0 = T0,
      Tf = Tf,
      cooling_rate = C,
      L = L,
      sds = sds,
      verbose = FALSE
    )
    t1 <- Sys.time()
    runtime <- as.numeric(difftime(t1, t0, units = "secs"))
    est_params <- sa_res$best
    ll_true <- loglik_weibull3(true_params, x)
    ll_est <- sa_res$best_ll
    results_list[[as.character(n)]] <- list(
      n = n,
      est = est_params,
      ll_true = ll_true,
      ll_est = ll_est,
      runtime = runtime,
      history = sa_res$history,
      x = x
    )

    # Plot history similarly to article: iterations vs log-likelihood
    df_h <- sa_res$history
    p <- ggplot(df_h, aes(x = eval, y = ll)) +
      geom_line() +
      geom_hline(yintercept = ll_true, linetype = "dashed", alpha = 0.8) +
      labs(
        title = paste0(plot_prefix, ": n = ", n),
        x = "Avaliações (iteração)",
        y = "Log-verossimilhança"
      ) +
      theme_minimal()
    plots[[as.character(n)]] <- p

    message(sprintf(
      "Done: n=%d | est=(%.4f, %.4f, %.4f) | ll_true=%.4f | ll_est=%.4f | time=%.2fs",
      n,
      est_params[1],
      est_params[2],
      est_params[3],
      ll_true,
      ll_est,
      runtime
    ))
  }

  # assemble table like artigo
  tab <- do.call(
    rbind,
    lapply(results_list, function(z) {
      c(
        n = z$n,
        est_b = z$est[1],
        est_g = z$est[2],
        est_c = z$est[3],
        ll_true = z$ll_true,
        ll_est = z$ll_est,
        runtime_s = z$runtime
      )
    })
  )
  tab <- as.data.frame(tab)
  tab$n <- as.integer(tab$n)
  rownames(tab) <- NULL

  # order rows descending by n to match article presentation
  tab <- tab %>% arrange(desc(n))

  return(list(results = results_list, table = tab, plots = plots))
}

# --- Executa para os três exemplos do artigo -------------------------------

examples <- list(
  ex1 = c(2, 2, 2),
  ex2 = c(2, 3, 4),
  ex3 = c(3, 2, 5)
)

sample_sizes <- c(2500, 1000, 500, 100)
all_results <- list()
all_tables <- list()

# Adjust sds for neighbor proposals (tuned empirically; você pode ajustar se quiser)
sds_default <- c(0.05, 0.05, 0.05) # menor para maior estabilidade

for (i in seq_along(examples)) {
  ex_name <- names(examples)[i]
  cat(
    "\n== Executando",
    ex_name,
    "com parâmetros verdadeiros:",
    examples[[i]],
    "==\n"
  )
  res <- run_experiment(
    examples[[i]],
    sample_sizes = sample_sizes,
    T0 = 100,
    Tf = 0.001,
    C = 0.99,
    L = 5,
    sds = sds_default,
    seed_base = 20252 + i * 10,
    plot_prefix = paste0("Example ", i)
  )
  all_results[[ex_name]] <- res$results
  all_tables[[ex_name]] <- res$table

  # Salva gráfico composto (4 subplots) semelhante às Figuras do artigo
  # colocamos as 4 curvas lado a lado (100, 500, 1000, 2500)
  glist <- lapply(as.character(sample_sizes), function(ss) res$plots[[ss]])
  # reordenar para mostrar (100,500,1000,2500) left-to-right as in paper they show (a)(b)(c)(d)
  # but paper arranged perhaps 100,500,1000,2500 as (a)-(d); aqui vamos organizar small->large
  png(filename = paste0("SA_example_", i, ".png"), width = 1400, height = 600)
  grid.arrange(grobs = glist, ncol = 4)
  dev.off()

  # Print table to console and save CSV
  print(res$table)
  write.csv(
    res$table,
    file = paste0("table_example_", i, ".csv"),
    row.names = FALSE
  )
}

# Salva todas as tabelas em um único arquivo
all_tables_df <- bind_rows(lapply(names(all_tables), function(nm) {
  df <- all_tables[[nm]]
  df$example <- nm
  df
}))
write.csv(all_tables_df, "all_tables_summary.csv", row.names = FALSE)

cat("\nExecução completa. Arquivos gerados:\n")
cat("- SA_example_1.png, SA_example_2.png, SA_example_3.png\n")
cat("- table_example_1.csv, table_example_2.csv, table_example_3.csv\n")
cat("- all_tables_summary.csv\n")
