# --- 1. CONFIGURAÇÃO E BIBLIOTECAS ---
# Instale os pacotes se necessário:
# install.packages(c("ggplot2", "dplyr", "tidyr", "gridExtra", "xtable"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(xtable)

# Criação da Pasta para Figuras
if (!dir.exists("figuras")) {
  dir.create("figuras")
}
print("Pasta 'figuras/' verificada/criada.")
set.seed(42) # Para reprodutibilidade

# --- 2. FUNÇÕES BASE (WEIBULL 3P E SA) ---

# Função de Log-Verossimilhança Weibull 3P
loglik_weibull3 <- function(params, x) {
  b <- params[1] # beta (forma)
  g <- params[2] # eta (escala)
  c <- params[3] # gamma (locação)

  if (b <= 0 || g <= 0 || any(x <= c)) {
    return(-Inf)
  }

  z <- (x - c) / g
  # Fórmula Log-Verossimilhança (Eq 7 do artigo)
  ll <- length(x) * (log(b) - log(g)) + sum((b - 1) * log(z) - (z^b))
  return(ll)
}

# Função para Gerar Vizinho (Proposta)
propose_neighbor <- function(curr, x, sds = c(0.1, 0.1, 0.1)) {
  attempt <- 0
  min_x <- min(x)

  repeat {
    attempt <- attempt + 1

    b1 <- curr[1] + rnorm(1, mean = 0, sd = sds[1])
    g1 <- curr[2] + rnorm(1, mean = 0, sd = sds[2])
    c1 <- curr[3] + rnorm(1, mean = 0, sd = sds[3])

    if (b1 <= 0) {
      b1 <- abs(b1) + 1e-6
    }
    if (g1 <= 0) {
      g1 <- abs(g1) + 1e-6
    }

    if (c1 < min_x - 1e-8) {
      return(c(b1, g1, c1))
    }

    if (attempt > 50) {
      sds[3] <- sds[3] * 1.5
    }

    if (attempt > 1000) {
      return(c(b1, g1, min_x - 1e-6))
    }
  }
}

# Algoritmo Simulated Annealing para Maximização
sa_weibull3 <- function(
  x,
  init = NULL,
  T0 = 100,
  Tf = 0.001,
  cooling_rate = 0.99,
  L = 5,
  sds = c(0.1, 0.1, 0.1)
) {
  n <- length(x)

  # Lógica de Inicialização (Estimativa Rápida por Regressão Log-Log)
  if (is.null(init)) {
    c0 <- min(x) - 0.1 * sd(x)
    if (c0 >= min(x)) {
      c0 <- min(x) - 1e-3
    }

    y <- x - c0
    y[y <= 0] <- min(y[y > 0]) * 0.1

    y_sorted <- sort(y)
    prob_rank <- (1:n) / (n + 1)

    ln_y <- log(y_sorted)
    ln_ln <- log(-log(1 - prob_rank))

    fit <- try(lm(ln_ln ~ ln_y), silent = TRUE)

    if (inherits(fit, "try-error")) {
      init <- c(1.5, sd(x), min(x) - 0.1)
    } else {
      b0 <- fit$coefficients[2]
      g0 <- exp(-fit$coefficients[1] / b0)
      init <- c(b0, g0, c0)

      if (
        !is.finite(b0) || !is.finite(g0) || !is.finite(c0) || b0 <= 0 || g0 <= 0
      ) {
        init <- c(1.5, sd(x), min(x) - 0.1)
      }
    }
  }

  curr <- init
  curr_ll <- loglik_weibull3(curr, x)
  best <- curr
  best_ll <- curr_ll
  T <- T0

  history <- data.frame(eval = 1, ll = curr_ll)
  eval_count <- 1

  # Loop SA
  while (T > Tf) {
    for (i in 1:L) {
      prop <- propose_neighbor(curr, x, sds = sds)
      prop_ll <- loglik_weibull3(prop, x)

      if (prop_ll > curr_ll) {
        curr <- prop
        curr_ll <- prop_ll
      } else {
        delta_ll <- prop_ll - curr_ll

        if (runif(1) < exp(delta_ll / T)) {
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
    }
    T <- cooling_rate * T
  }

  return(list(
    best_params = best,
    best_ll = best_ll,
    history = history,
    evals = eval_count
  ))
}

# --- 3. FUNÇÃO DE EXECUÇÃO E GERAÇÃO DE SAÍDA (COMPLETA) ---

run_and_save_example <- function(example_num, real_params) {
  sample_sizes <- c(2500, 1000, 500, 100)
  T0_val <- 100
  Tf_val <- 0.001
  cooling_rate_val <- 0.99
  L_val <- 5

  results_table <- data.frame(
    "N" = integer(),
    "beta_est" = numeric(),
    "eta_est" = numeric(),
    "gamma_est" = numeric(),
    "LL_real" = numeric(),
    "LL_est" = numeric(),
    "Time_s" = numeric()
  )
  plot_data <- list()

  # Loop para rodar o SA para cada tamanho de amostra
  for (N in sample_sizes) {
    # 1. Geração da Amostra (Weibull 3P)
    data_x <- rweibull(N, shape = real_params[1], scale = real_params[2]) +
      real_params[3]

    # 2. Log-Verossimilhança com os Parâmetros Reais
    ll_real <- loglik_weibull3(real_params, data_x)

    # 3. Execução do Simulated Annealing
    start_time <- Sys.time()
    sa_result <- sa_weibull3(
      x = data_x,
      T0 = T0_val,
      Tf = Tf_val,
      cooling_rate = cooling_rate_val,
      L = L_val
    )
    end_time <- Sys.time()
    run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # 4. Armazenamento dos Resultados
    best_p <- sa_result$best_params
    best_ll <- sa_result$best_ll

    new_row <- data.frame(
      "N" = N,
      "beta_est" = best_p[1],
      "eta_est" = best_p[2],
      "gamma_est" = best_p[3],
      "LL_real" = ll_real,
      "LL_est" = best_ll,
      "Time_s" = run_time
    )
    results_table <- rbind(results_table, new_row)

    # 5. Armazenamento do Histórico para o Gráfico
    plot_data[[as.character(N)]] <- sa_result$history %>%
      mutate(N_sample = N)
  }

  # --- Geração do Gráfico (PNG) ---

  full_history_df <- do.call(rbind, plot_data)

  panel_map <- data.frame(
    N_sample = c(100, 500, 1000, 2500),
    panel = c("(a)", "(b)", "(c)", "(d)")
  )

  plot_df <- full_history_df %>%
    left_join(panel_map, by = "N_sample") %>%
    mutate(facet_label = paste(panel, "Sample size of", N_sample))

  plot_convergence <- ggplot(plot_df, aes(x = eval, y = ll)) +
    geom_line(color = "darkblue", size = 0.2) +
    facet_wrap(~facet_label, scales = "free", ncol = 2) +
    labs(
      title = paste0(
        "Convergência SA para MLE (Exemplo ",
        example_num,
        ": (",
        paste(real_params, collapse = ", "),
        "))"
      ),
      x = "Número de Avaliações (Iterações)",
      y = "Log-Verossimilhança (ll)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold")
    )

  filename_png <- paste0(
    "figuras/sa_weibull_convergencia_exemplo",
    example_num,
    ".png"
  )
  ggsave(
    filename = filename_png,
    plot = plot_convergence,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300
  )

  print(paste("Gráfico salvo em:", filename_png))

  # --- Geração da Tabela LaTeX (TEX) ---

  latex_table_df <- results_table %>%
    mutate(
      Estimated_Parameters = sprintf(
        "(%.4f, %.4f, %.4f)",
        beta_est,
        eta_est,
        gamma_est
      ),
      LL_real_f = sprintf("%.4f", LL_real),
      LL_est_f = sprintf("%.4f", LL_est),
      Time_s_f = sprintf("%.4f", Time_s)
    ) %>%
    select(N, Estimated_Parameters, LL_real_f, LL_est_f, Time_s_f) %>%
    t() %>%
    as.data.frame()

  colnames(latex_table_df) <- latex_table_df[1, ]
  latex_table_df <- latex_table_df[-1, ]

  latex_table_df <- data.frame(
    "Measure" = c(
      "Estimated parameters (\\(\\beta\\), \\(\\eta\\), \\(\\gamma\\))",
      "Likelihood function at the real values",
      "Likelihood function at the estimated value",
      "Run time (s)"
    ),
    latex_table_df
  )

  caption_text <- paste0(
    "Tabela de Resultados para o Exemplo ",
    example_num,
    ": Weibull (",
    paste(real_params, collapse = ", "),
    ") via Simulated Annealing"
  )

  # Usando xtable para gerar o código LaTeX
  latex_code <- print(
    xtable(
      latex_table_df,
      caption = caption_text,
      align = c("l", "l", rep("r", length(sample_sizes))),
      label = paste0("tab:sa_weibull_results_ex", example_num)
    ),
    type = "latex",
    include.rownames = FALSE,
    floating.environment = "table",
    hline.after = c(-1, 0, nrow(latex_table_df)),
    booktabs = FALSE,
    add.to.row = list(
      pos = list(0, 1),
      command = c(
        paste0(
          "\\toprule\n\\multicolumn{1}{c}{Weibull parameters} & \\multicolumn{4}{c}{",
          paste0("(", paste(real_params, collapse = ", "), ")"),
          "}\\\\\n\\hline\nMeasure & "
        ),
        " \\midrule\n"
      )
    )
  )

  filename_tex <- paste0(
    "figuras/sa_weibull_results_exemplo",
    example_num,
    ".tex"
  )
  writeLines(latex_code, filename_tex)

  print(paste("Tabela LaTeX salva em:", filename_tex))

  return(list(
    png_file = filename_png,
    tex_file = filename_tex
  ))
}

# --- 4. EXECUÇÃO DOS EXEMPLOS (RODAR ESTE BLOCO) ---

# Exemplo 1: theta = (2, 2, 2)
real_params_ex1 <- c(2, 2, 2)
print(paste(
  "--- Executando Exemplo 1:",
  paste(real_params_ex1, collapse = ", "),
  "---"
))
results_ex1 <- run_and_save_example(1, real_params_ex1)

# Exemplo 2: theta = (2, 3, 4)
real_params_ex2 <- c(2, 3, 4)
print(paste(
  "--- Executando Exemplo 2:",
  paste(real_params_ex2, collapse = ", "),
  "---"
))
results_ex2 <- run_and_save_example(2, real_params_ex2)

# Exemplo 3: theta = (3, 2, 5)
real_params_ex3 <- c(3, 2, 5)
print(paste(
  "--- Executando Exemplo 3:",
  paste(real_params_ex3, collapse = ", "),
  "---"
))
results_ex3 <- run_and_save_example(3, real_params_ex3)

print(
  "Processamento concluído para os 3 exemplos. Verifique a pasta 'figuras/'"
)
