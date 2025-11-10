# =======================================================
# Simulated Annealing para Mistura de Duas Normais
# Estrutura padronizada (igual à Weibull 3P)
# =======================================================

set.seed(20252)

# --- 1. CONFIGURAÇÃO E BIBLIOTECAS ---

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(xtable)

if (!dir.exists("figuras")) dir.create("figuras")

# --- 2. FUNÇÕES BASE (GMM e SA) ---

# Função para gerar dados da mistura
rmixed_norm <- function(n, mu, sigma2, p) {
  sigma1 <- sqrt(sigma2[1])
  sigma2 <- sqrt(sigma2[2])
  s <- rbinom(n, 1, p)
  x <- s * rnorm(n, mu[1], sigma1) + (1 - s) * rnorm(n, mu[2], sigma2)
  return(x)
}

# Log-verossimilhança (negativa, pois SA minimiza)
loglik_gmm <- function(x, params) {
  delta <- params[1]
  mu1 <- params[2]
  mu2 <- params[3]
  sigma2_1 <- params[4]
  sigma2_2 <- params[5]
  
  if (delta < 0 || delta > 1 || sigma2_1 <= 0 || sigma2_2 <= 0)
    return(1e+100)
  
  f1 <- dnorm(x, mean = mu1, sd = sqrt(sigma2_1))
  f2 <- dnorm(x, mean = mu2, sd = sqrt(sigma2_2))
  L <- sum(log(delta * f1 + (1 - delta) * f2))
  
  return(-L) # minimiza
}

# Função de vizinhança (perturbação controlada)
perturbar <- function(params, step_size = 0.5) {
  params_new <- params + rnorm(length(params), 0, step_size)
  params_new[1] <- max(0.001, min(0.999, params_new[1]))
  params_new[4] <- max(0.01, params_new[4])
  params_new[5] <- max(0.01, params_new[5])
  return(params_new)
}

# Simulated Annealing para GMM
sa_gmm <- function(x, T0 = 10, Tf = 0.001, alpha = 0.99, step_init = 0.5, L = 5) {
  mu_mean <- mean(x)
  mu_sd <- sd(x)
  
  theta <- c(
    runif(1, 0.2, 0.8),
    rnorm(1, mu_mean - 0.5 * mu_sd, 0.5),
    rnorm(1, mu_mean + 0.5 * mu_sd, 0.5),
    runif(1, 0.5, 2),
    runif(1, 0.5, 2)
  )
  
  E_current <- loglik_gmm(x, theta)
  E_best <- E_current
  theta_best <- theta
  T <- T0
  step_size <- step_init
  
  history <- data.frame(eval = 1, energy = E_current)
  eval_count <- 1
  
  while (T > Tf) {
    for (i in 1:L) {
      theta_new <- perturbar(theta, step_size)
      E_new <- loglik_gmm(x, theta_new)
      delta_E <- E_new - E_current
      
      if (delta_E <= 0 || runif(1) < exp(-delta_E / T)) {
        theta <- theta_new
        E_current <- E_new
      }
      
      if (E_current < E_best) {
        E_best <- E_current
        theta_best <- theta
      }
      
      eval_count <- eval_count + 1
      history <- rbind(history, data.frame(eval = eval_count, energy = E_current))
    }
    
    T <- alpha * T
    step_size <- step_size * 0.999
  }
  
  return(list(best_params = theta_best, best_energy = E_best, history = history))
}

# --- 3. FUNÇÃO DE EXECUÇÃO (E GERAÇÃO DE GRÁFICOS/TABELAS) ---

run_and_save_example_gmm <- function(example_num, real_params) {
  sample_sizes <- c(100, 500, 1000, 2500)
  results_table <- data.frame()
  plot_data <- list()
  
  for (N in sample_sizes) {
    x <- rmixed_norm(N, real_params$mu, real_params$sigma2, real_params$p)
    ll_real <- -loglik_gmm(x, c(real_params$p, real_params$mu, real_params$sigma2))
    
    start_time <- Sys.time()
    sa_result <- sa_gmm(x)
    end_time <- Sys.time()
    
    runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
    best_p <- sa_result$best_params
    best_ll <- -sa_result$best_energy
    
    new_row <- data.frame(
      N = N,
      p_est = best_p[1],
      mu1_est = best_p[2],
      mu2_est = best_p[3],
      sigma1_est = best_p[4],
      sigma2_est = best_p[5],
      LL_real = ll_real,
      LL_est = best_ll,
      Time_s = runtime
    )
    results_table <- rbind(results_table, new_row)
    
    plot_data[[as.character(N)]] <- sa_result$history %>%
      mutate(N_sample = N)
  }
  
  # --- GRÁFICO DE CONVERGÊNCIA ---
  full_hist <- do.call(rbind, plot_data)
  panel_map <- data.frame(N_sample = c(100, 500, 1000, 2500),
                          panel = c("(a)", "(b)", "(c)", "(d)"))
  
  plot_df <- left_join(full_hist, panel_map, by = "N_sample") %>%
    mutate(facet_label = paste(panel, "Sample size:", N_sample))
  
  plot_convergence <- ggplot(plot_df, aes(x = eval, y = -energy)) +
    geom_line(color = "darkred", size = 0.2) +
    facet_wrap(~facet_label, scales = "free", ncol = 2) +
    labs(
      title = paste0("Convergência SA - Mistura Normal (Exemplo ", example_num, ")"),
      x = "Iterações",
      y = "Log-Verossimilhança"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), strip.text = element_text(face = "bold"))
  
  filename_png <- paste0("figuras/sa_gmm_convergencia_exemplo", example_num, ".png")
  ggsave(filename_png, plot_convergence, width = 10, height = 8, dpi = 300)
  print(paste("Gráfico salvo em:", filename_png))
  
  # --- TABELA LaTeX ---
  latex_df <- results_table %>%
    mutate(
      Est_Parametros = sprintf("(%.3f, %.3f, %.3f, %.3f, %.3f)", p_est, mu1_est, mu2_est, sigma1_est, sigma2_est),
      LL_real = sprintf("%.3f", LL_real),
      LL_est = sprintf("%.3f", LL_est),
      Tempo = sprintf("%.3f", Time_s)
    ) %>%
    dplyr::select(N, Est_Parametros, LL_real, LL_est, Tempo) %>%
    t() %>% as.data.frame()
  
  colnames(latex_df) <- latex_df[1, ]
  latex_df <- latex_df[-1, ]
  
  latex_df <- data.frame(
    "Medida" = c(
      "Parâmetros Estimados (p, μ₁, μ₂, σ₁², σ₂²)",
      "Log-verossimilhança real",
      "Log-verossimilhança estimada",
      "Tempo (s)"
    ),
    latex_df
  )
  
  caption_text <- paste0(
    "Resultados SA - Mistura Normal (Exemplo ", example_num, ") ",
    "com parâmetros reais: p=", real_params$p,
    ", μ=(", paste(real_params$mu, collapse=", "),
    "), σ²=(", paste(real_params$sigma2, collapse=", "), ")"
  )
  
  latex_code <- print(
    xtable(
      latex_df,
      caption = caption_text,
      align = c("l", "l", rep("r", length(sample_sizes))),
      label = paste0("tab:sa_gmm_results_ex", example_num)
    ),
    type = "latex",
    include.rownames = FALSE
  )
  
  filename_tex <- paste0("figuras/sa_gmm_results_exemplo", example_num, ".tex")
  writeLines(latex_code, filename_tex)
  print(paste("Tabela LaTeX salva em:", filename_tex))
  
  return(list(png_file = filename_png, tex_file = filename_tex))
}

# --- 4. EXECUÇÃO DOS EXEMPLOS ---

ex1 <- list(p = 0.7, mu = c(0, 5), sigma2 = c(1, 2))
print("Executando Exemplo 1 (Mistura Normal)")
res1 <- run_and_save_example_gmm(1, ex1)

ex2 <- list(p = 0.5, mu = c(1, 6), sigma2 = c(2, 4))
print("Executando Exemplo 2 (Mistura Normal)")
res2 <- run_and_save_example_gmm(2, ex2)

print("Processamento concluído. Verifique a pasta 'figuras/' para os resultados.")
