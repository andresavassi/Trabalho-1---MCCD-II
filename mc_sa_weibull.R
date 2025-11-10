run_example <- function(real_params) {
  sample_sizes <- 2500
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
  
  # Loop para rodar o SA para cada tamanho de amostra
  for (N in sample_sizes) {
    # 1. Geração da Amostra (Weibull 3P)
    data_x <- rweibull(N, shape = real_params[1], scale = real_params[2]) +
      real_params[3]
    
    # 2. Log-Verossimilhança com os Parâmetros Reais
    ll_real <- loglik_weibull3(real_params, data_x)
    
    # 3. Execução do Simulated Annealing
    start_time <- Sys.time()
    
    # Inicialização (Passada explicitamente para run_and_save_example)
    init_params <- c(1.5, sd(data_x), min(data_x) - 0.1)
    
    sa_result <- sa_weibull3(
      x = data_x,
      init = init_params,
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
    return(results_table)
  }}



estudo_mc <- function(n = 2500, real_params, m = 100) {
  
  beta_est <- numeric(m)
  eta_est <- numeric(m)
  gamma_est <- numeric(m)
  
  for(i in 1:m){
    
    
    # Rodar EM
    est <- run_example(real_params)
    
    
    beta_est[i] <- tail(est$beta, 1)
    eta_est[i] <- tail(est$eta, 1)
    gamma_est[i] <- tail(est$gamma, 1)
  }
  
  return(data.frame(beta_est, eta_est, gamma_est))
}


av1 <- estudo_mc(real_params = c(2,2,2))
av2 <- estudo_mc(real_params = c(3,5,7))
av3 <- estudo_mc(real_params = c(8,4,6))



avaliacao <- function(av, real_params){
  data.frame(
    beta_verdadeiro = real_params[1],
    beta_medio = mean(av$beta_est),
    beta_sd = sd(av$beta_est),
    eta_verdadeiro = real_params[2],
    eta_medio = mean(av$eta_est),
    eta_sd = sd(av$eta_est),
    gamma_verdadeiro = real_params[3],
    gamma_medio = mean(av$gamma_est),
    gamma_sd = sd(av$gamma_est)
  )
}


tab_resultados <- rbind(
  avaliacao(av1, c(2,2,2)),
  avaliacao(av2, c(3,5,7)),
  avaliacao(av3, c(8,4,6))
)



colnames(tab_resultados) <- c("\beta Verdadeiro", "\beta Médio", "\beta Desvio Padrão",
                              "\eta Verdadeiro", "\eta Médio", "\eta Desvio Padrão",
                              "\gamma Verdadeiro", "\gamma Médio", "\gamma Desvio Padrão")

kable(tab_resultados,
      digits = 4,
      align = "c",
      caption = "Avaliação dos Estimadores - Método de Monte Carlo") %>%
  kable_styling(bootstrap_options = "striped")
