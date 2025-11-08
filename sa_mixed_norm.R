set.seed(20252) # Semente para reprodutibilidade

# =======================================================
# 1. FUNÇÕES AUXILIARES E ALGORITMO SA
# =======================================================

# Função para gerar dados (fornecida no prompt)
rmixed_norm <- function(n, mu, sigma2, p) {
  sigma1 <- sqrt(sigma2[1])
  sigma2 <- sqrt(sigma2[2])
  s <- rbinom(n, 1, p)
  x <- s * rnorm(n, mu[1], sigma1) + (1 - s) * rnorm(n, mu[2], sigma2)
  return(x)
}

# Função de Energia (Negativo da Log-Verossimilhança)
loglik_gmm <- function(x, params) {
  delta <- params[1]
  mu1 <- params[2]
  mu2 <- params[3]
  sigma2_1 <- params[4]
  sigma2_2 <- params[5]

  # Penalidade alta para parâmetros inválidos
  if (delta < 0 || delta > 1 || sigma2_1 <= 0 || sigma2_2 <= 0) {
    return(1e+100)
  }

  f1 <- dnorm(x, mean = mu1, sd = sqrt(sigma2_1))
  f2 <- dnorm(x, mean = mu2, sd = sqrt(sigma2_2))

  L <- sum(log(delta * f1 + (1 - delta) * f2))

  # SA busca minimizar: Retorna o negativo
  return(-L)
}

# Função de Perturbação (Vizinhança)
perturbar <- function(params, step_size = 0.5) {
  params_new <- params + rnorm(length(params), 0, step_size)

  # Garantir restrições:
  params_new[1] <- max(0.001, min(0.999, params_new[1])) # Delta (0, 1)
  params_new[4] <- max(0.01, params_new[4]) # sigma2_1 (> 0)
  params_new[5] <- max(0.01, params_new[5]) # sigma2_2 (> 0)

  return(params_new)
}

# Algoritmo Simulated Annealing
SA_GMM <- function(
  x,
  T_init = 10,
  alpha = 0.99,
  step_init = 0.5,
  max_iter = 5000
) {
  # Inicialização
  mu_mean <- mean(x)
  mu_sd <- sd(x)

  # Chutes iniciais aleatórios (delta, mu1, mu2, sigma2_1, sigma2_2)
  theta <- c(
    runif(1, 0.2, 0.8),
    rnorm(1, mu_mean - 0.5 * mu_sd, 0.5),
    rnorm(1, mu_mean + 0.5 * mu_sd, 0.5),
    runif(1, 0.5, 2),
    runif(1, 0.5, 2)
  )

  T <- T_init
  E_current <- loglik_gmm(x, theta)
  theta_best <- theta
  E_best <- E_current
  step_size <- step_init

  for (t in 1:max_iter) {
    theta_new <- perturbar(theta, step_size = step_size)
    E_new <- loglik_gmm(x, theta_new)

    delta_E <- E_new - E_current

    # Critério de Metropolis
    if (delta_E <= 0 || runif(1) < exp(-delta_E / T)) {
      theta <- theta_new
      E_current <- E_new

      if (E_new < E_best) {
        theta_best <- theta_new
        E_best <- E_new
      }
    }

    # Resfriamento
    T <- alpha * T
    # Redução do passo de perturbação
    step_size <- step_size * 0.999
  }

  return(theta_best)
}

# =======================================================
# 2. ESTUDO DE MONTE CARLO (M=200 RÉPLICAS)
# =======================================================

estudo_mc_sa <- function(n, mu_true, sigma2_true, p_true, m = 200) {
  # Inicializa vetores de armazenamento
  p_est <- mu1_est <- mu2_est <- sigma1_est <- sigma2_est <- numeric(m)

  for (i in 1:m) {
    # Para garantir que cada réplica use uma nova semente no SA
    set.seed(20252 + i)

    x <- rmixed_norm(n, mu_true, sigma2_true, p_true)
    est <- SA_GMM(x)

    p_final <- est[1]
    mu1_final <- est[2]
    mu2_final <- est[3]
    sigma1_final <- est[4]
    sigma2_final <- est[5]

    # Reordenar para garantir que mu1_final < mu2_final (igual à lógica do EM na Questão 6)
    if (mu1_final > mu2_final) {
      # Troca Mu
      temp_mu <- mu1_final
      mu1_final <- mu2_final
      mu2_final <- temp_mu

      # Troca Sigma^2
      temp_sigma <- sigma1_final
      sigma1_final <- sigma2_final
      sigma2_final <- temp_sigma

      # Troca o peso p
      p_final <- 1 - p_final
    }

    # Armazena os resultados reordenados
    p_est[i] <- p_final
    mu1_est[i] <- mu1_final
    mu2_est[i] <- mu2_final
    sigma1_est[i] <- sigma1_final
    sigma2_est[i] <- sigma2_final
  }

  return(data.frame(p_est, mu1_est, mu2_est, sigma1_est, sigma2_est))
}

# Função de Avaliação
avaliacao_sa <- function(av, p_true, mu_true, sigma2_true) {
  data.frame(
    p_verdadeiro = p_true,
    p_medio = mean(av$p_est),
    p_sd = sd(av$p_est),
    mu1_verdadeiro = mu_true[1],
    mu1_medio = mean(av$mu1_est),
    mu1_sd = sd(av$mu1_est),
    mu2_verdadeiro = mu_true[2],
    mu2_medio = mean(av$mu2_est),
    mu2_sd = sd(av$mu2_est),
    sigma1_verdadeiro = sigma2_true[1],
    sigma1_medio = mean(av$sigma1_est),
    sigma1_sd = sd(av$sigma1_est),
    sigma2_verdadeiro = sigma2_true[2],
    sigma2_medio = mean(av$sigma2_est),
    sigma2_sd = sd(av$sigma2_est)
  )
}

# =======================================================
# 3. EXECUÇÃO PARA OS 4 CENÁRIOS DO PDF
# =======================================================

# Cenários (iguais à Questão 6 do seu PDF)
av1 <- estudo_mc_sa(
  n = 300,
  mu_true = c(0.5, 7),
  sigma2_true = c(3, 1),
  p_true = 0.8
)
av2 <- estudo_mc_sa(
  n = 300,
  mu_true = c(0.7, 3),
  sigma2_true = c(6, 4),
  p_true = 0.5
)
av3 <- estudo_mc_sa(
  n = 4000,
  mu_true = c(0.6, 4),
  sigma2_true = c(1.5, 6),
  p_true = 0.4
)
av4 <- estudo_mc_sa(
  n = 4000,
  mu_true = c(0.3, 8),
  sigma2_true = c(2, 2),
  p_true = 0.7
)

# Combinação e formatação da tabela final
tab_resultados_sa <- rbind(
  avaliacao_sa(av1, 0.8, c(0.5, 7), c(3, 1)),
  avaliacao_sa(av2, 0.5, c(0.7, 3), c(6, 4)),
  avaliacao_sa(av3, 0.4, c(0.6, 4), c(1.5, 6)),
  avaliacao_sa(av4, 0.7, c(0.3, 8), c(2, 2))
)

tab_resultados_sa <- round(tab_resultados_sa, 3)

# Renomeação das colunas para melhor visualização (usando LaTeX para os símbolos)
# Você pode usar a função kable() do R para renderizar com cabeçalhos bonitos.
# Usando print para saída simples.
print(tab_resultados_sa)
