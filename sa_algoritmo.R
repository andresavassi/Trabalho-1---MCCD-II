# 1. Definição da Matriz de Distância (Dados do Problema)
dist_matrix <- matrix(
  c(
    0,
    10,
    15,
    20,
    25,
    10,
    0,
    35,
    25,
    12,
    15,
    35,
    0,
    30,
    18,
    20,
    25,
    30,
    0,
    8,
    25,
    12,
    18,
    8,
    0
  ),
  nrow = 5,
  byrow = TRUE
)

# 2. Função Objetivo (f) - Calcula o custo da rota (Distância Total)
f <- function(S) {
  # Adiciona a volta para a cidade inicial (fechando o ciclo)
  S_ciclo <- c(S, S[1])
  total_distancia <- 0

  # Soma a distância entre pares de cidades consecutivas na rota
  for (i in 1:(length(S_ciclo) - 1)) {
    cidade_atual <- S_ciclo[i]
    cidade_proxima <- S_ciclo[i + 1]
    total_distancia <- total_distancia +
      dist_matrix[cidade_atual, cidade_proxima]
  }
  return(total_distancia)
}

# 3. Função de Geração de Vizinho (gerar_vizinho) - Muda a rota trocando 2 cidades
gerar_vizinho <- function(S) {
  n_cidades <- length(S)

  # Escolhe dois índices aleatórios diferentes
  indices <- sample(1:n_cidades, 2, replace = FALSE)
  i <- indices[1]
  j <- indices[2]

  # Troca as cidades nas posições i e j
  S_n <- S
  S_n[i] <- S[j]
  S_n[j] <- S[i]

  return(S_n)
}


# O código 'simulated_annealing' completo:
simulated_annealing <- function(f, S0, T0, L, alpha, T_min) {
  T <- T0
  S <- S0
  S_star <- S0
  f_S <- f(S)
  f_S_star <- f_S

  while (T > T_min) {
    for (n in 1:L) {
      S_n <- gerar_vizinho(S)
      f_S_n <- f(S_n)

      Delta <- f_S_n - f_S

      if (Delta <= 0) {
        S <- S_n
        f_S <- f_S_n
      } else {
        p <- exp(-Delta / T)
        if (runif(1) < p) {
          S <- S_n
          f_S <- f_S_n
        }
      }
    }

    if (f_S < f_S_star) {
      S_star <- S
      f_S_star <- f_S
    }

    T <- T * alpha
  }
  return(list(best_solution = S_star, best_cost = f_S_star))
}

# --- Parâmetros de SA ---
S0 <- 1:5 # Solução inicial (cidades 1, 2, 3, 4, 5)
T0 <- 100 # Temperatura inicial
L <- 100 # Iterações por nível de temperatura
alpha <- 0.99 # Fator de resfriamento
T_min <- 0.001 # Temperatura mínima (critério de parada)

# --- Execução ---
set.seed(42) # Para resultados reproduzíveis
resultado <- simulated_annealing(f, S0, T0, L, alpha, T_min)
print(paste("Custo inicial:", f(S0)))
print("--- Resultado do SA ---")
print(paste("Melhor Custo:", resultado$best_cost))
print(paste(
  "Melhor Rota (Cidades):",
  paste(resultado$best_solution, collapse = " -> ")
))
