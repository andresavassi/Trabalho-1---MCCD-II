# Geração da Binomial Negativa via Transformação Inversa

# --- Parâmetros ---
r_param <- 4      # número de sucessos
p_param <- 0.5    # Probabilidade de sucesso
n_amostra <- 1000 

# --- Função auxiliar para gerar um único valor ---
gerar_bn_inversa_aux <- function(r, p) {
  # Gera um número da Uniforme(0,1)
  u <- runif(1) 
  
  i <- 0
  pr <- p^r      # P(X=0)
  Fx <- pr       # FDA no ponto 0, F(0)
  
  # Loop para encontrar o valor k tal que F(k-1) <= u < F(k) 
  while (u >= Fx) {
    # Relação de recorrência para a Binomial Negativa
    pr <- pr * (i + r) / (i + 1) * (1 - p)
    Fx <- Fx + pr
    i <- i + 1
  }
  return(i)
}

# --- Função principal para gerar n valores ---
gerar_bn_inversa <- function(n, r, p) {
  # 'replicate' executa a função auxiliar 'n' vezes para criar a amostra
  return(replicate(n, gerar_bn_inversa_aux(r, p)))
}

# --- Execução e visualização do resultado ---
set.seed(123)
amostra_inversa <- gerar_bn_inversa(n = n_amostra, r = r_param, p = p_param)

print("Amostra gerada pelo Método da Inversa (primeiros 10 valores):")
print(head(amostra_inversa, 10))

# Histograma para validar
hist(amostra_inversa,
     breaks = seq(-0.5, max(amostra_inversa) + 0.5, 1),
     main = "Método da Transformação Inversa",
     xlab = "Valor Gerado", ylab = "Densidade", freq = FALSE)
# Adiciona a PMF teórica para comparação
points(0:max(amostra_inversa), dnbinom(0:max(amostra_inversa), size = r_param, prob = p_param),
       col = "red", pch = 19)
legend("topright", legend = "Prob. Teórica", col = "red", pch = 19)