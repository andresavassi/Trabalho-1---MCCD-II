# Geração da Binomial Negativa via Aceitação e Rejeição

# --- Parâmetros ---
r_param <- 4      # número de sucessos
p_param <- 0.5    # Probabilidade de sucesso
n_amostra <- 1000 

# --- Função principal para gerar n valores ---
gerar_bn_rejeicao <- function(n, r, p) {
  
  # Define a distribuição proposta (Geométrica) com a mesma média da alvo (BN)
  media_bn <- r * (1 - p) / p
  p_g <- 1 / (1 + media_bn) # Parâmetro da geométrica
  
  # Calcula a constante c, tal que c >= max( p(k) / q(k) ) 
  # (calculada numericamente para simplicidade)
  k_vals <- 0:100 
  razao_pq <- dnbinom(k_vals, size = r, prob = p) / dgeom(k_vals, prob = p_g)
  c <- max(razao_pq, na.rm = TRUE) * 1.01 # fator de segurança
  
  amostra <- numeric(n)
  contador_aceitos <- 0
  
  while (contador_aceitos < n) {
    # Passo 1: Gerar candidato 'w' da proposta Geométrica
    w <- rgeom(1, prob = p_g)
    
    # Passo 2: Gerar u ~ Uniforme(0,1) 
    u <- runif(1)
    
    # Probabilidades no ponto candidato 'w'
    p_w <- dnbinom(w, size = r, prob = p)
    q_w <- dgeom(w, prob = p_g)
    
    # Passo 3: Critério de aceitação 
    if (u < p_w / (c * q_w)) {
      contador_aceitos <- contador_aceitos + 1
      amostra[contador_aceitos] <- w
    }
  }
  return(amostra)
}

# --- Execução e visualização do resultado ---
set.seed(456)
amostra_rejeicao <- gerar_bn_rejeicao(n = n_amostra, r = r_param, p = p_param)

print("Amostra gerada pelo Método da Rejeição (primeiros 10 valores):")
print(head(amostra_rejeicao, 10))

# Histograma para validar
dir.create("figuras", showWarnings = FALSE)
png("figuras/graf_bn_ar.png", width = 800, height = 600, res = 120)
hist(amostra_rejeicao,
     col = "darkblue",
     breaks = seq(-0.5, max(amostra_rejeicao) + 0.5, 1),
     main = "Método da Aceitação e Rejeição",
     xlab = "Valor Gerado", ylab = "Densidade", freq = FALSE)
# Adiciona a PMF teórica para comparação
points(0:max(amostra_rejeicao), dnbinom(0:max(amostra_rejeicao), size = r_param, prob = p_param),
       col = "red", pch = 19)
legend("topright", legend = "Prob. Teórica", col = "red", pch = 19)

dev.off()