# Geração da Binomial Negativa via Amostragem por Importância (SIR)

# --- Parâmetros ---
r_param <- 4      # número de sucessos
p_param <- 0.5    # Probabilidade de sucesso
n_amostra <- 1000 # Tamanho da amostra a ser gerada
m_candidatos <- 20000 

# --- Função principal para gerar n valores ---
gerar_bn_sir <- function(n, r, p, m) {
  
  # Define a distribuição proposta (Geométrica) com a mesma média da alvo (BN)
  media_bn <- r * (1 - p) / p
  p_g <- 1 / (1 + media_bn)
  
  # Passo 1: Gerar 'm' candidatos da distribuição proposta
  candidatos <- rgeom(m, prob = p_g)
  
  # Passo 2: Calcular os pesos de importância 
  # f(y) = probabilidade do alvo (BN)
  # g(y) = probabilidade da proposta (Geométrica)
  pesos_alvo <- dnbinom(candidatos, size = r, prob = p)
  pesos_proposta <- dgeom(candidatos, prob = p_g)
  
  pesos_importancia <- pesos_alvo / pesos_proposta
  
  # Passo 3: Reamostrar dos candidatos com base nos pesos 
  amostra_final <- sample(
    x = candidatos, 
    size = n, 
    replace = TRUE, 
    prob = pesos_importancia
  )
  
  return(amostra_final)
}

# --- Execução e visualização do resultado ---
set.seed(789)
amostra_sir <- gerar_bn_sir(n = n_amostra, r = r_param, p = p_param, m = m_candidatos)

print("Amostra gerada pelo Método SIR (primeiros 10 valores):")
print(head(amostra_sir, 10))

# Histograma para validar
hist(amostra_sir,
     breaks = seq(-0.5, max(amostra_sir) + 0.5, 1),
     main = "Amostragem por Importância (SIR)",
     xlab = "Valor Gerado", ylab = "Densidade", freq = FALSE)
# Adiciona a PMF teórica para comparação
points(0:max(amostra_sir), dnbinom(0:max(amostra_sir), size = r_param, prob = p_param),
       col = "darkgreen", pch = 19)
legend("topright", legend = "Prob. Teórica", col = "darkgreen", pch = 19)