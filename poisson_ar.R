set.seed(20252)

# Função que gera um valor da Poisson(λ) via Aceitação e Rejeição
rpois_ar <- function(lambda) {
  
  # Candidato: distribuição geométrica (P(X = k) = p*(1-p)^k), k = 0,1,...
  # Escolhemos p de modo que E[X] = (1-p)/p ≈ λ  =>  p = 1/(1+λ)
  p <- 1 / (1 + lambda)
  
  # Função densidade da Poisson (não normalizada)
  f <- function(k) exp(-lambda) * lambda^k / factorial(k)
  
  # Função densidade da Geométrica
  g <- function(k) p * (1 - p)^k
  
  # Estimativa grosseira de M
  # M precisa ser >= max_k f(k)/g(k)
  # Calculamos para os primeiros valores até ter f/g pequeno
  ks <- 0: (lambda*10 + 10)
  M <- max(f(ks) / g(ks))
  
  repeat {
    # Gera X da Geométrica
    X <- rgeom(1, prob = p)
    # Gera U uniforme
    U <- runif(1)
    # Teste de aceitação
    if (U <= f(X) / (M * g(X))) {
      return(X)
    }
  }
}


# Gerar 10.000 valores para ver se o histograma parece correto
amostras <- replicate(10000, rpois_ar(4))

dir.create("figuras", showWarnings = FALSE)
png("figuras/graf_poisson_ar.png", width = 800, height = 600, res = 120)
hist(amostras, breaks = seq(-0.5, max(amostras)+0.5, 1),
     freq = FALSE, col = "darkblue", border = "white",
     main = "Distribuição Poisson(4) via Aceitação e Rejeição",
     xlab = "x")

# Comparar com a densidade teórica
x <- 0:max(amostras)
lines(x, dpois(x, 4), col = "red", lwd = 2)
legend("topright", legend = c("Simulado", "Teórico"), 
       fill = c("darkblue", NA), border = c("white", NA),
       lty = c(NA, 1), col = c(NA, "red"))


dev.off()
