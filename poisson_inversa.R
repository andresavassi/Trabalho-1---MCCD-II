set.seed(20252)
library(dplyr)
library(kableExtra)
## Implementação do algoritmo para gerar um valor da Poisson
rpois_aux <- function(lambda) {
  u <- runif(1, 0, 1)
  i <- 0
  pr <- exp(-lambda)
  Fx = pr
  while (u >= Fx) {
    pr <- pr * lambda / (i + 1)
    Fx <- Fx + pr
    i <- i + 1
  }
  return(i)
}

## Replicando para n valores
rpoisson <- function(n, lambda) {
  replicate(n, expr = rpois_aux(lambda), simplify = TRUE)
}
## Gerando uma amostra
lambda <- 5
amostra <- rpoisson(10000, lambda)
## Frequencias relativas na amostra gerada e prob teorica
tab <- rbind(
  round(prop.table(table(amostra)), 3),
  round(dpois(0:max(amostra), lambda), 3)
)
colnames(tab) <- 0:max(amostra)
rownames(tab) <- c("freq", "prob")


## imprimindo a tabela formatada com as 10 primeras colunas
tab <- tab[, 1:10] %>%
  kbl(
    align = c(rep("c", 5)),
    caption = "Comparação de frequências relativas com probabilidades teóricas"
  )
tab %>%
  kable_styling(
    latex_options = "hold_position",
    bootstrap_options = c("bordered"),
    font_size = 11
  )


## Graficamente
freq_rel <- prop.table(table(amostra))
par(mar = c(4, 4, 1, 1))
plot(freq_rel, main = "", xlab = "x", ylab = "Frequência relativa (prob)")
points(0:max(amostra), dpois(0:max(amostra), lambda), col = 2, cex = 1.5)
