# Importado de bellreg
dbell <- function(x, theta, log = FALSE) {
  Bx <- c()
  for (i in 1:length(x)) {
    Bx[i] <- numbers::bell(x[i])
  }
  lf <- x * log(theta) - exp(theta) + 1 + log(Bx) - lgamma(x + 1)
  if (log == TRUE) {
    return(lf)
  } else {
    return(exp(lf))
  }
}


# Por metodo da inversa

rbell_aux <- function(theta) {
  sapply(theta, function(t) {
    u <- runif(1, 0, 1)
    i <- 0
    pr <- dbell(0, t)
    Fx <- pr
    while (u >= Fx) {
      i <- i + 1
      pr <- dbell(i, t)
      Fx <- Fx + pr
    }
    return(i)
  })
}

rbell <- function(n, theta) {
  if (length(theta) == 1) {
    return(replicate(n, expr = rbell_aux(theta), simplify = TRUE))
  } else {
    return(rbell_aux(theta))
  }
}

#Conferindo resultados
#theta <- 1
#amostra <- rbell(100000, theta)

#tab <- rbind(
#  round(prop.table(table(amostra)), 3),
#  round(dbell(0:max(amostra), theta), 3)
#)
#colnames(tab) <- 0:max(amostra)
#rownames(tab) <- c("freq", "prob")

#tab <- tab[,1:10] %>%
#kbl(align = c(rep("c", 5)),
#caption = "Comparação de frequências relativas com probabilidades teóricas")
#tab %>%
#kable_styling(latex_options = "hold_position",
#bootstrap_options = c("bordered"), font_size = 8)



set.seed(20252)
theta <- 1
n <- 100000
amostra <- rbell(n, theta)

dir.create("figuras", showWarnings = FALSE)
png("figuras/graf_bell.png", width = 800, height = 600, res = 120)

# --- Histograma das frequências relativas ---
hist(amostra,
     breaks = seq(-0.5, max(amostra) + 0.5, 1),
     freq = FALSE,
     col = "darkblue",
     border = "white",
     main = sprintf("Distribuição Bell: Teórico vs Estimado (θ = %.2f, n = %d)", theta, n),
     xlab = "x",
     ylab = "Probabilidade / Frequência relativa"
)

# --- Probabilidades teóricas ---
x_vals <- 0:max(amostra)
prob_teo <- dbell(x_vals, theta)

# --- Sobreposição das probabilidades teóricas ---
points(x_vals, prob_teo, col = "red", pch = 19)
lines(x_vals, prob_teo, col = "red", lwd = 2)

# --- Legenda ---
legend("topright",
       legend = c("Frequência empírica", "Probabilidade teórica"),
       col = c("lightblue", "red"),
       pch = c(15, 19),
       pt.cex = c(1.5, 1.2),
       bty = "n")

dev.off()
