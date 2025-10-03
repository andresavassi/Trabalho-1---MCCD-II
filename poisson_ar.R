poisson_ar_1 <- function(lambda, kmax = 30) {
  f_pois <- function(k, lambda) {
    return(exp(-lambda) * lambda^k / factorial(k))
  }	
  fmax <- f_pois(floor(lambda), lambda) 
  repeat {	
    y <- sample(0:kmax, 1)
    fy <- f_pois(y, lambda)		
    gy <- 1 / (kmax + 1)				
    u <- runif(1)
    if (u < fy / (fmax * gy)) {
      return(y)
    }
  }
}
poisson_ar <- function(lambda, n = 1000, kmax = 30) {
  replicate(n, poisson_ar_1(lambda, kmax))
}