#' Solve single effect regression.
#'
#' @param X explanatory variable, n*p matrix.
#' @param Y response variable, n*1 matrix.
#' @param sigma2 variance of noise, scalar.
#' @param sigma0_2 variance of the coefficient, scalar.
#' @return The posterior parameters.
SER = function(X, Y, sigma2, sigma0_2){
  n = dim(X)[1]
  p = dim(X)[2]
  alpha = numeric(p)
  mu1 = numeric(p)
  sigma1_2 = numeric(p)

  # define the prior as simplest case
  Pi = rep(1/p, p)

  for(j in 1:p){
    re_j = BSLR(X[, j], Y, sigma2, sigma0_2)
    alpha[j] = Pi[j]*exp(re_j$logBF)
    mu1[j] = re_j$mu1
    sigma1_2[j] = re_j$sigma1_2
  }

  alpha = alpha/sum(alpha)

  return(list(alpha = alpha, mu1 = mu1, sigma1_2 = sigma1_2))
}
