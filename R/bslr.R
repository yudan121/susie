#' Solve Bayesian simple linear regression.
#'
#' @param x explanatory variable, n*1 vector.
#' @param y response variable, n*1 vector.
#' @param sigma2 variance of noise, scalar.
#' @param sigma0_2 variance of the coefficient, scalar.
#' @return The posterior mean and variance, and bayes factor.
BSLR = function(x, y, sigma2, sigma0_2){
  xTx = sum(x^2)
  xTy = sum(x*y)
  b_hat = xTy/xTx
  s2 = sigma2/xTx
  z = b_hat/sqrt(s2)

  # b|y ~ N(mu1, sigma1_2)
  sigma1_2 = 1/(1/s2 + 1/sigma0_2)
  mu1 = sigma1_2*b_hat/s2

  # Bayes Factor
  logBF = log(sqrt(s2/(sigma0_2 + s2))) + (z^2/2*sigma0_2/(sigma0_2 + s2))

  return(list(mu1 = mu1, sigma1_2 = sigma1_2, logBF = logBF))
}
