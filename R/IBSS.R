#' Iterative Bayesian Stepwise Selection.
#'
#' @param X explanatory variable, n*p matrix.
#' @param Y response variable, n*1 matrix.
#' @param sigma2 variance of the noise, scalar.
#' @param sigma0_2 variance of the coefficient, scalar.
#' @param L number of effect variables
#' @param iter.max maximum number of iterations
#' @return The estimated posterior parameters.
IBSS = function(X, Y, sigma2, sigma0_2, L, iter.max = 500){
  n = dim(X)[1]
  p = dim(X)[2]
  eps = 1e-3

  nb_bar = matrix(0, ncol = L, nrow = p)
  nalpha_mat = matrix(0, ncol = L, nrow = p)
  nmu1_mat = matrix(0, ncol = L, nrow = p)
  nsigma1_2_mat = matrix(0, ncol = L, nrow = p)

  b_bar = nb_bar
  alpha_mat = nalpha_mat
  mu1_mat = nmu1_mat
  sigma1_2_mat = nsigma1_2_mat

  for(iter in 1:iter.max){
    for(l in 1:L){
      # expected residuals without l-th single effect
      rl_bar = Y - rowSums(X%*%nb_bar) + X%*%nb_bar[, l]
      # fit SER to expected residuals
      result_l = SER(X, rl_bar, sigma2, sigma0_2)

      # get new parameters
      nalpha_mat[, l] = result_l$alpha
      nmu1_mat[, l] = result_l$mu1
      nsigma1_2_mat[, l] = result_l$sigma1_2
      nb_bar[, l] = nalpha_mat[, l]*nmu1_mat[, l]
    }

    # convergence criterion
    if(sum((alpha_mat - nalpha_mat)^2) + sum((mu1_mat - nmu1_mat)^2) + sum((sigma1_2_mat - nsigma1_2_mat)^2) < eps){
      break
    }

    b_bar = nb_bar
    alpha_mat = nalpha_mat
    mu1_mat = nmu1_mat
    sigma1_2_mat = nsigma1_2_mat
  }

  return(list(alpha_mat = alpha_mat, mu1_mat = mu1_mat, sigma1_2_mat = sigma1_2_mat))
}
