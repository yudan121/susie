#' Compute posterior inclusion probabilities.
#'
#' @param alpha_mat matrix containing posterior means of gammas.
#' @return Posterior inclusion probabilities.
PIP = function(alpha_mat){
  pip = 1 - apply(1 - alpha_mat, 1, prod)
  return(pip)
}
