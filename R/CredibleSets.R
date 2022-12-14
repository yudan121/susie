#' Compute credible sets.
#'
#' @param alpha_mat matrix containing posterior means of gammas.
#' @param rho minimal probability that credible sets should reach.
#' @return The required credible sets.
CredibleSet = function(alpha_mat, rho){
  L = dim(alpha_mat)[2]
  cs = list()
  for(l in 1:L){
    k0 = sum(cumsum(sort(alpha_mat[, l], decreasing = TRUE)) < rho) + 1
    cs = append(cs, list(order(alpha_mat[, l], decreasing = TRUE)[1:k0]))
  }
  return(cs)
}

#' Get purified credible sets
#'
#' @param X explanatory variable, n*p matrix.
#' @param cs credible sets
#' @param purity minimal purity as defined in the paper.
#' @param abs logical value, indicating whether the criterion should be abs of correlation.
#' @return Purified credible sets.
purify_cs = function(X, cs, purity = 0.5, abs = FALSE){
  L = length(cs)
  pure_cs = list()
  for(l in 1:L){
    if(length(cs[[l]]) > 1){
      if(abs == FALSE){
        if(min(cor(X[, cs[[l]]])) >= purity){
          pure_cs = c(pure_cs, list(cs[[l]]))
        }
      }else{
        if(min(abs(cor(X[, cs[[l]]]))) >= purity){
          pure_cs = c(pure_cs, list(cs[[l]]))
        }
      }
    }else{
      pure_cs = c(pure_cs, list(cs[[l]]))
    }
  }
  return(pure_cs)
}
