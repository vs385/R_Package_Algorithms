
#' algo_leverage function
#' @description
#'
#' The function implements algorithmic leveraging or linear regression using leverage score and uniform-based row subsampling.
#' It replicates the simulation setting of Fig 1 in Ma et Al's in https://onlinelibrary.wiley.com/doi/full/10.1002/wics.1324.
#'
#'
#'
#' @param y The response vector with n entries
#' @param X The covariates n by p matrix
#' @param r Percentage of subsamples to be drwan; default set to 0.1
#' @param sample_size Number of subsamples to consider for regression
#' @param method The type of subsampling: 'unif' for uniform-based, 'leverage' for leverage-score based
#'
#' @return The coefficient vector \code{beta}
#' @export
#'
#' @examples algo_leverage(rnorm(500), matrix(rnorm(1000),nrow=500), method='unif')
#'
algo_leverage <- function(y, X, r=0.1, sample_size = floor(r*length(y)), method='leverage') {
  n =length(y)

  unif_prob = rep(1, n)
  lev_prob = diag(tcrossprod(X, X))

  hii = diag(rep(1, sample_size))

  if (method=='unif') {
    sampled = sample(n, size=sample_size, prob = unif_prob)
  } else {
    sampled = sample(n, size=sample_size, prob=lev_prob)
    hii = diag(1/(lev_prob[sampled]))

    sampledX = X[sampled, ]
    sampledY= y[sampled]

    beta_hat = crossprod(solve(crossprod(sampledX, hii %*% sampledX)),
                         crossprod(sampledX, hii %*% as.matrix(sampledY)))

    return (beta_hat)
  }

}
