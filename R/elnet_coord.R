#' Elastic Net Regression
#' Fits elastic net to data using coordinate descent algorithm using simulation setting described in HW1
#' @usage elnet_coord()
#'
#' @return The function returns the matrix of regression coefficients for each \eqn{\lambda} given as column vectors, the simulated design matrix \eqn{X}, and the simulated response matrix \eqn{Y}
#' @export
#'

elnet_coord <- function() {
  library("MASS")

  # Generate the data
  # Generating the siumlated dataset using mvrnorm from library MASS

  p=20
  n=20
  alpha=0

  beta = c(2, 0, -2, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

  sigma = matrix (0, p, p)

  lambda_grid=seq(from=3,to=0,length.out = 100)

  sigma[1, 2] = 0.8
  sigma[2, 1] = 0.8
  sigma[5, 6] = 0.8
  sigma[6, 5] = 0.8
  for (i in 1:20) {
    sigma[i, i] = 1
  }
  mu = rep(0, p)
  X = mvrnorm(n=n, mu=mu, Sigma=sigma)

  # Normalization step
  for (j in 1:p) {
    X[, j] = sqrt(n)*X[, j]/sqrt(sum(X[, j]^2))
  }

  E = as.matrix(rnorm(n, 0, 1))
  y = X%*%beta + E

  # Coordinate descent algorithm
  soft_threshold <- function(mu, lmbda) {
    if (lmbda >= abs(mu)) {
      sol = 0
    } else if (mu > 0) {
      sol = mu - lmbda
    } else {
      sol = mu + lmbda
    }
    return (sol)
  }

  coord_des_step <- function(X, y, beta, j, lmbda) {
    residual = (1/n) * t(as.matrix(X[, j])) %*% (y - X[, -j] %*% beta[-j])
    beta[j] = soft_threshold(res, lmbda*alpha) / (1 + lmbda - lmbda*alpha)
    return (beta)

  }

  sim_func <- function(X, y, lmbda, alpha) {
    dif = 1
    beta_0 = matrix(nrow = p, ncol = 100)
    while (dif > 1e-4) {
      beta_0_old = beta_0
      for (i in 1: p) {
        beta_0 = coord_des_step(X, y, beta_0, i, lmbda)
      }
      dif = sqrt(sum(beta_0 - beta_0_old)^2) # Used RMSE to check convergence
    }
    return (beta_0)
  }

  beta_c0 = matrix(nrow = p, ncol = 100)
  beta_c0_l1 = matrix(nrow = 1, ncol = 100)
  beta_c0_l0 = matrix(nrow = 1, ncol = 100)

  for (j in 1:100) {
    beta_c0[, j] = sim_func(X=X, y=y, lmbda = lambda_grid[j], alpha= alpha)
    beta_c0_l1[, j] = sum(abs(beta_c0[, j]))
    beta_c0_l0[, j] = sum(beta_c0[, j] != 0)
  }

  matplot(x = t(beta_c0_l1), y = t(beta_c0), type = "l", lty = 1,
          xlab = "L1 Norm", ylab = "coefficients", col = 1:20)
  return (list(beta = beta_c0, X = X, Y = y))
}



