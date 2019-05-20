
#' solve_ols function
#'
#' @param A the (square) matrix of coefficients of the linear system
#' @param b vector of values for each of the linear equation in the system
#' @param n_cores specifies the number of cores for parallel computing framework
#' @param method specifies the method to solve the linear system, "GS" for Gauss-Seidel, "Jacobi" for Jacobi
#' @param n_iter states the number of iterations to convergence/stop updates. Default is 100
#'
#' @return The solution vector \code{x} for linear system \code{Ax = b}
#' @export
#'
#' @examples solve_ols(diag(rep(1, 3)), rep(1,3), n_cores = 2, method='GS', n_iter = 100)

solve_ols <- function(A, b, n_cores=2, method, n_iter=100) {

  gs_seq <- function(G_S, DL, x, n_iter, b) {
    x_init = x
    b_init = DL %*% b
    for (i in 1:n_iter) {
      x_new = G_S %*% x_init + b_init
      x_init = x_new
    }
  }


  jacobi_seq <- function(J_C, D, x, n_iter, b) {
    x_init = x
    b_init = D %*% b
    for (i in 1:n_iter) {
      x_new = J_C %*% x_init + b_init
      x_init = x_new
    }
  }

  jacobi_par <- function(A, x, n_iter, b, ncores) {
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", ncores))
    x_init = x
    for (i in 1:n_iter) {
      update_x <- function(j) {
        x_new = (b[j] - crossprod(A[j, -j], x_init[-j]))/A[j, j]
      }
      x_j = parallel::parLapply(cl=cl, X = 1:length(x_init), fun= update_x)
      x_init = x_j
    }
    parallel::stopCluster(cl)
    rm(cl)
    x_j
  }

  n = length(b)

  v = c()
  for (i in 1:50) {
    v[2*i] = 0
    v[2*i-1] = 1
  }

  x0 = array(data=0, dim=n)
  D = diag(diag(A))
  U = matrix(data=0, nrow=n, ncol=n)
  L = matrix(data=0, nrow=n, ncol=n)
  L[lower.tri(L)] = A[lower.tri(A)]
  U[upper.tri(U)] = A[upper.tri(A)]

  if (method=='GS') {
    DL = solve(D+L)
    G_S = -DL %*% U
    X = gs_seq(G_S, DL, x0, n_iter, b)
  } else if (method=='Jacobi') {
    if (n_cores > 1) {
      X = jacobi_par(A, x0, n_iter, b, n_cores)
    } else {
      D = solve(D)
      J = -D %*% (L+U)
      X = jacobi_seq(J, D, x0, n_iter, b)
    }
  }
}


