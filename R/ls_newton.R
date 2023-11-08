#' Least Squares Estimation using Newton's Method
#'
#' @param y response variable
#' @param X design matrix with intercept included
#' @param tolerance tolerance level for convergence
#' @param max_iterations maximum number of iterations
#' @return an object of class 'lsnewton' with the estimated coefficients
#' @export
#' @examples
#' y <- c(1, 2, 3)
#' X <- matrix(c(1, 1, 1, 1, 2, 3), ncol = 2)
#' result <- ls_newton(y, X)
#' print(result)
ls_newton <- function(y, X, tolerance = 1e-06, max_iterations = 100) {
  # assuming X contains the intercept
  p <- ncol(X)
  beta <- rep(0, p)
  for (iteration in 1:max_iterations) {
    r <-  y - X %*% beta
    step <- solve(t(X) %*% X) %*% t(X) %*% r
    beta <- beta + step
    if (max(abs(step)) < tolerance) {
      break
    }
  }
  class(beta) <- "lsnewton"
  return(beta)
}

#' Print method for lsnewton objects
#'
#' @param x an object of class 'lsnewton'
#' @param ... additional arguments affecting the summary produced
#' @export
print.lsnewton <- function(x, ...) {
  cat("Coefficients:\n")
  print(unclass(x))
}
