#' @title \code{\link{hydra}} with additional stress minimization
#'
#' @description  Runs the \code{\link{hydra}} method and then performs a further optimization step by minimizing the stress of the embedding and optimizing hyperbolic curvature
#'
#' @details See \url{https://arxiv.org/abs/1903.08977} for more details.
#'
#' @inheritParams hydra
#' @param curvature.bias Modify curvature before stress minimization by multiplying with \code{curvature.bias}
#' @param curvature.freeze Freeze the curvature returned by \code{\link{hydra}}. If TRUE then no optimization of curvature is attempted in the second stage of the algorithm. If FALSE then curvature is optimized in the second stage
#' @param curvature.max Upper bound for the curvature. If NULL, a defulat bound is used
#' @param maxit Maximal number of iterations. This parameter is passed to the optimization routine \code{\link[stats]{optim}}
#' @param ... Additional parameters are passed to \code{\link[stats]{optim}}, which performs the underlying stress minimization
#'
#' @inherit hydra return
#'
#' @author Martin Keller-Ressel <martin.keller-ressel@tu-dresden.de>
#'
#' @examples
#' data(karate)
#' embedding <- hydraPlus(karate$distance)
#' plot(embedding,labels=karate$label,node.col=karate$group,graph.adj=karate$adjacency)
#'
#' @importFrom stats optim
#'
#' @export

hydraPlus <- function(D,dim=2,curvature=1,alpha=1.10,equi.adj=0.5,control=list(), curvature.bias = 1,curvature.freeze = TRUE,curvature.max=NULL,maxit=1000,...) {

  n <- nrow(D)
  # calculate initial embedding
  control$polar <- FALSE
  h.embed <- hydra(D,dim,curvature=curvature,alpha=alpha,equi.adj=equi.adj,control)
  x.init <- as.vector(poincare.to.hyper(h.embed$r,h.embed$directional)) # convert to hyperbolic coordinatees
  x0 <- jitter(x.init) #add jitter to avoid duplicate points
  if(!curvature.freeze) x0 <- c(x0,h.embed$curvature) # should curvature be kept constant?

  # set optimization bounds
  lower <- rep(-Inf,dim*n)
  upper <- rep(Inf,dim*n)
  if(!curvature.freeze) {
    if(is.null(curvature.max)) curvature.max <- (24/max(D))^2 # heuristic for upper curvature bound
    lower <- c(lower,0.0001)
    upper <- c(upper,curvature.max)
  }
  if(curvature.freeze)
    # optimization with frozen curvature
    opt <- optim(x0,stress.objective,stress.gradient,nrows=n,ncols=dim,dist=D,curvature=curvature.bias*h.embed$curvature,...,method="L-BFGS-B",lower=lower,upper=upper,control=list(trace=1,maxit=maxit))
  else
    # optimization with variable curvature
    opt <- optim(x0,stress.objective,stress.gradient,nrows=n,ncols=dim,dist=D,curvature=NULL,...,method="L-BFGS-B",lower=lower,upper=upper,control=list(trace=1,maxit=maxit))

  # convert back to poincare coordintates
  X <- matrix(opt$par[1:(n*dim)], nrow = n, ncol = dim)
  out <- hyper.to.poincare(X)

  if(!curvature.freeze) out$curvature <- opt$par[n*dim+1]
  else out$curvature <- h.embed$curvature

  out$curvature.max <- curvature.max
  out$convergence.code <- opt$convergence
  out$stress <- get.stress(out$r,out$directional, out$curvature,D)

  class(out) <- c("hydra","list")
  return(out)
}

stress.gradient <- function(x, nrows, ncols, dist, weight=FALSE,curvature=NULL) {
  # This function calculates the gradient for stress-minimzation
  # x is the vectorization of the coordinate matrix X, which has dimensions nrows x ncols.
  # The rows of X are the embedded points and the columns the reduced hyperbolic coordinates
  # if curvature is not NULL, then curvature is appended to x as last element

  if(!is.null(curvature)) {
    n = length(x)
    c.grad <- FALSE
  } else {
    n = length(x) - 1
    curvature = x[n+1]
    c.grad <- TRUE
  }
  x = matrix(x[1:n], nrow = nrows, ncol = ncols)
  X = x %*% t(x)
  u_tilde = matrix(0, nrows, 1)
  for (i in 1:nrows) {
    u_tilde[i] = sqrt(X[i,i]+1)
  }
  H = X - (u_tilde %*% t(u_tilde))
  H <- pmin(H,-(1 + .Machine$double.eps))
  D = 1/sqrt(curvature)*acosh(-H)
  diag(D) = 0
  A = (D - dist) * (1 / sqrt(curvature*(H^2 - 1)))
  if(weight) A <- A / dist
  diag(A) <- 0.0
  B = (1 / u_tilde) %*% t(u_tilde)
  G = 2*(matrix(rowSums(A*B),nrow=nrows,ncol=ncols) * x - A %*% x)
  if(weight) {
    D.diff <- (D - dist) * D/dist
    diag(D.diff) <- 0.0
    g = -0.5 * sum(curvature^(-1) * (D.diff)) ## korrekt
  } else {
    g = -0.5 * sum(curvature^(-1) * ((D - dist) * D)) ## korrekt
  }
  z = matrix(G, ncol=1)
  if(c.grad)
    z <-  c(z, g)
  return(z)
}

stress.objective <- function(x, nrows, ncols, dist,weight=FALSE,curvature=NULL) {
  # This function evaluated the objective function for stress minimization
  if(!is.null(curvature))
    n = length(x)
  else {
    n = length(x) - 1
    curvature = x[n+1]
  }
  x = matrix(x[1:n], nrow = nrows, ncol = ncols)
  X = x %*% t(x)
  u_tilde = matrix(0, nrows, 1)
  for (i in 1:nrows) {
    u_tilde[i] = sqrt(X[i,i]+1)
  }
  H = X - (u_tilde %*% t(u_tilde))
  D = 1/sqrt(curvature)*acosh(pmax(-H,1))
  diag(D) = 0
  if(weight) {
    D.diff <- (D - dist)^2 / dist
    diag(D.diff) <- 0.0
    y = 0.5 * sum(D.diff)
  } else {
    y = 0.5 * sum((D - dist)^2)
  }
  return(y)
}





