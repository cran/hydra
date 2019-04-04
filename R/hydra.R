#' @title Calculate hyperbolic embedding of distance data
#'
#' @description Implements the HYDRA (hyperbolic distance recovery and approximation) method for embedding high-dimensional data points (represented by their distance matrix \code{D}) into low-dimensional hyperbolic space.
#'
#' @details See \url{https://arxiv.org/abs/1903.08977} for more details.
#'
#' @param D a square symmetric matrix of distances (or dissimiliarities) to be embdedded, can also be a \code{\link[stats]{dist}} object
#' @param dim embedding dimension
#' @param curvature embedding curvature; if this argument is NULL, hydra tries to find the optimal curvature
#' @param alpha real number greater one; adjusts the hyperbolic curvature. Values larger than one yield a more distorted embedding where points are pushed
#'          to the outer boundary (i.e. the ideal points) of hyperblic space. The interaction between \code{curvature} and \code{alpha} is non-linear.
#' @param equi.adj equi-angular adjustment; must be a real number between zero and one; only used if \code{dim} is 2. Value 0 means no ajustment, 1 adjusts
#'              embedded data points such that their angular coordinates in the Poincare disc are uniformly distributed. Other values interpolate between the two extremes. Setting the parameter to non-zero values can make the embedding result look more harmoniuous in plots.
#' @param control a list which may contain the following boolean flags:
#'            \itemize{
#'            \item polar - return polar coordinates in dimension 2 (default: TRUE if \code{dim} is 2. This flag is ignored in higher dimension)
#'            \item isotropic.adj - perform isotropic adjustment, ignoring Eigenvalues (default: TRUE if \code{dim} is 2, FALSE else)
#'            \item return.lorentz - return raw Lorentz coordinates (before projection to hyperbolic space) (default: FALSE)
#'            \item return.stress - return embedding stress (default: TRUE)
#'            \item return.dist - return hyperbolic distance matrix of embedded points (default: FALSE)
#'            \item use.eigs  - use \code{\link[RSpectra]{eigs}} function from \pkg{RSpectra} and \code{\link[Matrix]{norm}} function from \pkg{Matrix} to speed up computation (default: FALSE)
#'            }
#'
#' @return A `hydra' object, which is a list with all or some of the following components:
#'           \describe{
#'           \item{r}{a vector containing the radial coordinates of the embedded points}
#'           \item{directional}{a matrix with \code{dim} columns containing as rows the directional coordinates of the embedded points}
#'           \item{theta}{a vector containing the angular coordinates of the embedded points (only returned if \code{dim} is 2 and \code{polar} flag  is TRUE)}
#'           \item{curvature}{the curvature used for the returned embedding}
#'           \item{dim}{the dimension used for the returned embedding}
#'           \item{stress}{the stress (i.e. the mean-square difference) between distances supplied in \code{D} and the hyperbolic distance matrix of the returned embedding}
#'           \item{dist}{the hyperbolic distance matrix of the returned embedding (only returned if flag \code{return.dist} is true. Computation may be time- and memory-intensive.)}
#'           \item{x0}{a vector containing the 'time-like' coordinate of the raw Lorentz embedding (only returned if flag \code{return.lorentz} is true)}
#'           \item{X}{a matrix with \code{dim} columns containing as rows the 'space-like' coordinate of the raw Lorentz embedding (only returned if flag \code{return.lorentz} is true)}
#'           }
#'
#' @author Martin Keller-Ressel <martin.keller-ressel@tu-dresden.de>
#'
#' @examples
#' data(karate)
#' embedding <- hydra(karate$distance)
#' plot(embedding,labels=karate$label,lab.col=karate$group,graph.adj=karate$adjacency)
#'
#' ## Compare with Multidimensional scaling (MDS):
#' mds <- cmdscale(karate$distance) # Compute Euclidean embedding with MDS
#' mds.stress <- sqrt(sum((as.matrix(dist(mds)) - karate$distance)^2)) # Calculate embedding stress
#' c(embedding$stress,mds.stress) # Compare hyperbolic with Euclidean stress
#'
#' @importFrom stats optimize
#'
#' @export

hydra <- function(D,dim=2,curvature=1,alpha=1.1,equi.adj=0.5,control=list()) {
## This is the wrapper function for different variants of the hydra method

  if(!is.null(curvature)) { # Has a (non-null) curvature parameter been supplied?
    return(hydra.fixed.curvature(D,dim,curvature,alpha,equi.adj,control))  # run hydra with fixed curvature
  } else {
    # setup control parameter for curvature optimization
    control.inner <- control
    control.inner$return.stress <- TRUE
    control.inner$return.lorentz <- FALSE

    k.bounds <- c(.Machine$double.eps,(8/max(D))^2) ## a priori bounds for curvature
    k.optimal <- optimize(function(k) hydra.fixed.curvature(D,dim,k,alpha,equi.adj,control.inner)$stress, interval=k.bounds) # find curvature with lowest stress
    k1.objective <- hydra.fixed.curvature(D,dim,1,alpha,equi.adj,control.inner)$stress # stress for unit curvature
    if(k1.objective < k.optimal$objective) k.optimal$minimum <- 1 # make sure that returned result is never worse than unit curvature
    out <- hydra.fixed.curvature(D,dim,k.optimal$minimum,alpha,equi.adj,control) # calculate embedding using the optimal curvature parameter
    out$curvature <- k.optimal$minimum
    return(out)
  }
}

hydra.fixed.curvature <- function(D, dim=2, curvature = 1, alpha=1, equi.adj=0, control=list()) {
## Hydra method with fixed curvature. This is the workhorse of the hydra method

## sanitize/check input
  if(any(diag(D) != 0)) { # non-zero diagonal elements are set to zero
    diag(D) <- 0
    warning("Diagonal of input matrix D has been set to zero")
  }
  if(!isSymmetric(D)) warning("Input matrix D is not symmetric. Lower triangle part is used.")

  ## replace null parameters with default values
  if(is.null(control$return.lorentz)) control$return.lorentz <- FALSE
  if(is.null(control$return.dist)) control$return.dist <- FALSE
  if(is.null(control$return.stress)) control$return.stress <- TRUE
  if(is.null(control$use.eigs)) control$use.eigs <- FALSE

  if(dim == 2) {
    # set default values in dimension 2
    if(is.null(control$isotropic.adj)) control$isotropic.adj <- TRUE
    if(is.null(control$polar)) control$polar <- TRUE
  } else {
    # set default values in dimension > 2
    if(is.null(control$isotropic.adj)) control$isotropic.adj <- FALSE
    if(!is.null(control$polar) && control$polar) warning("Polar coordinates only valid in dimension two")
    control$polar <- FALSE
    if(equi.adj != 0.0) warning("Equiangular adjustment only possible in dimension two.")
  }

  A <- cosh(sqrt(curvature)*D) # convert distance matrix to 'hyperbolic Gram matrix'
  n <- ncol(A)

  ## check for large/infinite values
  A.max = max(A)
  if(A.max > 1e8) warning("Gram Matrix contains values > 1e8. Rerun with smaller curvature parameter or rescaled distances.")
  if(is.infinite(A.max)) stop("Gram matrix contains infinite values. Rerun with smaller curvature parameter or rescaled distances.")

  ## Compute Eigendecomposition of A
  if(!control$use.eigs) { # Use ordinary eigendecomposition
    spec <- eigen(A,symmetric=TRUE)

    ## Extract leading Eigenvalue and Eigenvector
    lambda0 <- spec$values[1]
    x0 <- spec$vectors[,1]

    ## Extract lower tail of spectrum
    X <- spec$vectors[,(n-dim+1):n] # Last dim Eigenvectors
    spec.tail <- spec$values[(n-dim+1):n] # Last dim Eigenvalues
    A.frob <- sqrt(sum(spec$values^2)) # Frobenius norm of A

  } else { # Use reduced Eigendecomposition from APACK ('eigs_sym')
    requireNamespace("RSpectra") || stop("Package RSpectra needed for reduced Eigendecomposition.")
    requireNamespace("Matrix") || stop("Package Matrix needed for computation of Frobenius norm.")

    ## Extract leading Eigenvalue and Eigenvector
    perron <- RSpectra::eigs_sym(A,1,which="LA")
    lambda0 <- perron$values[1]
    x0 <- perron$vectors[,1]

    ## Extract lower tail of spectrum
    spec <- RSpectra::eigs_sym(A,dim,which="SA")
    X <- spec$vectors
    spec.tail <- spec$values
    A.frob <- Matrix::norm(A,"F") # Frobenius norm of A
  }

  x0 <- x0*sqrt(lambda0) # scale by Eigenvalue
  if(x0[1]<0) x0 <- -x0 # Flip sign if first element negative
  x.min <- min(x0) # find minimum

  if(!control$isotropic.adj) { # no isotropic adjustment: rescale Eigenvectors by Eigenvalues
    if(any(spec.tail > 0)) warning("Spectral Values have been truncated to zero. Try to use lower embedding dimension")
    X <- X %*% diag(sqrt(pmax(-spec.tail,0)))
  }

  directional <- diag(apply(X,1,function(x) 1/sqrt(sum(x^2)))) %*% X # convert to directional coordinates

  out <- list() # Allocate output list

  ## Calculate radial coordinate
  r <- sqrt((alpha*x0 -x.min)/(alpha*x0 + x.min)) ## multiplicative adjustment (scaling)
  out$r <- r

  ## Calculate polar coordinates if dimension is 2
  if(dim==2) {
    ## calculate polar angle
    theta <- atan2(X[,2],X[,1])

    ## Equiangular adjustment
    if(equi.adj > 0.0) {
      delta <- 2*pi/n
      angles <- seq(-pi,pi - delta, length=n)
      theta.equi <- angles[rank(theta,ties.method="first")] # Equi-spaced angles
      theta <- (1-equi.adj)*theta + equi.adj*theta.equi # convex combination of original and equi-spaced angles
      directional <- cbind(cos(theta),sin(theta)) # update directional coordinate
    }
    out$theta <- theta
  }

  out$directional <- directional

  ## Set Additional return values

  if(control$return.lorentz) {
    out$x0 <- x0
    out$X <- X
  }

  #out$strain.R2 <- (lambda0^2 + sum(pmax(-spec.tail,0)^2)) / A.frob^2 # R^2 of residual strain

  if(control$return.dist) { ## Calculate hyperbolic distances (can be time consuming)
    out$dist <- get.distance(r,directional,curvature,D)
    out$stress <- attr(out$dist,"stress")
    attr(out$dist,"stress") <- NULL
  }

  if(control$return.stress && !control$return.dist) out$stress <- get.stress(r,directional,curvature,D) ## Calculate stress only

  out$curvature <- curvature
  out$dim <- dim
  class(out) <- c("hydra","list")
  return(out)
}

hyperbolic.distance <- function(r1,r2,directional1,directional2,curvature) {
  # compute the hyperbolic distance of two points given in radial/directional coordinates in the Poincare ball
  iprod <- min(max(-1.0,sum(directional1*directional2)),1.0) # force between numerical -1.0 and 1.0 to eliminate rounding errors
  acosh.arg <- 1.0 + max(0.0,2*(r1^2 + r2^2 - 2*r1*r2*iprod)/((1 - r1^2)*(1 - r2^2))) # hyperbolic 'angle'; force numerical >= 1.0
  distance <- 1/sqrt(curvature) * acosh(acosh.arg) # hyperbolic distance between points i and j
  return(distance)
}

get.distance <- function(r, directional, curvature, D = NULL) {
  ## Calculate hyperbolic distance matrix and stress from radial and directional embedding coordinates
  n <- length(r) # number of embeded points
  hyper.dist <- matrix(0,n,n)
  stress.sq <- 0.0 # allocate squared stress
  for(i in 1:n) {
    for(j in 1:n) {
      if(i==j) {
        hyper.dist[i,j] <- 0.0
      } else {
        dist.ij <- hyperbolic.distance(r[i],r[j], directional[i,], directional[j,],curvature)
        hyper.dist[i,j] <- dist.ij
        if(!is.null(D)) stress.sq <- stress.sq + (dist.ij - D[i,j])^2 # update stress
      }
    }
  }
  if(!is.null(D)) attr(hyper.dist,"stress") <- sqrt(stress.sq)
  return(hyper.dist)
}

get.stress <- function(r, directional, curvature, D) {
  ## Calculate stress of embedding from radial/directional coordinate
  n <- length(r) # number of embeded points
  stress.sq <- 0.0 # allocate squared stress
  for(i in 1:n) {
    for(j in 1:n) {
      if(i!=j){
        dist.ij <- hyperbolic.distance(r[i],r[j], directional[i,], directional[j,],curvature)
        stress.sq <- stress.sq + (dist.ij - D[i,j])^2
      }
    }
  }
  return(stress=sqrt(stress.sq))
}

poincare.to.hyper <- function(r,directional) {
  # convert coorindates in the Poincare ball to reduced hyperbolic coordinates
  return(diag(2*r/(1-r^2),nrow=length(r)) %*% directional)

}

hyper.to.poincare <- function(X) {
  # convert reduced hyperbolic coordinates to  coordinates in the Poincare ball

  l2.norms <- apply(X,1,function(x) sqrt(sum(x^2)))
  directional <- diag(1/l2.norms) %*% X
  directional[l2.norms==0.0] <- 0.0 # define directional coordinate as zero for center point
  r <- l2.norms / (1 + sqrt(1 + l2.norms^2))
  if(ncol(X) == 2) {
    theta <- apply(directional,1,function(x) atan2(x[2],x[1]))
    theta[l2.norms==0.0] <- 0.0
  }
  else
    theta <- NULL
  return(list(r=r,directional=directional,theta=theta))

}
