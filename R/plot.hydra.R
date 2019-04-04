#' @title Plot a hyperbolic embedding
#'
#' @description  Plot a two-dimensional hyperbolic embedding as returned by \code{\link{hydra}} in the Poincare disc
#'
#' @param x a hydra object as returned by \code{\link{hydra}} with dimension \code{dim} equal 2
#' @param labels character labels for the embedded points, supplied as a vector. NULL triggers default values
#' @param node.col colors for the labels and/or points, supplied as a vector. NULL triggers default values. See `Color Specification' in \code{\link[graphics]{par}} for details
#' @param graph.adj a graph adjacency matrix that is used to plot links between the embedded points (links are drawn for all non-zero elements of \code{graph.adj})
#' @param pch plotting 'characters' for the embedded points. supplied as a vector. NULL triggers default values. See \code{\link[graphics]{points}} for details
#' @param crop.disc should the Poincare disc be cropped or fully shown? Defaults to TRUE
#' @param shrink.disc if true, the Poincare disc is shrunk to tightly fit all plotted points. Defaults to FALSE
#' @param disc.col color of the Poincare disc. Set to "white" to hide disc
#' @param rotation rotate points by this angle (specified in degrees) around the center of the Poincare disc
#' @param mark.center Should a cross be placed at the center of the disc? If 0, nothing is drawn. Other values specify the relative size of the cross mark.
#' @param mark.angles Should the angular coordinates of points be marked at the boundary of the disc? If 0, nothing is drawn. Other values specify the relative size of the angle marks.
#' @param mildify large values reduce the curvature of links. Values around 3 are visually most appealing. Setting \code{mildify} to 1 shows the true hyperbolic curvature
#' @param cex character expansion for labels and points, supplied as a numerical vector. See also \code{\link[graphics]{points}}
#' @param ... all other parameters are passed on as additional graphical parameters (see \code{\link[graphics]{par}})
#'
#'
#' @author Martin Keller-Ressel <martin.keller-ressel@tu-dresden.de>
#'
#' @examples
#' data(karate)
#' embedding <- hydra(karate$distance)
#' plot(embedding,labels=karate$label,node.col=karate$group,graph.adj=karate$adjacency)
#'
#' # plot points instead of labels, hide Poincare disc and rotate by 90 degrees:
#' plot(embedding,pch=karate$group, node.col=karate$group,graph.adj=karate$adjacency, disc.col="white",
#'      rotation=90)
#'
#' # do not crop the Poincare disc, mark the center and mark angles:
#' plot(embedding,labels=karate$label, node.col=karate$group,graph.adj=karate$adjacency,
#'      crop.disc=FALSE, mark.center=0.05, mark.angles=0.025)
#'
#' @importFrom graphics lines plot points segments symbols text
#'
#' @export

plot.hydra <- function(x,labels=NULL,node.col=1,pch=NULL,graph.adj=NULL,crop.disc=TRUE,shrink.disc=FALSE, disc.col = "grey90", rotation = 0, mark.center=0, mark.angles=0, mildify = 3,cex=1.0,...) {

  hydra <- x
  if(shrink.disc) c.radius <- max(hydra$r) else c.radius <- 1.0
  c.radius <- c.radius + 0.03*cex

  rotation.rad <- rotation/180 * pi

  # compute x and y coordinate of points
  x.nodes <- hydra$r*cos(hydra$theta + rotation.rad)
  y.nodes <- hydra$r*sin(hydra$theta + rotation.rad)

  # plot points invisibly to set frame
  if(crop.disc) plot(x.nodes,y.nodes,type="n",asp=1, xlab="",ylab="",axes=FALSE,...)
  else plot(x.nodes,y.nodes,type="n",asp=1, xlab="",ylab="",axes=FALSE,xlim=c(-1,1)*c.radius,ylim=c(-1,1)*c.radius,...)

  # draw Poincare disc
  symbols(0,0,circles=c.radius,fg="white",bg=disc.col,inches=FALSE, add=TRUE)
  if(mark.center != 0.0) {
    lines(mark.center*c.radius*c(-1,1),c(0,0),...)
    lines(c(0,0),mark.center*c.radius*c(-1,1),...)
  }

  # draw graph links as hyperbolic geodesics
  if(!is.null(graph.adj)) {
    edgelist <- which(graph.adj != 0,arr.ind=TRUE)
    t <- seq(from=0,to=1,length=100)
    for(i in 1:nrow(edgelist)) {
      r <- hydra$r[edgelist[i,]] / mildify
      dir <- hydra$directional[edgelist[i,],]
      X <- poincare.to.hyper(r,dir)
      x0 <- apply(X,1,function(x) {1 + sum(x^2)})
      l.prod <- x0[1]*x0[2] - sum(X[1,] * X[2,])
      dist <- acosh(max(1,l.prod))
      aux <- (X[2,] - X[1,] * l.prod) / sqrt(l.prod^2 - 1)
      geodesic <- hyper.to.poincare(outer(cosh(t*dist),X[1,]) + outer(sinh(t*dist),aux))
      lines(geodesic$r* cos(geodesic$theta + rotation.rad) * mildify ,geodesic$r* sin(geodesic$theta + rotation.rad) * mildify,lty="dotted",col="black",...)
    }
  }

  # erase links around points/labels
  symbols(x.nodes,y.nodes,circles=rep(0.03*cex,length(hydra$r)),fg=disc.col,bg=disc.col,inches=FALSE, add=TRUE)
  # draw labels
  if(!is.null(labels)) text(x.nodes,y.nodes,labels,col=node.col,cex=cex,...)
  # draw points
  if(is.null(labels) && is.null(pch)) pch <- 21 # default pch
  if(!is.null(pch)) points(x.nodes,y.nodes,col=node.col,pch=pch,cex=cex,...)
  # draw angle marks
  if(mark.angles != 0) segments((1-mark.angles)*c.radius*cos(hydra$theta + rotation.rad),(1-mark.angles)*c.radius*sin(hydra$theta + rotation.rad),(1+mark.angles)*c.radius*cos(hydra$theta + rotation.rad),(1+mark.angles)*c.radius*sin(hydra$theta + rotation.rad),col=node.col,...)
}
