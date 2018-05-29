#' Precomputation of design matrix quantities.
#'
#' Computation of a variety of quantities from the design array
#' to be used in MCMC model fitting algorithms.
#' 
#' @usage precomputeX(X) 
#' @param X a three way array, the design array for an AME model 
#' @return the same three-way array but with derived quantities 
#' as attributes. 
#' @author Peter Hoff
#' @export precomputeX
precomputeX<-function(X)
{
  ### marginal means and regression sums of squares
  Xr<-apply(X,c(1,3),sum)            # row sum
  Xc<-apply(X,c(2,3),sum)            # col sum
  mX<- apply(X,3,c)                  # design matrix
  mXt<-apply(aperm(X,c(2,1,3)),3,c)  # dyad-transposed design matrix
  XX<-t(mX)%*%mX                     # regression sums of squares
  XXt<-t(mX)%*%mXt                   # crossproduct sums of squares

  attr(X,"Xr")<-Xr 
  attr(X,"Xc")<-Xc
  attr(X,"mX")<-mX
  attr(X,"mXt")<-mXt
  attr(X,"XX")<-XX
  attr(X,"XXt")<-XXt

  X
}


