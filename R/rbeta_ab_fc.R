#' Conditional simulation of additive effects and regression coefficients
#' 
#' Simulates from the joint full conditional distribution of (beta,a,b)
#' in a social relations regression model
#' 
#' @param Z n X n normal relational matrix
#' @param Sab row and column covariance
#' @param rho dyadic correlation
#' @param X n x n x p covariate array
#' @param s2 dyadic variance
#' @param offset a matrix of the same dimension as Z. It is assumed that 
#' Z-offset follows a SRRM, so the offset should contain any multiplicative 
#' effects (such as \code{U\%*\% t(V) } )
#' @param iV0 prior precision matrix for regression parameters
#' @param m0 prior mean vector for regression parameters 
#' @param g prior variance scale for g-prior when iV0 is unspecified 
#' 
#' @return \item{beta}{regression coefficients} \item{a}{additive row effects}
#' \item{b}{additive column effects}
#' @author Peter Hoff
#' @export rbeta_ab_fc
rbeta_ab_fc<-
function(Z,Sab,rho,X=NULL,s2=1,offset=0,iV0=NULL,m0=NULL,g=length(Z))
{
  Z<-Z-offset

  ### make fake design matrix if none provided
  if(is.null(X)){ X<-design_array(intercept=FALSE,n=nrow(Z)) }
  ###

  ### calculate statistics of X 
  if(is.null(attributes(X)$XX))
  { 
    X<-precomputeX(X) 
    warning("Summary statistics of X are not precomputed. ", 
             "Run X<-precomputeX(X) to speed up calculations.") 
  }
  ### 

  ### assign statistics of X
  Xr<-attributes(X)$Xr 
  Xc<-attributes(X)$Xc
  mX<-attributes(X)$mX
  mXt<-attributes(X)$mXt
  XX<-attributes(X)$XX
  XXt<-attributes(X)$XXt
  ###

  ### 
  p<-dim(X)[3]
  n<-nrow(Z)
  ###

  ### set priors 
  if(p>0 & is.null(iV0))
  { 
    # g-prior plus small ridge in case XX is singular
    iV0<-XX/g + diag(diag(XX),nrow=nrow(XX))/g^2  

    # now flatten prior on intercept
    if(all(mX[,1]==1))
    {  
      V0<-solve(iV0) 
      V0[1,1]<-V0[1,1] + sqrt(g) - g/n^2 
      iV0<-solve(V0)  
    } 
  } 

  if(is.null(m0)){ m0<-rep(0,dim(X)[3]) } 
  ###

  ### decorrelation
  Se<-matrix(c(1,rho,rho,1),2,2)*s2
  iSe2<-mhalf(solve(Se))
  td<-iSe2[1,1] ; to<-iSe2[1,2]
  Sabs<-iSe2%*%Sab%*%iSe2
  tmp<-eigen(Sabs)
  k<-sum(zapsmall(tmp$val)>0 )

  mXs<-td*mX+to*mXt                  # matricized transformed X
  XXs<-(to^2+td^2)*XX + 2*to*td*XXt  # sum of squares for transformed X
  Zs<-td*Z+to*t(Z)
  zr<-rowSums(Zs) ; zc<-colSums(Zs) ; zs<-sum(zc) 
  ###

  ### dyadic and prior contributions  
  if(p>0)
  {
    lb<- crossprod(mXs,c(Zs)) + iV0%*%m0
    Qb<- XXs + iV0 
  }
  ###

  ### row and column reduction
  ab<-matrix(0,nrow(Z),2)
  if(k>0)
  {
    G<-tmp$vec[,1:k] %*% sqrt(diag(tmp$val[1:k],nrow=k))
    K<-matrix(c(0,1,1,0),2,2)
    A<-n*t(G)%*%G + diag(k)
    B<-t(G)%*%K%*%G
    iA0<-solve(A)
    C0<- -solve(A+ n*B)%*%B%*%iA0

    iA<-G%*%iA0%*%t(G)
    C<-G%*%C0%*%t(G)

    if(p>0) 
    {    
    Xsr<-td*Xr + to*Xc  # row sums for transformed X
    Xsc<-td*Xc + to*Xr  # col sums for transformed X
    Xss<-colSums(Xsc)  
    lb<- lb - (iA[1,1]*crossprod(Xsr,zr) + iA[2,2]*crossprod(Xsc,zc) +
               iA[1,2]*(crossprod(Xsr,zc) + crossprod(Xsc,zr)) +
               sum(C)*Xss*zs )

    tmp<-crossprod(Xsr,Xsc)
    Qb<- Qb - (iA[1,1]*crossprod(Xsr,Xsr) + iA[2,2]*crossprod(Xsc,Xsc) +
               iA[2,1]*(tmp+t(tmp)) +  sum(C)*Xss%*%t(Xss) ) 
    }
  }
  ###

  ### if no covariates
  if(dim(X)[3]==0){ beta<-numeric(0) }
  ###

  ### if covariates 
  if(p>0) 
  { 
    V<-solve(Qb)
    m<-V%*%(lb)
    beta<-c(rmvnorm(1,m,V)) 
  }
  ###
 
  ### simulate a, b 
  if(k>0) 
  {
    E<- Zs-Xbeta(td*X+to*aperm(X,c(2,1,3)),beta)
    er<-rowSums(E) ; ec<-colSums(E) ; es<-sum(ec) 
    m<-t(t(crossprod(rbind(er,ec),t(iA0%*%t(G)))) + rowSums(es*C0%*%t(G)) )
    hiA0<-mhalf(iA0)
    e<-matrix(rnorm(n*k),n,k) 
    w<-m+ t( t(e%*%hiA0) - c(((hiA0-mhalf(iA0+n*C0))/n)%*% colSums(e) ) )
    ab<- w%*%t(G)%*%solve(iSe2) 
  }

list(beta=beta,a=ab[,1],b=ab[,2] )  
}
