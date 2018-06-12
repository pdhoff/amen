#' SRM log likelihood evaluated on a grid of rho-values
#'
#' Calculation of the SRM log-likelihood over a grid of rho-values
#'
#' @param Y sociomatrix assumed to follow a mean-zero SRM distribution
#' @param Sab covariance matrix of additive effects
#' @param rhos vector of rho-values at which to calculate the log-likelihood
#' @param s2 current value of s2
#' @return a vector of log-likelihood values
#' @author Peter Hoff
#' @export llsrmRho
llsrmRho<-function(Y,Sab,rhos,s2=1)
{
  n<-nrow(Y) 
  c<-(1/sqrt(1+rhos) + 1/sqrt(1-rhos))/(2*sqrt(s2))
  d<-(1/sqrt(1+rhos) - 1/sqrt(1-rhos))/(2*sqrt(s2)) 
  ihSe<-array(dim=c(2,2,length(rhos)))
  ihSe[1,1,]<-ihSe[2,2,]<-c 
  ihSe[1,2,]<-ihSe[2,1,]<-d
  Sabt<-array(apply(ihSe,3,function(x){x%*%Sab%*%x}),dim=c(2,2,length(rhos)))
  iSabt<-array(dim=dim(Sabt))
  iSabt[1,1,]<-Sabt[2,2,] ; iSabt[2,2,]<-Sabt[1,1,] 
  iSabt[1,2,]<-iSabt[2,1,]<- -Sabt[1,2,] 
  iSabt<-sweep(iSabt,3, Sabt[1,1,]*Sabt[2,2,] -  Sabt[2,1,]*Sabt[1,2,],"/") 
  G<-array(apply(iSabt,3,function(x){solve(x+n*diag(2))}),
           dim=c(2,2,length(rhos)))
  pH<-array(apply(iSabt,3,function(x){ -solve(x+n*matrix(1,2,2))}),
            dim=c(2,2,length(rhos)))
  H<-array(dim=dim(pH))
  H[1,1,]<-pH[1,1,]*G[2,1,] + pH[1,2,]*G[1,1,] 
  H[1,2,]<-pH[1,1,]*G[2,2,] + pH[1,2,]*G[1,2,]
  H[2,1,]<-pH[2,1,]*G[2,1,] + pH[2,2,]*G[1,1,]
  H[2,2,]<-pH[2,1,]*G[2,2,] + pH[2,2,]*G[1,2,]
  sumH<-apply(H,3,sum)  
  ## G, H look ok

  ldV<- n^2*log(s2)+choose(n+1,2)*log(1+rhos)+choose(n,2)*log(1-rhos) +
        (n-1)*apply(Sabt,3,function(x){ log(det(diag(2)+n*x))})  + 
        apply(Sabt,3,function(x){log(det(diag(2)+ n*x%*%matrix(1,2,2)))}) 
  # ldV looks correct

  sumZ<- (ihSe[1,1,] + ihSe[1,2,])*sum(Y)
  sumZ2<-(ihSe[1,1,]^2+ihSe[1,2,]^2)*sum(Y^2) + 
         2*ihSe[1,1,]*ihSe[1,2,]*sum(Y*t(Y)) 
  # sumZ and sumZ2 look correct

  rsy<-rowSums(Y) ; csy<-colSums(Y)
  rr<-sum(rsy^2) ; rc<-sum(rsy*csy) ; cc<-sum(csy^2)  
  c<-ihSe[1,1,] ; d<-  ihSe[1,2,]
 
  B<- G[1,1,]*( c^2*rr + d^2*cc + 2*c*d*rc )  +  
      2*G[1,2,]*( (c^2+d^2)*rc + c*d*(rr+cc)   )  + 
      G[2,2,]*( c^2*cc +d^2*rr + 2*c*d*rc )

  yPy<- sumZ2 - ( B + sumZ^2*sumH )

  -.5*(ldV+yPy)
}



#' Griddy Gibbs update for dyadic correlation
#'
#' Simulation of dyadic correlation from its approximate full conditional 
#' distribution using griddy Gibbs sampling
#'
#' @param Z n X n normal relational matrix
#' @param Sab covariance of additive effects
#' @param s2 residual variance
#' @param offset matrix of the same dimension as Z. It is assumed that 
#' Z-offset follows an SRM distribution, so the offset should contain any 
#' regression terms and multiplicative effects (such as 
#' \code{Xbeta(X,beta+ U\%*\%t(V) }   )
#' @param ngp the number of points for an unevenly-spaced 
#' grid on which to approximate the full conditional distribution
#' 
#' @return a value of rho
#' @author Peter Hoff
#' @export rrho_fc
rrho_fc<-function(Z,Sab,s2=1,offset=0,ngp=100)
{ 
  E<-Z-offset

  ## first obtain rough estimate of rho and its sd
  m<-mean(E,na.rm=TRUE)
  a<-apply(E,1,mean,na.rm=TRUE) - m
  b<-apply(E,2,mean,na.rm=TRUE) - m
  R<-E-(m+outer(a,b,"+"))
  RM<-cbind(R[upper.tri(R)],t(R)[upper.tri(R)] )/sqrt(s2)
  emcp<-sum(RM[,1]*RM[,2])
  emss<-sum(RM^2)
  nd<- nrow(RM)
  rho0<-cor(RM[,1],RM[,2])
  srho0<- 2*(1-rho0^2)/sqrt(nd)

  ## grid concentrated near likely regions of high probability 
  rhos<-qnorm( (1:ngp)/(ngp+1), rho0,2*srho0 )
  rhos<-rhos[ -1<rhos & rhos<1 ]

  ## griddy Gibbs
  ll<-llsrmRho(E,Sab,rhos,s2)
  prho<-exp( ll-max(ll) - .5*log(1-rhos^2) )
  Frho<-c(0,cumsum(prho)/sum(prho),1)
  rhos<-c(0,rhos,1)
  f<-runif(1)
  k<-max(which(Frho<f))
  rhos[k] + (rhos[k+1]-rhos[k])*(f-Frho[k])/(Frho[k+1]-Frho[k])
}




