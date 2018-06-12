#' Gibbs sampling of U and V
#' 
#' A Gibbs sampler for updating the multiplicative effect matrices U and V
#' 
#' @param Z n X n normal relational matrix
#' @param U current value of U
#' @param V current value of V 
#' @param Suv covariance of (U V) 
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @param offset a matrix of the same dimension as Z. It is assumed that   
#' Z-offset is equal to the multiplicative effects plus dyadic noise, so the 
#' offset should contain any additive effects (such as \code{Xbeta(X,beta+ 
#' outer(a,b,"+")  }  )
#'
#' @return \item{U}{a new value of U} \item{V}{a new value of V}
#' @author Peter Hoff
#' @export rUV_fc
rUV_fc <-
function(Z,U,V,Suv,rho,s2=1,offset=0)
{ 
  E<-Z-offset

  R<-ncol(U) 

  Se<-matrix(c(1,rho,rho,1),2,2)*s2
  iSe2<-mhalf(solve(Se))
  g<-iSe2[1,1] ; d<-iSe2[1,2]

  for(r in sample(1:R))
  {
    Er<- E - U[,-r]%*%t(V[,-r]) ;  Es<- (g^2+d^2)*Er+2*g*d*t(Er)

    ## update Ur
    vr<-V[,r]
    b0<- c(Suv[r,-r]%*%solve(Suv[-r,-r]))
    v0<- c(Suv[r,r] - b0%*%Suv[-r,r])
    m0<- cbind(U[,-r],V)%*%b0
    ssv<-max(sum(vr^2),1e-6)
    a<- (g^2+d^2)*ssv+1/v0 ; c<- -2*g*d/(a^2+a*2*g*d* ssv) 
    Esv<-Es%*%vr 
    m1<- Esv/a + c*vr*sum((Esv+m0/v0)*vr)  + m0/(a*v0)  
    #Vh<-sqrt(1/a)*diag(nrow(E))+(vr%*%t(vr))*(sqrt(1/a+ssv*c)-sqrt(1/a))/ssv 
    ah<-sqrt(1/a) ; bh<-(sqrt(1/a+ ssv*c)- sqrt(1/a) )/ssv ; e<-rnorm(nrow(E)) 
    U[,r]<- m1 + ah*e + bh*vr*sum(vr*e) 
    ##

    ## update Vr 
    ur<-U[,r] 
    rv<-R+r
    b0<- c(Suv[rv,-rv]%*%solve(Suv[-rv,-rv]))
    v0<- c(Suv[rv,rv] - b0%*%Suv[-rv,rv])
    m0<- cbind(U,V[,-r])%*%b0
    ssu<-max(sum(ur^2),1e-6)
    a<- (g^2+d^2)*ssu+1/v0 ; c<- -2*g*d/(a^2+a*2*g*d* ssu)
    tEsu<-t(Es)%*%ur
    m1<- tEsu/a + c*ur*sum((tEsu+m0/v0)*ur)  + m0/(a*v0)  
    #Vh<-sqrt(1/a)*diag(nrow(E))+(ur%*%t(ur))*(sqrt(1/a+ssu*c)-sqrt(1/a))/ssu
    ah<-sqrt(1/a) ; bh<-(sqrt(1/a+ ssu*c)- sqrt(1/a) )/ssu ; e<-rnorm(nrow(E)) 
    V[,r]<- m1 + ah*e + bh*ur*sum(ur*e) 
    ###
  }

list(U=U,V=V)
}
