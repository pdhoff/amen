#' Summary of an AME object
#' 
#' Summary method for an AME object
#' 
#' @param object the result of fitting an AME model
#' @param ... additional parameters (not used)
#' @return a summary of parameter estimates and confidence intervals for an AME
#' fit
#' @author Peter Hoff
#' @examples
#' 
#' ## The function is currently defined as
#' function(object, ...) 
#' {
#'   fit<-object
#' require(amen)
#'     tmp <- cbind(apply(fit$BETA, 2, mean), apply(fit$BETA, 2, 
#'         sd), apply(fit$BETA, 2, mean)/apply(fit$BETA, 2, sd), 
#'         2 * (1 - pnorm(abs(apply(fit$BETA, 2, mean)/apply(fit$BETA, 
#'             2, sd)))))
#'     colnames(tmp) <- c("pmean", "psd", "z-stat", "p-val")
#'     cat("\nbeta:\n")
#'     print(round(tmp, 3))
#'     tmp <- apply(fit$SABR, 2, mean)
#'     Sab <- matrix(tmp[c(1, 2, 2, 3)], 2, 2)
#'     dimnames(Sab) <- list(c("a", "b"), c("a", "b"))
#'     cat("\nSigma_ab pmean:\n")
#'     print(round(Sab, 3))
#'     cat("\nrho pmean:\n", round(tmp[4], 3), "\n")
#'   }
#' 
#' @S3method summary ame
summary.ame <-
function(object, ...)
{ 
  fit<-object
  require(amen)
  tmp<-cbind(apply(fit$BETA,2,mean), apply(fit$BETA,2,sd) ,
       apply(fit$BETA,2,mean)/apply(fit$BETA,2,sd) , 
       2*(1-pnorm( abs(apply(fit$BETA,2,mean)/apply(fit$BETA,2,sd)))))
  colnames(tmp)<-c("pmean","psd","z-stat","p-val") 
  cat("\nRegression coefficients:\n")
  print(round(tmp,3))


  tmp<-cbind(apply(fit$VC,2,mean), apply(fit$VC,2,sd) )
  colnames(tmp)<-c("pmean","psd") 
  cat("\nVariance parameters:\n")
  print(round(tmp,3))


}
