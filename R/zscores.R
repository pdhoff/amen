#' rank-based z-scores
#' 
#' Computes the normal scores corresponding to the ranks of a data vector
#' 
#' 
#' @usage zscores(y,ties.method="average")
#' @param y a numeric vector
#' @param ties.method method for dealing with ties
#' @return a numeric vector
#' @author Peter Hoff
#' @export zscores
zscores<-function(y,ties.method="average")
{
 qnorm( rank(y,na.last="keep",ties.method=ties.method)/(1+sum(!is.na(y))) )
}

