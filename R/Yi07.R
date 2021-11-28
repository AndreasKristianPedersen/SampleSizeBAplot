#' Sample size for the agreement test proposed in Yi's 2008 paper
#'
#' Computes the necessary sample size to achieve a given power by way of the formula presented in 'Reliability analysis for continuous measurements: Equivalence test for agreement' by Yi et al.
#'
#' @param sigma The standard deviation of the differences.
#'
#' @param delta The clinical acceptable variance of the differences.
#'
#' @param alpha Significance level. Default is set to 0.05.
#'
#' @param beta Type II error. Default is set to 0.2.
#'
#' @return Factor variable
#'
#' @examples
#' See *my paper*
#'@export
Yi08<-function (delta,sigma,alpha=0.05,beta=0.2){
  Yi98search<-function(n){
    v<-qchisq(1-beta,n)/qchisq(alpha,n)-delta^2/sigma^2
  }
  n<-2
  v1<-Yi98search(2)
  if (v1>0){
  while (v1>0 & n<100000){
    n<-n+1
    v1<-Yi98search(n)
    nfinal<-n
  }
  return(min(nfinal))
  }
}




