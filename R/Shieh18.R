#' Sample size to acquire a given width of the confidence interval for the limits of agreement proposed in Shieh's 2018 paper
#'
#' Computes the necessary sample size to achieve a given precision for the LOA by way of the confidence interval presented in 'The appropriateness of Bland-Altman's approximate confidence intervals for limits of agreement' by Shieh.
#'
#' @param sigma The standard deviation of the differences.
#'
#' @param ew The clinical acceptable width of the confidence interval.
#'
#' @param alpha alpha level of the confidence interval. Default is set to 0.05.
#'
#' @param gamma Coverage of the LOA. Default is set to 0.95.
#'
#' @param mu Mean of the differences.
#'
#' @return Factor variable
#'
#' @examples
#' See *my paper*
#'@export
Shieh18expectedwidth<-function(sigma,ew,gamma=0.95,alpha=0.05,mu=0){
  zp<-qnorm(1-(1-gamma)/2)
  sigsq<-sigma^2
  theta<-mu+zp*sigma
  coverp<-1-alpha
  n<-4
  ewe<-100
  while(ewe>ew){
    n<-n+1
    df<-n-1
    logc<-log(sqrt(df/2))+lgamma(df/2)-lgamma(n/2)
    c<-exp(logc)
    tl<-qt(alpha/2,df,zp*sqrt(n))
    tu<-qt(1-alpha/2,df,zp*sqrt(n))
    td<-tu-tl
    ewe<-td*sigma/(c*sqrt(n))
  }
  return(n)
}
