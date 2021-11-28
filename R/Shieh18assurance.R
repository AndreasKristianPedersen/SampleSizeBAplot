#' Sample size to acquire a given width with a given assurance by way of the confidence interval for the limits of agreement proposed in Shieh's 2018 paper
#'
#' Computes the necessary sample size to achieve a given precision with a given assurance probability for the LOA by way of the confidence interval presented in 'The appropriateness of Bland-Altman's approximate confidence intervals for limits of agreement' by Shieh.
#'
#' @param sigma The standard deviation of the differences.
#'
#' @param ew The clinical acceptable precision of the confidence interval.
#'
#' @param ap Assurance probability.
#'
#' @param me Mean differences. Default is set to 0.
#'
#' @param alpha alpha level of the confidence interval. Default is set to 0.05.
#'
#' @param gamma Coverage of the LOA. Default is set to 0.95.
#'
#' @return Factor variable
#'
#' @examples
#' See *my paper*
#'@export
Shieh18assurance<-function(sigma,ew,ap,alpha=0.05,me=0,gamma=0.95){
  zp<-qnorm(1-(1-gamma)/2)
  sigsq<-sigma^2
  theta<-mu+zp*sigma
  n<-4
  ape<-0
  while(ape<ap){
    n<-n+1
    df<-n-1
    logc<-log(sqrt(df/2))+lgamma(df/2)-lgamma(n/2)
    c<-exp(logc)
    tl<-qt((1-alpha)/2,df,zp*sqrt(n))
    tu<-qt(1-(1-alpha)/2,df,zp*sqrt(n))
    td<-tu-tl
    qe<-(n*df*ew^2)/(td^2*sigsq)
    ape<-pchisq(qe,df)
  }
  return(c(ap,n))
}
