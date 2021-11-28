#' Sample size to acquire a given width by way of the confidence interval for the limits of agreement proposed in Jan and Shiehs 2018 paper
#'
#' Computes the necessary sample size to achieve a given precision for the LOA by way of the confidence interval presented in 'The appropriateness of Bland-Altman's approximate confidence intervals for limits of agreement' by Jan and Shieh.
#'
#' @param sigma The standard deviation of the differences.
#'
#' @param ew The clinical acceptable precision of the confidence interval.
#'
#' @param alpha alpha level of the confidence interval. Default is set to 0.05.
#'
#' @param mu Mean of the differences. Default is set to 0.
#'
#' @param gamma Coverage of the LOA. Default is set to 0.95.
#'
#' @return Factor variable
#'
#' @examples
#' See *my paper*
#'@export
Jan18<-function(sigma,ew,alpha=0.05,mu=0,gamma=0.95){
  zp<-qnorm(1-(1-gamma)/2)
  sigsq<-sigma^2
  thetal<-mu-zp*sigma
  thetau<-mu+zp*sigma
  coverp<-1-alpha
  numint<-1000
  coevec<-c(1,rep(c(4,2),numint/2-1),4,1)
  gfun<-function () {
    cql<-10e-8
    cqu<-qchisq(1-cql,df)
    int<-cqu-cql
    intl<-int/numint
    cvec<-cql+intl*(0:numint)
    wcpdf<-(intl/3)*coevec*dchisq(cvec,df)
    gl<-0
    gu<-100;
    dd<-1
    while (abs(dd)>10e-09 | dd<0) {
      g<-(gl+gu)/2
      b<-sqrt(n)*(-zp+g*sqrt(cvec/df))
      cpt<-sum(wcpdf*((2*pnorm(b)-1)*(cvec>df*(zp/g)^2)))
      if (cpt>coverp) gu<-g
      else gl<-g
      dd<-cpt-coverp
    }
    return(g)
  }
  n<-4
  ehe<-1000
  while (ehe>ew & n<10000) {
    n<-n+1
    df<-n-1
    logc<-log(sqrt(df/2))+lgamma(df/2)-lgamma(n/2)
    c<-exp(logc)
    g<-gfun()
    ehe<-(g/c)*sigma
  }
  return(n)
}
