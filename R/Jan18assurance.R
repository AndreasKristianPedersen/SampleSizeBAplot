#' Sample size for the agreement test proposed in Jan and Shieh's 2018 paper
#'
#' Computes the necessary sample size to achieve a given width by way of the formula presented in 'Assessing Agreement Between Two Methods of Quantitative Measurements: Exact Test Procedure and Sample Size Calculation' by Jan and Shieh.
#'
#' @param sigma The standard deviation of the differences.
#'
#' @param delta The clinical acceptable width of the LOA.
#'
#' @param ap The assurance probability.
#'
#' @param me The mean of the differences. Default is set to 0.
#'
#' @param alpha Significance level. Default is set to 0.05.
#'
#' @param gamma Coverage of the LOA. Default is set to 0.95.
#'
#' @return Factor variable
#'
#' @examples
#' See *my paper*
#'@export
Jan18assurance<-function(sigma,delta,ap,me=0,alpha=0.05,gamma=0.95) {
  zp<-qnorm(1-(1-gamma)/2)
  sigsq<-sigma^2
  thetal<-me-zp*sigma
  thetau<-me+zp*sigma
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
  ape<-0
  while (ape<ap & n<20000) {
    n<-n+1
    df<-n-1
    g<-gfun()
    ape<-pchisq((df/g^2)*(delta/sigma)^2,df)
  }
  return(c(ap,n))
}
