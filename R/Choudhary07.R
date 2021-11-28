#' Sample size for the agreement test proposed in Choudhary and Nagarajs 2007 paper
#'
#' Computes the necessary sample size to achieve a given power by way of the formula presented in 'Tests for assessment of agreement using probability criteria' by Choudhary and Nagaraja.
#'
#' @param delta The clinical acceptable width of the LOA.
#'
#' @param sigma The standard deviation of the differences.
#'
#' @param me The mean of the differences. Default is set to 0.
#'
#' @param alpha Significance level. Default is set to 0.05.
#'
#' @param beta Type II error. Default is set to 0.2.
#'
#' @param gamma Coverage of the LOA. Default is set to 0.95.
#'
#' @return Factor variable
#'
#' @examples
#' See *my paper*
#'@export
Choudhary07<-function(delta, sigma,me=0,alpha=0.05,beta=0.2,gamma=0.95) {
  zp0<-qnorm(1-gamma)
  zalpha<-qnorm(1-alpha)
  zbeta<-qnorm(1-beta)
  p1<-pnorm(delta,mean=me,sd=sigma)
  zp1<-qnorm(1-p1)
  k<-(zbeta*zp0+zalpha*zp1)/(zbeta+zalpha)
  nChaudary07<-(1+0.5*k^2)*((zalpha+zbeta)/(zp0-zp1))^2
  return(nChaudary07)
}




