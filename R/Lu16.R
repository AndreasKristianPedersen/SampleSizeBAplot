#' Sample size for the agreement test proposed in Lu's 2016 paper
#'
#' Computes the necessary sample size to achieve a given power by way of the formula presented in 'Sample Size for Assessing Agreement between Two Methods of Measurement by Bland-Altman Method' by Lu et al.
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
Lu16<-function(delta, sigma,me=0,alpha=0.05,beta=0.2,gamma=0.95){
  zgamma<-qnorm(1-(1-gamma)/2)
  Lu16search<-function(n){
    tau1<-(delta-me-zgamma*sigma)/(sigma*sqrt(1/n+zgamma^2/(2*(n-1))))
    tau2<-(delta+me-zgamma*sigma)/(sigma*sqrt(1/n+zgamma^2/(2*(n-1))))
    quantilet<-qt(1-alpha/2,n-1)
    beta1<-pt(quantilet,n-1,tau1,lower.tail = FALSE)
    beta2<-pt(quantilet,n-1,tau2, lower.tail = FALSE)
    power<-(beta1+beta2)-1
    power
   }
  styrke<-1-beta
  pwr<-Lu16search(3)
  n<-3
  nfinal<-3
  while (pwr<styrke){
    n<-n+1
    pwr<-Lu16search(n)
    nfinal<-n
  }
  return(nfinal)
}


