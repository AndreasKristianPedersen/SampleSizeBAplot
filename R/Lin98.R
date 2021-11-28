#' Sample size for the agreement test proposed in Lin's 1998 paper
#'
#' Computes the necessary sample size to achieve a given power by way of the formula presented in 'Evaluation of statistical equivalence using limits of agreement and associated sample size calculation' by Lin et al.
#'
#' @param delta The clinical acceptable width of the LOA.
#'
#' @param sigma The standard deviation of the differences.
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
#'
#' @export
Lin98<-function (delta,sigma,alpha=0.05,beta=0.2,gamma=0.95) {
  zgamma<-qnorm(1-(1-gamma)/2)
  zalpha<-qnorm(1-alpha)
  zbeta<-qnorm(1-beta)
  Lin98n<-(sqrt((1+zgamma^2)/2)*zalpha+sqrt(zgamma^2/2)*zbeta)^2/(delta-zgamma*sigma)^2*sigma^2
  return(Lin98n)
}


