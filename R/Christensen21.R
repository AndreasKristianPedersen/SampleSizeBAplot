#' Sample size to acquire a given width by way of the confidence interval for the limits of agreement proposed in Christensen et al 2021 paper
#'
#' Computes the necessary sample size to achieve a given precision for the LOA by way of the confidence interval presented in 'The appropriateness of Bland-Altman's approximate confidence intervals for limits of agreement' by Christensen et al.
#'
#' @param sigmab The standard deviation between individuals.
#'
#' @param sigmaw The standard deviation within individuals.
#'
#' @param ew The clinical acceptable width of the confidence interval.
#'
#' @param gamma Coverage of the LOA. Default is set to 0.95.
#'
#' @param alpha alpha level of the confidence interval. Default is set to 0.05.
#'
#' @return Factor variable
#'
#' @examples
#' See *my paper*
#'@export
Christensen<-function(sigmab,sigmaw,ew,gamma=0.95, alpha=0.05){
Christensensearch<-function(n){
  hb<-1/qchisq((1-gamma)/2, 1, ncp = 0, lower.tail = TRUE, log.p = FALSE)-1
  lb<-1-1/qchisq(1-(1-gamma)/2, 1, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  hw<-1/qchisq((1-gamma)/2, n-1, ncp = 0, lower.tail = TRUE, log.p = FALSE)-1
  lw<-1-1/qchisq(1-(1-gamma)/2, n-1, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  H0<-sqrt(hb^2*(n*sigmab^2+sigmaw^2)^2+hw^2*(n-1)^2*(sigmaw^2)^2)
  L0<-sqrt(lb^2*(n*sigmab^2+sigmaw^2)^2+lw^2*(n-1)^2*(sigmaw^2)^2)
  diff<-qnorm(1-(1-alpha)/2)/sqrt(2*n)*(sqrt((n*sigmab^2+sigmaw^2)+(n-1)*sigmaw^2+H0)-sqrt((n*sigmab^2+sigmaw^2)+(n-1)*sigmaw^2-L0))-ew
  diff
}
n<-2
v1<-Christensensearch(2)
if (v1<0){
 return(n)
}
else{
while (v1>0 & n<10000){
  n<-n+1
  v1<-Christensensearch(n)
  nfinal<-n
}
  if (nfinal<10000){
return(min(nfinal))
  }
  else{
    return("Width is unachievable or demands over 10000 observations")
  }
}
}




