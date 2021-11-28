#' Sample size for the agreement test proposed in Shieh's 2020 paper
#'
#' Computes the necessary sample size to achieve a given power by way of the formula presented in 'Assessing Agreement Between Two Methods of Quantitative Measurements: Exact Test Procedure and Sample Size Calculation' by Shieh.
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
Shieh20<-function(delta, sigma,me=0,alpha=0.05,beta=0.2,gamma=0.95){
    prop0<-gamma #NULL CENTRAL PROPORTION
    Shieh20search<-function(n){
    #END OF SPECIFICATION
    pct<-1-(1-prop0)/2
    zp<-qnorm(pct)
    l<--delta
    u<-delta
    prop1<-pnorm((u-me)/sigma)-pnorm((l-me)/sigma)
    sigsq<-sigma^2
    thetal<-me-zp*sigma
    thetau<-me+zp*sigma
    df<-n-1
    std<-sqrt(sigsq/n)
    numint<-1000
    coevec<-c(1,rep(c(4,2),numint/2-1),4,1)
    cl<-1e-6
    cu<-qchisq(1-cl,df)
    int<-cu-cl
    intl<-int/numint
    cvec<-cl+intl*(0:numint)
    wcpdf<-(intl/3)*coevec*dchisq(cvec,df)
    fungam<-function () {
      gaml<-0
      gamu<-100
      dalpha<-1
      while(abs(dalpha)>1e-8 | dalpha<0){
        gam<-(gaml+gamu)/2
        h<-zp*sqrt(n)-gam*sqrt(cvec/df)
        ht<-h*(cvec<n*df*(zp/gam)^2)
        4
        alphat<-sum(wcpdf*(2*pnorm(ht)-1))
        if (alphat>alpha) gaml<-gam else gamu<-gam
        dalpha<-alphat-alpha
      }
      return(gam)
    }
    gam<-fungam()
    hel<-(l-me)/std+gam*sqrt(cvec/df)
    heu<-(u-me)/std-gam*sqrt(cvec/df)
    ke<-(n*df*(u-l)^2)/(4*gam^2*sigsq)
    helt<-hel*(cvec<ke)
    heut<-heu*(cvec<ke)
    gpower<-sum(wcpdf*(pnorm(heut)-pnorm(helt)))
    gpower
  }

  styrke<-1-beta
  pwr<-Shieh20search(3)
  n<-3
  nfinal<-n
  while (pwr<styrke & n<10000){
    n<-n+1
    pwr<-Shieh20search(n)
    nfinal<-n
  }
 return(nfinal)
}



