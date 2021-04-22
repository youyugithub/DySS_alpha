##' Monitor Longitudinal Data for Survival Outcomes
##'
##' @title Monitor Longitudinal Data for Survival Outcomes
##' @param yyijk_new an array of longitudinal observations. \cr
##' \code{yyijk_new[i,j,k]} is the jth observation of the kth dimension of the ith subject.
##' @param ttij_new a matrix of observation times. \cr
##' \code{ttij_new[i,j]} is the jth observation time of the ith subject. \cr
##' \code{yyijk_new[i,j,]} is observed at \code{ttij[i,j]}.
##' @param nobs_new an integer vector for number of observations. \cr
##' \code{nobs_new[i]} is the number of observations for the ith subject.
##' @param pattern the estimated longitudinal and survival pattern
##' @param method a string \cr
##' If \code{method="risk"}
##' If \code{method="joint"}
##' @param parameter a numeric value. \cr
##' The weighting parameter in the modified EWMA charts.
##' @param CL a numeric value. \cr
##' The control limit
mnt_long_surv<-function(
  yyijk_new,ttij_new,nobs_new,
  pattern,method,
  parameter=0.5,
  CL=Inf)
{
  if(any(dim(yyijk_new)[1:2]!=dim(ttij_new)))stop("Dimensions of 'yyijk_new' and 'ttij_new' don't match.")
  if(dim(yyijk_new)[1]!=length(nobs_new))stop("Dimensions of 'yyijk_new' and 'nobs_new' don't match.")
  
  nind<-dim(yyijk_new)[1]
  nmaxobs<-dim(yyijk_new)[2]
  ndim<-dim(yyijk_new)[3]
  
  ttij_new<-clean_matij_by_nobs(ttij_new,nobs_new,"ttij_new")
  for(kk in 1:ndim){
    yyijk_new[,,kk]<-clean_matij_by_nobs(yyijk_new[,,kk],nobs_new,"yyijk_new")
  }
  
  ttmax<-pattern$ttmax
  ttmin<-pattern$ttmin
  ntimepoints<-pattern$ntimepoints
  if(ttmin>min(ttij_new,na.rm=TRUE))stop("Value of 'ttmin' is too large.")
  if(ttmax<min(ttij_new,na.rm=TRUE))stop("Value of 'ttmax' is too small.")
  omega<-(ttmax-ttmin)/(ntimepoints-1)
  ttij_int_new<-round((ttij_new-ttmin)/omega)+1
  
  if(method=="risk"){
    if(parameter>1.0|parameter<0.0)
      stop("For method 'risk' parameter should be in the interval [0,1].")
  }else{
    stop("Error in the argument 'method'.")
  }
  
  if(method=="risk"){
    
    beta_est<-pattern$beta_est
    mean_rr_est<-pattern$mean_rr_est
    var_rr_est<-pattern$var_rr_est
    delta_bar<-pattern$delta_bar
    sd_rr_est<-sqrt(var_rr_est)
    
    rrij_new<-matrix(NA,nind,nmaxobs)
    eeij_new<-matrix(NA,nind,nmaxobs)
    for(ii in 1:nind){
      rrij_new[ii,1:nobs_new[ii]]<-yyijk_new[ii,1:nobs_new[ii],]%*%beta_est
      eeij_new[ii,1:nobs_new[ii]]<-
        (rrij_new[ii,1:nobs_new[ii]]-mean_rr_est[ttij_int_new[ii,1:nobs_new[ii]]])/
        sd_rr_est[ttij_int_new[ii,1:nobs_new[ii]]]
    }
    
    chartij_new<-chart_risk(eeij_new,ttij_int_new,nobs_new,parameter,
                            delta_bar,CL)
    
    result<-list(chart=chartij_new,
                 ttij_int=ttij_int_new,
                 standardized_values=eeij_new)
    
  }
}







