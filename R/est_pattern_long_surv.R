##' Estimate the Pattern of Longitudinal and Survival Data
##'
##' @title Estimate the Pattern of Longitudinal and Survival Data
##' @param yyijk a 3d-array of longitudinal observations. \cr
##' \code{yyijk[i,j,k]} is the jth observation of the kth dimension of the ith subject.
##' @param ttij a matrix of observation times. \cr
##' \code{ttij[i,j]} is the jth observation time of the ith subject. \cr
##' \code{yyijk[i,j,]} is observed at \code{ttij[i,j]}.
##' @param nobs an integer vector for number of observations. \cr
##' \code{nobs[i]} is the number of observations for the ith subject.
##' @param starttime a vector of survival times \cr
##' \code{starttime[i]} is the survival time of the ith subject.
##' @param survtime a vector of survival times \cr
##' \code{survtime[i]} is the survival time of the ith subject.
##' @param survevent a logical vector of survival events \cr
##' If \code{survevents[i]==TRUE}, then a survival event is observed at \code{survtime[i]}. \cr
##' If \code{survevents[i]==FALSE}, then no survival event is observed at \code{survtime[i]}.
##' @param ttmin,ttmax two numeric values. \cr
##' [\code{ttmin},\code{ttmax}] is the design time interval.
##' @param ntimepoints an integer value. \cr
##' In the analysis, the design interval and all observed times 
##' are discretized into multiples of basic time units.
##' \code{ntimepoints} is the number of basic time units in the design time interval. \cr
##' @param method a string. \cr
##' \code{method} is the number of basic timepoints in the design time interval. \cr
##' If \code{method="risk"}, apply the risk monitoring method (c.f., You and Qiu 2020). \cr
##' (Currently only the method "risk" is available.)
##' @param smoothing a string. \cr
##' If \code{smoothing="local constant"}, apply local constant approximation. \cr
##' If \code{smoothing="local linear"}, apply local linear approximation. \cr
##' @param hh_beta an integer value. \cr
##' The bandwidth parameter for estimating the regression coefficients beta in the Cox model.
##' @param hh_mean an integer value. \cr
##' The bandwidth parameter for estimating mean function.
##' @param hh_var an integer value. \cr
##' The bandwidth parameter for estimating variance function.
##' @param hh_cov an integer value. \cr
##' The bandwidth parameter for estimating covariance function.
##' @param hh_t an integer value. \cr
##' The bandwidth parameter in time axis for estimating distribution function.
##' @param hh_y a numeric value. \cr
##' The bandwidth parameter in y-axis for estimating distribution function.
est_pattern_long_surv<-function(
  yyijk,ttij,nobs,starttime,survtime,survevent,
  ttmin,ttmax,ntimepoints,
  method="risk",smoothing="local linear",
  hh_beta,hh_mean,hh_var,hh)
{
  if(any(dim(yyijk)[1:2]!=dim(ttij)))stop("Dimensions of 'yyijk' and 'ttij' don't match.")
  if(dim(yyijk)[1]!=length(nobs))stop("Dimensions of 'yyijk' and 'nobs' don't match.")
  
  nind<-dim(yyijk)[1]
  nmaxobs<-dim(yyijk)[2]
  ndim<-dim(yyijk)[3]
  
  if(length(starttime)!=nind)stop("Lengths of 'starttime' and 'nobs' don't match.")
  if(length(survtime)!=nind)stop("Lengths of 'survtime' and 'nobs' don't match.")
  if(length(survevent)!=nind)stop("Lengths of 'survevent' and 'nobs' don't match.")
  
  ttij<-clean_matij_by_nobs(ttij,nobs,"ttij")
  for(kk in 1:ndim){
    yyijk[,,kk]<-clean_matij_by_nobs(yyijk[,,kk],nobs,"yyijk")
  }
  
  if(missing(ttmin))ttmin<-min(ttij,na.rm=TRUE)
  if(missing(ttmax))ttmax<-max(ttij,na.rm=TRUE)
  if(ttmin>min(ttij,na.rm=TRUE))stop("Value of 'ttmin' is too large.")
  if(ttmax<min(ttij,na.rm=TRUE))stop("Value of 'ttmax' is too small.")
  omega<-(ttmax-ttmin)/(ntimepoints-1)
  ttij_int<-round((ttij-ttmin)/omega)+1
  starttime_int<-round((starttime-ttmin)/omega)+1
  survtime_int<-round((survtime-ttmin)/omega)+1
  
  if(is.numeric(survevent))survevent<-ifelse(survevent==1,TRUE,FALSE)
  if(!is.logical(survevent))stop("Values of 'survevent' should be logical or 0/1.")
  
  sumt<-0
  sumn<-0
  for(ii in 1:nind){
    sumt<-sumt+(ttij_int[ii,nobs[ii]]-ttij_int[ii,1])
    sumn<-sumn+(nobs[ii]-1)
  }
  delta_bar<-sumt/sumn
  
  if(method=="risk"){
    if(missing(hh_beta)|missing(hh_mean)|missing(hh_var))
      stop("Method 'risk' requires arguments 'hh_beta', 'hh_mean', 'hh_var'.")
  }else{
    stop("Error in the argument 'method'.")
  }
  
  if(method=="risk"){
    yfij<-matrix(NA,nind,ndim)
    for(ii in 1:nind){
      tt<-survtime_int[ii]
      idx<-which(ttij_int[ii,]<=tt)
      if(length(idx)>0){
        yfij[ii,]<-yyijk[ii,length(idx),]
      }else{
        yfij[ii,]<-yyijk[ii,1,]
      }
    }
    
    beta_est<-risk_estimate_beta(
      yyijk,ttij_int,nobs,
      starttime_int,survtime_int,survevent,
      yfij,hh_beta)
    rrij<-matrix(NA,nind,nmaxobs)
    for(ii in 1:nind){
      rrij[ii,1:nobs[ii]]<-yyijk[ii,1:nobs[ii],]%*%beta_est
    }
    if(smoothing=="local linear"){
      mean_rr_est<-local_linear_mean_est_faster(rrij,ttij_int,nobs,1:ntimepoints,hh_beta)
    }else if(smoothing=="local constant"){
      mean_rr_est<-local_const_mean_est_faster(rrij,ttij_int,nobs,1:ntimepoints,hh_beta)
    }
    
    epsij<-matrix(NA,nind,nmaxobs)
    eps2ij<-matrix(NA,nind,nmaxobs)
    for(ii in 1:nind){
      epsij[ii,1:nobs[ii]]<-rrij[ii,1:nobs[ii]]-mean_rr_est[ttij_int[ii,1:nobs[ii]]]
      eps2ij[ii,1:nobs[ii]]<-epsij[ii,1:nobs[ii]]^2
    }
    if(smoothing=="local linear"){
      var_rr_est<-local_linear_mean_est_faster(eps2ij,ttij_int,nobs,1:ntimepoints,hh_var)
    }else if(smoothing=="local constant"){
      var_rr_est<-local_const_mean_est_faster(eps2ij,ttij_int,nobs,1:ntimepoints,hh_var)
    }
    
    result<-list(
      grid=ttmin+omega*(0:(ntimepoints-1)),
      beta_est=beta_est,
      mean_rr_est=mean_rr_est,
      var_rr_est=var_rr_est,
      yyijk=yyijk,
      ttij_int=ttij_int,
      starttime_int=starttime_int,
      survtime_int=survtime_int,
      survevent=survevent,
      nobs=nobs,
      delta_bar=delta_bar,
      ttmin=ttmin,
      ttmax=ttmax,
      omega=omega,
      ntimepoints=ntimepoints,
      smoothing=smoothing,
      method=method,
      hh_beta=hh_beta,
      hh_mean=hh_mean,
      hh_var=hh_var)
    return(result)
  }
}

