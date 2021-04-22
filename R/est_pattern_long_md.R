##' Estimate the Pattern of Multivariate Longitudinal Data
##'
##' @title Estimate the Pattern of Multivariate Longitudinal Data
##' @param yyijk a 3d-array of longitudinal observations. \cr
##' \code{yyijk[i,j,k]} is the jth observation of the kth dimension of the ith subject.
##' @param ttij a matrix of observation times. \cr
##' \code{ttij[i,j]} is the jth observation time of the ith subject. \cr
##' \code{yyijk[i,j,]} is observed at \code{ttij[i,j]}.
##' @param nobs an integer vector for number of observations. \cr
##' \code{nobs[i]} is the number of observations for the ith subject.
##' @param ttmin,ttmax two numeric values. \cr
##' [\code{ttmin},\code{ttmax}] is the design time interval.
##' @param ntimepoints an integer value. \cr
##' In the analysis, the design interval and all observed times 
##' are descretized into multiples of basic time units.
##' \code{ntimepoints} is the number of basic time units in the design time interval. \cr
##' @param method a string. \cr
##' If \code{method="meanvar"}, the function will estimate the 
##' mean function (\eqn{\mathrm{E}[\mathbf{y}(t)]}), and 
##' variance function (\eqn{\mathrm{Var}(\mathbf{y}(t))}). 
##' Parameters \code{hh_mean} and \code{hh_var} are needed. \cr
##' If \code{method="meanvarcov"}, the function will estimate the 
##' mean function (\eqn{\mathrm{E}[\mathbf{y}(t)]}), 
##' variance function (\eqn{\mathrm{Var}(\mathbf{y}(t))}), and 
##' covariance function (\eqn{\mathrm{Cov}(\mathbf{y}(s),\mathbf{y}(t))}). 
##' Parameters \code{hh_mean}, \code{hh_var} and \code{hh_cov}.
##' @param hh_mean an integer value. \cr
##' The bandwidth parameter for estimating mean function.
##' @param hh_var an integer value. \cr
##' The bandwidth parameter for estimating variance function.
##' @param hh_cov an integer value. \cr
##' The bandwidth parameter for estimating covariance function.
##' @return the estimated longitudinal pattern. \cr
##' If \code{method="meanvar"}, return an object of class \code{patn_long_mult_meanvar}. \cr
##' If \code{method="meanvarcov"}, return an object of class \code{patn_long_mult_meanvarcov}.
##' @author LY
##' @export
##' @examples 
est_pattern_long_md<-
  function(yyijk,ttij,nobs,ttmin,ttmax,ntimepoints,
           method,hh_mean,hh_var,hh_cov){
    
    if(any(dim(yyijk)[1:2]!=dim(ttij)))stop("Dimensions of 'yyijk' and 'ttij' don't match.")
    if(dim(yyijk)[1]!=length(nobs))stop("Dimensions of 'yyijk' and 'nobs' don't match.")
    
    nind<-dim(yyijk)[1]
    nmaxobs<-dim(yyijk)[2]
    ndim<-dim(yyijk)[3]
    
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
    
    if(method=="meanvar"){
      if(missing(hh_mean)|missing(hh_var))
        stop("Method 'meanvar' requires arguments hh_mean, hh_var.")
    }else if(method=="meanvarcov"){
      if(missing(hh_mean)|missing(hh_var)|missing(hh_cov))
        stop("Method 'meanvarcov' requires arguments hh_mean, hh_var, hh_cov.")
    }else{
      stop("Error in the argument 'method'.")
    }
    
    mean_est<-matrix(0,ntimepoints,ndim)
    mean_est<-f90_local_const_mean_est_mult_wrap(
      yyijk,ttij_int,nobs,nind,nmaxobs,ndim,ntimepoints,hh_mean,mean_est)
    
    epsijk<-array(NA,c(nind,nmaxobs,ndim))
    for(ii in 1:nind){
      epsijk[ii,1:nobs[ii],]<-yyijk[ii,1:nobs[ii],]-mean_est[ttij_int[ii,1:nobs[ii]],]
    }
    
    var_est<-array(0,c(ntimepoints,ndim,ndim))
    var_est<-f90_local_const_var_est_mult_wrap(
      epsijk,ttij_int,nobs,nind,nmaxobs,ndim,ntimepoints,hh_var,var_est)
    
    if(method=="meanvar"){
      pattern<-list(
        grid=ttmin+omega*(0:(ntimepoints-1)),
        mean_est=mean_est,
        var_est=var_est,
        yyijk=yyijk,
        ttij_int=ttij_int,
        nobs=nobs,
        ttmin=ttmin,
        ttmax=ttmax,
        omega=omega,
        ntimepoints=ntimepoints,
        method=method)
      class(pattern)<-c(class(pattern),"patn_long_mult_meanvar")
      return(pattern)
    }
    
    cov_est<-array(0,c(ntimepoints,ntimepoints,ndim,ndim))
    cov_est<-f90_local_const_cov_est_mult_wrap(
      epsijk,ttij_int,nobs,nind,nmaxobs,ndim,ntimepoints,hh_cov,cov_est)
    for(kk in 1:ndim){
      diag(cov_est[,,kk,kk])<-var_est[,kk,kk]
    }
    
    if(method=="meanvarcov"){
      pattern<-list(
        grid=ttmin+omega*(0:(ntimepoints-1)),
        mean_est=mean_est,
        var_est=var_est,
        cov_est=cov_est,
        yyijk=yyijk,
        ttij_int=ttij_int,
        nobs=nobs,
        ttmin=ttmin,
        ttmax=ttmax,
        omega=omega,
        ntimepoints=ntimepoints,
        method=method)
      class(pattern)<-c(class(pattern),"patn_long_mult_meanvarcov")
      return(pattern)
    }
  }
