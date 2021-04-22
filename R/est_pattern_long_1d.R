##' Estimate the Pattern of Univariate Longitudinal Data
##'
##' @title Estimate the Pattern of Univariate Longitudinal Data
##' @param yyij a matrix of longitudinal observations. \cr
##' \code{yyij[i,j]} is the jth observation of the kth dimension of the ith subject.
##' @param ttij a matrix of observation times. \cr
##' \code{ttij[i,j]} is the jth observation time of the ith subject. \cr
##' \code{yyij[i,j]} is observed at \code{ttij[i,j]}.
##' @param nobs an integer vector for number of observations. \cr
##' \code{nobs[i]} is the number of observations for the ith subject.
##' @param ttmin,ttmax two numeric values. \cr
##' [\code{ttmin},\code{ttmax}] is the design time interval.
##' @param ntimepoints an integer value. \cr
##' In the analysis, the design interval and all observed times 
##' are descretized into multiples of basic time units.
##' \code{ntimepoints} is the number of basic time units in the design time interval. \cr
##' @param method a string. \cr
##' If \code{method="meanvar"}, the function will estimate 
##' the mean and variance functions using local smoothing
##' (c.f., Qiu and Xiang, 2014). 
##' Parameters \code{hh_mean} and \code{hh_var} are needed. \cr
##' If \code{method="meanvarcov"}, the function will estimate 
##' the mean, variance and covariance functions using local smoothing 
##' (c.f., Li and Qiu, 2016). 
##' Parameters \code{hh_mean}, \code{hh_var} and \code{hh_cov} are needed. \cr
##' If \code{method="meanvarcovmean"}, the function will estimate
##' the mean, variance and covariance functions (c.f., Li and Qiu, 2016).
##' In the last step, the mean function will be updated using the
##' covariance function. 
##' Parameters \code{hh_mean}, \code{hh_var} and \code{hh_cov} are needed. \cr
##' If \code{method="distribution"}, the function will estimate the 
##' distribution function (c.f., You and Qiu, 2020).
##' Parameters \code{hh_t} and \code{hh_y} are needed. \cr
##' If \code{method="distributionvarcov"}, the function will estimate the 
##' distribution function and the covariance function of standardized values
##' (c.f., You and Qiu 2020).
##' Parameters \code{hh_cov}, \code{hh_t} and \code{hh_y} are needed.
##' @param smoothing a string. \cr
##' If \code{smoothing="local constant"}, apply local constant approximation. \cr
##' If \code{smoothing="local linear"}, apply local linear approximation.
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
est_pattern_long_1d<-
  function(yyij,ttij,nobs,ttmin,ttmax,ntimepoints,
           method,smoothing="local linear",
           hh_mean,hh_var,hh_cov,hh_t,hh_y){
    
    if(any(dim(yyij)!=dim(ttij)))stop("Dimensions of 'yyij' and 'ttij' don't match.")
    if(dim(yyij)[1]!=length(nobs))stop("Dimensions of 'yyij' and 'nobs' don't match.")
    
    nind<-dim(yyij)[1]
    nmaxobs<-dim(yyij)[2]
    yyij<-clean_matij_by_nobs(yyij,nobs,"yyij")
    ttij<-clean_matij_by_nobs(ttij,nobs,"ttij")
    
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
    }else if(method=="meanvarcovmean"){
      if(missing(hh_mean)|missing(hh_var)|missing(hh_cov))
        stop("Method 'meanvarcovmean' requires arguments hh_mean, hh_var, hh_cov.")
    }else if(method=="distribution"){
      if(missing(hh_t)|missing(hh_y))
        stop("Method 'distribution' requires arguments hh_t, hh_y.")
    }else if(method=="disttributionvarcov"){
      if(missing(hh_t)|missing(hh_y)|missing(hh_cov))
        stop("Method 'disttributionvarcov' requires arguments hh_t, hh_y, hh_cov.")
    }else{
      stop("Error in the argument 'method'.")
    }
    
    if(!smoothing%in%c("local constant","local linear"))
      stop("Error in the argument 'smoothing'.")
    
    # ttij<-ttij_IC
    # yyij<-yyij_IC
    # nobs<-nobs_IC
    # omega<-4
    
    if(method%in%c("meanvar","meanvarcov","meanvarcovmean")){
      
      if(smoothing=="local linear"){
        mean_est<-local_linear_mean_est_faster(yyij,ttij_int,nobs,1:ntimepoints,hh_mean)
      }else if(smoothing=="local constant"){
        mean_est<-local_const_mean_est_faster(yyij,ttij_int,nobs,1:ntimepoints,hh_mean)
      }
      
      epsij<-matrix(NA,nind,nmaxobs)
      eps2ij<-matrix(NA,nind,nmaxobs)
      for(ii in 1:nind){
        epsij[ii,1:nobs[ii]]<-yyij[ii,1:nobs[ii]]-mean_est[ttij_int[ii,1:nobs[ii]]]
        eps2ij[ii,1:nobs[ii]]<-epsij[ii,1:nobs[ii]]^2
      }
      
      if(smoothing=="local linear"){
        var_est<-local_linear_mean_est_faster(eps2ij,ttij_int,nobs,1:ntimepoints,hh_var)
      }else if(smoothing=="local constant"){
        var_est<-local_const_mean_est_faster(eps2ij,ttij_int,nobs,1:ntimepoints,hh_var)
      }
      
      if(method=="meanvar"){
        pattern<-list(
          grid=ttmin+omega*(0:(ntimepoints-1)),
          mean_est=mean_est,
          var_est=var_est,
          yyij=yyij,
          ttij_int=ttij_int,
          nobs=nobs,
          ttmin=ttmin,
          ttmax=ttmax,
          omega=omega,
          ntimepoints=ntimepoints,
          smoothing=smoothing,
          method=method,
          hh_mean=hh_mean,
          hh_var=hh_var)
        class(pattern)<-c(class(pattern),"patn_long_univ_meanvar")
        return(pattern)
      }
      
      if(smoothing=="local linear"){
        cov_est<-local_linear_cov_est_faster(epsij,ttij_int,nobs,1:ntimepoints,hh_cov)
      }else if(smoothing=="local constant"){
        cov_est<-local_const_cov_est_faster(epsij,ttij_int,nobs,1:ntimepoints,hh_cov)
      }
      diag(cov_est)<-var_est
      
      if(method=="meanvarcov"){
        pattern<-list(
          grid=ttmin+omega*(0:(ntimepoints-1)),
          mean_est=mean_est,
          var_est=var_est,
          cov_est=cov_est,
          yyij=yyij,
          ttij_int=ttij_int,
          nobs=nobs,
          ttmin=ttmin,
          ttmax=ttmax,
          omega=omega,
          ntimepoints=ntimepoints,
          smoothing=smoothing,
          method=method,
          hh_mean=hh_mean,
          hh_var=hh_var,
          hh_cov=hh_cov)
        class(pattern)<-c(class(pattern),"patn_long_univ_meanvarcov")
        return(pattern)
      }
      
      if(smoothing=="local constant"){
        mean2_est<-local_const_mean_est_update_faster(yyij,ttij_int,nobs,1:ntimepoints,hh_mean,cov_est)
        mean2_est<-c(mean2_est)
      }else if(smoothing=="local linear"){
        mean2_est<-local_linear_mean_est_update_faster(yyij,ttij_int,nobs,1:ntimepoints,hh_mean,cov_est)
        mean2_est<-c(mean2_est)
      }
      
      if(method=="meanvarcovmean"){
        pattern<-list(
          grid=ttmin+omega*(0:(ntimepoints-1)),
          mean_est=mean2_est,
          var_est=var_est,
          cov_est=cov_est,
          mean0_est=mean_est,
          yyij=yyij,
          ttij_int=ttij_int,
          nobs=nobs,
          ttmin=ttmin,
          ttmax=ttmax,
          omega=omega,
          ntimepoints=ntimepoints,
          smoothing=smoothing,
          method=method,
          hh_mean=hh_mean,
          hh_var=hh_var,
          hh_cov=hh_cov)
        class(pattern)<-c(class(pattern),"patn_long_univ_meanvarcov")
        return(pattern)
      }
    }
    
    if(method%in%c("distribution","distributionvarcov")){
      yy_raw<-c()
      tt_raw<-c()
      for(ii in 1:nind){
        yy_raw<-c(yy_raw,yyij[ii,1:nobs[ii]])
        tt_raw<-c(tt_raw,ttij_int[ii,1:nobs[ii]])
      }
      tt_order<-order(tt_raw)
      tt_ref<-tt_raw[tt_order]
      yy_ref<-yy_raw[tt_order]
      rm(yy_raw)
      rm(tt_raw)
      
      starting_idx<-rep(NA,ntimepoints)
      ending_idx<-rep(NA,ntimepoints)
      upper_line<-rep(NA,ntimepoints)
      for(tt in 1:ntimepoints){
        starting_idx[tt]<-min(which(abs(tt_ref-tt)<hh_t))
        ending_idx[tt]<-max(which(abs(tt_ref-tt)<hh_t))
        upper_line[tt]<-quantile(yy_ref[starting_idx[tt]:ending_idx[tt]],probs=0.8,na.rm=TRUE)
      }
      
      zzij<-local_const_percentile_est_faster(
        yyij,ttij_int,nobs,yy_ref,tt_ref,
        starting_idx,ending_idx,upper_line,ntimepoints,hh_t,hh_y)
      
      if(method=="distribution"){
        pattern<-list(grid=ttmin+omega*(0:(ntimepoints-1)),
                      yy_ref=yy_ref,
                      tt_ref=tt_ref,
                      starting_idx=starting_idx,
                      ending_idx=ending_idx,
                      upper_line=upper_line,
                      yyij=yyij,
                      zzij=zzij,
                      ttij_int=ttij_int,
                      nobs=nobs,
                      ttmin=ttmin,
                      ttmax=ttmax,
                      omega=omega,
                      ntimepoints=ntimepoints,
                      smoothing=smoothing,
                      method=method,
                      hh_t=hh_t,
                      hh_y=hh_y)
        class(pattern)<-c(class(pattern),"patn_long_univ_distribution")
        return(pattern)
      }
      
      cov_est<-local_const_cov_est_faster(zzij,ttij_int,nobs,1:ntimepoints,hh_cov)
      diag(cov_est)<-1.0
      
      if(method=="distributionvarcov"){
        pattern<-list(grid=ttmin+omega*(0:(ntimepoints-1)),
                      yy_ref=yy_ref,
                      tt_ref=tt_ref,
                      starting_idx=starting_idx,
                      ending_idx=ending_idx,
                      upper_line=upper_line,
                      cov_est=cov_est,
                      yyij=yyij,
                      zzij=zzij,
                      ttij_int=ttij_int,
                      nobs=nobs,
                      ttmin=ttmin,
                      ttmax=ttmax,
                      omega=omega,
                      ntimepoints=ntimepoints,
                      smoothing=smoothing,
                      method=method,
                      hh_t=hh_t,
                      hh_y=hh_y,
                      hh_cov=hh_cov)
        class(pattern)<-c(class(pattern),"patn_long_univ_distribution")
        return(pattern)
      }
    }
  }
