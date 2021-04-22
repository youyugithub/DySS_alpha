##' Monitor Multivariate Longitudinal Data
##'
##' @title Monitor Multivariate Longitudinal Data
##' @param yyijk_new an array of longitudinal observations. \cr
##' \code{yyijk_new[i,j,k]} is the jth observation of the kth dimension of the ith subject.
##' @param ttij_new a matrix of observation times. \cr
##' \code{ttij_new[i,j]} is the jth observation time of the ith subject. \cr
##' \code{yyijk_new[i,j,]} is observed at \code{ttij[i,j]}.
##' @param nobs_new an integer vector for number of observations. \cr
##' \code{nobs_new[i]} is the number of observations for the ith subject.
##' @param pattern the estimated regular longitudinal pattern
##' @param side a string \cr
##' If \code{side="upward"}, control charts aim to detect upward shifts.\cr
##' If \code{side="downward"}, control charts aim to detect downward shifts.\cr
##' If \code{side="both"}, control charts aim to detect shifts in both sides.
##' @param method a string \cr
##' If \code{method="simultaneous CUSUM"}, 
##' apply simultaneous CUSUM charts. (See SIMUL in You et al, 2020.)\cr
##' If \code{method="simultaneous EWMA"}
##' apply simultaneous EWMA charts. (See SIMUL in You et al, 2020.)\cr
##' If \code{method="multivariate CUSUM"}
##' apply multivariate CUSUM charts.\cr
##' If \code{method="multivariate EWMA"}
##' apply multivariate EWMA charts. (See Qiu and Xiang, 2015 or QX-1S/QS-2S in You et al, 2020.)\cr
##' If \code{method="decorrelation CUSUM"}
##' apply decorrelation CUSUM charts. (See Li and Qiu, 2017 or LQ-1S/LQ-2S in You et al, 2020) \cr
##' If \code{method="decorrelation EWMA"}
##' apply decorrelation EWMA charts. (See Li and Qiu, 2017 or LQ-1S/LQ-2S in You et al, 2020) \cr
##' If \code{method="nonparametric CUSUM"} \cr
##' If \code{method="nonparametric EWMA"} \cr
##' @param parameter a numeric value. \cr
##' \code{parameter} is the allowance constant if \code{method} is a CUSUM chart. \cr
##' \code{parameter} is the weighting parameter if \code{method} is an EWMA chart.
##' @param CL a numeric value \cr
##' \code{CL} is the control limit. 
##' A signal will be given if charting statistics are larger than the control limit. 
##' (Note: in this package, signs of charting statistics may be reversed such that 
##' larger values of charting statistics indicate worse performance of processes.)
##' After the signal is given, the algorithm stops calculating the charting statistics for the remaining observation times.
##' The default value of control limit is infinity, which means we will calculate the charting statistics for all observation times.
##' @return a list \cr
##' \code{$chart} a numeric matrix, \code{$chart[i,j]} is the jth charting statistic of the ith subject calculated at time \code{ttij_new[i,j]}. \cr
##' \code{$ttij_int} an integer matrix, the observation times rounded to the multiples of basic time unit. \cr
##' (Supplementary) \code{$SSijk}, \code{$SSijk_upward}, \code{$SSijk_downward} numeric arrays, 
##' the multivariate statistics used in the calculation of (downward/upward/two-sided) control charts.
##' \code{$SSijk[i,j,]} is the multivariate statistic at time \code{ttij_new[i,j]}.\cr
##' (Supplementary) \code{$standardized_values}, \code{$standardized_values_upward}, \code{$standardized_values_downward} numeric matrices, 
##' the standardized values used in the calculation of (downward/upward/two-sided) control charts.
mnt_long_md<-function(
  yyijk_new,ttij_new,nobs_new,pattern,
  side="both",method="multivariate EWMA",
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
  
  if(!side%in%c("upward","downward","both"))stop("Error in argument 'side'.")
  
  all_methods<-c("simultaneous CUSUM",
                 "simultaneous EWMA",
                 "multivariate CUSUM",
                 "multivariate EWMA",
                 "decorrelation CUSUM",
                 "decorrelation EWMA",
                 "nonparametric CUSUM",
                 "nonparametric EWMA")
  if(!method%in%all_methods)stop("Error in argument 'method'.")
  
  ######################
  # simultaneous CUSUM #
  ######################
  
  if(method=="simultaneous CUSUM"){
    if(!inherits(pattern,c("patn_long_mult_meanvar","patn_long_mult_meanvarcov")))
      stop("Pattern does not match method.")
    mean_est<-pattern$mean_est
    var_est<-pattern$var_est
    sd_est<-matrix(NA,ntimepoints,ndim)
    for(kk in 1:ndim)sd_est[,kk]<-sqrt(var_est[,kk,kk])
    eeijk_new<-array(NA,c(nind,nmaxobs,ndim))
    for(ii in 1:nind){
      eeijk_new[ii,,]<-(yyijk_new[ii,,]-mean_est[ttij_int_new[ii,],])/
        sd_est[ttij_int_new[ii,],]
    }
    if(side=="both"){
      CCij_new<-matrix(0.0,nind,nmaxobs)
      SSijk_upward_new<-array(0.0,c(nind,nmaxobs,ndim))
      SSijk_downward_new<-array(0.0,c(nind,nmaxobs,ndim))
      chart_output<-f90_mchart_simultaneous_CUSUM_both_wrap(
        eeijk_new,nobs_new,nind,nmaxobs,ndim,parameter,CL,
        CCij_new,SSijk_upward_new,SSijk_downward_new)
      chartij_new<-chart_output[[1]]
      SSijk_upward_new<-chart_output[[2]]
      SSijk_downward_new<-chart_output[[3]]
      result<-list(chart=chartij_new,
                   ttij_int=ttij_int_new,
                   SSijk_upward=SSijk_upward_new,
                   SSijk_downward=SSijk_downward_new,
                   standardized_values=eeijk_new)
      return(result)
    }else if(side=="upward"){
      CCij_new<-matrix(0.0,nind,nmaxobs)
      SSijk_upward_new<-array(0.0,c(nind,nmaxobs,ndim))
      chart_output<-f90_mchart_simultaneous_CUSUM_upward_wrap(
        eeijk_new,nobs_new,nind,nmaxobs,ndim,parameter,CL,
        CCij_new,SSijk_upward_new)
      chartij_new<-chart_output[[1]]
      SSijk_upward_new<-chart_output[[2]]
      result<-list(chart=chartij_new,
                   ttij_int=ttij_int_new,
                   SSijk_upward=SSijk_upward_new,
                   standardized_values=eeijk_new)
      return(result)
    }else if(side=="downward"){
      CCij_new<-matrix(0.0,nind,nmaxobs)
      SSijk_upward_new<-array(0.0,c(nind,nmaxobs,ndim))
      chart_output<-f90_mchart_simultaneous_CUSUM_upward_wrap(
        -eeijk_new,nobs_new,nind,nmaxobs,ndim,parameter,CL,
        CCij_new,SSijk_upward_new)
      chartij_new<-chart_output[[1]]
      SSijk_upward_new<-chart_output[[2]]
      result<-list(chart=chartij_new,
                   ttij_int=ttij_int_new,
                   SSijk_upward=SSijk_upward_new,
                   standardized_values=eeijk_new)
      return(result)
    }
  }
  
  #####################
  # simultaneous EWMA #
  #####################
  
  if(method=="simultaneous EWMA"){
    if(!inherits(pattern,c("patn_long_mult_meanvar","patn_long_mult_meanvarcov")))
      stop("Pattern does not match method.")
    mean_est<-pattern$mean_est
    var_est<-pattern$var_est
    sd_est<-matrix(NA,ntimepoints,ndim)
    for(kk in 1:ndim)sd_est[,kk]<-sqrt(var_est[,kk,kk])
    eeijk_new<-array(NA,c(nind,nmaxobs,ndim))
    for(ii in 1:nind){
      eeijk_new[ii,,]<-(yyijk_new[ii,,]-mean_est[ttij_int_new[ii,],])/
        sd_est[ttij_int_new[ii,],]
    }
    if(side=="both"){
      CCij_new<-matrix(0.0,nind,nmaxobs)
      SSijk_new<-array(0.0,c(nind,nmaxobs,ndim))
      chart_output<-f90_mchart_simultaneous_EWMA_both_wrap(
        eeijk_new,nobs_new,nind,nmaxobs,ndim,parameter,CL,
        CCij_new,SSijk_new)
      chartij_new<-chart_output[[1]]
      SSijk<-chart_output[[2]]
      result<-list(chart=chartij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new,
                   standardized_values=eeijk_new)
      return(result)
    }else if(side=="upward"){
      CCij_new<-matrix(0.0,nind,nmaxobs)
      SSijk_new<-array(0.0,c(nind,nmaxobs,ndim))
      chart_output<-f90_mchart_simultaneous_EWMA_upward_wrap(
        eeijk_new,nobs_new,nind,nmaxobs,ndim,parameter,CL,
        CCij_new,SSijk_new)
      chartij_new<-chart_output[[1]]
      SSijk<-chart_output[[2]]
      result<-list(chart=chartij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new,
                   standardized_values=eeijk_new)
      return(result)
    }else if(side=="downward"){
      CCij_new<-matrix(0.0,nind,nmaxobs)
      SSijk_new<-array(0.0,c(nind,nmaxobs,ndim))
      chart_output<-f90_mchart_simultaneous_EWMA_upward_wrap(
        -eeijk_new,nobs_new,nind,nmaxobs,ndim,parameter,CL,
        CCij_new,SSijk_new)
      chartij_new<-chart_output[[1]]
      SSijk<-chart_output[[2]]
      result<-list(chart=chartij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new,
                   standardized_values=eeijk_new)
      return(result)
    }
  }
  
  ######################
  # multivariate CUSUM #
  ######################
  
  if(method=="multivariate CUSUM"){
    if(!inherits(pattern,c("patn_long_mult_meanvar","patn_long_mult_meanvarcov")))
      stop("Pattern does not match method.")
    mean_est<-pattern$mean_est
    var_est<-pattern$var_est
    epsijk_new<-array(NA,c(nind,nmaxobs,ndim))
    for(ii in 1:nind){
      epsijk_new[ii,,]<-yyijk_new[ii,,]-mean_est[ttij_int_new[ii,],]
    }
    if(side=="both"){
      CCij_new<-matrix(NA,nind,nmaxobs)
      SSijk_new<-array(NA,c(nind,nmaxobs,ndim))
      for(ii in 1:nind){
        eps_ii<-t(epsijk_new[ii,,])
        var_cube_ii<-var_est[ttij_int_new[ii,],,]
        chart_output<-mchart1_multivariate_CUSUM_multivariate_both(
          eps_ii,nobs_new[ii],var_cube_ii,ndim,parameter,CL)
        CCij_new[ii,1:nobs_new[ii]]<-chart_output[[1]]
        SSijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[2]])
      }
      result<-list(chart=CCij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new)
      return(result)
    }else if(side=="upward"){
      CCij_new<-matrix(NA,nind,nmaxobs)
      SSijk_new<-array(NA,c(nind,nmaxobs,ndim))
      for(ii in 1:nind){
        eps_ii<-t(epsijk_new[ii,,])
        var_cube_ii<-var_est[ttij_int_new[ii,],,]
        chart_output<-mchart1_multivariate_CUSUM_multivariate_upward(
          eps_ii,nobs_new[ii],var_cube_ii,ndim,parameter,CL)
        CCij_new[ii,1:nobs_new[ii]]<-chart_output[[1]]
        SSijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[2]])
      }
      result<-list(chart=CCij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new)
      return(result)
    }else if(side=="downward"){
      CCij_new<-matrix(NA,nind,nmaxobs)
      SSijk_new<-array(NA,c(nind,nmaxobs,ndim))
      for(ii in 1:nind){
        eps_ii<-t(epsijk_new[ii,,])
        var_cube_ii<-var_est[ttij_int_new[ii,],,]
        chart_output<-mchart1_multivariate_CUSUM_multivariate_upward(
          -eps_ii,nobs_new[ii],var_cube_ii,ndim,parameter,CL)
        CCij_new[ii,1:nobs_new[ii]]<-chart_output[[1]]
        SSijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[2]])
      }
      result<-list(chart=CCij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new)
      return(result)
    }
  }
  
  #####################
  # multivariate EWMA #
  #####################
  
  if(method=="multivariate EWMA"){
    if(!inherits(pattern,c("patn_long_mult_meanvar","patn_long_mult_meanvarcov")))
      stop("Pattern does not match method.")
    mean_est<-pattern$mean_est
    var_est<-pattern$var_est
    epsijk_new<-array(NA,c(nind,nmaxobs,ndim))
    for(ii in 1:nind){
      epsijk_new[ii,,]<-yyijk_new[ii,,]-mean_est[ttij_int_new[ii,],]
    }
    if(side=="both"){
      CCij_new<-matrix(NA,nind,nmaxobs)
      SSijk_new<-array(NA,c(nind,nmaxobs,ndim))
      eeijk_new<-array(NA,c(nind,nmaxobs,ndim))
      for(ii in 1:nind){
        eps_ii<-t(epsijk_new[ii,,])
        var_cube_ii<-var_est[ttij_int_new[ii,],,]
        chart_output<-mchart1_multivariate_EWMA_multivariate_both(
          eps_ii,nobs_new[ii],var_cube_ii,ndim,parameter,CL)
        CCij_new[ii,1:nobs_new[ii]]<-chart_output[[1]]
        SSijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[2]])
        eeijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[3]])
      }
      result<-list(chart=CCij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new,
                   standardized_values=eeijk_new)
      return(result)
    }else if(side=="upward"){
      CCij_new<-matrix(NA,nind,nmaxobs)
      SSijk_new<-array(NA,c(nind,nmaxobs,ndim))
      eeijk_new<-array(NA,c(nind,nmaxobs,ndim))
      for(ii in 1:nind){
        eps_ii<-t(epsijk_new[ii,,])
        var_cube_ii<-var_est[ttij_int_new[ii,],,]
        chart_output<-mchart1_multivariate_EWMA_multivariate_upward(
          eps_ii,nobs_new[ii],var_cube_ii,ndim,parameter,CL)
        CCij_new[ii,1:nobs_new[ii]]<-chart_output[[1]]
        SSijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[2]])
        eeijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[3]])
      }
      result<-list(chart=CCij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new,
                   standardized_values=eeijk_new)
      return(result)
    }else if(side=="downward"){
      CCij_new<-matrix(NA,nind,nmaxobs)
      SSijk_new<-array(NA,c(nind,nmaxobs,ndim))
      eeijk_new<-array(NA,c(nind,nmaxobs,ndim))
      for(ii in 1:nind){
        eps_ii<-t(epsijk_new[ii,,])
        var_cube_ii<-var_est[ttij_int_new[ii,],,]
        chart_output<-mchart1_multivariate_CUSUM_multivariate_upward(
          -eps_ii,nobs_new[ii],var_cube_ii,ndim,parameter,CL)
        CCij_new[ii,1:nobs_new[ii]]<-chart_output[[1]]
        SSijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[2]])
        eeijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[3]])
      }
      result<-list(chart=CCij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new,
                   standardized_values=eeijk_new)
      return(result)
    }
  }
  
  #######################
  # decorrelation CUSUM #
  #######################
  
  if(method=="decorrelation CUSUM"){
    if(!inherits(pattern,c("patn_long_mult_meanvarcov")))
      stop("Pattern does not match method.")
    mean_est<-pattern$mean_est
    var_est<-pattern$var_est
    cov_est<-pattern$cov_est
    cov_mat_est<-matrix(aperm(cov_est,c(3,1,4,2)),ntimepoints*ndim,ntimepoints*ndim)
    epsijk_new<-array(NA,c(nind,nmaxobs,ndim))
    for(ii in 1:nind){
      epsijk_new[ii,,]<-yyijk_new[ii,,]-mean_est[ttij_int_new[ii,],]
    }
    if(side=="both"){
      CCij_new<-matrix(NA,nind,nmaxobs)
      SSijk_new<-array(NA,c(nind,nmaxobs,ndim))
      eeijk_new<-array(NA,c(nind,nmaxobs,ndim))
      for(ii in 1:nind){
        idx_ii<-rep(ttij_int_new[ii,]-1,each=ndim)*ndim+1:ndim
        cov_ii<-cov_mat_est[idx_ii,idx_ii]
        eps_ii<-c(t(epsijk_new[ii,1:nobs_new[ii],]))
        chart_output<-mchart1_decorrelation_CUSUM_multivariate_both(
          eps_ii,cov_ii,nobs_new[ii],ndim,parameter,CL)
        CCij_new[ii,1:nobs_new[ii]]<-chart_output[[1]]
        SSijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[2]])
        eeijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[3]])
      }
      result<-list(chart=CCij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new,
                   standardized_values=eeijk_new)
      return(result)
    }else if(side=="upward"){
      CCij_new<-matrix(NA,nind,nmaxobs)
      SSijk_new<-array(NA,c(nind,nmaxobs,ndim))
      eeijk_new<-array(NA,c(nind,nmaxobs,ndim))
      for(ii in 1:nind){
        idx_ii<-rep(ttij_int_new[ii,]-1,each=ndim)*ndim+1:ndim
        cov_ii<-cov_mat_est[idx_ii,idx_ii]
        eps_ii<-c(t(epsijk_new[ii,1:nobs_new[ii],]))
        chart_output<-mchart1_decorrelation_CUSUM_multivariate_upward(
          eps_ii,cov_ii,nobs_new[ii],ndim,parameter,CL)
        CCij_new[ii,1:nobs_new[ii]]<-chart_output[[1]]
        SSijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[2]])
        eeijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[3]])
      }
      result<-list(chart=CCij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new,
                   standardized_values=eeijk_new)
      return(result)
    }else if(side=="downward"){
      CCij_new<-matrix(NA,nind,nmaxobs)
      SSijk_new<-array(NA,c(nind,nmaxobs,ndim))
      eeijk_new<-array(NA,c(nind,nmaxobs,ndim))
      for(ii in 1:nind){
        idx_ii<-rep(ttij_int_new[ii,]-1,each=ndim)*ndim+1:ndim
        cov_ii<-cov_mat_est[idx_ii,idx_ii]
        eps_ii<-c(t(epsijk_new[ii,1:nobs_new[ii],]))
        chart_output<-mchart1_decorrelation_CUSUM_multivariate_upward(
          -eps_ii,cov_ii,nobs_new[ii],ndim,parameter,CL)
        CCij_new[ii,1:nobs_new[ii]]<-chart_output[[1]]
        SSijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[2]])
        eeijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[3]])
      }
      result<-list(chart=CCij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new,
                   standardized_values=-eeijk_new)
      return(result)
    }
  }
  
  ######################
  # decorrelation EWMA #
  ######################
  
  if(method=="decorrelation EWMA"){
    if(!inherits(pattern,c("patn_long_mult_meanvarcov")))
      stop("Pattern does not match method.")
    mean_est<-pattern$mean_est
    var_est<-pattern$var_est
    cov_est<-pattern$cov_est
    cov_mat_est<-matrix(aperm(cov_est,c(3,1,4,2)),ntimepoints*ndim,ntimepoints*ndim)
    epsijk_new<-array(NA,c(nind,nmaxobs,ndim))
    for(ii in 1:nind){
      epsijk_new[ii,,]<-yyijk_new[ii,,]-mean_est[ttij_int_new[ii,],]
    }
    if(side=="both"){
      CCij_new<-matrix(NA,nind,nmaxobs)
      SSijk_new<-array(NA,c(nind,nmaxobs,ndim))
      eeijk_new<-array(NA,c(nind,nmaxobs,ndim))
      for(ii in 1:nind){
        idx_ii<-rep(ttij_int_new[ii,]-1,each=ndim)*ndim+1:ndim
        cov_ii<-cov_mat_est[idx_ii,idx_ii]
        eps_ii<-c(t(epsijk_new[ii,1:nobs_new[ii],]))
        chart_output<-mchart1_decorrelation_EWMA_multivariate_both(
          eps_ii,cov_ii,nobs_new[ii],ndim,parameter,CL)
        CCij_new[ii,1:nobs_new[ii]]<-chart_output[[1]]
        SSijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[2]])
        eeijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[3]])
      }
      result<-list(chart=CCij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new,
                   standardized_values=eeijk_new)
      return(result)
    }else if(side=="upward"){
      CCij_new<-matrix(NA,nind,nmaxobs)
      SSijk_new<-array(NA,c(nind,nmaxobs,ndim))
      eeijk_new<-array(NA,c(nind,nmaxobs,ndim))
      for(ii in 1:nind){
        idx_ii<-rep(ttij_int_new[ii,]-1,each=ndim)*ndim+1:ndim
        cov_ii<-cov_mat_est[idx_ii,idx_ii]
        eps_ii<-c(t(epsijk_new[ii,1:nobs_new[ii],]))
        chart_output<-mchart1_decorrelation_EWMA_multivariate_upward(
          eps_ii,cov_ii,nobs_new[ii],ndim,parameter,CL)
        CCij_new[ii,1:nobs_new[ii]]<-chart_output[[1]]
        SSijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[2]])
        eeijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[3]])
      }
      result<-list(chart=CCij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new,
                   standardized_values=eeijk_new)
      return(result)
    }else if(side=="downward"){
      CCij_new<-matrix(NA,nind,nmaxobs)
      SSijk_new<-array(NA,c(nind,nmaxobs,ndim))
      eeijk_new<-array(NA,c(nind,nmaxobs,ndim))
      for(ii in 1:nind){
        idx_ii<-rep(ttij_int_new[ii,]-1,each=ndim)*ndim+1:ndim
        cov_ii<-cov_mat_est[idx_ii,idx_ii]
        eps_ii<-c(t(epsijk_new[ii,1:nobs_new[ii],]))
        chart_output<-mchart1_decorrelation_EWMA_multivariate_downward(
          -eps_ii,cov_ii,nobs_new[ii],ndim,parameter,CL)
        CCij_new[ii,1:nobs_new[ii]]<-chart_output[[1]]
        SSijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[2]])
        eeijk_new[ii,1:nobs_new[ii],]<-t(chart_output[[3]])
      }
      result<-list(chart=CCij_new,
                   ttij_int=ttij_int_new,
                   SSijk=SSijk_new,
                   standardized_values=-eeijk_new)
      return(result)
    }
  }
  
  ##############################
  # nonparametric and standard #
  ##############################
  
  if(method=="nonparametric CUSUM"){
    stop("Method 'nonparametric CUSUM' currently not supported")
  }
  
  ##############################
  # nonparametric and standard #
  ##############################
  
  if(method=="nonparametric EWMA"){
    stop("Method 'nonparametric EWMA' currently not supported")
  }
}