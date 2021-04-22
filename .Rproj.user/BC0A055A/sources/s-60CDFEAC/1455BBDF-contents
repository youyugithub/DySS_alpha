##' Monitor Univariate Longitudinal Data
##'
##' @title Monitor Univariate Longitudinal Data
##' @param yyij_new a matrix of longitudinal observations \cr
##' \code{yyij_new[i,j]} is the jth observation of the ith subject.
##' @param ttij_new a matrix of observation times \cr
##' \code{ttij_new[i,j]} is the jth observation time of the ith subject. \cr
##' \code{yyij_new[i,j]} is observed at \code{ttij_new[i,j]}.
##' @param nobs_new an integer vector for number of observations. \cr
##' \code{nobs_new[i]} is the number of observations for the ith subject.
##' @param pattern the estimated regular longitudinal pattern
##' @param side a string \cr
##' If \code{side="upward"}, apply control charts that aim to detect upward shifts.\cr
##' If \code{side="downward"}, apply control charts that aim to detect downward shifts.\cr
##' If \code{side="both"}, apply control charts that aim to detect shifts in both sides.
##' @param chart a string \cr
##' If \code{chart="CUSUM"}, apply CUSUM charts.\cr
##' If \code{chart="EWMA"}, apply EWMA charts.
##' @param method a string \cr
##' If \code{method="standard"}, standardize observations by mean and variance (cf., Qiu and Xiang, 2014).\cr
##' If \code{method="decorrelation"}, standardize and decorrelate observations by mean and covariance (cf., Li and Qiu, 2016).\cr
##' If \code{method="sprint"}, standardize and decorrelate observations within sprint length by mean and covariance (cf., You and Qiu 2018).\cr
##' If \code{method="distribution and standard"}, standardize observations by distribution (cf., You and Qiu, 2020).\cr
##' If \code{method="distribution and decorrelation"}, standardize observations by distribution and covariance (cf., You and Qiu, 2020).\cr
##' If \code{method="distribution and sprint"},standardize and decorrelate observations within sprint length by distribution and covariance (cf., You and Qiu, 2020).\cr
##' \code{method="nonparametric and standard"} currently not supported.\cr
##' \code{method="nonparametric and decorrelation"} currently not supported
##' @param parameter a numeric value \cr
##' If \code{chart="CUSUM"}, \code{parameter} is the allowance constant in the control chart.\cr
##' If \code{chart="EWMA"}, \code{parameter} is the weighting in the control chart.
##' @param CL a numeric value \cr
##' \code{CL} is the control limit. 
##' A signal will be given if charting statistics are larger than the control limit. 
##' (Note: in this package, signs of charting statistics may be reversed such that 
##' larger values of charting statistics indicate worse performance of processes.)
##' After the signal is given, the algorithm stops calculating the charting statistics for the remaining observation times.
##' The default value of control limit is infinity, which means we will calculate the charting statistics for all observation times.
##' @return a list \cr
##' \code{$chart} a numeric matrix, \code{$chart[i,j]} is the jth charting statistic of the ith subject. \cr
##' \code{$ttij_int} an integer matrix, the observation times rounded to the multiples of basic time unit. \cr
##' (Supplementary) \code{$chart_upward}, \code{$chart_downward} numeric matrices, the upward/downward charting statistics. \cr
##' (Supplementary) \code{$standardized_values}, \code{$standardized_values_upward}, \code{$standardized_values_downward} numeric matrices, 
##' the standardized values used in the calculation of (downward/upward/two-sided) control charts.
mnt_long_1d<-function(
  yyij_new,ttij_new,nobs_new,pattern,
  side="upward",chart="CUSUM",
  method="standard",
  parameter=0.5,
  CL=Inf)
{
  if(any(dim(yyij_new)!=dim(ttij_new)))stop("Dimensions of 'yyij_new' and 'ttij_new' don't match.")
  if(dim(yyij_new)[1]!=length(nobs_new))stop("Dimensions of 'yyij_new' and 'nobs_new' don't match.")
  
  nind<-dim(yyij_new)[1]
  nmaxobs<-dim(yyij_new)[2]
  yyij_new<-clean_matij_by_nobs(yyij_new,nobs_new,"yyij")
  ttij_new<-clean_matij_by_nobs(ttij_new,nobs_new,"ttij")
  
  ttmax<-pattern$ttmax
  ttmin<-pattern$ttmin
  ntimepoints<-pattern$ntimepoints
  if(ttmin>min(ttij_new,na.rm=TRUE))stop("Value of 'ttmin' is too large. Please re-estimate the pattern with a smaller ttmin.")
  if(ttmax<min(ttij_new,na.rm=TRUE))stop("Value of 'ttmax' is too small. Please re-estimate the pattern with a larger ttmax.")
  omega<-(ttmax-ttmin)/(ntimepoints-1)
  ttij_int_new<-round((ttij_new-ttmin)/omega)+1
  
  if(!side%in%c("upward","downward","both"))stop("Error in argument 'side'.")
  if(chart=="CUSUM"){
    if(parameter<0.0)
      stop("For chart 'CUSUM', parameter should be greater than 0.")
  }else if(chart=="EWMA"){
    if(parameter>1.0|parameter<0.0)
      stop("For chart 'EWMA', parameter should be in the interval [0,1].")
  }else{
    stop("Error in argument 'chart'.")
  }
  
  all_methods<-c("standard",
                 "decorrelation",
                 "sprint",
                 "distribution and standard",
                 "distribution and decorrelation",
                 "distribution and sprint",
                 "nonparametric and standard",
                 "nonparametric and decorrelation")
  if(!method%in%all_methods)stop("Error in argument 'method'.")
  
  ############
  # standard #
  ############
  
  if(method=="standard"){
    if(!inherits(pattern,c("patn_long_univ_meanvar","patn_long_univ_meanvarcov")))
      stop("Pattern does not match method.")
    mean_est<-pattern$mean_est
    var_est<-pattern$var_est
    sd_est<-sqrt(var_est)
    eeij_new<-matrix(NA,nind,nmaxobs)
    for(ii in 1:nind){
      eeij_new[ii,1:nobs_new[ii]]<-
        (yyij_new[ii,1:nobs_new[ii]]-mean_est[ttij_int_new[ii,1:nobs_new[ii]]])/
        sd_est[ttij_int_new[ii,1:nobs_new[ii]]]
    }
    if(side=="upward"){
      if(chart=="CUSUM"){
        chartij_new<-chart_CUSUM_univariate_std(eeij_new,ttij_int_new,nobs_new,parameter,CL)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }else if(chart=="EWMA"){
        chartij_new<-chart_EWMA_univariate_std(eeij_new,ttij_int_new,nobs_new,parameter,CL)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }
    }else if(side=="downward"){
      if(chart=="CUSUM"){
        chartij_new<-chart_CUSUM_univariate_std(-eeij_new,ttij_int_new,nobs_new,parameter,CL)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }else if(chart=="EWMA"){
        chartij_new<-chart_EWMA_univariate_std(-eeij_new,ttij_int_new,nobs_new,parameter,CL)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }
    }else if(side=="both"){
      if(chart=="CUSUM"){
        chartij_upward_new<-chart_CUSUM_univariate_std(eeij_new,ttij_int_new,nobs_new,parameter,CL)
        chartij_downward_new<-chart_CUSUM_univariate_std(-eeij_new,ttij_int_new,nobs_new,parameter,CL)
        chartij_new<-pmax(chartij_upward_new,chartij_downward_new)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     chart_upward=chartij_upward_new,
                     chart_downward=chartij_downward_new,
                     standardized_values=eeij_new)
        return(result)
      }else if(chart=="EWMA"){
        chartij_upward_new<-chart_EWMA_univariate_std(eeij_new,ttij_int_new,nobs_new,parameter,CL)
        chartij_downward_new<-chart_EWMA_univariate_std(-eeij_new,ttij_int_new,nobs_new,parameter,CL)
        chartij_new<-pmax(chartij_upward_new,chartij_downward_new)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     chart_upward=chartij_upward_new,
                     chart_downward=chartij_downward_new,
                     standardized_values=eeij_new)
        return(result)
      }
    }
  }
  
  #################
  # decorrelation #
  #################
  
  if(method=="decorrelation"){
    if(!inherits(pattern,"patn_long_univ_meanvarcov"))
      stop("Pattern does not match method.")
    mean_est<-pattern$mean_est
    var_est<-pattern$var_est
    cov_est<-pattern$cov_est
    epsij_new<-matrix(NA,nind,nmaxobs)
    for(ii in 1:nind){
      epsij_new[ii,1:nobs_new[ii]]<-yyij_new[ii,1:nobs_new[ii]]-mean_est[ttij_int_new[ii,1:nobs_new[ii]]]
    }
    if(side=="upward"){
      if(chart=="CUSUM"){
        chart_output<-chart_CUSUM_univariate_dec(epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<-chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }else if(chart=="EWMA"){
        chart_output<-chart_EWMA_univariate_dec(epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<-chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }
    }else if(side=="downward"){
      if(chart=="CUSUM"){
        chart_output<-chart_CUSUM_univariate_dec(-epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<- -chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }else if(chart=="EWMA"){
        chart_output<-chart_EWMA_univariate_dec(-epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<- -chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }
    }else if(side=="both"){
      if(chart=="CUSUM"){
        chart_upward_output<-chart_CUSUM_univariate_dec(epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chart_downward_output<-chart_CUSUM_univariate_dec(-epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_upward_new<-chart_upward_output[[1]]
        chartij_downward_new<- -chart_downward_output[[1]]
        eeij_upward_new<-chart_upward_output[[2]]
        eeij_downward_new<- -chart_downward_output[[2]]
        chartij_new<-pmax(chartij_upward_new,-chartij_downward_new)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     chart_upward=chartij_upward_new,
                     chart_downward=chartij_downward_new,
                     standardized_values_upward=eeij_upward_new,
                     standardized_values_downward=eeij_downward_new)
        return(result)
      }else if(chart=="EWMA"){
        chart_upward_output<-chart_EWMA_univariate_dec(epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chart_downward_output<-chart_EWMA_univariate_dec(-epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_upward_new<-chart_upward_output[[1]]
        chartij_downward_new<- -chart_downward_output[[1]]
        eeij_upward_new<-chart_upward_output[[2]]
        eeij_downward_new<- -chart_downward_output[[2]]
        chartij_new<-pmax(chartij_upward_new,-chartij_downward_new)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     chart_upward=chartij_upward_new,
                     chart_downward=chartij_downward_new,
                     standardized_values_upward=eeij_upward_new,
                     standardized_values_downward=eeij_downward_new)
        return(result)
      }
    }
  }
  
  ##########
  # sprint #
  ##########
  
  if(method=="sprint"){
    if(!inherits(pattern,"patn_long_univ_meanvarcov"))
      stop("Pattern does not match method.")
    mean_est<-pattern$mean_est
    var_est<-pattern$var_est
    cov_est<-pattern$cov_est
    epsij_new<-matrix(NA,nind,nmaxobs)
    for(ii in 1:nind){
      epsij_new[ii,1:nobs_new[ii]]<-yyij_new[ii,1:nobs_new[ii]]-mean_est[ttij_int_new[ii,1:nobs_new[ii]]]
    }
    if(side=="upward"){
      if(chart=="CUSUM"){
        chart_output<-chart_CUSUM_univariate_sprint(epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<-chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }else if(chart=="EWMA"){
        chart_output<-chart_EWMA_univariate_sprint(epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<-chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }
    }else if(side=="downward"){
      if(chart=="CUSUM"){
        chart_output<-chart_CUSUM_univariate_sprint(-epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<- -chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }else if(chart=="EWMA"){
        chart_output<-chart_EWMA_univariate_sprint(-epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<- -chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }
    }else if(side=="both"){
      if(chart=="CUSUM"){
        chart_upward_output<-chart_CUSUM_univariate_sprint(epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chart_downward_output<-chart_CUSUM_univariate_sprint(-epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_upward_new<-chart_upward_output[[1]]
        chartij_downward_new<- -chart_downward_output[[1]]
        eeij_upward_new<-chart_upward_output[[2]]
        eeij_downward_new<- -chart_downward_output[[2]]
        chartij_new<-pmax(chartij_upward_new,-chartij_downward_new)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     chart_upward=chartij_upward_new,
                     chart_downward=chartij_downward_new,
                     standardized_values_upward=eeij_upward_new,
                     standardized_values_downward=eeij_downward_new)
        return(result)
      }else if(chart=="EWMA"){
        chart_upward_output<-chart_EWMA_univariate_sprint(epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chart_downward_output<-chart_EWMA_univariate_sprint(-epsij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_upward_new<-chart_upward_output[[1]]
        chartij_downward_new<- -chart_downward_output[[1]]
        eeij_upward_new<-chart_upward_output[[2]]
        eeij_downward_new<- -chart_downward_output[[2]]
        chartij_new<-pmax(chartij_upward_new,-chartij_downward_new)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     chart_upward=chartij_upward_new,
                     chart_downward=chartij_downward_new,
                     standardized_values_upward=eeij_upward_new,
                     standardized_values_downward=eeij_downward_new)
        return(result)
      }
    }
  }
  
  #############################
  # distribution and standard #
  #############################
  
  if(method=="distribution and standard"){
    if(!inherits(pattern,c("patn_long_univ_dist","patn_long_univ_distvarcov")))
      stop("Pattern does not match method.")
    yy_ref<-pattern$yy_ref
    tt_ref<-pattern$tt_ref
    starting_idx<-pattern$starting_idx
    ending_idx<-pattern$ending_idx
    upper_line<-pattern$upper_line
    hh_t<-pattern$hh_t
    hh_y<-pattern$hh_y
    zzij_new<-local_const_percentile_est_faster(
      yyij_new,ttij_int_new,nobs_new,yyref,tt_ref,
      starting_idx,ending_idx,upper_line,ntimepoints,hh_t,hh_y)
    if(side=="upward"){
      if(chart=="CUSUM"){
        chartij_new<-chart_CUSUM_univariate_std(zzij_new,ttij_int_new,nobs_new,parameter,CL)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=zzij_new)
        return(result)
      }else if(chart=="EWMA"){
        chartij_new<-chart_EWMA_univariate_std(zzij_new,ttij_int_new,nobs_new,parameter,CL)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=zzij_new)
        return(result)
      }
    }else if(side=="downward"){
      if(chart=="CUSUM"){
        chartij_new<-chart_CUSUM_univariate_std(-zzij_new,ttij_int_new,nobs_new,parameter,CL)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=zzij_new)
        return(result)
      }else if(chart=="EWMA"){
        chartij_new<-chart_EWMA_univariate_std(-zzij_new,ttij_int_new,nobs_new,parameter,CL)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=zzij_new)
        return(result)
      }
    }else if(side=="both"){
      if(chart=="CUSUM"){
        chart_upward_output<-chart_CUSUM_univariate_std(zzij_new,ttij_int_new,nobs_new,parameter,CL)
        chart_downward_output<-chart_CUSUM_univariate_std(zzij_new,ttij_int_new,nobs_new,parameter,CL)
        chartij_upward_new<-chart_upward_output[[1]]
        chartij_downward_new<- -chart_downward_output[[1]]
        eeij_upward_new<-chart_upward_output[[2]]
        eeij_downward_new<- -chart_downward_output[[2]]
        chartij_new<-pmax(chartij_upward_new,-chartij_downward_new)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     chart_upward=chartij_upward_new,
                     chart_downward=chartij_downward_new,
                     standardized_values_upward=eeij_upward_new,
                     standardized_values_downward=eeij_downward_new)
        return(result)
      }else if(chart=="EWMA"){
        chart_upward_output<-chart_EWMA_univariate_std(zzij_new,ttij_int_new,nobs_new,parameter,CL)
        chart_downward_output<-chart_EWMA_univariate_std(zzij_new,ttij_int_new,nobs_new,parameter,CL)
        chartij_upward_new<-chart_upward_output[[1]]
        chartij_downward_new<- -chart_downward_output[[1]]
        eeij_upward_new<-chart_upward_output[[2]]
        eeij_downward_new<- -chart_downward_output[[2]]
        chartij_new<-pmax(chartij_upward_new,-chartij_downward_new)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     chart_upward=chartij_upward_new,
                     chart_downward=chartij_downward_new,
                     standardized_values_upward=eeij_upward_new,
                     standardized_values_downward=eeij_downward_new)
        return(chartij_new)
      }
    }
  }
  
  ##################################
  # distribution and decorrelation #
  ##################################
  
  if(method=="distribution and decorrelation"){
    if(!inherits(pattern,"patn_long_univ_distvarcov"))
      stop("Pattern does not match method.")
    yy_ref<-pattern$yy_ref
    tt_ref<-pattern$tt_ref
    starting_idx<-pattern$starting_idx
    ending_idx<-pattern$ending_idx
    upper_line<-pattern$upper_line
    cov_est<-pattern$cov_est
    hh_t<-pattern$hh_t
    hh_y<-pattern$hh_y
    hh_cov<-pattern$hh_cov
    zzij_new<-local_const_percentile_est_faster(
      yyij_new,ttij_int_new,nobs_new,yyref,tt_ref,
      starting_idx,ending_idx,upper_line,ntimepoints,hh_t,hh_y)
    if(side=="upward"){
      if(chart=="CUSUM"){
        chart_output<-chart_CUSUM_univariate_dec(zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<-chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }else if(chart=="EWMA"){
        chart_output<-chart_EWMA_univariate_dec(zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<-chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }
    }else if(side=="downward"){
      if(chart=="CUSUM"){
        chart_output<-chart_CUSUM_univariate_dec(-zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<- -chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }else if(chart=="EWMA"){
        chart_output<-chart_EWMA_univariate_dec(-zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<- -chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }
    }else if(side=="both"){
      if(chart=="CUSUM"){
        chart_upward_output<-chart_CUSUM_univariate_dec(zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chart_downward_output<-chart_CUSUM_univariate_dec(-zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_upward_new<-chart_upward_output[[1]]
        chartij_downward_new<- -chart_downward_output[[1]]
        eeij_upward_new<-chart_upward_output[[2]]
        eeij_downward_new<- -chart_downward_output[[2]]
        chartij_new<-pmax(chartij_upward_new,-chartij_downward_new)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     chart_upward=chartij_upward_new,
                     chart_downward=chartij_downward_new,
                     standardized_values_upward=eeij_upward_new,
                     standardized_values_downward=eeij_downward_new)
        return(result)
      }else if(chart=="EWMA"){
        chart_upward_output<-chart_EWMA_univariate_dec(zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chart_downward_output<-chart_EWMA_univariate_dec(-zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_upward_new<-chart_upward_output[[1]]
        chartij_downward_new<- -chart_downward_output[[1]]
        eeij_upward_new<-chart_upward_output[[2]]
        eeij_downward_new<- -chart_downward_output[[2]]
        chartij_new<-pmax(chartij_upward_new,-chartij_downward_new)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     chart_upward=chartij_upward_new,
                     chart_downward=chartij_downward_new,
                     standardized_values_upward=eeij_upward_new,
                     standardized_values_downward=eeij_downward_new)
        return(result)
      }
    }
  }
  
  ###########################
  # distribution and sprint #
  ###########################
  
  if(method=="distribution and sprint"){
    if(!inherits(pattern,"patn_long_univ_distvarcov"))
      stop("Pattern does not match method.")
    yy_ref<-pattern$yy_ref
    tt_ref<-pattern$tt_ref
    starting_idx<-pattern$starting_idx
    ending_idx<-pattern$ending_idx
    upper_line<-pattern$upper_line
    cov_est<-pattern$cov_est
    hh_t<-pattern$hh_t
    hh_y<-pattern$hh_y
    hh_cov<-pattern$hh_cov
    zzij_new<-local_const_percentile_est_faster(
      yyij_new,ttij_int_new,nobs_new,yyref,tt_ref,
      starting_idx,ending_idx,upper_line,ntimepoints,hh_t,hh_y)
    if(side=="upward"){
      if(chart=="CUSUM"){
        chart_output<-chart_CUSUM_univariate_sprint(zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<-chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }else if(chart=="EWMA"){
        chart_output<-chart_EWMA_univariate_sprint(zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<-chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }
    }else if(side=="downward"){
      if(chart=="CUSUM"){
        chart_output<-chart_CUSUM_univariate_sprint(-zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<- -chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }else if(chart=="EWMA"){
        chart_output<-chart_EWMA_univariate_sprint(-zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_new<-chart_output[[1]]
        eeij_new<- -chart_output[[2]]
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     standardized_values=eeij_new)
        return(result)
      }
    }else if(side=="both"){
      if(chart=="CUSUM"){
        chart_upward_output<-chart_CUSUM_univariate_sprint(zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chart_downward_output<-chart_CUSUM_univariate_sprint(-zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_upward_new<-chart_upward_output[[1]]
        chartij_downward_new<- -chart_downward_output[[1]]
        eeij_upward_new<-chart_upward_output[[2]]
        eeij_downward_new<- -chart_downward_output[[2]]
        chartij_new<-pmax(chartij_upward_new,-chartij_downward_new)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     chart_upward=chartij_upward_new,
                     chart_downward=chartij_downward_new,
                     standardized_values_upward=eeij_upward_new,
                     standardized_values_downward=eeij_downward_new)
        return(result)
      }else if(chart=="EWMA"){
        chart_upward_output<-chart_EWMA_univariate_sprint(zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chart_downward_output<-chart_EWMA_univariate_sprint(-zzij_new,ttij_int_new,nobs_new,cov_est,parameter,CL)
        chartij_upward_new<-chart_upward_output[[1]]
        chartij_downward_new<- -chart_downward_output[[1]]
        eeij_upward_new<-chart_upward_output[[2]]
        eeij_downward_new<- -chart_downward_output[[2]]
        chartij_new<-pmax(chartij_upward_new,-chartij_downward_new)
        result<-list(chart=chartij_new,
                     ttij_int=ttij_int_new,
                     chart_upward=chartij_upward_new,
                     chart_downward=chartij_downward_new,
                     standardized_values_upward=eeij_upward_new,
                     standardized_values_downward=eeij_downward_new)
        return(result)
      }
    }
  }
  
  ##############################
  # nonparametric and standard #
  ##############################
  
  if(method=="nonparametric and standard"){
    stop("Method 'nonparametric and standard' currently not supported")
  }
  
  ###################################
  # nonparametric and decorrelation #
  ###################################
  
  if(method=="nonparametric and decorrelation"){
    stop("Method 'nonparametric and decorrelation' currently not supported")
  }
}
