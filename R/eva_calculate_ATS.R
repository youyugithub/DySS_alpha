##' Calculate ATS
##'
##' @title Calculate ATS
##' @param chartij a matrix of charting statistics. \cr
##' \code{chartij[i,j]} is the jth charting statistic of the ith subject.
##' @param ttij a matrix of observation times. \cr
##' \code{ttij[i,j]} is the jth observation time of the ith subject.
##' @param nobs an integer vector for number of observations. \cr
##' \code{nobs[i]} is the number of observations for the ith subject.
##' @param starttime a vector of times from the start of monitoring. \cr
##' \code{starttime[i]} is the time that the ith subject starts to be monitored.
##' @param endtime a vector of times from the start of monitoring. \cr
##' \code{endtime[i]} is the time that the ith subject is lost to be monitored.
##' @param ttmin,ttmax two numeric values. \cr
##' [\code{ttmin},\code{ttmax}] is the design time interval.
##' @param ntimepoints an integer value. \cr
##' \code{ntimepoints} is the number of basic timepoints in the design time interval.
##' @param CL a numeric value. \cr
##' \code{CL} is the control limit, signals will be given if charting statistics are greater than the control limit. 
##' @return ATS, a numeric value.
##' @author LY
##' @export
##' @examples 
eva_calculate_ATS<-function(
  chartij,ttij,nobs,starttime,endtime,ttmin,ttmax,ntimepoints,
  CL,no_signal_action="omit")
{
  if(any(dim(chartij)!=dim(ttij)))stop("Dimensions of 'chartij' and 'ttij' don't match.")
  if(dim(chartij)[1]!=length(nobs))stop("Dimensions of 'chartij' and 'nobs' don't match.")
  
  nind<-dim(chartij)[1]
  nmaxobs<-dim(chartij)[2]
  chartij<-clean_matij_by_nobs(chartij,nobs,"chartij")
  ttij<-clean_matij_by_nobs(ttij,nobs,"ttij")
  
  if(missing(starttime))ttij<-ttij[,1]
  if(missing(endtime))endtime<-ttij[cbind(1:nind,nobs)]
  if(missing(ttmin))ttmin<-min(ttij,na.rm=TRUE)
  if(missing(ttmax))ttmax<-max(ttij,na.rm=TRUE)
  if(missing(ntimepoints))stop("Error in argument 'ntimepoints'.")
  if(missing(CL))stop("Error in argument 'CL'.")
  
  if(!no_signal_action%in%c("omit","maxtime","centime"))
    stop("Error in argument 'no_signal_action'.")
  
  omega<-(ttmax-ttmin)/(ntimepoints-1)
  ttij_int<-round((ttij-ttmin)/omega)+1
  starttime_int<-round((starttime-ttmin)/omega)+1
  endtime_int<-round((endtime-ttmin)/omega)+1
  ttij_int_shifted<-sweep(ttij_int,1,starttime_int)
  
  if(no_signal_action=="omit"){
    ATS<-eva_calculate_ATS_omit(
      chartij,ttij_int_shifted,nobs,endtime_int,CL)
    return(ATS)
  }
  
  if(no_signal_action=="maxtime"){
    imputetime<-ntimepoints-starttime_int
    ATS<-eva_calculate_ATS_impute(
      chartij,ttij_int_shifted,nobs,endtime_int,CL)
    return(ATS)
  }
  
  if(no_signal_action=="endtime"){
    imputetime<-endtime_int-starttime_int
    ATS<-eva_calculate_ATS_impute(
      chartij,ttij_int_shifted,nobs,endtime_int,CL)
    return(ATS)
  }
}