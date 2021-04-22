##' Search Control Limit
##'
##' @title Search Control Limit
##' @param chartij a matrix of charting statistics. \cr
##' \code{chartij[i,j]} is the jth charting statistic of the ith subject.
##' @param ttij a matrix of observation times.  \cr
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
##' @param ATS_nominal a numeric value. \cr
##' \code{ATS_nominal} is the targeted/nominal ATS one wished to achieve. 
##' @param CL_lower,CL_step,CL_upper three numeric values. \cr
##' The control limit will be searched within the interval [\code{CL_lower},\code{CL_upper}]. \cr
##' When applying grid search, the algorithm will use a step size of \code{CL_step}. \cr
##' (Namely, the algorithm will start with \code{CL_lower}, 
##' and search through the sequences \code{CL_lower}, \code{CL_lower+CL_step}, \code{CL_lower+2*CL_step}, ... until \code{CL_upper}.)
##' @param ATS_tol a numeric value. \cr
##' Error tolerence for ATS.
##' @param CL_tol a numeric value. \cr
##' Error tolerence for control limit.
##' @return control limit, a numeric value.
##' @author LY
##' @export
##' @examples 
search_CL<-function(
  chartij,ttij,nobs,starttime,endtime,ttmin,ttmax,ntimepoints,ATS_nominal,
  CL_lower,CL_step,CL_upper,
  no_signal_action="omit",ATS_tol,CL_tol)
{
  nind<-dim(chartij)[1]
  nmaxobs<-dim(chartij)[2]
  if(missing(starttime))ttij<-ttij[,1]
  if(missing(endtime))endtime<-ttij[cbind(1:nind,nobs)]
  if(missing(ttmin))ttmin<-min(ttij,na.rm=TRUE)
  if(missing(ttmax))ttmax<-max(ttij,na.rm=TRUE)
  if(missing(ntimepoints))stop("Error in argument 'ntimepoints'.")
  if(missing(ATS_nominal))stop("Error in argument 'ATS_nominal'.")
  if(missing(CL_lower))CL_lower<-min(chartij,na.rm=TRUE)
  if(missing(CL_upper))CL_upper<-max(chartij,na.rm=TRUE)
  if(missing(CL_step))CL_step<-(CL_upper-CL_lower)/20
  if(missing(ATS_tol))ATS_tol<-ATS_nominal*1e-7
  if(missing(CL_tol))CL_tol<-(CL_upper-CL_lower)*1e-7
  
  if(!no_signal_action%in%c("omit","maxtime","centime"))stop("Error in argument 'no_signal_action'.")
  
  omega<-(ttmax-ttmin)/(ntimepoints-1)
  ttij_int<-round((ttij-ttmin)/omega)+1
  starttime_int<-round((starttime-ttmin)/omega)+1
  endtime_int<-round((endtime-ttmin)/omega)+1
  ttij_int_shifted<-sweep(ttij,1,starttime_int)
  
  if(no_signal_action=="omit"){
    CL<-add_search_control_limit_omit(
      chartij,ttij_int_shifted,nobs,endtime,ATS_nominal,ATS_tol,
      CL_lower,CL_step,CL_upper,CL_tol)
    if(CL-CL_lower<CL_tol|CL_upper-CL<CL_tol)
      stop("Control limit not found in the range: ","[",CL_lower,",",CL_upper,"]")
    return(CL)
  }
  
  if(no_signal_action=="maxtime"){
    imputetime<-ntimepoints-starttime_int
    CL<-add_search_control_limit_impute(
      chartij,ttij_int_shifted,nobs,imputetime,ATS_nominal,ATS_tol,
      CL_lower,CL_step,CL_upper,CL_tol)
    if(CL-CL_lower<CL_tol|CL_upper-CL<CL_tol)
      stop("Control limit not found in the range","[",CL_lower,",",CL_upper,"]")
    return(CL)
  }
  
  if(no_signal_action=="endtime"){
    imputetime<-endtime_int-starttime_int
    CL<-add_search_control_limit_impute(
      chartij,ttij_int_shifted,nobs,imputetime,ATS_nominal,ATS_tol,
      CL_lower,CL_step,CL_upper,CL_tol)
    if(CL-CL_lower<CL_tol|CL_upper-CL<CL_tol)
      stop("Control limit not found in the range","[",CL_lower,",",CL_upper,"]")
    return(CL)
  }
}
