##' Evaluate Control Charts
##'
##' @title Evaluate Control Charts
##' @param chartij a matrix of charting statistics. \cr
##' \code{chartij[i,j]} is the jth charting statistic of the ith subject.
##' @param ttij a matrix of observation times. \cr
##' \code{ttij[i,j]} is the jth observation time of the ith subject. \cr
##' \code{chartij[i,j,]} is the charting statistic of the ith subject at \code{ttij[i,j]}.
##' @param nobs an integer vector for number of observations. \cr
##' \code{nobs[i]} is the number of observations for the ith subject.
##' @param starttime a numeric vector.
##' \code{starttime[i]} is the time when monitoring starts for ith subject.
##' @param endtime a numeric vector, times when monitoring end.
##' \code{endtime[i]} is the time when monitoring ends for ith subject.
##' @param ttmin,ttmax two numeric values.
##' [\code{ttmin},\code{ttmax}] is the design time interval.
##' @param ntimepoints an integer value.
##' \code{ntimepoints} is the number of basic timepoints in the design time interval.
eva_control_chart_one<-function(
  chartij,ttij,nobs,
  starttime,endtime,event,
  ttmin,ttmax,ntimepoints,
  no_signal_action="maxtime"){
  
  if(any(dim(chartij)!=dim(ttij)))stop("Dimensions of 'chartij' and 'ttij' don't match.")
  if(!all.equal(dim(chartij)[1],length(starttime),length(endtime),length(event)))
    stop("Dimensions of 'starttime', 'endtime' and 'event' don't match.")
  
  chartij<-clean_matij_by_nobs(chartij,nobs,"chartij")
  ttij<-clean_matij_by_nobs(ttij,nobs,"ttij")
  
  add_thres<-min(chartij,na.rm=TRUE)-0.1
  
  omega<-(ttmax-ttmin)/(ntimepoints-1)
  ttij_int<-round((ttij-ttmin)/omega)+1
  starttime_int<-round((starttime-ttmin)/omega)+1
  endtime_int<-round((endtime-ttmin)/omega)+1
  ttij_int_shifted<-sweep(ttij_int,1,starttime_int)
  
  if(no_signal_action=="omit"){
    four_elements_output<-eva_four_elements_omit(
      chartij[!event,],ttij_int_shifted[!event,],nobs[!event],
      chartij[event,],ttij_int_shifted[event,],nobs[event],
      add_thres=add_thres)
    
    result<-list(
      thres=four_elements_output$thres,
      nFP=four_elements_output$nFP,
      nTP=four_elements_output$nTP,
      FPR=four_elements_output$FPR,
      TPR=four_elements_output$TPR,
      sumtime_ctrl=four_elements_output$sumtime_ctrl,
      sumtime_case=four_elements_output$sumtime_case,
      ATS0=four_elements_output$sumtime_ctrl/four_elements_output$nFP,
      ATS1=four_elements_output$sumtime_case/four_elements_output$nTP,
      ttmin=ttmin,
      ttmax=ttmax,
      ntimepoints=ntimepoints,
      no_signal_action=no_signal_action)
    
    class(result)<-"eva_control_chart"
    return(result)
  }
  
  if(no_signal_action=="maxtime"){
    
    imputetime<-ntimepoints-starttime_int
    
    four_elements_output<-eva_four_elements_impute(
      chartij[!event,],ttij_int_shifted[!event,],nobs[!event],imputetime[!event],
      chartij[event,],ttij_int_shifted[event,],nobs[event],imputetime[event],
      add_thres=add_thres)
    
    result<-list(
      thres=four_elements_output$thres,
      nFP=four_elements_output$nFP,
      nTP=four_elements_output$nTP,
      FPR=four_elements_output$FPR,
      TPR=four_elements_output$TPR,
      sumtime_ctrl=four_elements_output$sumtime_ctrl,
      sumtime_case=four_elements_output$sumtime_case,
      ATS0=four_elements_output$sumtime_ctrl/sum(!event),
      ATS1=four_elements_output$sumtime_case/sum(event),
      ttmin=ttmin,
      ttmax=ttmax,
      ntimepoints=ntimepoints,
      no_signal_action=no_signal_action)
    
    class(result)<-"eva_control_chart"
    return(result)
  }
  
  if(no_signal_action=="endtime"){
    
    imputetime<-endtime_int-starttime_int
    
    four_elements_output<-eva_four_elements_impute(
      chartij[!event,],ttij_int_shifted[!event,],nobs[!event],imputetime[!event],
      chartij[event,],ttij_int_shifted[event,],nobs[event],imputetime[event],
      add_thres=add_thres)
    
    result<-list(
      thres=four_elements_output$thres,
      nFP=four_elements_output$nFP,
      nTP=four_elements_output$nTP,
      FPR=four_elements_output$FPR,
      TPR=four_elements_output$TPR,
      sumtime_ctrl=four_elements_output$sumtime_ctrl,
      sumtime_case=four_elements_output$sumtime_case,
      ATS0=four_elements_output$sumtime_ctrl/sum(!event),
      ATS1=four_elements_output$sumtime_case/sum(event),
      ttmin=ttmin,
      ttmax=ttmax,
      ntimepoints=ntimepoints,
      no_signal_action=no_signal_action)
    
    class(result)<-"eva_control_chart"
    return(result)
  }
}

##' Evaluate Control Charts
##'
##' @title Evaluate Control Charts
##' @param chartij_IC,chartij_OC matrices of charting statistics. \cr
##' \code{chartij_IC[i,j]} is the jth charting statistic of the ith IC subject. \cr
##' \code{chartij_OC[i,j]} is the jth charting statistic of the ith OC subject.
##' @param ttij_IC,ttij_OC matrices of observation times. \cr
##' \code{ttij_IC[i,j]} is the jth observation time of the ith IC subject. \cr
##' \code{ttij_OC[i,j]} is the jth observation time of the ith OC subject. \cr
##' \code{chartij_IC[i,j]} is the charting statistic of the ith IC subject at \code{ttij[i,j]}. \cr
##' \code{chartij_OC[i,j]} is the charting statistic of the ith OC subject at \code{ttij[i,j]}.
##' @param nobs_IC,nobs_OC integer vectors for number of observations. \cr
##' \code{nobs_IC[i]} is the number of observations for the ith subject. \cr
##' \code{nobs_OC[i]} is the number of observations for the ith subject.
##' @param starttime_IC,starttime_OC a numeric vector. \cr
##' The times when monitoring start.
##' @param endtime_IC,endtime_OC a numeric vector, times when monitoring end. \cr
##' The times when monitoring end.
##' @param ttmin,ttmax two numeric values.
##' [\code{ttmin},\code{ttmax}] is the design time interval.
##' @param ntimepoints an integer value.
##' \code{ntimepoints} is the number of basic timepoints in the design time interval.
eva_control_chart_two<-function(
  chartij_IC,ttij_IC,nobs_IC,
  starttime_IC,endtime_IC,
  chartij_OC,ttij_OC,nobs_OC,
  starttime_OC,endtime_OC,
  ttmin,ttmax,ntimepoints,
  no_signal_action="maxtime"){
  
  if(any(dim(chartij_IC)!=dim(ttij_IC)))stop("Dimensions of 'chartij_IC' and 'ttij_IC' don't match.")
  if(any(dim(chartij_OC)!=dim(ttij_OC)))stop("Dimensions of 'chartij_OC' and 'ttij_OC' don't match.")
  if(!all.equal(dim(chartij_IC)[1],length(starttime_IC),length(endtime_IC)))
    stop("Dimensions of 'starttime_IC', 'endtime_IC' and 'event_IC' don't match.")
  if(!all.equal(dim(chartij_OC)[1],length(starttime_OC),length(endtime_OC)))
    stop("Dimensions of 'starttime_OC', 'endtime_OC' and 'event_OC' don't match.")
  
  chartij_IC<-clean_matij_by_nobs(chartij_IC,nobs_IC,"chartij_IC")
  chartij_OC<-clean_matij_by_nobs(chartij_OC,nobs_OC,"chartij_OC")
  ttij_IC<-clean_matij_by_nobs(ttij_IC,nobs_IC,"ttij_IC")
  ttij_OC<-clean_matij_by_nobs(ttij_OC,nobs_OC,"ttij_OC")
  nind_IC<-nrow(ttij_IC)
  nind_OC<-nrow(ttij_OC)
  
  add_thres<-min(c(chartij_IC,chartij_OC),na.rm=TRUE)-0.1
  
  omega<-(ttmax-ttmin)/(ntimepoints-1)
  ttij_IC_int<-round((ttij_IC-ttmin)/omega)+1
  ttij_OC_int<-round((ttij_OC-ttmin)/omega)+1
  starttime_IC_int<-round((starttime_IC-ttmin)/omega)+1
  starttime_OC_int<-round((starttime_OC-ttmin)/omega)+1
  endtime_IC_int<-round((endtime_IC-ttmin)/omega)+1
  endtime_OC_int<-round((endtime_OC-ttmin)/omega)+1
  ttij_IC_int_shifted<-sweep(ttij_IC_int,1,starttime_IC_int)
  ttij_OC_int_shifted<-sweep(ttij_OC_int,1,starttime_OC_int)
  
  if(no_signal_action=="omit"){
    four_elements_output<-eva_four_elements_omit(
      chartij_IC,ttij_IC_int,nobs_IC,
      chartij_OC,ttij_OC_int,nobs_OC,
      add_thres=add_thres)
    
    result<-list(
      thres=four_elements_output$thres,
      nFP=four_elements_output$nFP,
      nTP=four_elements_output$nTP,
      FPR=four_elements_output$FPR,
      TPR=four_elements_output$TPR,
      sumtime_ctrl=four_elements_output$sumtime_ctrl,
      sumtime_case=four_elements_output$sumtime_case,
      ATS0=four_elements_output$sumtime_ctrl/four_elements_output$nFP,
      ATS1=four_elements_output$sumtime_case/four_elements_output$nTP,
      ttmin=ttmin,
      ttmax=ttmax,
      ntimepoints=ntimepoints,
      no_signal_action=no_signal_action)
    
    class(result)<-"eva_control_chart"
    return(result)
  }
  
  if(no_signal_action=="maxtime"){
    
    imputetime_IC<-ntimepoints-starttime_IC_int
    imputetime_OC<-ntimepoints-starttime_OC_int
    
    four_elements_output<-eva_four_elements_impute(
      chartij_IC,ttij_IC_int_shifted,nobs_IC,imputetime_IC,
      chartij_OC,ttij_OC_int_shifted,nobs_OC,imputetime_OC,
      add_thres=add_thres)
    
    result<-list(
      thres=four_elements_output$thres,
      nFP=four_elements_output$nFP,
      nTP=four_elements_output$nTP,
      FPR=four_elements_output$FPR,
      TPR=four_elements_output$TPR,
      sumtime_ctrl=four_elements_output$sumtime_ctrl,
      sumtime_case=four_elements_output$sumtime_case,
      ATS0=four_elements_output$sumtime_ctrl/nind_IC,
      ATS1=four_elements_output$sumtime_case/nind_OC,
      ttmin=ttmin,
      ttmax=ttmax,
      ntimepoints=ntimepoints,
      no_signal_action=no_signal_action)
    
    class(result)<-"eva_control_chart"
    return(result)
  }
  
  if(no_signal_action=="endtime"){
    
    imputetime_IC<-endtime_IC_int-starttime_IC_int
    imputetime_OC<-endtime_OC_int-starttime_OC_int
    
    four_elements_output<-eva_four_elements_impute(
      chartij_IC,ttij_IC_int_shifted,nobs_IC,imputetime_IC,
      chartij_OC,ttij_OC_int_shifted,nobs_OC,imputetime_OC,
      add_thres=add_thres)
    
    result<-list(
      thres=four_elements_output$thres,
      nFP=four_elements_output$nFP,
      nTP=four_elements_output$nTP,
      FPR=four_elements_output$FPR,
      TPR=four_elements_output$TPR,
      sumtime_ctrl=four_elements_output$sumtime_ctrl,
      sumtime_case=four_elements_output$sumtime_case,
      ATS0=four_elements_output$sumtime_ctrl/nind_IC,
      ATS1=four_elements_output$sumtime_case/nind_OC,
      ttmin=ttmin,
      ttmax=ttmax,
      ntimepoints=ntimepoints,
      no_signal_action=no_signal_action)
    
    class(result)<-"eva_control_chart"
    return(result)
  }
}

##' Evaluate and Visualize Control Charts
##'
##' @title Evaluate and Visualize Control Charts
##' @param eva_control_chart an object of class eva_control_chart. \cr
##' \code{eva_control_chart} is an output from \code{eva_control_chart_one} or \code{eva_control_chart_two}.
plot.eva_control_chart<-function(eva_control_chart){
  if(eva_control_chart$no_signal_action=="omit"){
    old_par<-par(mfrow=c(1,1))
    plot(eva_control_chart$FPR,eva_control_chart$TPR,type="l",
         xlim=c(0,1),ylim=c(0,1),
         xlab="FPR",ylab="TPR",
         main="Positive Rates")
    par(old_par)
  }else{
    old_par<-par(mfrow=c(1,2))
    plot(eva_control_chart$FPR,eva_control_chart$TPR,type="l",
         xlim=c(0,1),ylim=c(0,1),
         xlab="FPR",ylab="TPR",
         main="Positive Rates")
    plot(eva_control_chart$ATS0,eva_control_chart$ATS1,type="l",
         xlim=c(0,eva_control_chart$ntimepoints),ylim=c(0,eva_control_chart$ntimepoints),
         xlab="ATS0",ylab="ATS1",
         main="Average Time to Signal")
    par(old_par)
  }
}

##' Evaluate and Visualize Control Charts
##'
##' @title Evaluate and Visualize Control Charts
##' @param eva_control_chart an object of class eva_control_chart. \cr
##' \code{eva_control_chart} is an output from \code{eva_control_chart_one} or \code{eva_control_chart_two}.
plot_PMROC<-function(eva_control_chart){
  if(eva_control_chart$no_signal_action=="omit"){
    message('plot_PMROC is not applicable when no_signal_action=="omit"')
  }else if(eva_control_chart$no_signal_action=="maxtime"){
    old_par<-par(mfrow=c(1,1),pty="s")
    VATS0<-with(eva_control_chart,(ATS0-min(ATS0))/(max(ATS0)-min(ATS0)))
    VATS1<-with(eva_control_chart,(ATS1-min(ATS1))/(max(ATS1)-min(ATS1)))
    DFPR<-(1-VATS0)*eva_control_chart$FPR
    TFPR<-(1-VATS1)*eva_control_chart$TPR
    plot(DFPR,TFPR,type="l",
         xlim=c(0,1),ylim=c(0,1),
         xlab="DFPR",ylab="TFPR",
         main="PM-ROC")
    par(old_par)
  }else if(eva_control_chart$no_signal_action=="endtime"){
    old_par<-par(mfrow=c(1,1),pty="s")
    VATS0<-with(eva_control_chart,(ATS0-min(ATS0))/(max(ATS0)-min(ATS0)))
    VATS1<-with(eva_control_chart,(ATS1-min(ATS1))/(max(ATS1)-min(ATS1)))
    DFPR<-(1-VATS0)*eva_control_chart$FPR
    TFPR<-(1-VATS1)*eva_control_chart$TPR
    plot(DFPR,TFPR,type="l",
         xlim=c(0,1),ylim=c(0,1),
         xlab="DFPR",ylab="TFPR",
         main="PM-ROC")
    par(old_par)
  }
}
