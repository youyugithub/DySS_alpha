##' Ancillary Function for Cleaning Matrices
##'
##' @title Ancillary Function for Cleaning Matrices
##' @param matij a matrix to be cleaned by \code{nobs}. \cr
##' \code{matij[i,j]} is the jth observation of the ith subject.
##' @param nobs an integer vector of observation times. \cr
##' \code{nobs[i]} is the number of observations for the ith subject.
##' @param object_name name of the matrix object.
##' @return a cleaned matrix of \code{matij}. \cr
##' \code{matij[i,j]} is set to NA if \code{j>nobs[i]}
##' @author LY
##' @export
##' @examples 
clean_matij_by_nobs<-function(matij,nobs,object_name="matij"){
  nind<-dim(matij)[1]
  nmaxobs<-dim(matij)[2]
  if(any(nobs%%1!=0))stop("Error in 'nobs': 'nobs' should be a positive integer vector.")
  if(any(nobs<=0))stop("Error in 'nobs': 'nobs' should be a positive integer vector.")
  if(nmaxobs<max(nobs))stop(paste0("Error in dimensions of '",object_name,"' and 'nobs'."))
  if(nrow(matij)!=length(nobs))stop(paste0("Error in dimensions of '",object_name,"' and 'nobs'."))
  for(ii in 1:nind){
    if(any(!is.numeric(matij[ii,1:nobs[ii]])))stop(paste0("'",object_name,"' contains non-numeric values."))
    matij[ii,1:nmaxobs>nobs[ii]]<-NA
  }
  return(matij)
}
