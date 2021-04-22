#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List chart_joint(arma::cube yyij,
                 arma::icube ttij,
                 arma::imat nobs,
                 int ntimepoints,
                 arma::imat ttij_eval,
                 arma::ivec nobs_eval,
                 arma::mat mum_est,
                 arma::cube Sigmam_est,
                 arma::mat Sigmae_est,
                 double lambda,
                 double delta_bar,
                 int hh,
                 int niter=200){
  // #define ARMA_NO_DEBUG
  // arma::mat all_iter_beta(ndim,niter);
  
  List result_temp(10);
  
  arma::wall_clock timer;
  double timing;
  timer.tic();
  
  const int nind=yyij.n_rows;
  const int ndim=yyij.n_slices;
  const int ndim3=ndim*3;
  const arma::uword I0=0,I1=ndim-1,II0=ndim,II1=2*ndim-1,III0=2*ndim,III1=3*ndim-1;
  int ii,jj,kk,ll,tt,ttdiff,iter,iter1,iitt,iijj;//,kk,iter2
  int nmax_tt_eval=ttij_eval.n_cols;
  double omega=1.0/double(ntimepoints);
  double ttdiff_double;
  
  const double tol=1e-5*double(ndim);
  double epsilon;
  
  ///--- Kernel/Timepoints vectors
  arma::vec allttdiff(hh*2+1,arma::fill::zeros);
  arma::vec ker(hh*2+1,arma::fill::zeros),tker(hh*2+1,arma::fill::zeros),ttker(hh*2+1,arma::fill::zeros);
  arma::vec tttker(hh*2+1,arma::fill::zeros),ttttker(hh*2+1,arma::fill::zeros);
  double tt_double,hh_double;
  hh_double=double(hh)*omega;
  for(tt=-hh;tt<=hh;tt++){
    tt_double=double(tt)*omega;
    allttdiff(tt+hh)=tt_double;
    ker(tt+hh)=std::pow(1.0-lambda,std::abs(double(tt)/delta_bar));
  }
  timing=timer.toc();
  Rcpp::Rcout<<"1: "<<timing<<"s"<<std::endl;
  
  ///--- Initialize parameters
  timer.tic();
  arma::field<arma::vec> F_mum(ntimepoints);
  arma::field<arma::mat> F_Sigmam(ntimepoints);
  arma::field<arma::mat> F_Sigmam_inv(ntimepoints);
  arma::field<arma::vec> F_Sigmae_vec(ntimepoints);
  arma::field<arma::vec> F_Sigmae_inv_vec(ntimepoints);
  
  arma::field<arma::vec> F_Emi(nmax_tt_eval,nind);
  // arma::field<arma::mat> F_Vmi(nmax_tt_eval,nind);
  
  arma::vec mum(ndim,arma::fill::zeros),mum_old(ndim,arma::fill::zeros);
  arma::mat Sigmam(ndim,ndim,arma::fill::eye);
  arma::mat Sigmam_inv(ndim,ndim,arma::fill::eye);
  Sigmam=10000.0*Sigmam;
  Sigmam_inv=arma::inv_sympd(Sigmam);
  arma::vec Sigmae_vec(ndim,arma::fill::ones);
  arma::vec Sigmae_inv_vec(ndim,arma::fill::ones);
  
  arma::vec Emi(ndim,arma::fill::zeros);
  arma::mat Vmi(ndim,ndim,arma::fill::eye);
  arma::mat Vmi_half(ndim,ndim,arma::fill::eye);
  
  for(iitt=0;iitt<ntimepoints;iitt++){
    mum=mum_est.col(iitt);
    Sigmam=Sigmam_est.slice(iitt);
    Sigmae_vec=Sigmae_est.col(iitt);
    
    F_mum(iitt)=mum;
    F_Sigmam(iitt)=Sigmam;
    F_Sigmam_inv(iitt)=arma::inv_sympd(Sigmam);
    F_Sigmae_vec(iitt)=Sigmae_vec;
    F_Sigmae_inv_vec(iitt)=1.0/Sigmae_vec;
  }
  
  // arma::vec SS(ndim);
  // arma::mat II(ndim,ndim);
  // double denominator,expterm;
  // arma::vec numerator_vec(ndim);
  // arma::mat numerator_mat(ndim,ndim);
  // arma::vec avec(ndim3);
  // arma::mat amat(ndim3,ndim3);
  
  double X0,X1,X2,X3,X4;
  double XY0,XY1,XY2;
  double SUMX0,SUMX1,SUMX2,SUMX3,SUMX4;
  arma::vec XKSKY(ndim);
  arma::mat XKSKX(ndim,ndim);
  double temp_double;
  int temp_int;
  arma::vec temp_ndim(ndim);
  arma::vec temp_ndim3(ndim3);
  arma::mat temp_ndim_ndim(ndim,ndim);
  arma::mat temp_ndim3_ndim3(ndim3,ndim3);
  arma::mat temp_sum_eee_across_jj(ndim,ndim);
  arma::vec hatmu(ndim),hatee(ndim);
  arma::vec3 temp_3;
  arma::mat33 temp_33;
  
  ///--- Create F_XKX
  arma::field<arma::mat> F_XKX(nmax_tt_eval,nind);
  arma::field<arma::vec> F_XKY(nmax_tt_eval,nind);
  arma::mat XKX(ndim,ndim);
  arma::vec XKY(ndim);
  // arma::vec F_nobs_in_ker(nmax_tt_eval);
  // arma::vec F_nind_at_risk(nmax_tt_eval);
  arma::mat F_sum_of_ker(nmax_tt_eval,ndim,arma::fill::zeros);
  // double nobs_in_ker,nind_at_risk,sum_of_ker;
  for(ii=0;ii<nind;ii++){
    for(iitt=0;iitt<nobs_eval(ii);iitt++){
      tt=ttij_eval(ii,iitt);
      
      XKX.fill(0.0);
      XKY.fill(0.0);
      for(kk=0;kk<ndim;kk++){
        // X0=0.0;X1=0.0;X2=0.0;X3=0.0;X4=0.0;
        // XY0=0.0;XY1=0.0;XY2=0.0;
        for(jj=0;jj<nobs(ii,kk);jj++){
          ttdiff=ttij(ii,jj,kk)-tt;
          if(std::abs(ttdiff)>=hh)continue;
          if(ttij(ii,jj,kk)>tt)continue;
          ttdiff_double=double(ttdiff)*omega;
          
          XKX(I0+kk,I0+kk)=XKX(I0+kk,I0+kk)+ker(ttdiff+hh);
          // XKX(II0+kk,I0+kk)=XKX(II0+kk,I0+kk)+ttdiff_double*ker(ttdiff+hh);
          // XKX(III0+kk,I0+kk)=XKX(III0+kk,I0+kk)+ttdiff_double*ttdiff_double*ker(ttdiff+hh);
          // XKX(I0+kk,II0+kk)=XKX(I0+kk,II0+kk)+ttdiff_double*ker(ttdiff+hh);
          // XKX(II0+kk,II0+kk)=XKX(II0+kk,II0+kk)+ttdiff_double*ttdiff_double*ker(ttdiff+hh);
          // XKX(III0+kk,II0+kk)=XKX(III0+kk,II0+kk)+ttdiff_double*ttdiff_double*ttdiff_double*ker(ttdiff+hh);
          // XKX(I0+kk,III0+kk)=XKX(I0+kk,III0+kk)+ttdiff_double*ttdiff_double*ker(ttdiff+hh);
          // XKX(II0+kk,III0+kk)=XKX(II0+kk,III0+kk)+ttdiff_double*ttdiff_double*ttdiff_double*ker(ttdiff+hh);
          // XKX(III0+kk,III0+kk)=XKX(III0+kk,III0+kk)+ttdiff_double*ttdiff_double*ttdiff_double*ttdiff_double*ker(ttdiff+hh);
          
          XKY(I0+kk)=XKY(I0+kk)+ker(ttdiff+hh)*yyij(ii,jj,kk);
          // XKY(II0+kk)=XKY(II0+kk)+ttdiff_double*ker(ttdiff+hh)*yyij(ii,jj,kk);
          // XKY(III0+kk)=XKY(III0+kk)+ttdiff_double*ttdiff_double*ker(ttdiff+hh)*yyij(ii,jj,kk);
          
          F_sum_of_ker(iitt,kk)=F_sum_of_ker(iitt,kk)+ker(ttdiff+hh);
        }
      }
      F_XKX(iitt,ii)=XKX;
      F_XKY(iitt,ii)=XKY;
    }
    // F_nobs_in_ker(iitt)=nobs_in_ker;
    // F_nind_at_risk(iitt)=nind_at_risk;
  }
  
  Rcout<<"2: "<<std::endl;
  
  double sum_for_nind;
  double ee;
  arma::vec sum_for_Emi(ndim);
  arma::mat sum_for_Sigmam(ndim,ndim);
  arma::vec sum_for_Sigmae(ndim);
  arma::vec temp_sum_XVX_across_jj(ndim);
  arma::vec temp_sum_eeeeK_across_jj(ndim);
  
  for(ii=0;ii<nind;ii++){
    for(iitt=0;iitt<nobs_eval(ii);iitt++){
      tt=ttij_eval(ii,iitt);
      mum=F_mum(tt-1);
      Sigmam=F_Sigmam(tt-1);
      Sigmam_inv=F_Sigmam_inv(tt-1);
      Sigmae_vec=F_Sigmae_vec(tt-1);
      Sigmae_inv_vec=F_Sigmae_inv_vec(tt-1);
      
      XKSKY=F_XKY(iitt,ii);
      XKSKX=F_XKX(iitt,ii);
      for(kk=0;kk<ndim;kk++){
        XKSKX.row(I0+kk)=XKSKX.row(I0+kk)*Sigmae_inv_vec(kk);
        // XKSKX.row(II0+kk)=XKSKX.row(II0+kk)*Sigmae_inv_vec(kk);
        // XKSKX.row(III0+kk)=XKSKX.row(III0+kk)*Sigmae_inv_vec(kk);
        XKSKY(I0+kk)=XKSKY(I0+kk)*Sigmae_inv_vec(kk);
        // XKSKY(II0+kk)=XKSKY(II0+kk)*Sigmae_inv_vec(kk);
        // XKSKY(III0+kk)=XKSKY(III0+kk)*Sigmae_inv_vec(kk);
      }
      
      Vmi=arma::inv_sympd(Sigmam_inv+XKSKX);
      Emi=mum+Vmi*(XKSKY-XKSKX*mum);// Emi=mum+Sigmam*(XKSKY-XKSKX*Vmi*XKSKY);
      // F_Vmi(iitt,ii)=Vmi;
      F_Emi(iitt,ii)=Emi;
    }
  }
  
  timing=timer.toc();
  Rcpp::Rcout<<"5: "<<timing<<"s"<<std::endl;
  
  List result(1);
  
  result(0)=F_Emi;
  
  return(result);
}
