#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]

List joint_local_const(arma::cube yyijk,
                       arma::icube ttijk,
                       arma::imat nobs,
                       arma::ivec st,
                       arma::ivec ot,
                       LogicalVector delta,
                       arma::ivec alltimepoints,
                       int hh,
                       int niter=200){
  
  const int nind=yyijk.n_rows;
  const int ndim=yyijk.n_slices;
  
  arma::mat all_iter_beta(ndim,niter);
  
  List result_temp(10);
  
  const int ndim1=ndim;
  const arma::uword I0=0,I1=ndim-1;
  int ii,jj,kk,ll,tt,ttdiff,iter,iter1,iitt;//,kk,iter2
  int ntimepoints=alltimepoints.n_elem;
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
    ker(tt+hh)=std::max(0.0,1.0-(tt_double/hh_double)*(tt_double/hh_double));
    tker(tt+hh)=ker(tt+hh)*tt_double;
    ttker(tt+hh)=tker(tt+hh)*tt_double;
    tttker(tt+hh)=ttker(tt+hh)*tt_double;
    ttttker(tt+hh)=tttker(tt+hh)*tt_double;
  }

  ///--- Initialize parameters
  arma::field<arma::vec> F_mub(ntimepoints);
  arma::field<arma::mat> F_Sigmab(ntimepoints);
  arma::field<arma::mat> F_Sigmab_inv(ntimepoints);
  arma::field<arma::mat> F_Sigmae(ntimepoints);
  arma::field<arma::mat> F_Sigmae_inv(ntimepoints);
  arma::field<arma::vec> F_Sigmae_vec(ntimepoints);
  arma::field<arma::vec> F_Sigmae_inv_vec(ntimepoints);
  
  arma::field<arma::vec> F_Ebbi(ntimepoints,nind);
  arma::field<arma::mat> F_Ebbibbit(ntimepoints,nind);
  arma::field<arma::mat> F_Eeeeeeet(ntimepoints,nind);
  
  arma::field<double> F_Eexp(ntimepoints,nind);
  arma::field<arma::vec> F_Eyexp(ntimepoints,nind);
  arma::field<arma::mat> F_Eyyexp(ntimepoints,nind);
  
  arma::field<arma::vec> F_Ebi(ntimepoints,nind);
  arma::field<arma::mat> F_Vbi(ntimepoints,nind);
  arma::field<arma::mat> F_Vbi_half(ntimepoints,nind);
  
  arma::vec beta(ndim,arma::fill::zeros);
  arma::vec beta_old(ndim,arma::fill::zeros);
  arma::vec mub(ndim1,arma::fill::zeros),mub_old(ndim1,arma::fill::zeros);
  arma::mat Sigmab(ndim1,ndim1,arma::fill::eye);
  arma::mat Sigmab_inv(ndim1,ndim1,arma::fill::eye);
  Sigmab=10000.0*Sigmab;
  Sigmab_inv=arma::inv_sympd(Sigmab);
  arma::mat Sigmae(ndim,ndim,arma::fill::eye);//Sigmae=Sigmae*0.25;
  arma::mat Sigmae_inv(ndim,ndim,arma::fill::eye);//Sigmae_inv=Sigmae_inv*4.0;
  arma::vec Sigmae_vec(ndim,arma::fill::ones);
  arma::vec Sigmae_inv_vec(ndim,arma::fill::ones);
  
  arma::vec Ebbi(ndim1,arma::fill::zeros);
  arma::mat Ebbibbit(ndim1,ndim1,arma::fill::eye);
  arma::mat Eeeeeeet(ndim,ndim,arma::fill::eye);
  
  arma::vec Eyexp(ndim,arma::fill::zeros);
  arma::mat Eyyexp(ndim,ndim,arma::fill::zeros);
  
  arma::vec Ebi(ndim1,arma::fill::zeros);
  arma::mat Vbi(ndim1,ndim1,arma::fill::eye);
  arma::mat Vbi_half(ndim1,ndim1,arma::fill::eye);
  
  for(iitt=0;iitt<ntimepoints;iitt++){
    F_mub(iitt)=mub;
    F_Sigmab(iitt)=Sigmab;
    F_Sigmab_inv(iitt)=Sigmab_inv;
    F_Sigmae(iitt)=Sigmae;
    F_Sigmae_inv(iitt)=Sigmae_inv;
    F_Sigmae_vec(iitt)=Sigmae_vec;
    F_Sigmae_inv_vec(iitt)=Sigmae_inv_vec;
  }
  
  arma::vec SS(ndim);
  arma::mat II(ndim,ndim);
  double denominator,expterm;
  arma::vec numerator_vec(ndim);
  arma::mat numerator_mat(ndim,ndim);
  arma::vec avec(ndim1);
  arma::mat amat(ndim1,ndim1);
  
  double X0;//,X1,X2,X3,X4
  double XY0;//,XY1,XY2
  double SUMX0;//,SUMX1,SUMX2,SUMX3,SUMX4
  arma::vec XKSKY(ndim1);
  arma::mat XKSKX(ndim1,ndim1);
  double temp_double;
  int temp_int;
  arma::vec temp_ndim(ndim);
  arma::vec temp_ndim1(ndim1);
  arma::mat temp_ndim_ndim(ndim,ndim);
  arma::mat temp_ndim1_ndim1(ndim1,ndim1);
  arma::mat temp_sum_eee_across_jj(ndim,ndim);
  arma::vec hatmu(ndim),hatee(ndim);
  
  ///--- Create F_XKX
  arma::field<arma::mat> F_XKX(ntimepoints,nind);
  arma::field<arma::vec> F_XKY(ntimepoints,nind);
  arma::mat XKX(ndim1,ndim1);
  arma::vec XKY(ndim1);
  arma::vec F_nobs_in_ker(ntimepoints);
  arma::vec F_nind_at_risk(ntimepoints);
  arma::mat F_sum_of_ker(ntimepoints,ndim,arma::fill::zeros);
  double nobs_in_ker,nind_at_risk,sum_of_ker;
  for(iitt=0;iitt<ntimepoints;iitt++){
    tt=alltimepoints(iitt);
    nobs_in_ker=0.0;
    nind_at_risk=0.0;
    for(ii=0;ii<nind;ii++){
      if(st(ii)>tt)continue;
      if(ot(ii)<tt)continue;
      nind_at_risk=nind_at_risk+1.0;
      XKX.fill(0.0);
      XKY.fill(0.0);
      for(kk=0;kk<ndim;kk++){
        for(jj=0;jj<nobs(ii,kk);jj++){
          ttdiff=ttijk(ii,jj,kk)-tt;
          if(std::abs(ttdiff)>=hh)continue;
          ttdiff_double=double(ttdiff)*omega;
          
          XKX(I0+kk,I0+kk)=XKX(I0+kk,I0+kk)+ker(ttdiff+hh);
          // XKX(II0+kk,I0+kk)=XKX(II0+kk,I0+kk)+tker(ttdiff+hh);
          // XKX(I0+kk,II0+kk)=XKX(I0+kk,II0+kk)+tker(ttdiff+hh);
          // XKX(II0+kk,II0+kk)=XKX(II0+kk,II0+kk)+ttker(ttdiff+hh);
          
          XKY(I0+kk)=XKY(I0+kk)+ker(ttdiff+hh)*yyijk(ii,jj,kk);
          // XKY(II0+kk)=XKY(II0+kk)+tker(ttdiff+hh)*yyijk(ii,jj,kk);
          
          F_sum_of_ker(iitt,kk)=F_sum_of_ker(iitt,kk)+ker(ttdiff+hh);
        }
      }
      F_XKX(iitt,ii)=XKX;
      F_XKY(iitt,ii)=XKY;
    }
    F_nobs_in_ker(iitt)=nobs_in_ker;
    F_nind_at_risk(iitt)=nind_at_risk;
  }
  
  // Rcout<<"2: "<<std::endl;
  
  double sum_for_nind;
  double ee;
  arma::vec sum_for_Ebi(ndim1);
  arma::mat sum_for_Sigmab(ndim1,ndim1);
  arma::vec sum_for_Sigmae(ndim);
  arma::vec temp_sum_XVX_across_jj(ndim);
  arma::vec temp_sum_eeeeK_across_jj(ndim);
  
  for(iter=0;iter<niter;iter++){
    for(iitt=0;iitt<ntimepoints;iitt++){
      tt=alltimepoints(iitt);
      mub=F_mub(iitt);
      Sigmab=F_Sigmab(iitt);
      Sigmab_inv=F_Sigmab_inv(iitt);
      Sigmae=F_Sigmae(iitt);
      Sigmae_inv=F_Sigmae_inv(iitt);
      Sigmae_vec=F_Sigmae_vec(iitt);
      Sigmae_inv_vec=F_Sigmae_inv_vec(iitt);
      
      sum_for_nind=0.0;
      sum_for_Ebi.fill(0.0);
      sum_for_Sigmab.fill(0.0);
      sum_for_Sigmae.fill(0.0);
      
      for(ii=0;ii<nind;ii++){
        if(st(ii)>tt)continue;
        if(ot(ii)<tt)continue;
        
        ///--- Calculate Vbi Ebi
        
        XKSKY=F_XKY(iitt,ii);
        XKSKX=F_XKX(iitt,ii);
        for(kk=0;kk<ndim;kk++){
          XKSKX.row(I0+kk)=XKSKX.row(I0+kk)*Sigmae_inv_vec(kk);
          // XKSKX.row(II0+kk)=XKSKX.row(II0+kk)*Sigmae_inv_vec(kk);
          
          XKSKY(I0+kk)=XKSKY(I0+kk)*Sigmae_inv_vec(kk);
          // XKSKY(II0+kk)=XKSKY(II0+kk)*Sigmae_inv_vec(kk);
        }
        
        Vbi=arma::inv_sympd(Sigmab_inv+XKSKX);
        Ebi=mub+Vbi*(XKSKY-XKSKX*mub);// Ebi=mub+Sigmab*(XKSKY-XKSKX*Vbi*XKSKY);
        F_Vbi(iitt,ii)=Vbi;
        F_Ebi(iitt,ii)=Ebi;
        
        sum_for_nind=sum_for_nind+1.0;
        sum_for_Ebi=sum_for_Ebi+Ebi;
        sum_for_Sigmab=sum_for_Sigmab+Vbi+(Ebi-mub)*(Ebi-mub).t();
        
        temp_sum_XVX_across_jj.fill(0.0);
        temp_sum_eeeeK_across_jj.fill(0.0);
        
        for(kk=0;kk<ndim;kk++){
          for(jj=0;jj<nobs(ii,kk);jj++){
            ttdiff=ttijk(ii,jj,kk)-tt;
            if(std::abs(ttdiff)>=hh)continue;
            ttdiff_double=double(ttdiff)*omega;
            
            ee=yyijk(ii,jj,kk)-Ebi(I0+kk);
            // -ttdiff_double*ttdiff_double*Ebi(III0+kk);//ker(hh+ttdiff)*ee*ee
            // temp_ndim_ndim=ker(hh+ttdiff)*ee*ee.t();
            // ee=ker(hh+ttdiff)*ee*ee;
            temp_double=ker(hh+ttdiff)*ee*ee;
            temp_sum_eeeeK_across_jj(kk)=temp_sum_eeeeK_across_jj(kk)+temp_double;
            
            temp_double=0.0;
            temp_double=temp_double+Vbi(I0+kk,I0+kk)*ker(hh+ttdiff);
            // temp_double=temp_double+Vbi(II0+kk,I0+kk)*tker(hh+ttdiff);
            // temp_double=temp_double+Vbi(I0+kk,II0+kk)*tker(hh+ttdiff);
            // temp_double=temp_double+Vbi(II0+kk,II0+kk)*ttker(hh+ttdiff);
            
            temp_sum_XVX_across_jj(kk)=temp_sum_XVX_across_jj(kk)+temp_double;
          }
        }
        sum_for_Sigmae=sum_for_Sigmae+temp_sum_eeeeK_across_jj+temp_sum_XVX_across_jj;
        // Rcpp::Rcout<<"eeee"<<temp_sum_eeeeK_across_jj.t();
        // Rcpp::Rcout<<"sXVX"<<temp_sum_XVX_across_jj.t();
      }
      
      mub=sum_for_Ebi/F_nind_at_risk(iitt);
      Sigmab=sum_for_Sigmab/F_nind_at_risk(iitt);
      for(kk=0;kk<ndim;kk++){
        Sigmae_vec(kk)=sum_for_Sigmae(kk)/F_sum_of_ker(iitt,kk);
      }
      Sigmae_inv_vec=1.0/Sigmae_vec;
      
      F_mub(iitt)=mub;
      F_Sigmab(iitt)=Sigmab;
      F_Sigmab_inv(iitt)=arma::inv_sympd(Sigmab);
      F_Sigmae(iitt)=Sigmae;
      F_Sigmae_inv(iitt)=arma::inv_sympd(Sigmae);
      F_Sigmae_vec(iitt)=Sigmae_vec;
      F_Sigmae_inv_vec(iitt)=Sigmae_inv_vec;
    }
    
    SS.fill(0.0);
    II.fill(0.0);
    for(ii=0;ii<nind;ii++){
      if(!delta(ii))continue;
      tt=ot(ii);
      avec=F_Ebi(tt-1,ii).subvec(I0,I1);
      SS=SS+avec;
      denominator=0.0;
      numerator_vec.fill(0.0);
      numerator_mat.fill(0.0);
      for(ll=0;ll<nind;ll++){
        if(st(ll)>tt)continue;
        if(ot(ll)<tt)continue;
        avec=F_Ebi(tt-1,ll).subvec(I0,I1);
        expterm=std::exp(arma::dot(avec,beta));
        denominator=denominator+expterm; //F_Eexp(tt-1,ll)
        numerator_vec=numerator_vec+expterm*avec; //F_Eyexp(tt-1,ll);
        numerator_mat=numerator_mat+expterm*avec*avec.t();
      }
      SS=SS-numerator_vec/denominator;
      II=II-numerator_mat/denominator+(numerator_vec/denominator)*(numerator_vec/denominator).t();
    }
    beta_old=beta;
    beta=beta-arma::solve(II,SS);
    all_iter_beta.col(iter)=beta;

    epsilon=arma::sum(arma::abs(beta-beta_old));
    all_iter_beta.col(iter)=beta;
    Rcpp::Rcout<<iter<<": "<<"epsilon"<<epsilon<<"  beta:"<<beta.t();
    if(epsilon<tol)break;
  }
  
  List result(5);
  
  result(0)=F_mub;
  result(1)=F_Sigmab;
  result(2)=F_Sigmae_vec;
  result(3)=F_Ebi;
  result(4)=beta;
  
  return(result);
}
