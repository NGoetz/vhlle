#include <iostream>
#include <iomanip>
#include "eos.h"
#include "trancoeff.h"
#include "inc.h"


TransportCoeff::TransportCoeff(double _etaS, double _zetaS, double _etaS0, double _eps0, double _ahr, double _ah,double _al, double _rhodecay, double _T0, double _D, double _E, double _F, EoS *_eos, std::string outputdir, int _nx, int _ny, int _nz, double _eCrit ) {
 etaS = _etaS;
 zetaS0 = _zetaS;
 etaS0=_etaS0;
 eps0=_eps0;
 ah=_ah;
 al=_al;
 ahr=_ahr;
 rhodecay=_rhodecay;
 T0=_T0;
 eos = _eos;
 D=_D;
 E=_E;
 F=_F;
 sum_eta_s=0;
 sum_epsilon=0;
 sum_eta_s_current=0;
 sum_epsilon_current=0;
 tau=0;
 std::string outcells = outputdir;
 outcells.append("/cells.dat");
 fcells.open(outcells.c_str());
 //fcells.close();
 if(etaS0==0 && etaS>0)
  etaS0=etaS;
 nx=_nx;
 ny=_ny;
 nz=_nz;
 eCrit=_eCrit;
}

void TransportCoeff::printZetaT()
{
 std::cout << "------zeta/s(T):\n";
 for(double e=0.1; e<3.0; e+=0.1){
  double T, mub, muq, mus, p;
  eos->eos(e, 0., 0., 0., T, mub, muq, mus, p);
  std::cout << std::setw(14) << T << std::setw(14) << zetaS(e, T) << std::endl;
 }
 std::cout << "---------------:\n";
}

double TransportCoeff::zetaS(double e, double T)
{
 return zetaS0 * (1. / 3. - eos->cs2(e)) / (exp((0.16 - T) / 0.001) + 1.);
}

double TransportCoeff::etaSfun(double e,double rho, double T)
{
 double etaSreturner;
 if(etaS0<0){
   return etaS;
 }
 if(T0>0){
   etaSreturner=etaS0+  ((T>T0) ? ah*(T-T0) :  al*(T-T0));
 }else{
  etaSreturner=std::max(0.0,etaS0+((e>eps0) ? (1.0/(1+D*rho*std::tanh(e-eps0))*(ah*(e-eps0)+ahr*rho*(1+F*(e-eps0))) ): al*(e-eps0)+ahr*rho*(1+E*(e-eps0*rho*rho*rho))));
 }
 return etaSreturner;
}

void TransportCoeff::getEta(double e, double rho, double T, double &_etaS, double &_zetaS) {
 _etaS = etaSfun(e,rho, T);
 _zetaS = zetaS(e,T);
}

void TransportCoeff::saveEta(double e, double rho, double T,double muB, int ix, int iy, int iz, double tau_){
  //if(ix<nx/2.0 +1 && ix>nx/2.0 -2 && iy<ny/2.0 +1 && iy>ny/2.0 -2 &&iz<nz/2.0 +1 && iz>nz/2.0 -2 ){
  if(e>=eCrit ){
    double etaS=etaSfun(e,rho,T);
    sum_eta_s+=etaS*e;
    sum_epsilon+=e;
    if(tau==tau_){
      sum_eta_s_current+=etaS*e;
      sum_epsilon_current+=e;

    }else{
      tau=tau_;
      sum_eta_s_current=etaS*e;
      sum_epsilon_current=e;
    }
    outputCell(e,rho,T,muB,tau_,ix,iy,iz);

  }


}

void TransportCoeff::outputCell(double e, double rho, double T, double muB, double tau, int ix, int iy, int iz) {
    fcells.precision(15);
    fcells << tau <<","<<ix<<","<<iy<<","<<iz<<","<<e<<","<<rho<<","<<T<<","<<muB<<std::endl;
}

void TransportCoeff::getTau(double e, double rho, double T, double &_taupi, double &_tauPi) {
 if (T > 0.) {
  _taupi = 5. / 5.068 * etaSfun(e,rho, T) / T;
  _tauPi = 1. / 5.068 * zetaS(e,T) / (15. * pow(0.33333-eos->cs2(e),2) * T);
 } else {
  _taupi = _tauPi = 0.;
 }
}
