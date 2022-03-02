
#include <fstream>
class EoS;

// this class contains the information about the transport coefficients
// of the fluid: eta/s, zeta/s and the corresponding relaxation times,
// taupi (\tau_\pi) and tauPi (\tau_\Pi)
class TransportCoeff {
 double etaS, zetaS0, etaS0,eps0,ahr,ah,al,rhodecay, T0, D, E, F,zwidth, taupi, tauPi;
 EoS *eos;  // EoS instance is needed optionally for zeta/s parametrization,
            // which depends on the speed of sound
 double zetaSfun(double e, double T);
 double etaSfun(double e, double rho, double T);
public:

 TransportCoeff(double _etaS, double _zetaS, double etaS0, double _eps, double _ahr, double _ah,double _al, double _rhodecay, double _T0, double _D, double _E, double _F,double _zwidth, EoS *_eos, std::string outputdir,double eCrit);

 ~TransportCoeff(){fcells.close();};
 void printZetaT();
 // returns (optionally temperature dependent) eta/s and zeta/s
 void getEta(double e, double rho, double T, double &_etaS, double &_zetaS);
 void saveEta(double e, double rho, double T, double muB, int ix, int iy, int iz, double tau);
 // returns shear and bulk relaxation times
 void outputCell(double e, double rho,double T, double muB, double tau, int ix, int iy, int iz);
 void getTau(double e, double rho, double T, double &_taupi, double &_tauPi);
 // deltapipi, taupipi, lambdapiPi * divided by tau_pi * !
 void getOther(double e, double nb, double nq, double ns,
   double &deltapipi, double &taupipi, double &lambdapiPi, double &phi7) {
  deltapipi = 4./3.;  taupipi = 10./7.;  lambdapiPi = 6./5.;
  phi7 = 9./70./eos->p(e, nb, nq, ns);
  if(std::isinf(phi7)) phi7=0.0;
 }
 void getOtherBulk(double e, double nb, double nq, double ns,
   double &delPiPi, double &lamPipi) {
  delPiPi = 2./3.;  lamPipi = 8./5.*(1. / 3. - eos->cs2(e));
 }
 // isViscous tells whether the fluid is viscous or inviscid
 inline bool isViscous() {
  if (etaS > 0. || zetaS0 > 0.)
   return true;
  else
   return false;
 }
 long double sum_eta_s;
 long double sum_zeta_s;
 long double sum_epsilon;
 long double sum_eta_s_current;
 long double sum_zeta_s_current;
 long double sum_epsilon_current;
 double tau;
 double nx;
 double ny;
 double nz;
 double eCrit;
 std::string OutputDir;
 std::ofstream fcells;
 std::string outcells;
};
