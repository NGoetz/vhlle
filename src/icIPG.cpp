#include <cfloat>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstring>

#include "eos.h"
#include "eoChiral.h"
#include "fld.h"
#include "icIPG.h"
#include "rmn.h"


using namespace std;

IcIPG::IcIPG(Fluid* f, const char* filename, double _tau0, const char* setup, int _momentum_aniso) {
 cout << "loading IPG\n";
 nx = f->getNX();
 ny = f->getNY();
 nz = f->getNZ();
 dx = f->getDx();
 dy = f->getDy();
 dz = f->getDz();
 xmin = f->getX(0);
 xmax = f->getX(nx - 1);
 ymin = f->getY(0);
 ymax = f->getY(ny - 1);
 zmin = f->getZ(0);
 zmax = f->getZ(nz - 1);

 tau0 = _tau0;

 sNN = 200;
 eta0 = 1.5; // midrapidity plateau
 sigEta = 1.4; // diffuseness of rapidity profile
 etaM = 3.36;
 ybeam = 5.36; // beam rapidity
 alphaMix = 0.145; // WN/binary mixing
 Rg = 0.4; // Gaussian smearing in transverse dir
  //sNorm = 0.56; // normalization of initial entropy profile
 A = 0.0 ; // 5e-4; // initial shear flow
 nsigma = 0.6;
 neta0 = 1.4;
 cout << "IPG: setup for 200 GeV RHIC\n";
 momentum_aniso=_momentum_aniso;
 

 nsmoothx = (int)(3.0 * Rg / dx);  // smoothly distribute to +- this many cells
 nsmoothy = nsmoothx;

 rho = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  rho[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   rho[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    rho[ix][iy][iz] = 0.0;
   }
  }
 }
 utau = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  utau[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   utau[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    utau[ix][iy][iz] = 0.0;
   }
  }
 }
 ux = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  ux[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   ux[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    ux[ix][iy][iz] = 0.0;
   }
  }
 }
 uy = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  uy[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   uy[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    uy[ix][iy][iz] = 0.0;
   }
  }
 }
 ueta = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  ueta[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   ueta[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    ueta[ix][iy][iz] = 0.0;
   }
  }
 }

 nrho = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  nrho[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   nrho[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    nrho[ix][iy][iz] = 0.0;
   }
  }
 }

 // ---- read the events
 nevents = 0;
 ifstream fin(filename);
 if (!fin.good()) {
  cout << "I/O error with " << filename << endl;
  exit(1);
 }
 int npart = 0;  // number of participants
 string line;
 istringstream instream;
 while (!fin.eof()) {
  for (int i = 0; i < 3; i++) { // read first 3 lines - the 3rd line contains npart
   getline(fin, line);
  }
  if (line == "") break;
  instream.str(line);
  instream.seekg(9);
  instream >> npart; // read npart
  for (int i = 0; i < 5; i++) { // read the rest 5 lines from header
   getline(fin, line);
  }
  getline(fin, line); // read grid-step
  getline(fin, line); // read grid-nsteps
  instream.str(line);
  instream.seekg(15);
  instream >> n_grid;

  // allocate source array
  source = new double*[n_grid];
  for (int ix = 0; ix < n_grid; ix++) {
   source[ix] = new double[n_grid];
   for (int iy = 0; iy < n_grid; iy++) {
    source[ix][iy] = 0.0;
   }
  }
  source_utau = new double*[n_grid];
  for (int ix = 0; ix < n_grid; ix++) {
   source_utau[ix] = new double[n_grid];
   for (int iy = 0; iy < n_grid; iy++) {
    source_utau[ix][iy] = 0.0;
   }
  }
  source_ux = new double*[n_grid];
  for (int ix = 0; ix < n_grid; ix++) {
   source_ux[ix] = new double[n_grid];
   for (int iy = 0; iy < n_grid; iy++) {
    source_ux[ix][iy] = 0.0;
   }
  }
  source_uy = new double*[n_grid];
  for (int ix = 0; ix < n_grid; ix++) {
   source_uy[ix] = new double[n_grid];
   for (int iy = 0; iy < n_grid; iy++) {
    source_uy[ix][iy] = 0.0;
   }
  }
  source_ueta = new double*[n_grid];
  for (int ix = 0; ix < n_grid; ix++) {
   source_ueta[ix] = new double[n_grid];
   for (int iy = 0; iy < n_grid; iy++) {
    source_ueta[ix][iy] = 0.0;
   }
  }

  getline(fin, line); // read grid-max
  instream.str(line);
  instream.seekg(12);
  instream >> xmaxG;
  ymaxG = xmaxG;
  xminG = -xmaxG;
  yminG = -xmaxG;
  if (nevents == 0) cout << "IPG grid: x,ymaxG = " << xmaxG << "  n_grid = " << n_grid << endl;
  for (int iy = 0; iy < n_grid; iy++) {
   for (int ix = 0; ix < n_grid; ix++) {
    fin >> source[ix][iy];
   }
  }
  for (int iy = 0; iy < n_grid; iy++) {
   for (int ix = 0; ix < n_grid; ix++) {
    fin >> source_utau[ix][iy];
   }
  }
  for (int iy = 0; iy < n_grid; iy++) {
   for (int ix = 0; ix < n_grid; ix++) {
    fin >> source_ux[ix][iy];
   }
  }
  for (int iy = 0; iy < n_grid; iy++) {
   for (int ix = 0; ix < n_grid; ix++) {
    fin >> source_uy[ix][iy];
   }
  }
  for (int iy = 0; iy < n_grid; iy++) {
   for (int ix = 0; ix < n_grid; ix++) {
    fin >> source_ueta[ix][iy];
   }
  }
  nevents += 1;
  cout << "event " << nevents << "  npart = " << npart << "\n";
  makeSmoothTable(npart);

  // delete source array
  for (int ix = 0; ix < n_grid; ix++) {
   delete[] source[ix];
   delete[] source_utau[ix];
   delete[] source_ux[ix];
   delete[] source_uy[ix];
   delete[] source_ueta[ix];
  }
  delete[] source;
  delete[] source_utau;
  delete[] source_ux;
  delete[] source_uy;
  delete[] source_ueta;
  getline(fin, line);
 }

 // autocalculation of sNorm and nNorm
 sNorm = 1.0;
 double old_sNorm = 0.0;
 do {
   old_sNorm = sNorm;
   sNorm = pow(setNormalization(npart), 0.75)*old_sNorm;
 } while (abs(sNorm-old_sNorm) > 0.0001);
 cout << "sNorm set to " << sNorm << endl;
 double old_nNorm = 0.0;
}

IcIPG::~IcIPG() {
 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   delete[] rho[ix][iy];
   delete[] nrho[ix][iy];
   delete[] utau[ix][iy];
   delete[] ux[ix][iy];
   delete[] uy[ix][iy];
   delete[] ueta[ix][iy];
  }
  delete[] rho[ix];
  delete[] nrho[ix];
  delete[] utau[ix];
  delete[] ux[ix];
  delete[] uy[ix];
  delete[] ueta[ix];
 }
 delete[] rho;
 delete[] nrho;
 delete[] utau;
 delete[] ux;
 delete[] uy;
 delete[] ueta;
}

double IcIPG::interpolateGrid(double x, double y, double** grid) {
 const double dxG = (xmaxG - xminG) / (n_grid - 1);
 const double dyG = (ymaxG - yminG) / (n_grid - 1);
 int ix = (int)((x - xminG) / dxG);
 int iy = (int)((y - yminG) / dyG);
 if (ix < 0) ix = 0;
 if (iy < 0) iy = 0;
 if (ix > n_grid - 2) ix = n_grid - 2;
 if (iy > n_grid - 2) iy = n_grid - 2;
 const double xm = x - xminG - ix * dxG;
 const double ym = y - yminG - iy * dyG;
 double wx[2] = {1. - xm / dxG, xm / dxG};
 double wy[2] = {1. - ym / dyG, ym / dyG};
 double return_val = 0.;
 for (int jx = 0; jx < 2; jx++)
  for (int jy = 0; jy < 2; jy++) {
   return_val += wx[jx] * wy[jy] * grid[ix + jx][iy + jy];
  }
 return_val = std::max(return_val, 0.);
 return return_val;
}

void IcIPG::makeSmoothTable(int npart) {
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    // longidudinal profile here
    const double x = xmin + ix * dx;
    const double y = ymin + iy * dy;
    const double eta = zmin + iz * dz;
    double fEta = 0.;
    if(fabs(eta)<eta0) fEta = 1.0;
    else if (fabs(eta)<ybeam) fEta = exp(-0.5*pow((fabs(eta)-eta0)/sigEta,2));
    rho[ix][iy][iz] += interpolateGrid(x,y, source) * fEta;
    utau[ix][iy][iz] += interpolateGrid(x,y, source_utau) * fEta;
    ux[ix][iy][iz] += interpolateGrid(x,y, source_ux) * fEta;
    uy[ix][iy][iz] += interpolateGrid(x,y, source_uy) * fEta;
    ueta[ix][iy][iz] += interpolateGrid(x,y, source_ueta) * fEta;
 } // Z(eta) loop
}

void IcIPG::setIC(Fluid* f, EoS* eos) {
 double E = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0, Nb = 0.0, S = 0.0;
 double Jy0 = 0.0, Jint1 = 0.0, Jint3 = 0.0, Xcm = 0.0, Ycm = 0.0, Zcm = 0.0;
 double E_midrap = 0.0, Jy0_midrap = 0.0;  // same quantity at midrapidity
 double Tcm = 0.0;
 double e, p, nb,utau_val,ux_val,uy_val,ueta_val;
 double total_energy = 0.0;
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    e = sNorm * rho[ix][iy][iz] / nevents / dx / dy;
    nb = 0.;
    p = eos->p(e, 0., 0., 0.);
    Cell* c = f->getCell(ix, iy, iz);
    utau_val=sNorm * utau[ix][iy][iz]/ nevents / dx / dy;
    ux_val=sNorm * ux[ix][iy][iz]/ nevents / dx / dy;
    uy_val=sNorm * uy[ix][iy][iz]/ nevents / dx / dy;
    ueta_val=sNorm * ueta[ix][iy][iz]/ nevents / dx / dy;
    if(momentum_aniso==1){
         ueta_val=sqrt(ueta_val*ueta_val+ux_val*ux_val+uy_val*uy_val);
         uy_val=0;
         ux_val=0;
    }
    double u[4] = {utau_val,ux_val,uy_val,ueta_val};
    c->setPrimVar(eos, tau0, e, nb, 0.4*nb, 0.,u[1]/u[0], u[2]/u[0], u[3]/u[0]);
    if (e > 0.) c->setAllM(1.);
    double eta = zmin + iz * dz;
    double coshEta = cosh(eta);
    double sinhEta = sinh(eta);
    total_energy += tau0*e*dx*dy*dz*coshEta;
    double u0lab = u[0] * coshEta + u[3] * sinhEta;
    double uzlab = u[0] * sinhEta + u[3] * coshEta;
    double dE = tau0 * ((e + p) * u[0] * (u[0] * cosh(eta) + u[3] * sinh(eta)) -p * cosh(eta)) *dx * dy * dz;;
    // if(zmin+iz*dz>0.)
    E+=dE;
    Pz +=  tau0 * ((e + p) * u[0] * (u[0] * sinh(eta) + u[3] * cosh(eta)) -p * sinh(eta)) *dx * dy * dz;
    Px += tau0 * (e + p) * u[1] * u[0] * dx * dy * dz;
    Py += tau0 * (e + p) * u[2] * u[0] * dx * dy * dz;
    Nb += nb*tau0*dx*dy*dz;
    S += tau0 * eos->s(e, 0., 0., 0.) * u[0] * dx * dy * dz;
    // angular momentum calculation
    const double t = tau0 * coshEta;
    const double z = tau0 * sinhEta;
    const double x = xmin + ix * dx;
    const double y = ymin + iy * dy;
    Xcm += x * dE;
    Ycm += y * dE;
    Zcm += z * dE;
    Tcm += t * dE;
    Jy0 +=
        tau0 * (e + p) * u[0] * (z * u[1] - x * uzlab) * dx * dy * dz * gevtofm;
    Jint1 += tau0 * (e + p) * u[0] * u[1] * dx * dy * dz * gevtofm;
    Jint3 += tau0 * (e + p) * u[0] * uzlab * dx * dy * dz * gevtofm;
    if (iz > nz / 2 - 2 && iz < nz / 2 + 2) {
     E_midrap += dE;
     Jy0_midrap += tau0 * (e + p) * u[0] * (z * u[1] - x * uzlab) * dx * dy *
                   dz * gevtofm;
    }
   }
 Xcm = Xcm / E;
 Ycm = Ycm / E;
 Zcm = Zcm / E;
 Tcm = Tcm / E;
 double Jy = Jy0 - Zcm * Jint1 + Xcm * Jint3;
 cout << "hydrodynamic E = " << E << "  Pz = " << Pz << "  Nbar = " << Nb
      << endl
      << "  Px = " << Px << "  Py = " << Py << endl;
 cout << "initial_entropy S_ini = " << S << endl;
 cout << "Xcm: " << sqrt(Tcm * Tcm - Zcm * Zcm) << "  " << Xcm << "  " << Ycm
      << "  " << 0.5 * log((Tcm + Zcm) / (Tcm - Zcm)) << endl;
 cout << "initial/corrected J_y  " << Jy0 << " " << Jy << endl;
 cout << "J_to_analyze " << setw(14) << E << setw(14) << Nb << setw(14) << Jy
      << endl;
 cout << "1/tau*dE/dy_ini: : " << E_midrap/(3.0*dz*tau0) << endl;
 cout << "1/tau*dJ/dy_ini: " << Jy0_midrap/(3.0*dz*tau0) << endl;
 //exit(1);
}

double IcIPG::setNormalization(int npart) {
 double e;
 double total_energy = 0.0;
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    e = sNorm * rho[ix][iy][iz] / nevents / dx / dy;
    double eta = zmin + iz * dz;
    double coshEta = cosh(eta);
    total_energy += tau0*e*dx*dy*dz*coshEta;
   }
 return npart*0.5*sNN/total_energy;
}

