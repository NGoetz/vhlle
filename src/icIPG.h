
class Fluid;
class EoS;

class IcIPG {
private:
 double eta0; // midrapidity plateau
 double sigEta; // diffuseness of rapidity profile
 double ybeam; // beam rapidity
 double alphaMix; // WN/binary mixing
 double Rg; // Gaussian smearing in transverse dir
 double sNorm; // normalization of initial entropy profile
 double nNorm; // normalization of baryon density
 double sNN; // energy
 double nsigma; // width of gaussian for baryon density
 double neta0; // mean of gaussian for baryon density
 double etaM;
 double A; // initial shear flow constant
 int nx, ny, nz, nevents;
 int n_grid;
 double xminG, xmaxG, yminG, ymaxG; // grid sizes
 double xmin, xmax, ymin, ymax, zmin, zmax;
 double dx, dy, dz;
 double ***rho;
 double ***utau;
 double ***ux;
 double ***uy;
 double ***ueta;
 double ***nrho;
 double **source;
 double **source_utau;
 double **source_ux;
 double **source_uy;
 double **source_ueta;
 static const int NP = 10000;  // dimension for particle arrays
 // auxiliary particle arrays

 double tau0;
 int nsmoothx;  // smoothly distribute to +- this many cells
 int nsmoothy;
 int nsmoothz;
 void makeSmoothTable(int npart);

public:
 IcIPG(Fluid *f, const char *filename, double tau0, const char* setup);
 ~IcIPG();
 double interpolateGrid(double x, double y,double** grid);
 void setIC(Fluid *f, EoS *eos);
 double setNormalization(int npart);
};
