#ifndef DSIGMA_H
#define DSIGMA_H
#define PROTONS_PER_KTON 6.67e31
double dsigma(double e_nu,double e_e,double costheta);
float enu(float ee,double costheta);
void genedir(double emin, double emax, double *ep, double *cth, double *phi);
double dsigmasv_max(double ep);
float weight_enu(float enu, float ctheta, float low, float bwidth, int nbins, double *nuspectrum, double livetime = 2790.1, double fiducial = 22.5);
float weight_enu_interacted(float enu, float low, float bwidth, int nbins, double *nuspectrum);
float weight_ep(float ep, float low, float bwidth, int nbins, double *pspectrum, double rate = 0.86, double livetime = 2790.1, double fiducial = 22.5);
//void getspec_(double *emin, double *emax, double *ux, double *uy, double *uz, double *en, double *upx, double *upy, double *upz);
#endif
