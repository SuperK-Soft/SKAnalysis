#ifndef PAIRLIKELIHOOD
#define PAIRLIKELIHOOD
#include "timefit.h"
#include "fit_param.h"
#include "fitquality.h"
#include "bonsaifit.h"
#include "searchgrid.h"
#include "plato.h"

// *************************************************************
// * Defines BONSAI pairlikelihood and maximization                *
// *************************************************************
class pairlikelihood: public fit_param,public fitquality, public timefit
{
  hitsel       *prompt,*delayed;
  float        cang0,plusdang,minusdang,devp,devd; // direction fit constraint
  float        verfit[5],dirfit[10],like[4];   // vertex, direction, pairlikelihoods
                                               // like[0] prompt pure timing likelihood
                                               // like[1] prompt likelihood
                                               // like[2] delayed pure timing likelihood
                                               // like[3] delayed likelihood
  float        dt;                             // timefit uncertainty
  short int    nlike;                          // number of fit vertices
  dodecahedron dod;                            // defines how surrounding points
  axes         orientation;                    // are calculated

 public:
  inline pairlikelihood(float r,float z);            // construct timefit and pairlikelihood
  inline void set_hits(hitsel *p,hitsel *d);        // define hits to fit
         virtual float    quality(float *vertex);// calculate pairlikelihood
  // surround a point with test vertices
  inline virtual void check_around(float *vertex,float *result,
				   float radius,float *q,int &max);
  inline virtual char ncheck(void);              // return # of surrounding test vertices
  // interpolate best fit position
  inline virtual void interpolate(float *vertex,
				  float radius,float *q,float *inter);
  inline virtual int nresult(void);              // return size of result array
  inline virtual void get_result(float *r);      // get result array
  inline virtual void set_result(float *r);      // set result array
  // maximizing procedure
  inline void      maximize(bonsaifit *fit,searchgrid *grid, bool useAngle=true);
  inline void      maximize(bonsaifit *fit,float *point, bool useAngle=true);
  inline float     get_zero(int ev);             // get t0
  inline float     residual(short int hit);      // get time residual
  inline int       nwind(int ev,float *vertex,float tmin,float tmax);
  inline float     ntgood(int ev,float *vertex,float bgrate,float &guncor);
  inline float     get_ll(int ev);               // get pairlikelihood
  inline float     get_ll0(int ev);              // get pairlikelihood w/o angle constraint
  inline void      get_dir(float *dir);          // get direction fit
  inline int       nfit(void);                   // return # of test vertices
};

#include "pairlikelihood.inline"
#endif
