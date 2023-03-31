#include "pairlikelihood.h"

// *************************************************************
// * calculate likelihood with and without timing constraint   *
// *************************************************************
float pairlikelihood::quality(float *vertex)
{
  verfit[0]=vertex[0];
  verfit[1]=vertex[1];
  verfit[2]=vertex[2];
  // calculate likelihood for delayed event
  event_hits=delayed;
  like[2]=fittime(1,verfit,dirfit+5,dt);
  verfit[4]=verfit[3]; // save emission time for delayed event
  event_hits=prompt;
  // calculate likelihood for prompt event
  like[0]=fittime(1,verfit,dirfit,dt);
  // apply Cherenkov angle correction for both events
  devp=dirfit[3]-cang0;
  devd=dirfit[8]-cang0;
  if (devp>0)
    if (plusdang==FIT_PARAM_NONE)
      like[1]=like[0];
    else
      like[1]=like[0]-devp*devp*plusdang;
  else
    if (minusdang==FIT_PARAM_NONE)
      like[1]=like[0];
    else
      like[1]=like[0]-devp*devp*minusdang;
  if (devd>0)
    {
      if (plusdang==FIT_PARAM_NONE)
	like[3]=like[2];
      else
	like[3]=like[2]-devd*devd*plusdang;
    }
  else
    {
      if (minusdang==FIT_PARAM_NONE)
	like[3]=like[2];
      else
	like[3]=like[2]-devd*devd*minusdang;
    }
  if (nlike++==0) set_worst(like[1]+like[3]); else check_worst(like[1]+like[3]);
  //printf("%9.3lf\n",like);
  return(like[1]+like[3]);
}
