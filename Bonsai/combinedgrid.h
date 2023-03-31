#include "searchgrid.h"

#ifndef COMBINEDGRID
#define COMBINEDGRID

class combinedgrid: public fit_param,public searchgrid
{
 public:
  inline combinedgrid(double r,double z,
		      searchgrid *first,searchgrid *second);
};

inline combinedgrid::combinedgrid(double r,double z,
				  searchgrid *first,searchgrid *second):
		    searchgrid(r,z,dwall4hit())
{
  int nset,set,size1,size2,point;
  float *array;

  nset=first->nset();
  if (nset!=second->nset())
    {
      printf("combination grid impossible as the two grids have different number of sets: %d and %d\n",
	     nset,second->nset());
      return;
    }
  //printf("combining %d sets\n",nset);
  for(set=0; set<=nset; set++)
    {
      size1=first->size(set);
      size2=second->size(set);
      //printf("set %d: sizes %d %d\n",set,size1,size2);
      expand_size(size1+size2);
      if (size1>size2) array=new float[3*size1]; else array=new float[3*size2];
      first->copy_points(set,array,3,0,1,2);
      for(point=0; point<size1; point++)
	add_point(array[3*point],array[3*point+1],array[3*point+2]);
      second->copy_points(set,array,3,0,1,2);
      for(point=0; point<size2; point++)
	add_point(array[3*point],array[3*point+1],array[3*point+2]);
      delete [] array;
      if (set<nset) close();
    }
}

#endif
