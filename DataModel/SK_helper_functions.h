#include <stdint.h>
#include "fortran_routines.h"

class TVector3;

void set_rflist_zbs( int lun, const char *filename, bool write );

float GetDWallInDirection(TVector3 vtx, TVector3 dir);
float GetDWall(TVector3 vtx);
