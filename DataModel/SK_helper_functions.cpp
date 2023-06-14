#include "SK_helper_functions.h"
#include "Calculator.h"
#include "geotnkC.h"

// helper function
void set_rflist_zbs( int& lun, const char *filename, bool write ) {
	std::string file = filename;
	std::string s1 = "LOCAL";
	std::string s2 = " ";
	std::string s3 = "RED";
	std::string s4 = " ";
	std::string s5 = " ";
	std::string s6 = "recl=5670 status=old"; // XXX
	std::string s7 = " ";
	std::string s8 = " ";
	if( write ){
	s3 = "WRT";
	s6 = "recl=5670 status=new";
	}
	set_rflist_( &lun, file.data(),
	             s1.data(), s2.data(), s3.data(), s4.data(),
	             s5.data(), s6.data(), s7.data(), s8.data(),
	             file.length(),
	             s1.length(), s2.length(), s3.length(), s4.length(),
	             s5.length(), s6.length(), s7.length(), s8.length() );
}
/* for argument meanings:
  set_rflist_(lun, fname, unit, lgname, access,
              stg, multi, comment, hname, mode,
              strlen(fname),strlen(unit),strlen(lgname),
              strlen(access),strlen(stg),strlen(multi),
              strlen(comment),strlen(hname),strlen(mode));
*/

// s6 may also need to be "form=unformatted" or "form=formatted",
// optionally also with "status=new" or "status=old" or "status=unknown"

/*
   we can also export the environmental variable 'RFLIST' to the path of a file
   within which we can specify any number of files, along with their associated
   LUNs to assign and format strings, in the form of the set_rflist argument list
   e.g.
   export RFLIST=rflist.$$.`hostname`
   cat RFLIST:
   -------------
   10{{"/path/to/file_1",LOCAL,,RED,,,"form=unformatted",}}
   11{{"/path/to/file_2",LOCAL,,RED,,,"recl=5670 status=old",}}
   12{{"/path/to/file_3",LOCAL,,WRT,,,"recl=5670 status=new",}}
   -------------
   we then pass the corresponding index to skopenf:
     int LUN;            << this will be set according to the value in the RFLIST file
     int fileIndex = 1;  << index in RFLIST file
     char ftype='f';     << fortran file
     int ihndl=1;        << ?? not sure this is used for fortran files
     skopenf_(LUN, fileIndex, "f", get_ok, ihndl);
     if(LUN==-1) std::cerr<<"failed to open RFLIST file "<<fileIndex<<std::endl;
     else        std::cout<<"opened RFLIST file "<<fileIndex<<" with LUN "<<LUN<<std::endl;
*/


// These got moved here from Calculator.cpp because they depend on sk functions and constants.
// but that results in a situation where we can't build root dictionaries for the NTag classes
float GetDWall(TVector3 vtx)
{
    float vertex[3] = {(float)vtx.x(), (float)vtx.y(), (float)vtx.z()};
    return wallsk_(vertex);
}


float GetDWallInDirection(TVector3 vtx, TVector3 dir)
{
    dir = dir.Unit();

    float dot = vtx.Dot(dir) - vtx.z()*dir.z();
    float dirSq = dir.Perp2(); float vtxSq = vtx.Perp2();

    // Calculate distance to barrel and distance to top/bottom
    float distR = (-dot + sqrt(dot*dot + dirSq*(RINTK*RINTK - vtxSq))) / dirSq;
    float distZ = dir.z() > 0 ? (ZPINTK-vtx.z())/dir.z() : (ZMINTK-vtx.z())/dir.z();

    // Return the smaller
    return distR < distZ ? distR : distZ;
}
