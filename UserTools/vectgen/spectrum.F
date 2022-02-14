**********************************************************************
*     ----------------------------------------
      SUBROUTINE spectrum(iflag,idate,itime,pos,seed,emin, emax, ranpos)
*     ----------------------------------------
*     
*     ( purpose )
*       VECTOR GENERATION ( SRN NEUTRINO )
*
*     ( input ) 
*       IFLAG : see spcsunsl.F
*       POS(3)   : VERTEX POSITION
*       idate    : run date
*       itime    : run time
*       pos      : vertex position
*       seed     : random seed for spectrum generation
*       emin, emax: energy bin in which to generate events
*
*     ( output ) 
*       COMMON VCWORK,VCVRTX
*
*     ( creation date and author )
*       1995.11.16 ; First version by Y.Koshio
*       1998.11.05 ; modified by Y.Takeuchi
*                         without HEP (only B8) version
*
**********************************************************************
#include "vcwork.h"
#include "vcvrtx.h"

      integer iflag, j, nvct
      real emin, emax, me, rlu
      real pos(3,10),postmp(3)
      real ep, enu, eneutron, mneutron, mproton
      real ct,phi
      real dirpx, dirpy, dirpz
      integer*4 seed
      integer idate(3),itime(3),neutron,ip(3)
      real dirtmp(3), dir(3,3), pmass,pmomabs(3),pmom(3,10)
      real dwall,tvtx(3)
      character*200 specparamfile
      integer ranpos

      call kzeclr

      nvct=3
      dwall=150.0
      mproton=938.272
      mneutron=939.5654
      me = 0.511

      ip(1)= -12
      ip(2)= -11
      ip(3)= 2112
c     If position not preset
      if (ranpos.eq.1) then
          call rnpos(postmp,dwall)
      else
          do j = 1,3
              postmp(j) = pos(j, 1)
          enddo
      endif
      do j = 1,3
          pos(1,j)=postmp(1)
          pos(2,j)=postmp(2)
          pos(3,j)=postmp(3)
          tvtx(j) = 0.
      enddo
      call rndir(dirtmp)
      dir(1,1)=dirtmp(1)
      dir(2,1)=dirtmp(2)
      dir(3,1)=dirtmp(3)
c     Use Strumia-Vissani cross-section to get positron direction
c     and deduce the one for the neutron
      ep=rlu() * (emax - emin) + emin
      ct = rlu() * 2 - 1
      phi = rlu() * 2 * 3.1415926
      call getspec(ep,ct,phi,dir(1,1),dir(2,1),dir(3,1),enu,dir(1,2),dir(2,2),dir(3,2))
      pmomabs(1)=enu
      pmomabs(2)=sqrt(ep*ep-me*me)
      eneutron = enu + mproton - ep
      pmomabs(3)=sqrt(eneutron * eneutron - mneutron * mneutron)
      do j = 1,3
          pmom(j,1) = dir(j,1) * pmomabs(1);
          pmom(j,2) = dir(j,2) * pmomabs(2)
          pmom(j,3) = pmom(j,1) - pmom(j,2)
          dir(j,3) = pmom(j,3)/pmomabs(3);
      enddo

c     write(*,*) '**************************'
c     write(*,*) 'EV, Vertex::',i,pos(1,1),pos(2,1),pos(3,1)
c     write(*,*) '1st neutrino::',pmomabs(1),dir(1,1),dir(2,1),dir(3,1)
c     write(*,*) '1st positron::',pmomabs(2),dir(1,2),dir(2,2),dir(3,2)
c     write(*,*) '1st neutron::',pmomabs(3),dir(1,3),dir(2,3),dir(3,3)
c     write(30,*) pos(1,2), pos(2,2), pos(3,2)


      call VCMPMK(NVCT,POS,PMOM,TVTX,IP)
C
      CALL VCMKVC
      CALL VCMKVX
C
      RETURN
      END
