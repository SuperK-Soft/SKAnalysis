      subroutine bonsai_ini()

      implicit none

#include "skhead.h"
#include "skparm.h"
#include "sktq.h"
#include "skbadc.h"
#include "geopmt.h"
#include "skroot_lowe.h"
#include "skday.h"
#include "skwt.h"

      integer lun
      parameter (lun=20)

      real watert
      real dark, darklf, darks
      
c***  init water transparency table
      call skrunday
      call skwt
c***  in the case of MC, assign dark counts based on run 62428
! gain corr for SK4 only, and run must match.., 85000-88000 for SK6
c      nrunsk = 62428
      nrunsk = 86000
cc*** init. bonsai 
      CALL cfbsinit(MAXPM,xyzpm)
      call lfwater(relapse(nrunsk), watert)
      print *, "nrunsk/watert = ", nrunsk, watert

      return
 
      end


c------------------------------------------------------------------------------

      subroutine bonsai_fit(dt, tisksend, qisksend, cabsend, ntisksend,
     & tenergy, tvx, tvy, tvz, tgood)

      implicit none

#include "skhead.h"
#include "skparm.h"
#include "sktq.h"
#include "skbadc.h"
#include "geopmt.h"
#include "skroot_lowe.h"
#include "skday.h"
#include "skwt.h"

      integer lun

      real tisksend(1000), qisksend(1000)
      integer cabsend(1000)
      integer ntisksend, i, lfflag
      real dt,toff
      real x1,y1,z1, dis1, ttof(2000)
      real watert
      real effwallf, effwal
      real tenergy, tvx, tvy, tvz, tgood, g2d2

c*** number of hit limit for clusfit
      integer NHITCUT
      parameter (NHITCUT =   1000)

c*** branches to skip
      character*100 temp
      integer   nintoskip, nouttoskip
      parameter (nintoskip  = 6)
      parameter (nouttoskip = 11)
      character*100 in_branches_to_skip(nintoskip)
      character*100 out_branches_to_skip(nouttoskip)

      real pawc
      common/pawc/pawc(6000000)
      integer ii, cab, n10, nc, nback, nlow
      real*8  like
      real offset,total,avg,ntot
      integer nhit, nbf, bsnsel, bonsaifit, enelf_version, enelf_sk4
      real qarray(MAXPM),tarray(MAXPM)
      external bonsaifit

c*** loop
         watert = 12431.3
         
         total = 0
         do i = 1, ntisksend
            total = total + tisksend(cabsend(i))
         enddo
         avg = total/(ntisksend)
         offset = dt + 1000 - avg

         do i = 1, ntisksend 
            ihcab(i) = cabsend(i)
            tisk(ihcab(i)) = tisksend(ihcab(i)) - dt + offset
            qisk(ihcab(i)) = qisksend(ihcab(i))
c            print *, "tisk/qisk/cab : " , tisksend(ihcab(i)), tisk(ihcab(i)),qisk(ihcab(i)), ihcab(i)
         end do
         nqisk = ntisksend
         nhit = nqisk

         call lfclear_all
         call lbfset0(qarray,MAXPM)
         call lbfset0(tarray,MAXPM)
c         print *, "Set0"
         do i = 1, nqisk
            qarray(i) = qisk(ihcab(i))
            tarray(i) = tisk(ihcab(i))
         enddo
c         print *, "Ready"
 
         nbf=bonsaifit(bsvertex,bsresult,bsgood,bsnsel,
     $        nhit,ihcab,tarray,qarray)
         call lfneweff_sk6(bsvertex, bsdir, bseffhit, watert, 1)
         enelf_version = 6
c c*** save info
 
c         print *,"bseffhit ", bseffhit(1)
         bseffh=bseffhit(1)
c         tenergy = enelf_sk4(bseffh, enelf_version)
         tenergy = bseffh
         tvx = bsvertex(1)
         tvy = bsvertex(2)
         tvz = bsvertex(3)
         tgood = bsgood(1)
c         print *, "Bonsai Results: ",  tenergy, tvx, tvy, tvz, tgood
      return

 1002 continue
 1001 continue
 10   continue
      stop
      end


      subroutine bonsai_end()

      implicit none

#include "skhead.h"
#include "skparm.h"
#include "sktq.h"
#include "skbadc.h"
#include "geopmt.h"
#include "skroot_lowe.h"
#include "skday.h"
#include "skwt.h"

      call cfbsexit

      return
      end
