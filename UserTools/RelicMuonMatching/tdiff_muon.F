************************************************************************
*     -------------------------------------
      SUBROUTINE tdiff_muon(nevhwsk_tar, it0xsk_tar, mode, timediff)
*     -------------------------------------
*
*     (Purpose)
*       Calculate time difference to a given event
*
*     (Input)
*       nevhwsk_tar     
*       it0xsk_tar      
*       mode          0 = time diff. to previous event
*                     otherwise to following event 
*
*     (Output)
*       timediff
*
*     (Creation Date and Author)
*       2010.12.27 Haibing (modified from tdiff.F)
*
************************************************************************


      implicit none
#include "skhead.h"
#include "skparm.h"
      
      integer   nevhwsk_tar, it0xsk_tar, mode
      real*8    timediff
      integer*8 icurrent, itarget, itmp	


      integer*8 idiff
      integer*8 bitone
      parameter (bitone = 1)

      if (SK_GEOMETRY .lt. 4) then
         print *, 'tdiff_muon: only for SK-IV'
         stop
      endif

      icurrent = ishft(IT0XSK,-16)
      icurrent = ishft(icurrent,+16) 
      icurrent = ior(icurrent, iand(IT0XSK, Z'0000FFFF')) 
      itmp = ishft(nevhwsk,-17) ! +15-32 = -17, to make upper 32bit w/ 15bit shift
      icurrent = icurrent + ishft(itmp, +32) 
         
      itarget = ishft(it0xsk_tar, -16)
      itarget = ishft(itarget, +16)
      itarget = ior(itarget, iand(it0xsk_tar, Z'0000FFFF'))
      itmp = ishft(nevhwsk_tar, -17)
      itarget = itarget + ishft(itmp, +32) 
      
      idiff = icurrent - itarget ! normal case
c      if (mode.ne.0) idiff = idiff * -1
      if (mode.ne.0) idiff = idiff * (-1)
      if (idiff.lt.0) then 
         idiff = ishft( bitone, +47 ) + idiff ! rollover case
c*** for debugging **************************************
         print *, '***** rollover detected! *****'
         print *, 'current t0'
         print 100, 'it0xsk   = 0x', it0xsk
         icurrent = ishft(it0xsk,-16)
         print 100, 'icurrent = 0x', icurrent
         icurrent = ishft(icurrent, +16)
         print 100, 'icurrent = 0x', icurrent
         icurrent = ior(icurrent, iand(IT0XSK, Z'0000FFFF')) 
         print 100, 'icurrent = 0x', icurrent
         itmp = ishft(nevhwsk,-17) 
         print 100, 'itmp     = 0x', itmp
         icurrent = icurrent + ishft(itmp, +32)
         print 100, 'icurrent = 0x', icurrent

         print *, ''
         print *, 'target t0'
         print 100, 'it0xsk_tar = 0x', it0xsk_tar
         icurrent = ishft(it0xsk_tar,-16)
         print 100, 'itarget  = 0x', icurrent
         icurrent = ishft(icurrent, +16)
         print 100, 'itarget  = 0x', icurrent
         icurrent = ior(icurrent, iand(it0xsk_tar, Z'0000FFFF')) 
         print 100, 'itarget  = 0x', icurrent
         itmp = ishft(nevhwsk_tar,-17) 
         print 100, 'itmp     = 0x', itmp
         icurrent = icurrent + ishft(itmp, +32)
         print 100, 'itarget  = 0x', icurrent
 100     format(a, z16.16)
c********************************************************
      endif

c***  modify timediff
      timediff = dfloat(idiff) / count_per_nsec                
c      if (mode.ne.0) timediff = timediff * -1
      if (mode.ne.0) timediff = timediff * (-1)

 1000 continue
      return
      end
