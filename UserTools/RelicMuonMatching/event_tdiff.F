************************************************************************
*     -------------------------------------
      SUBROUTINE event_tdiff(nevhwsk_cur, it0xsk_cur, nevhwsk_tar, it0xsk_tar, mode, timediff, idiff)
*     -------------------------------------
*
*     (Purpose)
*       Calculate time difference between event 1 and event 2
*       use mode 0 if event 1 comes after event 2, or mode 1 if event 1 is before event 2
*       specifying this is required in case of rollover between the two events
*
*     (Input)
*       nevhwsk_cur   counter32 aka nevhwsk for event 1
*       it0xsk_cur    it0xsk aka t0 for event 1
*       nevhwsk_tar   counter32 aka nevhwsk for event 2
*       it0xsk_tar    it0xsk aka t0 for event 1
*       mode          0 = whether event 2 is before or after current event
*
*     (Output)
*       timediff in ns
*
*     (Creation Date and Author)
*       2023.09.29 MOflaher (modified from tdiff_muon.F)
*
************************************************************************


      implicit none
#include "skhead.h"
#include "skparm.h"
      
      ! default 'integer' is 4-bytes (32-bit)
      ! note that all fortran integers are signed
      integer   nevhwsk_tar, it0xsk_tar, nevhwsk_cur, it0xsk_cur, mode
      ! 8-byte variables (64-bit)
      real*8    timediff
      integer*8 icurrent, itarget, itmp, itmp2, itest
      
      integer*8 idiff
      integer*8 bitone
      parameter (bitone = 1)
      logical take2, notsame
      
      if (SK_GEOMETRY .lt. 4) then
         print *, 'event_tdiff: only for SK-IV'
         stop
      endif
      
 !    define a printing format for 32-bit hex numbers
 100     format(a, z32.32)
      
CM      print *,"current: ",nevhwsk_cur," , ",it0xsk_cur," , target: ",nevhwsk_tar," , ",it0xsk_tar
      
      !print 100, "start t it0xsk_cur: 0x",it0xsk_cur
      ! copy 32-bit it0sk value into a 64-bit variable
      ! we need to use bitwise operation as just using = with negative it0sk values
      ! will fill out the upper 32 bits with 1's (2's complement representation)
      icurrent = iand(Z'FFFFFFFF', it0xsk_cur)
      ! print binary to check
      ! write (*,fmt='(b64.64)') it0xsk_cur
      
      !print 100,"start nevhwsk: 0x",nevhwsk_cur
      ! nevhwsk is never negative so we don't need to do the bitwise copy
      itmp = iand(nevhwsk_cur, not(Z'1FFFF')) ! mask lower 17 bits
      itmp = ishft(itmp, 15)                  ! shift 15 bits up
      icurrent = icurrent + itmp              ! combine with it0xsk
      
      ! repeat for target
      itarget = iand(Z'FFFFFFFF', it0xsk_tar)
      itmp = iand(nevhwsk_tar, not(Z'1FFFF')) ! mask lower 17 bits
      itmp = ishft(itmp, 15)                  ! shift 15 bits up
      itarget = itarget + itmp                ! combine with it0xsk
      
CM      print 100, "combined icurrent: ",icurrent
CM      print *,"as bits:"
CM      write (*,fmt='(b64.64)') icurrent
      
CM      print 100, "combined itarget: ",itarget
CM      print *,"as bits:"
CM      write (*,fmt='(b64.64)') itarget
      
      idiff = icurrent - itarget ! normal case
      if (mode.ne.0) idiff = idiff * (-1)
      if (idiff.lt.0) then 
         idiff = ishft( bitone, +47 ) + idiff ! rollover case
      endif
CM      print 100, "idiff: ",idiff
CM      print *,"as bits:"
CM      write (*,fmt='(b64.64)') idiff
      
c***  convert ticks to ns and fix sign if necessary
      timediff = dfloat(idiff) / count_per_nsec
      if (mode.ne.0) timediff = timediff * (-1)
CM      print *, "converted to ns: ",timediff
CM      print *, "as bits"
CM      write (*,fmt='(b64.64)') timediff
      

 1000 continue
      return
      end
