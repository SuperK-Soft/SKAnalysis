C ok, this seems to be the same routine as kzinit (just with different idim,
C and i'd wager both values are just "arbitrarily large").
C So probably we should just use kzinit.
C
C N.B. There are multiple functions that use this in SKOFL_ROOT and ATMPD_ROOT,
C each of which has their own copy of the sourcefiles (not all the same),
C but none of which are built into any standard libraries we can link against.
C So if we really want to call this, we're going to need to our own copy....
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine zbsinit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none

      integer idim
      parameter (idim=4000000)
      real zbs(idim)
      common/kzbs/zbs

C      CALL KZINIT
      SK_FILE_FORMAT = 0
      call kzlimit(idim)
      return
      end
