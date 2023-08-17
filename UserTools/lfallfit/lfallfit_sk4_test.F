cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lfallfit_sk4_test(watert, nhitcut, flag_skip,
     $     flag_log, lfflag)
c
c     (input)
c        watert     water transparency for fit (will be stored in lwatert)
c        nhitcut    max. nqisk for clusfit
c        flag_skip  ==0  clear variables, then do everything
c                   ==1,others  just recalculate energy and related variables
c                   ==2  recalculate energy, clik, msg and related variables      
c        flag_log   LOGLV in findconst
c                   1) Lots of output
c                   2) only prints filenames
c                   3) only prints when not found
c                   4) do not print
c                   5) only prints when found
c     
c     (output)
c        lfflag   status flag (may not be used)
c
c     Note;
c        lfflag = 0  fit went through okay
c               = -1 too many hits in event for lowenergy
c               = -2 .le.0 hits in event
c
c  adapted from lfallfit_sk4_final_qe43.F for local debugging in ToolFramework
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer nhitcut, lfflag, flag_skip, flag_log
      real watert

#include "skhead.h"
#include "skparm.h"
#include "geopmt.h"
#include "sktq.h"
#include "skroot_lowe.h" 

      integer nvalid,listatm(3),nhitatm(3),atmchnum(3),mxhit8,nqrng(3)
      integer i 
      real atmratio(3), rqrng(3), enelf_sk4, wallsk, lfcal_r02

c*** for BONSAI 3.1
      integer clusfit, bonsaifit, ncf, nbf, clnsel, bsnsel 
      external clusfit, bonsaifit

      integer icbl(MAXPM), bshit(MAXPM), clhit(MAXPM)

      integer nhit, ierr
      real qarray(MAXPM),tarray(MAXPM), patliklf_sk4

      real effwallf !, poswal(3)
c*** for SLE
      integer IDATA(100), ndata, iloc, bsn20
      real rdata(100)
      equivalence (idata(1),rdata(1))

c*** for Ariadne
      INTEGER hit
cmi 20150112      REAL    av(4),hits(3,500) 
      REAL    av(3),hits(3,400) 

c*** enelf version
      integer enelf_version
      
      write(*,77) "nqisk is at ",LOC(nqisk)
77    FORMAT(A Z7)
      write(*,fmt="(a,Z7,a,Z7,a,Z7)",advance="yes")
     &         "ihcab is at ",  LOC(ihcab(1)),
     &         ", qisk is at ", LOC(qisk(ihcab(1))),
     &         ",tisk is at ",  LOC(tisk(ihcab(1)))
      print *,"printing tq common block hits: "
      do i = 1, nqisk
         print *,"cable: ",ihcab(i),", Q: ",qisk(ihcab(i)),", T: ",tisk(ihcab(i))
      enddo
      ! note that when accessing from c++ we need to us skq_.qisk[skchnl_.ihcab[i]-1],
      ! to correct for fortran-based indexing of the cable numbers

c*** check SK_GEOMETRY
      if (SK_GEOMETRY .lt. 4) then
         print *,'lfallfit_sk4_final: error SK_GEOM=',SK_GEOMETRY
         stop
      endif

c ----------------------------------------------------------------------
c ---- initialize variables: (both lowe & muon)
c ----------------------------------------------------------------------
      if (flag_skip .eq. 0) call lfclear_all
      nhit = 0
      lwatert = watert

c*** set enelf version
      enelf_version = 7

C-------------------------------------------------------------------
C general variables
C-------------------------------------------------------------------
      if (flag_skip .eq. 0) then

c*** ATM hit parameters
         call skatmmap2(3,nvalid,listatm,nhitatm,atmratio,atmchnum)
         latmnum = atmchnum(1)
         latmh   = nhitatm(1)

c*** mx24
         call lfneihit(mxqisk,mxhit8,lmx24)

c*** time to the previous event
         call tdiff(nt48sk,ltimediff)

c*** NS ratio
         call lfqhit(nqrng,rqrng)
         if (nqrng(1)+nqrng(2)+nqrng(3).gt.0)then
            lnsratio = float(nqrng(2))/float(nqrng(1)+nqrng(2)+nqrng(3))
         else
            lnsratio = 0.
         end if

c*** linfo=1: nflf
         call lfflasher(lmx24, qimxsk, nflf)
c*** linfo=2: lnahit
         call lfnhita(lnahit)
c*** linfo=3: lnqisk
         lnqisk = nqisk
         if (lninfo.lt.3) lninfo = 3 

c*** setup arrays for bonsai
         call lbfset0(qarray,MAXPM)
         call lbfset0(tarray,MAXPM)
         do i = 1, nqisk
            qarray(i) = qisk(ihcab(i))
            tarray(i) = tisk(ihcab(i))
         enddo

      endif

c***  check number of hit
      if (nqisk.gt.NHITCUT) then
         nhit = NHITCUT 
         print *,'=== nqisk.gt.NHITCUT skip clusfit === ',
     $        nrunsk,nevsk,nqisk,NHITCUT
         goto 100 
      else
         nhit = nqisk
      endif
         
C-------------------------------------------------------------------
C do the cluster finding algorithm
C-------------------------------------------------------------------
      ncf = 0
      if (flag_skip .eq. 0) then
         ncf=clusfit(clvertex,clresult,clgoodness,clnsel,
     $        nhit,ihcab,tarray,qarray)
         
         if (ncf .gt. 0) then
            lfflag = 0
            
c***  clusfit N50
            CALL lfnhit(50.,clvertex,cln50,clhit)
         
c*** linfo=9: clwallsk
            clwallsk = wallsk(clvertex)
            if (lninfo.lt.9) lninfo = 9 

            if (clwallsk.gt.0.0) then
c***  clusfit dirks
               CALL lfdir2(clvertex, cldir, cldirks) ! for Sk-I
               if(cldirks .ne. -1.)then
                  call lfdirks(clvertex, cldir, cldirks)
               else
                  cldirks = 2.0
               endif

c*** calc. cossun and direction of sun
               call slsundir3(ndaysk, ntimsk, cldir, lsdir, clcossun, dist)

c*** linfo=4: dist
cc             rlinfo(4) = dist
            endif
         else
            lfflag = -1
            if(nqisk.le.0)lfflag = -2
         endif
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c for non-zero flag_skip do here
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (flag_skip .ne. 0 .or. ncf .gt.0) then
c***clusfit energy 
         call lfneweff_sk4_final_qe43(clvertex, cldir, cleffhit, lwatert, flag_log, 0)
         cleffh = cleffhit(1)
         clenergy = enelf_sk4(cleffh, enelf_version) ! for SK-IV

c*** patlik
         clpatlik = patliklf_sk4(clvertex,cldir,clenergy,ierr) 
c*** linfo=5: cleffwal
         cleffwal = effwallf(1, clvertex, cldir, poswal) 
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         
c*** selection for bonsai fit for data
 100  continue
cc      if (nqisk.ge.18 .and. ltimediff.gt.50000.d0) then
      if (nqisk.ge.10) then

c*************************************************************
ccc store the last T/Q (for debug)
c*************************************************************
c         open(99,file='tq.dat',form='formatted')
c         write(99,*) nrunsk, nsubsk, nevsk, nqisk
c         do i = 1, nqisk
c            write(99,*) i, ihcab(i), tarray(i), qarray(i)
c         enddo
c         close(99)
c*************************************************************

         if (flag_skip .eq. 0) then
c*** apply bonsai fit
            nbf=bonsaifit(bsvertex,bsresult,bsgood,bsnsel,
     $           nhit,ihcab,tarray,qarray)

c*** original bonsai direction (obsolete???)     
c           bsll = bstimedirfit(0, bsvertex, bsdirect, bsdt)

c*** bonsai N50
            CALL lfnhit(50.,bsvertex,bsn50,bshit)

c*** linfo=10: bswallsk
            bswallsk = wallsk(bsvertex)
            if (lninfo.lt.10) lninfo = 10 
         endif

         if (flag_skip .eq. 0 .or. flag_skip .eq. 2) then
            if (bswallsk.gt.0.0) then

c*** bonsai direction from lfdir2() for SK-I
               call lfdir2(bsvertex, bsdir_lfdir2, bsdirks)

c*** linfo=23: bonsai N20raw
               call lfnhit_tisk(20., bsvertex, bsn20raw, icbl, bshit)

c*** linfo=24: bonsai R02
               call lfnhit2(20., bsvertex, bsn20, icbl, bshit)
               if (bsn20.eq.0) then
                  bsr02 = -1.0
               else
                  bsr02 = lfcal_r02(bsn20, icbl, bsvertex)
               endif

c*** linfo=27-33: Ariadne
               amsg=-1
               IF (bsn20 .LE. 390) THEN
                  DO i=1,bsn20
                     if (icbl(i).lt.1 .or. icbl(i).gt.11146) then
                        print *, 'strange icbl: ', icbl(i), i
                        stop
                     endif
                     hits(1,i)=xyzpm(1,icbl(i))
                     hits(2,i)=xyzpm(2,icbl(i))
                     hits(3,i)=xyzpm(3,icbl(i))
                  END DO
                  av(1)=bsvertex(1)
                  av(2)=bsvertex(2)
                  av(3)=bsvertex(3)
c                 print*,  'vertex:',av
                  CALL lfariadne(av,bsn20,hits,adir,amsg,aratio,anscat,acosscat)
               END IF
               if (lninfo.lt.33) lninfo = 33
            endif
         endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c for non-zero flag_skip do here (when nqisk>=10)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (bswallsk.gt.0.0) then

c***  bonsai energy based on lfdir2()
            call lfneweff_sk4_final_qe43(bsvertex, bsdir_lfdir2, bseffhit, lwatert, flag_log, 0)
            bseffh = bseffhit(1)
c*** linfo=14: bsenergy_lfdir2
            bsenergy_lfdir2 = enelf_sk4(bseffh, enelf_version) ! for SK-IV
            if (lninfo.lt.14) lninfo = 14

c*** bonsai direction from lfdir4() for SK-IV 
            call lfdir4(bsvertex, bsenergy_lfdir2, bsdir, bsdirks)

c*** bonsai dirks
            if (bsdirks .ne. -1.) then
               call lfdirks(bsvertex, bsdir, bsdirks)
            else
               bsdirks = 2.0
            endif

c*** calc. cossun and direction of sun
            call slsundir3(ndaysk, ntimsk, bsdir, lsdir, bscossun, dist)

c*** bonsai energy 
            call lfneweff_sk4_final_qe43(bsvertex, bsdir, bseffhit, lwatert, flag_log, 0)
            bseffh = bseffhit(1)
            bsenergy = enelf_sk4(bseffh, enelf_version) ! for SK-IV

c*** patlik
            bspatlik = patliklf_sk4(bsvertex,bsdir,bsenergy,ierr) 

c*** linfo=6: bseffwal
            bseffwal = effwallf(1, bsvertex, bsdir, poswal)
            if (lninfo.lt.36) lninfo = 36  ! for poswal(3)

c*** linfo=25: bonsai clik
            bsclik = float(bsn20raw) / bseffh * bsr02

c*** linfo=26: bonsai ovaq
            bsovaq = bsgood(2)**2-bsdirks**2
            if (lninfo.lt.26) lninfo = 26
         endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      endif

 999  continue

      return
      end

