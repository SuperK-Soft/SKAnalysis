************************************************************************
*     -------------------------------------
      SUBROUTINE SET_TIMING_GATE_M(it0)
*     -------------------------------------
*
*     (Purpose)
*        set timing gate for TISK, QISK, TASK, QASK, etc.
*        set IT0XSK
*
*     (Input)
*        t0:  time 0 in tiskz/taskz (count)
*
*     (Output)
*        written to common 
*       
*     (Return)
*        none:
*       
*     (Creation Date and Author)
*       2008.07.25 ; First version by Y.Takeuchi
*
************************************************************************
      implicit none

#include "skparm.h"
#include "skhead.h"
#include "sktq.h"

#include "zbs.h"

      integer it0, i, icab

c*** for ZBS
      INTEGER IPBANK, NDATA, IBUF(10), LDAT, NENT
      REAL RBUF(10)
      EQUIVALENCE(IBUF, RBUF)
 
      INTEGER ISZHDZ
      PARAMETER(ISZHDZ=7) 
      INTEGER IBUFZ(3,MAXPM*30)
      REAL    RBUFZ(3,MAXPM*30)
      EQUIVALENCE (IBUFZ(1,1),RBUFZ(1,1))

c*** -300ns  +1000ns gate
c*** -300 * 1.92 = -576 count
c*** +1000 * 1.92 = +1920 count

c*** set timing gates
      call  SET_TIMING_GATE_ID_M(it0)
      call  SET_TIMING_GATE_OD_M(it0)
      call  SET_TIMING_GATE_VETO_M(it0) ! should be called after ID

c*** set IT0XSK        
      IT0XSK = it0
      NTRIGSK= (it0xsk-it0sk)/count_per_nsec/10 ! 20090403 added by y.t.

c*** store results in tqreal bank/branch
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (SK_FILE_FORMAT .eq. 0 ) then ! zbs
cccccccccccccccccccccccccccccccccccccccccccccccccccc

c*** refill tqreal bank, if exist
         CALL KZBLOC('TQREAL',IPBANK)
         IF(IPBANK.LE.0) GOTO 100
         CALL KZGET1('TQREAL',0,NDATA,IBUF)
         if(ibuf(3).lt.2) then
            print *, 'set_timing_gate: cannot be used for TQREAL ver.1'
            stop
         endif

c*** replace IT0XSK         
         IBUF(7) = IT0XSK
         CALL KZREP1('TQREAL',0,ndata,ibuf)

c*** refill common (is this ok?? 20080728 y.t.)
         IF (NQISKZ.gt.0) THEN
            do i = 1, nqisk_raw
               IBUFZ(1,i) = icabbf_raw(i)
               RBUFZ(2,i) = tbuf_raw(i)
               RBUFZ(3,i) = qbuf_raw(i)
            enddo
            CALL KZLDAT('TQREAL',LDAT)
            nent = 3*nqisk_raw
            CALL LBFCOPY(IBUFZ,IZBS(LDAT+ISZHDZ+1),nent)
         END IF

 100     continue
c*** refill tqareal bank, if exist
         CALL KZBLOC('TQAREAL',IPBANK)
         IF(IPBANK.LE.0) GOTO 200
         CALL KZGET1('TQAREAL',0,NDATA,IBUF)
         if(ibuf(3).lt.2) then
            print *, 'set_timing_gate: cannot be used for TQAREAL ver.1'
            stop
         endif

c*** replace IT0XSK         
         IBUF(7) = IT0XSK
         CALL KZREP1('TQAREAL',0,ndata,ibuf)

c*** refill common (is this ok?? 20080728 y.t.)
         IF (NHITAZ_raw.gt.0) THEN
            do i = 1, NHITAZ_raw
               IBUFZ(1,i) = icabaz_raw(i)
               RBUFZ(2,i) = TASKZ_raw(i)
               RBUFZ(3,i) = QASKZ_raw(i)
            enddo
            CALL KZLDAT('TQAREAL',LDAT)
            nent = 3*nhitaz_raw
            CALL LBFCOPY(IBUFZ,IZBS(LDAT+ISZHDZ+1),nent)
         END IF

 200     continue

      endif

cc      print *, 'set_timing_gate: ',it0sk,it0,id_start,id_end,od_start,od_end

      return
      end

************************************************************************
*     -------------------------------------
      SUBROUTINE SET_TIMING_GATE_ID_M(it0)
*     -------------------------------------
*
*     (Purpose)
*        set timing gate for TISK, QISK, etc.
*
*     (Input)
*        t0:  time 0 in tiskz/taskz (count)
*
*     (Output)
*        written to common 
*       
*     (Return)
*        none:
*       
*     (Creation Date and Author)
*       2008.07.25 ; First version by Y.Takeuchi
*
************************************************************************
      implicit none

#include "skparm.h"
#include "skhead.h"
#include "sktq.h"

      integer it0, id_start, id_end, i
      integer p0,p1,p2,p3,p4,p5
      real tlowest, thighest

      p0=0
      p1=0
      p2=0
      p3=0
      p4=0
      p5=0
      tlowest=999999999
      thighest=-99999999
      
c*** For ID (moved from tqskz.F)
      if (nrunsk.lt.50791) then
         id_start =  1424 + it0 - it0sk
         id_end   =  3920 + it0 - it0sk
      else if (nrunsk.eq.52026) then
         id_start = -1960 + it0 - it0sk
         id_end   =   536 + it0 - it0sk
      else
         id_start = SKGATE_START_COUNT + it0 - it0sk
         id_end   = SKGATE_END_COUNT   + it0 - it0sk
      endif
      ! convert to nsec since we don't have itiskz
      id_start = id_start/COUNT_PER_NSEC
      id_end = id_end/COUNT_PER_NSEC
      print *, "gate: ",id_start," -> ",id_end

      if (nqisk_raw.gt.0) then
c*** clear iab(1:MAXPM) to choose only the earliest hits
         do i = 1,maxpm
            iab(i) = 0
         enddo

c*** bad channel mask will be removed         
         nqiskz = nqisk_raw
         do i = 1, nqisk_raw
            !itiskz(i)  = itiskz_raw(i)   ! not available
            !iqiskz(i)  = iqiskz_raw(i)   ! not available
            ! for processed files itiskz_raw is not available, but then
            ! itiskz is also not normally populated, so this seems unneeded
            tiskz(i)   = tbuf_raw(i) -(it0-it0sk)/count_per_nsec ! ok
            if(tiskz(i).lt.tlowest) tlowest = tiskz(i)
            if(tiskz(i).gt.thighest) thighest = tiskz(i)
            qiskz(i)   = qbuf_raw(i) ! ok
            ICABIZ(i)  = iand(icabbf_raw(i), 2**16-1) ! ok
            IHTIFLZ(i) = ishft(icabbf_raw(i), -16) ! ok
            if (icabiz(i).ge.1.and.icabiz(i).le.MAXPM) then
               IHTIFLZ(i) = iand(IHTIFLZ(i),2**16-2) ! make 0bit = 0
               p0 = p0+1
	           if (iand(ihtiflz(i),2).eq.2) then
	           		p1 = p1+1
	           		if(iab(icabiz(i)).eq.0) then
	           			p2 = p2+1
	           			if(tbuf_raw(i).ge. id_start) then
	           				p3 = p3+1
	           				if(tbuf_raw(i).le.id_end) then
	           					p4 = p4+1
	           				endif
	           			endif
	           		endif
	           endif
               if(iand(IHTIFLZ(i),2).eq.2 .and.
     $              iab(icabiz(i)).eq.0.and.
     $              tbuf_raw(i) .ge. id_start .and.   ! changed itiskz -> tbuf_raw
     $              tbuf_raw(i) .le. id_end) then     ! (tiskz is shifted by t0_sub)
c***  In 1.3usec & in-gate
                  IHTIFLZ(i) = IHTIFLZ(i) + 1 ! 0bit
                  p5 = p5+1
                  iab(icabiz(i)) = 1
               endif
               icabbf_raw(i) = ICABIZ(i) + (IHTIFLZ(i) * 2**16)
            endif
         enddo
      endif

      print *, 'set_timing_gate_id: ',it0sk,it0,id_start,id_end
      print *, 'tlowest,highest:',tlowest,thighest
      print *,'passing hits:',p0,p1,p2,p3,p4,p5

      return
      end

************************************************************************
*     -------------------------------------
      SUBROUTINE SET_TIMING_GATE_OD_M(it0)
*     -------------------------------------
*
*     (Purpose)
*        set timing gate for TASK, QASK, etc.
*
*     (Input)
*        t0:  time 0 in tiskz/taskz (count)
*
*     (Output)
*        written to common 
*       
*     (Return)
*        none:
*       
*     (Creation Date and Author)
*       2008.07.25 ; First version by Y.Takeuchi
*
************************************************************************
      implicit none

#include "skparm.h"
#include "skhead.h"
#include "sktq.h"
#include "skbadc.h"

      integer it0, od_start, od_end, i
      integer nhitaztmp, icaba, j
      INTEGER ISKBIT
      COMMON /SKOPT/ISKBIT

c*** For OD (moved from tqskaz.F)
c      od_start = -576 + it0 - it0sk
c      od_end   = 1920 + it0 - it0sk

c*** For OD (copied from ID)
      if (nrunsk.lt.50791) then
         od_start =  1424 + it0 - it0sk
         od_end   =  3920 + it0 - it0sk
      else if (nrunsk.eq.52026) then
         od_start = -1960 + it0 - it0sk
         od_end   =   536 + it0 - it0sk
      else
         od_start = SKGATE_START_COUNT + it0 - it0sk
         od_end   = SKGATE_END_COUNT   + it0 - it0sk
      endif
      ! convert to nsec since we don't have itiskz
      od_start = od_start/COUNT_PER_NSEC
      od_end = od_end/COUNT_PER_NSEC

c*** bad channel mask will be removed         
      NHITAZ = nhitaz_raw
      do i = 1, nhitaz_raw
         itaskz(i) = itaskz_raw(i)
         iqaskz(i) = iqaskz_raw(i)
         taskz(i)  = taskz_raw(i) -(it0-it0sk)/count_per_nsec
         qaskz(i)  = qaskz_raw(i)
         icabaz(i) = iand(icabaz_raw(i),2**16-1)
         IHTFLZ(I) = ishft(icabaz_raw(i), -16)

         IHTFLZ(i) = iand(IHTFLZ(i),2**16-2) ! make 0bit = 0
         if (iand(IHTFLZ(i),2).eq.2 .and.
     $        TASKZ(i).ge.od_start .and. 
     $        TASKZ(i).le.od_end) then
c***  In 1.3usec & in-gate
            IHTFLZ(i) = IHTFLZ(i) + 1 ! 0bit
c            nqask = nqask + 1   ! in gate & 1.3usec
c               print*,'in gate ',nqask,icabaz(i),itaskz(i),iqaskz(i)
         endif
         icabaz_raw(i) = ICABAZ(i) + (IHTFLZ(i) * 2**16)
      enddo

cc      print *, 'set_timing_gate_od: ',it0sk,it0,od_start,od_end

      return
      end

************************************************************************
*     -------------------------------------
      SUBROUTINE SET_TIMING_GATE_VETO_M(it0)
*     -------------------------------------
*
*     (Purpose)
*        set timing gate for TVSK, QVSK, etc.
*
*     (Input)
*        t0:  time 0 in tiskz/taskz (count)
*
*     (Output)
*        written to common 
*       
*     (Return)
*        none:
*       
*     (Creation Date and Author)
*       2008.12.09 ; First version by Y.Koshio
*
************************************************************************
      implicit none

#include "skparm.h"
#include "skhead.h"
#include "sktq.h"
#include "skveto.h"

      integer it0, id_start, id_end, i

c*** For ID (moved from tqskz.F)
      if (nrunsk.lt.50791) then
         id_start =  1424 + it0 - it0sk
         id_end   =  3920 + it0 - it0sk
      else if (nrunsk.eq.52026) then
         id_start = -1960 + it0 - it0sk
         id_end   =   536 + it0 - it0sk
      else
         id_start = SKGATE_START_COUNT + it0 - it0sk
         id_end   = SKGATE_END_COUNT   + it0 - it0sk
      endif
      ! convert to nsec since we don't have itiskz
      id_start = id_start/COUNT_PER_NSEC
      id_end = id_end/COUNT_PER_NSEC

      if (nqiskz.gt.0) then
c*** clear iabv(MAXPM+1,MAXPMV) to choose only the earliest hits
         do i = maxpm+1,maxpmv
            iabv(i) = 0
         enddo
         
         do i = 1, nqiskz
            if (icabiz(i).ge.MAXPM+1 .and. icabiz(i).le.MAXPMV) then
               IHTIFLZ(i) = iand(IHTIFLZ(i),2**16-2) ! make 0bit = 0
               if(iand(IHTIFLZ(i),2).eq.2 .and.
     $              iabv(icabiz(i)).eq.0.and.
     $              TISKZ(i).ge.id_start .and. 
     $              TISKZ(i).le.id_end) then
c***  In 1.3usec & in-gate
                  IHTIFLZ(i) = IHTIFLZ(i) + 1 ! 0bit
                  iabv(icabiz(i)) = 1
               endif
               icabbf_raw(i) = ICABIZ(i) + (IHTIFLZ(i) * 2**16)
            endif
         enddo
      endif

      return
      end
