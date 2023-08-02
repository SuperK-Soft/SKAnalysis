*******************************************************************
      FUNCTION EXPLKYMF_KIRK(QFLAG)
*
*     calcurate likelyfood function from
*     residual time of each PMTs
*
*     FUN = Σ1/σ*exp(-(t)**2/σ**2)
*           -------------------------
*                  Σ1/σ
*
*******************************************************************
      IMPLICIT NONE
      REAL*8 EXPLKYMF_KIRK
      INTEGER QFLAG

#include "skparm.h"
#include "skhead.h"
#include "sktq.h"
#include "skpdst.h"
#include "geotnk.h"
#include "geopmt.h"
#include "skbadc.h"
#include "mufit.h"

      INTEGER IPMT
      REAL SIGMA, WEIGHT,lq,lamdaq
      REAL*8 LIKE, NORM, TMP
      REAL*8 TDCAVE, TSUM, TMPAVE

      LIKE = 0.
      NORM = 0.

      TDCAVE = 0.
      TSUM = 0.
      DO 98 IPMT=1,MAXPM
         IF (MFTRES(IPMT).gt.5000.) GOTO 98
         Call SKTSIG(QISK(IPMT),SIGMA)
CC         SIGMA = SQRT(SIGMA)
         if(sigma.gt.10) then
CC            print *, 'sigma: ', sigma
         end if
         TDCAVE = TDCAVE + MFTRES(IPMT)/SIGMA
         TSUM = TSUM + 1.0/SIGMA
 98   ENDDO
      TDCAVE = TDCAVE/TSUM
      TMPAVE = TDCAVE
      TDCAVE = 0.
      TSUM = 0.
      DO 99 IPMT=1,MAXPM
CC         print *, 'IAB, QISK, MFTRES, TMPAVE='
CC     &        , IAB(IPMT), QISK(IPMT), MFTRES(IPMT), TMPAVE
         IF (MFTRES(IPMT).gt.(TMPAVE+100.)
     &        .or.MFTRES(IPMT).lt.(TMPAVE-100.)) GOTO 99
         Call SKTSIG(QISK(IPMT),SIGMA)
CC         SIGMA = SQRT(SIGMA)
         TDCAVE = TDCAVE + MFTRES(IPMT)/SIGMA
         TSUM = TSUM + 1.0/SIGMA
 99   ENDDO
      TDCAVE = TDCAVE/TSUM

      DO 100 IPMT=1,MAXPM
CC         print *, 'IAB, QISK, MFTRES, TDCAVE='
CC     &        , IAB(IPMT), QISK(IPMT), MFTRES(IPMT), TDCAVE
         IF (MFTRES(IPMT).gt.(TDCAVE+200.)
     &        .or.MFTRES(IPMT).lt.(TDCAVE-200.)) GOTO 100
         Call SKTSIG(QISK(IPMT),SIGMA)
CC         SIGMA = SQRT(SIGMA)
         if(sigma.gt.10) then
CC            print *, 'sigma: ', sigma
         end if


c         if((MFTRES(IPMT)-tdcave.lt.5).or.qisk(ipmt).gt.30) then
         TMP = -0.5D0*DBLE(MFTRES(IPMT)-TDCAVE)**2D0
     &        /(DBLE(SIGMA)**2D0)
c         else
c            lq = lamdaq(qisk(ipmt))/29.94
c            tmp = -(MFTRES(IPMT)-TDCAVE)/lq
c         end if


         WEIGHT = 1.0
CC         WEIGHT = 100.0/(EMDIS(IPMT)+100.)
         LIKE = LIKE + 1.0D0/DBLE(SIGMA)*DBLE(WEIGHT)*DEXP(TMP)
         NORM = NORM + 1.0D0/DBLE(SIGMA)*DBLE(WEIGHT)
 100  ENDDO

CC      print *, 'like,norm,nhitmf42,nhitmf:'
CC      print *, like, norm, nhitmf42,nhitmf
      IF (QFLAG.eq.0) THEN
         EXPLKYMF_KIRK = LIKE / NORM
      ELSE
         EXPLKYMF_KIRK = LIKE / NORM * 
     &                     (DBLE(NHITMF42)/DBLE(NHITMF))
CC         print *, 'LIKELYFOOD changed', LIKE/NORM, '*'
CC     &        , (FLOAT(NHITMF42)/FLOAT(NHITMF))**4, '=', EXPLKYMF
      ENDIF
CC      EXPLKYMF = LIKE / NORM * DBLE(QTOT42/QTOT)**5D0
CC      print *, 'LIKELYFOOD changed', LIKE/NORM, '*'
CC     &     , (QTOT42/QTOT)**5, '=', EXPLKYMF
CC      EXPLKYMF = LIKE / NORM

      RETURN
      END

