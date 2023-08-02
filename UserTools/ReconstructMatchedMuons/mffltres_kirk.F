************************************************************************
*     -------------------------------------
      subroutine MFFLTRES_KIRK(POSIN,POSOUT)
*     -------------------------------------
*
*     (Purpose)
*       deduce residual time for each PMTs and fill arrays.
*       dedicated for cosmic ray muon fitting.
*
************************************************************************
      IMPLICIT NONE
      REAL POSIN(3), POSOUT(3)

#include "skparm.h"
#include "skhead.h"
#include "sktq.h"
#include "sktabl.h"
#include "geotnk.h"
#include "geopmt.h"
#include "skbadc.h"
#include "mufit.h"

      REAL TRACK
      REAL WATERN, WATERNC, COS42
C      PARAMETER (WATERN=1.33)
      PARAMETER (WATERN=1.344)
      PARAMETER (WATERNC=1.4)
      REAL A, B, C, TMPDIR(3), D, X, MUDIR(3), DIST
      REAL R,tempr
      INTEGER IPMT,bin,i, threshold
      REAL QOUT,tv(3),phi,alpha,theta1,theta2

      COS42 = 1./WATERN

c      open(unit=1,file='fitres.txt')

c      if(qismsk.gt.100000) threshold = 10.0
c      if(qismsk.le.100000) threshold = 6.0

      if(qismsk.le.50000) then
         threshold = 1.6*sqrt(qismsk/1000.) + 1
      end if
      
      if(qismsk.gt.50000) then
         threshold = 1.8*sqrt(qismsk/1000.) + 1
      end if


      TRACK = SQRT((POSIN(1)-POSOUT(1))**2
     &     +(POSIN(2)-POSOUT(2))**2+(POSIN(3)-POSOUT(3))**2)
      MUDIR(1) = POSOUT(1) - POSIN(1)
      MUDIR(2) = POSOUT(2) - POSIN(2)
      MUDIR(3) = POSOUT(3) - POSIN(3)
      DIST = SQRT(MUDIR(1)**2+MUDIR(2)**2+MUDIR(3)**2)
      MUDIR(1) = MUDIR(1)/DIST
      MUDIR(2) = MUDIR(2)/DIST
      MUDIR(3) = MUDIR(3)/DIST
c      print *, mudir,dist,posin,posout
      QTOTMF42 = 0.
      QTOTMF   = 0.
      QOUT   = 0.
      NHITMF42 = 0
      NHITMF = 0

         theta2 = atan(mudir(1)/mudir(3))
         alpha = sin(theta2)*mudir(1)+cos(theta2)*mudir(3)
         theta1 = asin(mudir(2)/(mudir(2)**2+alpha**2))

         tv(1) = cos(theta2)*mudir(1)-sin(theta2)*mudir(3)
         tv(2) = cos(theta1)*mudir(2)+sin(theta1)*alpha
         tv(3) = sin(theta1)*mudir(2)-cos(theta1)*alpha

c         print *, '!!!!!!!!tv!!!!!!!!!'
c         print *, tv

      DO 100 IPMT = 1,MAXPM
         IF (IAB(IPMT).eq.0
     &        .or.QISK(IPMT).lt.threshold.or.QISK(IPMT).gt.400.0
     &        .or.TISK(IPMT).eq.0) then
c     &        .or.NATMCH(IPMT).gt.3) THEN
            EMDIRMF(1,IPMT) = 0.
            EMDIRMF(2,IPMT) = 0.
            EMDIRMF(3,IPMT) = 0.
            EMPOSMF(1,IPMT) = 0.
            EMPOSMF(2,IPMT) = 0.
            EMPOSMF(3,IPMT) = 0.
            MFTRES(IPMT) = 10000.
c            print *, ipmt,qisk(ipmt),tisk(ipmt),natmch(ipmt),iab(ipmt)
            GOTO 100
         ENDIF
c         if(qisk(ipmt).gt.20) print *, ipmt,qisk(ipmt),iab(ipmt)
         QTOTMF = QTOTMF + QISK(IPMT)
         NHITMF = NHITMF + 1
C find emitting point of cherenkov light
C
C                x entrance(COMPOS(XYZ))
C               /    | |
C              /     | |
C             x E( COMPOS + MUDIR x X )
C            /-\     | |
C           /42 \    | |
C          /     \   | |
C      mu /       \  | |
C        /         \ | |
C       /           (D-| 
C      /          PMT(XYZPM)   
C     /              | |
C    x exit (POSOUT(XYZ))
C=======================

         TMPDIR(1) = XYZPM(1,IPMT)-POSIN(1)
         TMPDIR(2) = XYZPM(2,IPMT)-POSIN(2)
         TMPDIR(3) = XYZPM(3,IPMT)-POSIN(3)
 
         A = 1.0-COS42**2
         B = -2.0*(1.0-COS42**2)
     &        *(TMPDIR(1)*MUDIR(1)+TMPDIR(2)*MUDIR(2)
     &        +TMPDIR(3)*MUDIR(3))
         C = (TMPDIR(1)*MUDIR(1)+TMPDIR(2)*MUDIR(2)
     &        +TMPDIR(3)*MUDIR(3))**2
     &        -COS42**2*(TMPDIR(1)**2+TMPDIR(2)**2+TMPDIR(3)**2)
         D = B**2 - 4.*A*C
         IF (D.lt.0.) THEN
c            print *, ' %trans-hanbetsu < 0'
c            print *, 'A,B,C=',A,B,C
c            print *, 'XYZPM(XYZ)='
c     &           ,XYZPM(1,IPMT),XYZPM(2,IPMT),XYZPM(3,IPMT)
c            print *, 'POSIN=',POSIN(1),POSIN(2),POSIN(3)
c            print *, 'MUDIR=',MUDIR(1),MUDIR(2),MUDIR(3)
            EMDIRMF(1,IPMT) = 0.
            EMDIRMF(2,IPMT) = 0.
            EMDIRMF(3,IPMT) = 0.
            EMPOSMF(1,IPMT) = 0.
            EMPOSMF(2,IPMT) = 0.
            EMPOSMF(3,IPMT) = 0.
            MFTRES(IPMT) = 10000.
            GOTO 100
         ENDIF
         X = (-1.0*B-SQRT(B**2-4.0*A*C))/(2.0*A)
         IF (X.lt.0.or.X.gt.TRACK) THEN
            EMDIRMF(1,IPMT) = 0.
            EMDIRMF(2,IPMT) = 0.
            EMDIRMF(3,IPMT) = 0.
            EMPOSMF(1,IPMT) = 0.
            EMPOSMF(2,IPMT) = 0.
            EMPOSMF(3,IPMT) = 0.
            MFTRES(IPMT) = 10000.
            QOUT = QOUT + QISK(IPMT)
            GOTO 100
         ENDIF

C emitting point of cherenkov light
         EMPOSMF(1,IPMT) = POSIN(1) + MUDIR(1)*X
         EMPOSMF(2,IPMT) = POSIN(2) + MUDIR(2)*X
         EMPOSMF(3,IPMT) = POSIN(3) + MUDIR(3)*X
C directon of emitted light
         EMDIRMF(1,IPMT) = XYZPM(1,IPMT)-EMPOSMF(1,IPMT)
         EMDIRMF(2,IPMT) = XYZPM(2,IPMT)-EMPOSMF(2,IPMT)
         EMDIRMF(3,IPMT) = XYZPM(3,IPMT)-EMPOSMF(3,IPMT)
C traveling length of cherenkov light 
         R = SQRT(EMDIRMF(1,IPMT)**2+EMDIRMF(2,IPMT)**2
     &        +EMDIRMF(3,IPMT)**2)
         EMDIRMF(1,IPMT) = EMDIRMF(1,IPMT)/R
         EMDIRMF(2,IPMT) = EMDIRMF(2,IPMT)/R
         EMDIRMF(3,IPMT) = EMDIRMF(3,IPMT)/R
C angle between PMT and light
C         COSTH(IPMT) = -EMDIRMF(1,IPMT)*DXYZPM(1,IPMT)
C     &        -EMDIRMF(2,IPMT)*DXYZPM(2,IPMT)
C     &        -EMDIRMF(3,IPMT)*DXYZPM(3,IPMT)

         EMDISMF(IPMT) = X
C         MFTRES(IPMT) = TISK(IPMT)-X/30.-R/30.*WATERN
         MFTRES(IPMT) = TISK(IPMT)-X/30.-R/30.*WATERNC
         QTOTMF42 = QTOTMF42+QISK(IPMT)
         NHITMF42 = NHITMF42+1

         if(sqrt(emposmf(1,ipmt)**2+emposmf(2,ipmt)**2).lt.1490.
     ,      .and.abs(emposmf(3,ipmt)).lt.1610.) then

         tv(1) = emdirmf(1,ipmt)*cos(theta2)-sin(theta2)*emdirmf(3,ipmt)
         tv(2) = cos(theta1)*emdirmf(2,ipmt)+sin(theta1)*
     ,            (sin(theta2)*emdirmf(1,ipmt)+cos(theta2)*emdirmf(3,ipmt))
         tv(3) = sin(theta1)*emdirmf(2,ipmt)-cos(theta1)*
     ,            (sin(theta2)*emdirmf(1,ipmt)+cos(theta2)*emdirmf(3,ipmt))

         tempr = mudir(1)*emdirmf(1,ipmt)+mudir(2)*emdirmf(2,ipmt)+
     ,           mudir(3)*emdirmf(3,ipmt)

         phi = atan2(tv(2),tv(1))

         phi = phi * 180./3.14159 + 180

c         write(1,*) phi, tempr, sqrt(tv(1)**2+tv(2)**2+tv(3)**2), 
c     ,            tv(1),tv(2),tv(3),qisk(ipmt)

         bin = phi/10 + 1

         end if

 100  ENDDO

      close(1)

      RETURN
      END


