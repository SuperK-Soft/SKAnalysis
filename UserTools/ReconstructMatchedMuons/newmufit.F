	subroutine newmufit(in,out,maxgood)
	implicit none

#include "skparm.h"
#include "skhead.h"
#include "sktq.h"
#include "skpdst.h"
#include "lfextinfo.h"
#include "geopmt.h"
#include "skwt.h"
#include "skday.h"
#include "skbadc.h"
#include "mufit.h"

	integer i,j,k,flag,total,spacing,l
	real in(3),out(3),dir(3),maxgood,ex,ey,ez
	real ox,oy,oz,cx,cy,cz,costheta,mod2
	real*8 gdn, explkymf_kirk
	real x(25000),y(25000),z(25000),mod1,qin,qall
	real tarray(12000),qarray(12000),vec(3),nin
	real rmod,qmax(25000)

	maxgood = 0
	flag = 0

c	open(unit=1,file='xyzpmt.txt')
c	open(unit=8,file='gridgood.txt')

	spacing = 120

        do i=1,nqisk
	   j=ihcab(i)
	   tarray(i)=tisk(j)
	   qarray(i)=qisk(j)
	end do

	call cylgrid(spacing,x,y,z,total)

	do i = 1,total
	   qmax(i) = 0
	   do j = 1,nqisk
	      k = ihcab(j)
	      rmod = sqrt((x(i)-xyzpm(1,k))**2+(y(i)-xyzpm(2,k))**2
     ,                    + (z(i)-xyzpm(3,k))**2)
	      if(rmod.le.200.and.qisk(k).gt.qmax(i)) qmax(i)=qisk(k)
	   end do
	end do

	print *, 'starting coarse fit:', total, qismsk

	do i = 1, total

c	   print *, i,x(i),y(i),z(i),qmax(i)

c	   if(mod(i,1000).eq.0) then
c	      print *, i, ' of ', total
c	   end if
	   if(qismsk.gt.100000) then
	      if(qmax(i).lt.20) GOTO 33
	   end if
	   if(qismsk.gt.250000) then
	      if(qmax(i).lt.26) GOTO 33
	   end if
	   if(qismsk.gt.400000) then
	      if(qmax(i).lt.33) GOTO 33
	   end if
	   if(qismsk.le.100000) then
	      if(qmax(i).lt.12) GOTO 33
	   end if
	   in(1) = x(i)
	   in(2) = y(i)
	   in(3) = z(i)
	   do j = 1,total
c	      print *, in,out,qmax(i),qmax(j)
	      if(qmax(j).lt.70) GOTO 88
c	      print *, in,out,qmax(i),qmax(j)
	      flag = 0
	      out(1) = x(j)
	      out(2) = y(j)
	      out(3) = z(j)
	      if(in(3).gt.1767.and.out(3).gt.1750) flag = 1
	      if(in(3).lt.-1767.and.out(3).lt.-1750) flag = 1
	      if(sqrt(in(1)**2+in(2)**2).eq.1690.and.in(1).eq.out(1)
     ,            .and.in(2).eq.out(2)) flag = 1
	      if(flag.eq.0) then
		 nin = 0
		 qin = 0
		 qall = 0
		 do k = 1,nqisk
		    costheta = 0
		    do l = 1,3
		       vec(l)=XYZPM(l,ihcab(k))-in(l)
		       dir(l) = out(l)-in(l)
		    end do
		    mod1 = sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
		    mod2 = sqrt(dir(1)**2+dir(2)**2+dir(3)**2)
		    do l = 1,3
		       vec(l) = vec(l)/mod1
		       dir(l) = dir(l)/mod2
		       costheta = costheta + vec(l)*dir(l)
		    end do
		    qall = qall + qarray(k)
		    if(costheta.gt.0.74) then
		       nin = nin + 1
		       qin = qin + qarray(k)
		    end if
		 end do	      

		 call MFFLTRES_KIRK(in,out)

c		 gdn = EXPLKYMF_KIRK(1)*exp((qin/qall)**2/(nin/nqisk)**2)
		 gdn = EXPLKYMF_KIRK(1)
c		 write(8,*) gdn, in, out
c		 print *,  gdn, in, out

		 if(gdn.gt.maxgood) then

c		    print *, 'new max good for coarse search;'
c		    print *, gdn, in, out

		    maxgood = gdn
		    
		    ex = in(1)
		    ey = in(2)
		    ez = in(3)
		    
		    ox = out(1)
		    oy = out(2)
		    oz = out(3)

		 end if
	      end if
 88	   end do
 33	end do

	in(1) = ex
	in(2) = ey
	in(3) = ez
	
	out(1) = ox
	out(2) = oy
	out(3) = oz

	print *, 'finished coarse fit, starting fine fit.', in, out

	call finesearch(in,out,maxgood)

	print *, 'finished fine search', in,out,maxgood

c	close(1)

	return
	END
	

#include "cylgrid.F"
#include "mffltres_kirk.F"
#include "explkymf_kirk.F"

	subroutine finesearch(in,out,good)

	real in(3),out(3),dir(3),tin(3,4),tout(3,4),good
        real dist,oldin(3),oldout(3), oangle,nangle, newangle
	real*8 gdn, explkymf_kirk, oldgdn, savegdn
	integer i,j,k

	do i = 1,3
	   oldin(i) = in(i)
	   oldout(i) = out(i)
	end do

	call MFFLTRES_KIRK(in,out)
	oldgdn = EXPLKYMF_KIRK(1)

	dist = 30
 10	if(dist.lt.2) GOTO 11

	print *, oldgdn, savegdn, dist, in, out

	if(oldgdn.gt.-100.and.oldgdn.lt.100) then
	   continue
	else
	   print *, 'WTF IS HAPPENING!?!?!?!'
	   print *, 'oldgdn,dist,in,out'
	   print *, oldgdn,dist,in,out
	   STOP
	end if

	if(abs(in(3)).lt.1810.) then
	   if(abs(sqrt(in(1)**2+in(2)**2)-1690).gt.0.01) then
	      print *, 'entrance not on detector!', in
	   end if
	end if

	if(abs(out(3)).lt.1810.) then
	   if(abs(sqrt(out(1)**2+out(2)**2)-1690).gt.0.01) then
	      print *, 'exit pos not on detector!', out
	   end if
	end if

	savegdn = oldgdn

	if(abs(in(3)).gt.1810.0) then
	   print *, 'entrance position out of bounds!'
	   STOP
	end if

	if(abs(out(3)).gt.1810.0) then
	   print *, 'exit position out of bounds!'
	   STOP
	end if

	do j = 1,3
	   do i = 1,4
	      tin(j,i) = in(j)
	      tout(j,i) = out(j)
	   end do
	end do

cccccc entrance position is on side ccccccccccc

	if(abs(in(3)).lt.1810.0) then

	      nangle = dist/1690. * 6.2832

	      oangle = atan2(in(2),in(1))
	      
	      if(abs(in(1)-1690.*cos(oangle)).gt.1) then
		 print *, 'whoa, angle isnt right!!!! 1'
		 print *,in(1),in(2),oangle,abs(in(1)-1690.*cos(oangle))
	      end if

	      tin(3,1) = tin(3,1) + dist
	      if(tin(3,1).gt.1810.) then
		 tin(3,1) = 1810.
		 if(tin(1,1).lt.0) then
		    tin(1,1) = tin(1,1)+1
		 else
		    tin(2,1) = tin(2,1)-1
		 end if
		 if(tin(2,1).lt.0) then
		    tin(2,1) = tin(2,1) +1
		 else
		    tin(2,1) = tin(2,1) -1
		 end if
	      end if

	      tin(3,2) = tin(3,2) - dist
	      if(tin(3,2).lt.-1810.) then
		 tin(3,2) = -1810.
		 if(tin(1,2).lt.0) then
		    tin(1,2) = tin(1,2)+1
		 else
		    tin(2,2) = tin(2,2)-1
		 end if
		 if(tin(2,2).lt.0) then
		    tin(2,2) = tin(2,2) +1
		 else
		    tin(2,2) = tin(2,2) -1
		 end if
	      end if
	      
	      tin(1,3) = 1690. * cos(oangle+newangle)
	      tin(2,3) = 1690. * sin(oangle+newangle)

	      tin(1,4) = 1690. * cos(oangle-newangle)
	      tin(2,4) = 1690. * sin(oangle-newangle)

cccc   exit position is on top or bottom 111111111111

	   if(abs(out(3)).eq.1810..and.sqrt(out(1)**2+out(2)**2)
     ,     .le.1690.) then

	      tout(1,1) = tout(1,1) + dist
	      if(sqrt(tout(1,1)**2+tout(2,1)**2).gt.1690) then
		 nangle = atan2(tout(2,1),tout(1,1))
		 tout(1,1) = 1690*cos(nangle)
		 tout(2,1) = 1690*sin(nangle)
		 if(out(3).lt.0) then
		    tout(3,1) = -1800.
		 else
		    tout(3,1) = 1800.
		 end if
	      end if

	      tout(1,2) = tout(1,2) - dist
	      if(sqrt(tout(1,2)**2+tout(2,2)**2).gt.1690) then
		 nangle = atan2(tout(2,2),tout(1,2))
		 tout(1,2) = 1690*cos(nangle)
		 tout(2,2) = 1690*sin(nangle)
		 if(out(3).lt.0) then
		    tout(3,2) = -1800.
		 else
		    tout(3,2) = 1800.
		 end if
	      end if

	      tout(2,3) = tout(2,3) + dist
	      if(sqrt(tout(1,3)**2+tout(2,3)**2).gt.1690) then
		 nangle = atan2(tout(2,3),tout(1,3))
		 tout(1,3) = 1690*cos(nangle)
		 tout(2,3) = 1690*sin(nangle)
		 if(out(3).lt.0) then
		    tout(3,3) = -1800.
		 else
		    tout(3,3) = 1800.
		 end if
	      end if

	      tout(2,4) = tout(2,4) - dist
	      if(sqrt(tout(1,4)**2+tout(2,4)**2).gt.1690) then
		 nangle = atan2(tout(2,4),tout(1,4))
		 tout(1,4) = 1690*cos(nangle)
		 tout(2,4) = 1690*sin(nangle)
		 if(out(3).lt.0) then
		    tout(3,4) = -1800.
		 else
		    tout(3,4) = 1800.
		 end if
	      end if

	      do i = 1,4
		 do j = 1,4		    
		    call MFFLTRES_KIRK(tin(1,i),tout(1,j))
		    gdn = EXPLKYMF_KIRK(1)
		    if(gdn.gt.oldgdn) then
		       oldgdn = gdn
		       do k = 1,3
			  in(k) = tin(k,i)
			  out(k) = tout(k,j)
		       end do
		    end if
		 end do
	      end do

	      print *, 'did 1, new in,out is:'
	      print *, in, out
	   
	   end if

cccc   exit position is on side 11111111111

	   if(abs(out(3)).lt.1810) then

	      nangle = dist/1690. * 6.2832

	      oangle = atan2(out(2),out(1))
	      
	      if(abs(out(1)-1690.*cos(oangle)).gt.1) then
		 print *, 'whoa, angle isnt right!!!! 2'
		 print *,out(1),out(2),oangle,abs(out(1)-1690.*cos(oangle))
	      end if

	      tout(3,1) = tout(3,1) + dist
	      if(tout(3,1).gt.1810.) then
		 tout(3,1) = 1810.
		 if(tout(1,1).lt.0) then
		    tout(1,1) = tout(1,1)+1
		 else
		    tout(2,1) = tout(2,1)-1
		 end if
		 if(tout(2,1).lt.0) then
		    tout(2,1) = tout(2,1) +1
		 else
		    tout(2,1) = tout(2,1) -1
		 end if
	      end if

	      tout(3,2) = tout(3,2) - dist
	      if(tout(3,2).lt.-1810.) then
		 tout(3,2) = -1810.
		 if(tout(1,2).lt.0) then
		    tout(1,2) = tout(1,2)+1
		 else
		    tout(2,2) = tout(2,2)-1
		 end if
		 if(tout(2,2).lt.0) then
		    tout(2,2) = tout(2,2) +1
		 else
		    tout(2,2) = tout(2,2) -1
		 end if
	      end if

	      tout(1,3) = 1690. * cos(oangle+newangle)
	      tout(2,3) = 1690. * sin(oangle+newangle)

	      tout(1,4) = 1690. * cos(oangle-newangle)
	      tout(2,4) = 1690. * sin(oangle-newangle)

	      do i = 1,4
		 do j = 1,4		    
		    call MFFLTRES_KIRK(tin(1,i),tout(1,j))
		    gdn = EXPLKYMF_KIRK(1)
		    if(gdn.gt.oldgdn) then
		       oldgdn = gdn
		       do k = 1,3
			  in(k) = tin(k,i)
			  out(k) = tout(k,j)
		       end do
		    end if
		 end do
	      end do

	      print *, 'did 2, new in,out is:'
	      print *, in,out

	   end if
	   
	   if(oldgdn.eq.savegdn) then
	      dist = dist/2.
	   end if

	   GOTO 10

	end if



cccccc entrance position is on top or bottom cccccccc

	if(abs(in(3)).eq.1810..and.sqrt(in(1)**2+in(2)**2)
     ,     .le.1690.01) then

	      tin(1,1) = tin(1,1) + dist
	      if(sqrt(tin(1,1)**2+tin(2,1)**2).gt.1690) then
		 nangle = atan2(tin(2,1),tin(1,1))
		 tin(1,1) = 1690*cos(nangle)
		 tin(2,1) = 1690*sin(nangle)
		 if(in(3).lt.0) then
		    tin(3,1) = -1800.
		 else
		    tin(3,1) = 1800.
		 end if
	      end if

	      tin(1,2) = tin(1,2) - dist
	      if(sqrt(tin(1,2)**2+tin(2,2)**2).gt.1690) then
		 nangle = atan2(tin(2,2),tin(1,2))
		 tin(1,2) = 1690*cos(nangle)
		 tin(2,2) = 1690*sin(nangle)
		 if(in(3).lt.0) then
		    tin(3,2) = -1800.
		 else
		    tin(3,2) = 1800.
		 end if
	      end if

	      tin(2,3) = tin(2,3) + dist
	      if(sqrt(tin(1,3)**2+tin(2,3)**2).gt.1690) then
		 nangle = atan2(tin(2,3),tin(1,3))
		 tin(1,3) = 1690*cos(nangle)
		 tin(2,3) = 1690*sin(nangle)
		 if(in(3).lt.0) then
		    tin(3,3) = -1800.
		 else
		    tin(3,3) = 1800.
		 end if
	      end if

	      tin(2,4) = tin(2,4) - dist
	      if(sqrt(tin(1,4)**2+tin(2,4)**2).gt.1690) then
		 nangle = atan2(tin(2,4),tin(1,4))
		 tin(1,4) = 1690*cos(nangle)
		 tin(2,4) = 1690*sin(nangle)
		 if(in(3).lt.0) then
		    tin(3,4) = -1800.
		 else
		    tin(3,4) = 1800.
		 end if
	      end if

cccc exit position is on top or bottom 222222222

	   if(abs(out(3)).eq.1810..and.sqrt(out(1)**2+out(2)**2)
     ,     .lt.1690.) then

	      tout(1,1) = tout(1,1) + dist
	      if(sqrt(tout(1,1)**2+tout(2,1)**2).gt.1690) then
		 nangle = atan2(tout(2,1),tout(1,1))
		 tout(1,1) = 1690*cos(nangle)
		 tout(2,1) = 1690*sin(nangle)
		 if(out(3).lt.0) then
		    tout(3,1) = -1800.
		 else
		    tout(3,1) = 1800.
		 end if
	      end if

	      tout(1,2) = tout(1,2) - dist
	      if(sqrt(tout(1,2)**2+tout(2,2)**2).gt.1690) then
		 nangle = atan2(tout(2,2),tout(1,2))
		 tout(1,2) = 1690*cos(nangle)
		 tout(2,2) = 1690*sin(nangle)
		 if(out(3).lt.0) then
		    tout(3,2) = -1800.
		 else
		    tout(3,2) = 1800.
		 end if
	      end if

	      tout(2,3) = tout(2,3) + dist
	      if(sqrt(tout(1,3)**2+tout(2,3)**2).gt.1690) then
		 nangle = atan2(tout(2,3),tout(1,3))
		 tout(1,3) = 1690*cos(nangle)
		 tout(2,3) = 1690*sin(nangle)
		 if(out(3).lt.0) then
		    tout(3,3) = -1800.
		 else
		    tout(3,3) = 1800.
		 end if
	      end if

	      tout(2,4) = tout(2,4) - dist
	      if(sqrt(tout(1,4)**2+tout(2,4)**2).gt.1690) then
		 nangle = atan2(tout(2,4),tout(1,4))
		 tout(1,4) = 1690*cos(nangle)
		 tout(2,4) = 1690*sin(nangle)
		 if(out(3).lt.0) then
		    tout(3,4) = -1800.
		 else
		    tout(3,4) = 1800.
		 end if
	      end if

	      do i = 1,4
		 do j = 1,4		    
		    call MFFLTRES_KIRK(tin(1,i),tout(1,j))
		    gdn = EXPLKYMF_KIRK(1)
c	            print *, 'trying:',tin(1,i),tin(2,i),tin(3,i),
c     ,                        tout(1,j),tout(2,j),tout(3,j),gdn

		    if(gdn.gt.oldgdn) then
		       oldgdn = gdn
		       do k = 1,3
			  in(k) = tin(k,i)
			  out(k) = tout(k,j)
		       end do
		    end if
		 end do
	      end do

	      print *, 'did 3.  New in,out is:'
	      print *, in,out

	   end if

ccccc	exit position is on side 222222222222   

	   if(abs(out(3)).lt.1810) then

	      nangle = dist/1690. * 6.2832

	      oangle = atan2(out(2),out(1))
	      
	      if(abs(out(1)-1690.*cos(oangle)).gt.1) then
		 print *, 'whoa, angle isnt right!!!! 3'
		 print *,out(1),out(2),oangle,abs(out(1)-1690.*cos(oangle))
	      end if

	      tout(3,1) = tout(3,1) + dist
	      if(tout(3,1).gt.1810.) then
		 tout(3,1) = 1810.
		 if(tout(1,1).lt.0) then
		    tout(1,1) = tout(1,1)+1
		 else
		    tout(2,1) = tout(2,1)-1
		 end if
		 if(tout(2,1).lt.0) then
		    tout(2,1) = tout(2,1) +1
		 else
		    tout(2,1) = tout(2,1) -1
		 end if
	      end if

	      tout(3,2) = tout(3,2) - dist
	      if(tout(3,2).lt.-1810.) then
		 tout(3,2) = -1810.
		 if(tout(1,2).lt.0) then
		    tout(1,2) = tout(1,2)+1
		 else
		    tout(2,2) = tout(2,2)-1
		 end if
		 if(tout(2,2).lt.0) then
		    tout(2,2) = tout(2,2) +1
		 else
		    tout(2,2) = tout(2,2) -1
		 end if
	      end if

	      tout(1,3) = 1690. * cos(oangle+newangle)
	      tout(2,3) = 1690. * sin(oangle+newangle)

	      tout(1,4) = 1690. * cos(oangle-newangle)
	      tout(2,4) = 1690. * sin(oangle-newangle)

	      do i = 1,4
		 do j = 1,4		    
		    call MFFLTRES_KIRK(tin(1,i),tout(1,j))
		    gdn = EXPLKYMF_KIRK(1)
c		    print *, '4, trying in out pos, good:'
c		    print *, tin(1,i),tout(1,j),gdn
		    if(gdn.gt.oldgdn) then
		       oldgdn = gdn
		       do k = 1,3
			  in(k) = tin(k,i)
			  out(k) = tout(k,j)
		       end do
		    end if
		 end do
	      end do

	      print *, 'did 4.  new in,out is:'
	      print *, in,out
	      
	   end if

	   if(oldgdn.eq.savegdn) then
	      dist = dist/2.
	   end if

	   GOTO 10

	end if

	print *, 'neither branch hit!!!!', in, out

 11	continue

	print *, dist, gdn, in, out

 	good = oldgdn

	return
	end

