      subroutine cylgrid(spacing,x,y,z,total)

      integer spacing, total,i,j,t1,t2,t3,t4
      real x(20000),y(20000),z(20000),r(3),angle

      total = 0

      r(3) = 1810

      do i = 1,20000
         x(i) = 0
         y(i) = 0
         z(i) = 0
      end do

      t1 = 1690/spacing
      t2 = 2*t1 + 1

      do i = 1, t2
      r(1) = 0 - t1*spacing + (i-1)*spacing
         do j = 1,t2
            r(2) = 0 - t1*spacing + (j-1)*spacing
            if(sqrt(r(1)**2+r(2)**2).le.1690) then
               total = total + 1
               x(total) = r(1)
               y(total) = r(2)
               z(total) = 1810.
            end if
         end do
      end do

      do i = 1, t2
      r(1) = 0 - t1*spacing + (i-1)*spacing
         do j = 1,t2
            r(2) = 0 - t1*spacing + (j-1)*spacing
            if(sqrt(r(1)**2+r(2)**2).le.1690) then
               total = total + 1
               x(total) = r(1)
               y(total) = r(2)
               z(total) = -1810.
            end if
         end do
      end do

c      print *, 'done with 2', total

      t1 = 1810/spacing
      t2 = t1*2+1

      t3 = (6.2832*1690)/spacing
      angle = 6.2832/t3

c      print *, 't3, angle = ', t3,angle

      do i = 1,t2
         r(3) = 0 - t1*spacing + (i-1)*spacing         
         do j = 1,t3
            r(1) = 1690*cos(j*angle)
            r(2) = 1690*sin(j*angle)
c            print *, j, j*angle, r(1), r(2)
            total = total + 1
            x(total) = r(1)
            y(total) = r(2)
            z(total) = r(3)
         end do
      end do

      return
      end
