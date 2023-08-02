           subroutine makededx(e,mdir,ihcab,qisk,tisk,xyzpm,nqisk,dedx)
           implicit none
   
           real length,e(4),mdir(3)
           external mujecttime, mujectangle
           integer mujecttime, mujectangle,mujresult,nqisk
           real kpmt(4),kentry(4),kdist(2),kdir(3),cosang(2)
           integer ti,stop_bin,index
           real dedx(200),sum,save,temp,percent,zerobins
           integer lbins,i,j,ihcab(11146)
           real qisk(11146),tisk(11146)
           real*4 xyzpm(3,11146)
   
           do i = 1,200
              dedx(i) = 0
           end do
             kentry(1) = e(1)
             kentry(2) = e(2)
             kentry(3) = e(3)
             kentry(4) = e(4)/29.94
   
             kdir(1) = mdir(1)
             kdir(2) = mdir(2)
             kdir(3) = mdir(3)
             
             do i=1,nqisk
               ti= ihcab(i)
   
               if(qisk(ti).ge.10) then
   
                  kdist(1) = 0
                  kdist(2) = 0
                  
                  kpmt(1) = xyzpm(1,ti)
                  kpmt(2) = xyzpm(2,ti)
                  kpmt(3) = xyzpm(3,ti)
                  kpmt(4) = tisk(ti)
   
             mujresult = mujecttime(kentry,kdir,kpmt,kdist,cosang)
             do j = 1,mujresult
             print *, "KDIST IS"
             print *, kdist(j)
             stop_bin=(kdist(j)/50.)+1.
             if((stop_bin.ge.1).and.(stop_bin.le.200)) then
               print *, "Setting dedx"
               dedx(stop_bin)=dedx(stop_bin)+qisk(ti)
             end if
             end do
             end if
             end do
   
             save = 0
             do i =1,191
                sum = 0
                do j = 1,10
                   sum = sum + dedx(i+j-1)
                end do
                if(sum.gt.save) then
                   index = i
                   save = sum
                end if
             end do
            lbins = (length /50.)+1
            j = lbins
            zerobins = 0
            do i = 1,j
               if(dedx(i).eq.0) zerobins=zerobins+1
            end do
   
            percent = zerobins/lbins
            temp = 50. *(index + 4)
            
            
            return
            end
   