
      subroutine fastatom

      include "common.f"

      dimension tmin(3), tmax(3)
      logical   skip

      if (ipka.eq.0) return

      if (numstep.eq.1) then
        open(unit=12,file='fastatom.txt',status='unknown')
      end if

      distmin=1001.0
      atomkemax=0.0
      num1=0
      num10=0

      do ind=1,3
        tmin(ind)=0.0
        tmax(ind)=0.0
      end do

      do i=ifix+1,natom
        skip=.false.
        do ind=1,3
          if (x(i,ind).lt.0.0) then
             tmin(ind)=min(tmin(ind),x(i,ind))
             tmax(ind)=max(tmax(ind),x(i,ind))
            skip=.true.
          end if
        end do
        if (xepka)  skip=.false.
        if (c60pka) skip=.false.

        if (.not.skip) then
          atomke=0.0
          do ind=1,3
            atomke=atomke + 0.5*fracmass(i)*vx(i,ind)**2
          end do

          if (atomke.gt.10.0) then
            do ind=1,3
              dxx= x(i,ind)-tmin(ind)
              dxx=dxx - box(ind)*dnint(dxx/box(ind))
              distmin=min(distmin,abs(dxx))

              dxx= x(i,ind)-tmax(ind)
              dxx=dxx - box(ind)*dnint(dxx/box(ind))
              distmin=min(distmin,abs(dxx))
            end do
            
            if ((itherm.ne.0).and.(distmin.lt.5.0)) then
               write(6,*) 'i,atomke,distmin=',i,atomke,distmin
               stop 'Fast atom too close to boundary'
            end if
          end if
  
          if (atomke.gt.10.0) num10=num10+1
          if (atomke.gt.1.0) num1=num1+1
          atomkemax=max(atomke,atomkemax)
        end if
      end do

      if (mod(numstep,nprint).ne.0) return

      tsim=timetau*(t-t0)
      if (distmin.gt.1000.0) distmin=0.0

      write(12,100) numstep,tsim,atomkemax,num10,num1,distmin
      flush(12)
 100  format(I9,1X,F10.5,1X,F10.1,2I6,1X,F7.1)

      end

