
      subroutine varystep

      include "common.f"

      save      forcemax

c +--------------------------------------------------------+
c | Only use variable timestep for cascades or if flag set |
c +--------------------------------------------------------+
      if (.not.variablestep) return

c +------------+
c | Initialize |
c +------------+
      if (numstep.eq.0) then
        forcemax=0.0
        open(unit=10,file='maxforce.txt',status='unknown')
      end if
   
c +-----------------------------------+
c | Compute maximum force on any atom |
c +-----------------------------------+
      forcemaxnow=0.0
      do i=ifix+1,natom
        forcenorm=0.0
        do ind=1,3
          forcenorm=forcenorm + fx(i,ind)**2
        end do
        forcenorm=sqrt(forcenorm)
        forcemaxnow=max(forcemaxnow,forcenorm)
      end do

      forcemax=max(forcemax,forcemaxnow)
      !write(10,*) numstep, h*forcemax

c +-----------------------------------------------------------+
c | Immediately halve h if instantaneous forcemax is too high |
c +-----------------------------------------------------------+
      if (h*forcemaxnow.gt.0.2) then
        h=h/2.0
        vflag=1.0
        write(6,*)
        write(6,*) 'Decreasing timestep to',h
        write(6,*)
      end if

c +--------------------------------------------+
c | Exceptions to testing if h should increase |
c +--------------------------------------------+
      if (numstep.eq.0)                     return  
      if ((ipka.ne.0).and.(numstep.le.200)) return
      if (mod(numstep,100).ne.0)            return

c +-----------------------------------------------------+
c | Calculate metric and increase h by 50% if warranted |
c +-----------------------------------------------------+
      write(10,*) numstep, h*forcemax
      if (h*forcemax.lt.0.04) then
        h=h*1.5
        vflag=1.0
        write(6,*)
        write(6,*) 'Increasing timestep to',h
        write(6,*)
      end if
      forcemax=0.0

      end

