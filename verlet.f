
      subroutine verlet

      include "common.f"

      dimension xnext(NMAX,3),xprev(NMAX,3), fx_prev(NMAX,3)
      save      xnext

c +---------------------------------------------------+
c |  Increment counters and store previous positions  |
c +---------------------------------------------------+
      t=t+h
      istep=istep+1
      numstep=numstep+1

c +-----------------------------+
c |  Steepest Descent Algorithm |
c +-----------------------------+
      if (steepestdescent) then
        call force
        do i=ifix+1,natom
          do ind=1,3
            x(i,ind)= x(i,ind) + h*fx(i,ind)
            vx(i,ind)= fx(i,ind)
          end do
        end do

        if (sphereconstrain) call spherecoord(x)

        call properties
        return
      end if

c +-----------------------------+
c |  Velocity-Verlet Integrator |
c +-----------------------------+
      if (velocityverlet) then
        ! Force already available from previous step (or init for first step)
        !
        ! x(t+h)=x(t) + h*v(t) + 0.5*h^2*a(t)
        do i=ifix+1,natom
          do ind=1,3
            x(i,ind)= x(i,ind) + h*vx(i,ind) + 0.5d0*h*h*fx(i,ind)
            fx_prev(i,ind)= fx(i,ind)
          end do
        end do

        if (sphereconstrain) call spherecoord(x)

        ! v(t+h)=v(t)=0.5*h*[a(t)+a(t+h)]
        call force
        do i=ifix+1,natom
          do ind=1,3
            vx(i,ind)= vx(i,ind) + 0.5d0*h*(fx(i,ind) + fx_prev(i,ind))
          end do
        end do

        call properties
        return
      end if


      do i=ifix+1,natom
	do ind=1,3
          xprev(i,ind) = x(i,ind)
	end do
      end do

c +----------------------------------+
c |  R(t+h)=R(t)+V(t)*h+F(t)*h**2/2  |
c +----------------------------------+
      if (vflag.gt.0.5) then
        do i=ifix+1,natom
	  do ind=1,3
            xnext(i,ind)=x(i,ind) + vx(i,ind)*h + fx(i,ind)*h*h*0.5d0
	  end do
        end do
        vflag=0.0d0 

        if (sphereconstrain) call spherecoord(xnext)
      end if

c +---------------------------------------+
c |  assign R(t) = R(t+h) from last step  |
c +---------------------------------------+
      do i=ifix+1,natom
	do ind=1,3
          x(i,ind)=xnext(i,ind)
	end do
      end do

c  +-------------------------------------------------------+
c  |  Calc F(t+h), and thus find R(t+2h) and V(t+h) using  |
c  |           R(t+h)=2R(t)-R(t-h)+F(t)*h**2               |
c  |           V(t)=(R(t+h)-R(t-h))/2h                     |
c  +-------------------------------------------------------+
c     write(*,*)'In Verlet...'
      call force

      do i=ifix+1,natom
	do ind=1,3
          xnext(i,ind)= 2.0d0*x(i,ind) - xprev(i,ind) + fx(i,ind)*h*h
          vx(i,ind)= (xnext(i,ind) - xprev(i,ind))/(2.0d0*h)
	end do
      end do

      if (sphereconstrain) call spherecoord(xnext)

      call stress(0)
      call properties
      tav=tav + tempk

      end

