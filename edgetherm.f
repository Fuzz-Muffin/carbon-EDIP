
      subroutine edgetherm

      include "common.f"

      logical edge(NMAX)


c +--------------------------------------+
c | Decide whether to do anything at all |
c +--------------------------------------+
      if (mod(numstep,50).ne.0) return

c +-----------------------------------------+
c | Identify atoms in the thermostat region |
c +-----------------------------------------+
      ntherm=0
      do i=1,natom
        edge(i)=.false.
	do ind=1,3
	  if (x(i,ind).lt.0.0) edge(i)=.true.
	end do
	if (edge(i)) ntherm=ntherm+1
      end do

c +------------------------------------------+
c | Compute temperature in thermostat region |
c +------------------------------------------+
      eke=0.0
      do i=1,natom
        if (edge(i)) then
	  do ind=1,3
	    eke= eke + 0.5*vx(i,ind)**2
	  end do
        end if
      end do
      tempk= tfac * eke/dfloat(ntherm)
      vfac=dsqrt(startt/tempk)

c +------------------------------------------------+
c | Rescale thermostat atoms to target temperature |
c +------------------------------------------------+
      call properties
      ekeold=eke
      vflag=1.0

      do i=1,natom
	if (edge(i)) then
          do ind=1,3
            vx(i,ind)=vx(i,ind)*vfac
          end do
	end if
      end do

c +-------------------------------------+
c | Find how much energy lost or gained |
c +-------------------------------------+
      call properties
      elosttherm=elosttherm + ekeold - eke
      nrescale=nrescale+1

      end

