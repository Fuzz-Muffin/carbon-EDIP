
      subroutine sphereforce

      include "common.f"

c +-----------------------------------------------------+
c | Assign sphereradius automatically if it is negative |
c +-----------------------------------------------------+
      if (sphereradius.lt.0.0) then
        radius=0.0
        do i=1,natom
          radius=radius + sqrt(x(i,1)**2 + x(i,2)**2 + x(i,3)**2)
        end do
        sphereradius=radius/real(natom)
        write(6,*) 'sphereradius=',sphereradius
      end if

c +--------------------------------------+
c | Remove normal component of the force |
c |      F'= F - rhat(F.rhat)            |
c +--------------------------------------+
      do i=1,natom
        dot=0.0
        do ind=1,3
          dot=dot + fx(i,ind)*x(i,ind)/sphereradius
        end do

        do ind=1,3
          fx(i,ind)=fx(i,ind) - dot*x(i,ind)/sphereradius
        end do
      end do

      end
 
