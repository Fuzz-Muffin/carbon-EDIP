
      subroutine spherecoord(xx)

      include "common.f"

      dimension xx(NMAX,3)

c +------------------------------------------------------+
c | Project coordinates onto sphere of prescribed radius |
c +------------------------------------------------------+
      do i=1,natom
        radius=sqrt(xx(i,1)**2 + xx(i,2)**2 + xx(i,3)**2)

        do ind=1,3
          xx(i,ind)=xx(i,ind)/radius * sphereradius
        end do
      end do


      end
 
