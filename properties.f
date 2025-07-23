
      subroutine properties

      include "common.f"
      dimension vsum(3), summ(3), deltav(3)

c +------------------------------------------------+
c | Calculate the Energy, Momentum and Temperature |
c +------------------------------------------------+
      mv=0.0
      eke=0.0

      do ind=1,3
        vxsum=0.0
        do i=ifix+1,natom
          vxsum=vxsum + fracmass(i)*vx(i,ind)
          eke=eke + 0.5*fracmass(i)*vx(i,ind)**2
        end do
        mv=mv + vxsum**2
      end do

      mv=dabs(mv)
      tempk=eke/(dfloat(natom-ifix))*tfac


      ! For Xe atom, compute temperature of carbon atoms with Z>0.01 relative to centre-of-mass motion
      if (xepka) then
        ncluster=0
        do ind=1,3
          vsum(ind)=0.0d0
          summ(ind)=0.0d0
        end do

        do i=1,natom
          if (z(i).gt.0.01d0) then
            ncluster=ncluster+1
          
            do ind=1,3
              vsum(ind)=vsum(ind) + fracmass(i)*vx(i,ind)
              summ(ind)=summ(ind) + fracmass(i)
            end do
          end if
        end do

        do ind=1,3
          deltav(ind)=vsum(ind)/summ(ind)
        end do

        ekecluster=0.0d0
        do i=1,natom
          if (z(i).gt.0.01d0) then
            do ind=1,3
              vnet= vx(i,ind) - deltav(ind)
              ekecluster=ekecluster + 0.5*fracmass(i)*vnet**2
            end do
          end if
        end do
        tempk=ekecluster/(dfloat(ncluster))*tfac
      end if

      end
 
