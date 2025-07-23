
      subroutine xecheckcell(i,ix,iy,iz,jstart)

      include "common.f"

      ! No wrap-around cell indices
      if (ix.eq.0) return
      if (iy.eq.0) return
      if (iz.eq.0) return

      if (ix.gt.numcells(1)) return
      if (iy.gt.numcells(2)) return
      if (iz.gt.numcells(3)) return

      ! Loop over the required atoms in the cell
      do jloop=jstart,ncell(ix,iy,iz)
        j=icell(ix,iy,iz,jloop)

        call checkrij(i,j)
      end do

      end
