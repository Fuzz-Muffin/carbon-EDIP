
      subroutine force

      include "common.f"

c +-----------------------+
c | Get Number of Threads |
c +-----------------------+
      ncpu=1
!$    ncpu=OMP_GET_MAX_THREADS()

c +------------------------+
c | Clear temporary arrays |
c +------------------------+
      do i=1,natom
        do icpu=1,ncpu
          u2i(i,icpu)=0.0d0
          u3i(i,icpu)=0.0d0
        end do
      end do

      do i=1,natom
      do ind=1,3
          do icpu=1,ncpu
            fxx(i,ind,icpu)=0.0d0
          end do
        end do
      end do

c +--------------------------+
c | Special case for Xe atom |
c +--------------------------+
      if (xepka) natom=natom-1

c +--------------------------+
c | Main loop over the atoms |
c +--------------------------+
      call neighbour
      call distance
      call cutoff

      if (zbl) then
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i)
        do i=1,natom
          call coordination(i)
          call pair_zbl(i)
          call triple(i)
        end do
!$OMP END PARALLEL DO

      else
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i)
        do i=1,natom
          call coordination(i)
          call pair(i)
          call triple(i)
        end do
!$OMP END PARALLEL DO
      end if

c +--------------------------+
c | Manually reduce energies |
c +--------------------------+
      u2=0.0d0
      u3=0.0d0
      do i=1,natom
        peatom(i)=0.0d0
        do icpu=1,ncpu
          u2=u2 + u2i(i,icpu)
          u3=u3 + u3i(i,icpu)
          peatom(i)=peatom(i) + u2i(i,icpu) + u3i(i,icpu)
        end do
      end do
      pe=u2+u3

c +------------------------+
c | Manually reduce forces |
c +------------------------+
      do i=1,natom
        do ind=1,3
          fx(i,ind)=0.0d0
          do icpu=1,ncpu
            fx(i,ind)=fx(i,ind) + fxx(i,ind,icpu)
          end do
        end do
      end do

c +------------------------+
c | Special Routine for Xe |
c +------------------------+
      if (xepka) then
        natom=natom+1
        call xeforce
      end if

c +--------------------------------+
c | Scale forces according to mass |
c +--------------------------------+
      do i=1,natom
        do ind=1,3
          fx(i,ind)=fx(i,ind)/fracmass(i)
        end do
      end do

      if (sphereconstrain) call sphereforce

      if (slabconstrain) then
        do i=1,natom
           fx(i,3)=0.0d0
        end do
      end if

      end
