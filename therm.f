
      subroutine therm

      include "common.f"

      save temp0, dtherm

c +---------------------------------+
c | Decide what do do - if anything |
c +---------------------------------+
      if (steepestdescent) return

      if (nloop.eq.1) then
        if (numstep.eq.1) temp0= tempk
        if (itherm.eq.2) dtherm= (startt-temp0)/nstep
        if (itherm.eq.3) dtherm= (startt/temp0)**(1.0/nstep)
      end if

      if (itherm.eq.0) temp0=tempk
      if (itherm.eq.1) temp0=startt
      if (itherm.eq.2) temp0=temp0 + dtherm
      if (itherm.eq.3) temp0=temp0 * dtherm
      if (itherm.eq.4) call edgetherm

      if (itherm.eq.0)    return
      if (itherm.eq.4)    return
c     if (nloop.eq.nstep) return

c +--------------------------------------+
c | Do the velocity rescale as specified |
c +--------------------------------------+
      call properties
      tempfrac= abs(tempk/temp0 - 1.0)
      if (tempfrac.lt.1.5/sqrt(dfloat(natom-ifix))) return
c     if (tempfrac.lt.1.0/sqrt(dfloat(natom))) return

      vflag=1.0
      ekeold=eke
      tempkold=tempk
      vfac=dsqrt(temp0/tempk)

      do i=ifix+1,natom
        do ind=1,3
          vx(i,ind)=vx(i,ind)*vfac
        end do
      end do

c +-------------------------------------+
c | Find how much energy lost or gained |
c +-------------------------------------+
      call properties
      elosttherm=elosttherm + ekeold - eke
      nrescale=nrescale+1

      end

