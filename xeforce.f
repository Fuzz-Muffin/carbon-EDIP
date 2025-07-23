      subroutine xeforce

      include "common.f"
      dimension dxx(3)
      double precision ljexp, ljsf

      include "xezblpar.f"

c +-------------------+
c | Get thread number |
c +-------------------+
c     icpu=1
!$    icpu=OMP_GET_THREAD_NUM()+1

! Need to activate parallel threads
! Or does it need to be serial ??
! If calculation done before energy/force reduction can split across CPUs.

C Need to manually clear out arrays for the Xe atom
      do ind=1,3
        fx(natom,ind)=0.0d0
      end do

      j=natom
      do 100 i=1,natom-1
        do ind=1,3
          delta= x(i,ind) - x(j,ind)
          dxx(ind)=delta - box(ind)*dnint(delta/box(ind))
        end do
        rij=dsqrt(dxx(1)**2 + dxx(2)**2 + dxx(3)**2)

        if (rij.gt.10.0d0) goto 100

c +-------------------+
c | The ZBL potential | (note: no factor of two here)
c +-------------------+
        const= (z_1*z_2*E_CHAR*1D+10)/(4d0*pi*EPSILON)
        z1= z_c1*exp(-z_d1*rij)
        z2= z_c2*exp(-z_d2*rij)
        z3= z_c3*exp(-z_d3*rij)
        z4= z_c4*exp(-z_d4*rij)
        zphi= const*(z1 + z2 + z3 + z4)
        dzphi= const*(-z_d1*z1 - z_d2*z2 - z_d3*z3 - z_d4*z4)
        
        uzbl= zphi/rij

c +------------------+
c | The LJ pair term |
c +------------------+
        fact6=(sigmaxe/rij)**6
        fact12=fact6*fact6
        ulj= 4.0d0*epsxe * (fact12 - fact6)
        dulj= 4.0d0*epsxe * (-12.0d0*fact12/rij + 6.0d0*fact6/rij)

c +-------------------------------+
c | The Fermi switching functions | (check all the parameters are right)
c +-------------------------------+
        zexp = exp((rij - zcut + zshift)/zskin)
        zblsf = 1.0d0/(zexp + 1.0d0)

        ljexp = exp((rij - zcut - zshift)/zskin)
        ljsf = 1.0d0-1.0d0/(ljexp + 1.0d0)

        dzblsf = -zexp/(zskin*(zexp+1.0d0)**2)
        dljsf  = ljexp/(zskin*(ljexp+1.0d0)**2)

c +--------------------------------------------+
c | Add potential energy to relevant variables |
c +--------------------------------------------+
        uxe= uzbl*zblsf + ulj*ljsf

        pe=pe + uxe
        peatom(i)=peatom(i) + 0.5d0*uxe
        peatom(j)=peatom(j) + 0.5d0*uxe

c +----------------------------------------------+
c | Calculate the ZBL force using the chain rule |
c +----------------------------------------------+
        fzbl=(dzphi/rij-zphi/(rij*rij))*zblsf+(dzblsf)*(zphi/rij)

c +----------------------------------------------+
c | Calculate the LJ force using the chain rule |
c +----------------------------------------------+
        flj= ulj*dljsf + dulj*ljsf

c Hard code just the LJ term
c       fzbl=0.0d0
c       flj=dulj
c       pe=ulj
c       write(6,*) 'rij,pe,flj=',rij,pe,flj

c Hard code just the ZBL term
c       flj=0.0d0
c       fzbl=dulj
c       fzbl=dzphi/rij-zphi/(rij*rij)
c       pe=uzbl

c +-----------------------------+
c | Force between atoms i and j |
c +-----------------------------+
        do ind=1,3
          dxr= dxx(ind)/rij
          totalf = fzbl*dxr + flj*dxr
!          write(6,*) 'i,j,ind,totalf=',i,j,ind,totalf

!          write(6,*) 'Before: fx(i,ind)=',i,fx(i,ind)
!          write(6,*) 'Before: fx(j,ind)=',j,fx(j,ind)

          fx(i,ind) = fx(i,ind) - totalf
          fx(j,ind) = fx(j,ind) + totalf

!          write(6,*) 'After: fx(i,ind)=',i,fx(i,ind)
!          write(6,*) 'After: fx(j,ind)=',j,fx(j,ind)
        end do

c Hard code energy & force for testing
c       do ind=2,3
c         fx(1,ind)=0.0d0
c         fx(2,ind)=0.0d0
c       end do
c       pe= 100.0d0/rij**2
c       fx(1,1)= -200d0/rij**3
c       fx(2,1)=  200d0/rij**3

        !write(6,*) 'rij,uzbl,ulj',rij,uxe,fx(1,1)
        !write(6,*) 'rij,uzbl,ulj',rij,ulj*ljsf,flj
        !write(6,*) 'rij,uzbl,ulj',rij,uzbl*zblsf,fzbl

 100  continue
              
      end
