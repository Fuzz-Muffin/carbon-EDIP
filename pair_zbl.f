      subroutine pair_zbl(i)

      include "common.f"
      include "zblpar.f"

c +-------------------+
c | Get thread number |
c +-------------------+
      icpu=1
!$    icpu=OMP_GET_THREAD_NUM()+1

c +-----------------------------------------------+
c | Loop over neighbours j to compute pair energy |
c +-----------------------------------------------+
      do jj=1, num(i)
      if (dr(i,jj).lt.(a1+a2*z(i)-0.001)) then
        j = near(i,jj)
        rij = dr(i,jj)

c +-------------------+
c | The ZBL potential |
c +-------------------+
        const= 0.5d0*(z_1*z_2*E_CHAR*1D+10)/(4d0*pi*EPSILON)
        z1= z_c1*exp(-z_d1*rij)
        z2= z_c2*exp(-z_d2*rij)
        z3= z_c3*exp(-z_d3*rij)
        z4= z_c4*exp(-z_d4*rij)
        zphi= const*(z1 + z2 + z3 + z4)
        dzphi= const*(-z_d1*z1 - z_d2*z2 - z_d3*z3 - z_d4*z4)
        
        uzbl= zphi/rij

c +--------------------+
c | The EDIP pair term |
c +--------------------+
        arg = 1.0d0/(rij-a1 - a2*z(i))
        bond = exp(-beta*z(i)*z(i))
        r4 = 1.0d0/(rij**4)
        r5 = r4/rij
        part1= bb*r4 - bond
        part2= aa*exp(sigma*arg)
        uedip = part1*part2

        ! Expressions needed for derivatives
        dedipp1 = -4.0d0*bb*r4*part2/rij 
        dedipp2 = -bb*r4*(sigma*arg**2)*part2
        dedipp3 = bond*(sigma*arg**2)*part2
        dedip = dedipp1 + dedipp2 + dedipp3

        dedzp1 = part2*(sigma*arg**2)*(a2)*part1
        dedzp2 = part2*2.0d0*beta*z(i)*bond
        dedipdz = dedzp1 + dedzp2

c +-------------------------------+
c | The Fermi switching functions |
c +-------------------------------+
        zexp = exp((rij - zcut + shift)/zskin)
        zblsf = 1.0d0/(zexp + 1.0d0)
        eexp = exp((rij - zcut - shift)/zskin)
        edipsf = 1.0d0-1.0d0/(eexp + 1.0d0)

        u2i(i,icpu)=u2i(i,icpu) + uedip*edipsf + uzbl*zblsf

        dzblsf = -zexp/(zskin*(zexp+1.0d0)**2)
        dedipsf = eexp/(zskin*(eexp+1.0d0)**2)

c +----------------------------------------------+
c | Calculate the ZBL force using the chain rule |
c +----------------------------------------------+
        fzbl=(dzphi/rij-zphi/(rij*rij))*zblsf+(dzblsf)*(zphi/rij)

c +-----------------------------+
c | Force between atoms i and j |
c +-----------------------------+
        do ind=1,3
          dxr= dxdr(i,jj,ind)

          part3= -4.0d0*bb*r5*dxr + 2.0d0*beta*z(i)*bond*dzdx(jj,ind)
          part4= -sigma*arg*arg * (dxr - a2*dzdx(jj,ind))
          f2= part2*(part3 + part1*part4)

          totalf = f2*edipsf + uedip*dedipsf*dxr + fzbl*dxr

          fxx(i,ind,icpu) = fxx(i,ind,icpu) - totalf
          fxx(j,ind,icpu) = fxx(j,ind,icpu) + totalf
          !if (calcstress) call stress(i,jj,f2,ind)
        end do

c +---------------------------------+
c | Forces due to neighbours k.ne.j |
c +---------------------------------+
        do kk=1,num(i)
        if (finiteforce(kk)) then
          k=near(i,kk)
          if (k.ne.j) then
            do ind=1,3
              f2= dedipdz*dzdx(kk,ind)*edipsf
              fxx(i,ind,icpu) = fxx(i,ind,icpu) - f2
              fxx(k,ind,icpu) = fxx(k,ind,icpu) + f2
              !if (calcstress) call stress(i,kk,f2,ind)
            end do
          end if

c +---------------------------------------------+
c | Forces due to neighbours-of-neighbours of i |
c +---------------------------------------------+
          do mm=1,num(k)
          if (finiteforce2(kk,mm)) then
            m=near(k,mm)
            do ind=1,3
              f2= dedipdz*dzdxx(kk,mm,ind)*edipsf
              fxx(k,ind,icpu) = fxx(k,ind,icpu) - f2
              fxx(m,ind,icpu) = fxx(m,ind,icpu) + f2
              !if (calcstress) call stress(k,mm,f2,ind)
            end do
          end if
          end do
  
        end if
        end do
              
      end if
      end do
              
      end
