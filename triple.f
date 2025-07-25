
      subroutine triple(i)

      include "common.f"

      dimension dcosj(3), dcosk(3)

c +-------------------+
c | Get thread number |
c +-------------------+
      icpu=1
!$    icpu=OMP_GET_THREAD_NUM()+1

c +-----------------------------------+
c | Ideal angle varies according to Z |
c +-----------------------------------+
      t1=-6.0d0*2.5d0
      t2= 6.0d0

      fac= 1.0d0/12.0d0 * (1.0d0 + tanh(t1+t2*z(i)))
      tau=1.0d0 - z(i)*fac
      dtau= -z(i)*t2/(12.0d0*cosh(t1+t2*z(i))**2) - fac

      do jj=1,num(i)-1
      if (dr(i,jj).lt.a1+a2*z(i)-0.001) then

        do kk=jj+1,num(i)
        if (dr(i,kk).lt.a1+a2*z(i)-0.001) then

c +----------------------------------------------+
c | Calculate cos(theta_jik) and its derivatives |
c +----------------------------------------------+
          cosi=0.0
	  do ind=1,3
	    cosi=cosi + dxdr(i,jj,ind)*dxdr(i,kk,ind)
	  end do

	  do ind=1,3
	    dcosj(ind)=(dxdr(i,kk,ind) - cosi*dxdr(i,jj,ind))/dr(i,jj)
	    dcosk(ind)=(dxdr(i,jj,ind) - cosi*dxdr(i,kk,ind))/dr(i,kk)
	  end do

c +---------------------------------------------------+
c | Terms for U3 energy and constants for force loops |
c +---------------------------------------------------+
          arg1=1.0d0/(dr(i,jj)-a1-a2*z(i))
          arg2=1.0d0/(dr(i,kk)-a1-a2*z(i))

 	  zmu= 3.0d0 + 0.069d0/xmu    ! 0.069=log(lam4/lam3)
          zexpgam= exp(gamma*(arg1+arg2) - xmu*(z(i)-4.0)*(z(i)-zmu))

 	  expcos= exp(-qq*(cosi+tau)**2)
          theta=  xlam/qq * (1.0d0 - expcos)
 	  dtheta= xlam * expcos * 2.0d0*(cosi+tau)

          arg11= zexpgam * theta * gamma*arg1*arg1
          arg22= zexpgam * theta * gamma*arg2*arg2
          argz= -zexpgam * theta * xmu * (2.0d0*z(i) - 4.0d0 - zmu)
          f33=   zexpgam * dtheta*dtau + a2*(arg11 + arg22) + argz

          u3i(i,icpu)=u3i(i,icpu) + zexpgam*theta

c +-------------------------------+
c | Forces due to neighbours of i |
c +-------------------------------+
          do mm=1,num(i)
	  if ((jj.eq.mm).or.(kk.eq.mm).or.finiteforce(mm)) then
	    m=near(i,mm)

            do ind=1,3
              dri=arg11*(a2*dzdx(mm,ind) - kron(jj,mm)*dxdr(i,jj,ind))
              drj=arg22*(a2*dzdx(mm,ind) - kron(kk,mm)*dxdr(i,kk,ind))
	      drijz= dri + drj + argz*dzdx(mm,ind)

	      cosjk= dcosj(ind)*kron(jj,mm) + dcosk(ind)*kron(kk,mm)
	      f3= drijz + zexpgam*dtheta*(cosjk + dtau*dzdx(mm,ind))

              fxx(i,ind,icpu) = fxx(i,ind,icpu) - f3
              fxx(m,ind,icpu) = fxx(m,ind,icpu) + f3
	      !if (calcstress) call stress(i,mm,f3,ind)
            end do

c +---------------------------------------------+
c | Forces due to neighbours-of-neighbours of i |
c +---------------------------------------------+
	    do nn=1,num(m)
	    if (finiteforce2(mm,nn)) then
              n=near(m,nn)
              do ind=1,3
                f3= f33*dzdxx(mm,nn,ind)
                fxx(m,ind,icpu) = fxx(m,ind,icpu) - f3
                fxx(n,ind,icpu) = fxx(n,ind,icpu) + f3
	        !if (calcstress) call stress(m,nn,f3,ind)
              end do
	    end if
	    end do

          end if
          end do
	end if
        end do
      end if
      end do

      end
   
