
      subroutine coordination(i)

      include 'common.f'


      z(i)=zz(i)

c +----------------------------+
c | Setup neighbours of atom i |
c +----------------------------+
      do jj=1,num(i)
	if ((dr(i,jj).gt.zlow).and.(dr(i,jj).lt.zhigh)) then
	  finiteforce(jj)=.true.
	else
	  finiteforce(jj)=.false.
	end if
        do ind=1,3
          dzdx(jj,ind)=dzz(i,jj)*dxdr(i,jj,ind)
        end do
      end do

c +-------------------------------------+
c | Setup neighbours-of-neighbours of i |
c +-------------------------------------+
      do jj=1,num(i)
      do kk=1,num(near(i,jj))
        finiteforce2(jj,kk)=.false.
        do ind=1,3
          dzdxx(jj,kk,ind)=0.0d0
        end do
      end do
      end do

c +-------------------------------+
c | Call sp/sp2 terms if required |
c +-------------------------------+
      call dihedral(i)
      call repulsion(i)
      call linear(i)

      end

