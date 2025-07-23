
      subroutine resetmsd

      include 'common.f'

      if (imsd.eq.1) return

      do i=1,natom
        do ind=1,3
	  x0(i,ind)=x(i,ind)
	end do
      end do

      end

