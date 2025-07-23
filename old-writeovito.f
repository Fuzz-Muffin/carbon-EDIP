
      subroutine writeovito

      include "common.f"

      dimension xx(NMAX,3)

c +-------------------------+
c | Decide whether to print |
c +-------------------------+
      if (nsnap.eq.0) return

      if (numstep.eq.0) goto 100
      if ((nsnap.lt.0).and.(mod(numstep,abs(nsnap)).eq.0)) goto 100
      if ((nsnap.gt.0).and.(mod(numstep,abs(nsnap)).eq.0)) goto 110
      return

c +-------------------+
c | Write header line |
c +-------------------+
 100  call system("/bin/rm -f in.xyz")

 110  continue
      open(unit=8,file='in.xyz',position='append')

      tnow=(t-t0)*timetau

      do i=1,natom
        do ind=1,3
          if (xbspbc) then
            xx(i,ind)=x(i,ind)-box(ind)*nint(x(i,ind)/box(ind))
          else
            xx(i,ind)=x(i,ind)
          end if
        end do
      end do

      write(8,*) natom
      write(8,120) pe,(box(ind),ind=1,3)

      do i=1,natom
        write(8,130) (xx(i,ind),ind=1,3)
        write(8,140) (fx(i,ind),ind=1,3)
      end do
      close(unit=8)

 120  format(F10.4,3(1X,F8.2))
 130  format($,'C',3(1X,F8.2))
 140  format(      3(1X,F8.2))

      end
