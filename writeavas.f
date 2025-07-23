
      subroutine writeavas

      include "common.f"


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
 100  call system("/bin/rm -f avas.xyz")
      open(unit=8,file='avas.xyz',status='unknown')
      write(8,115) (box(ind),ind=1,3)
      close(unit=8)

 110  continue
      open(unit=8,file='avas.xyz',position='append')

      tnow=(t-t0)*timetau

      write(8,*) natom
      write(8,120) tnow,tempk

      do i=1,natom
        atomke=0.0
        do ind=1,3
          atomke=atomke + 0.5*vx(i,ind)**2
        end do
        if (steepestdescent) then
          write(8,130) (x(i,ind),ind=1,3),atomke,peatom(i),zz(i)
        else
          write(8,140) (x(i,ind),ind=1,3),atomke,peatom(i)
        end if
      end do
      close(unit=8)

 115  format(3(F12.5,1X))
 120  format('time=',f8.4,' ps   temp=',f7.1,' K')
 130  format('C',3(1X,F7.2),1X,F7.2,1X,F8.3,1X,F3.1)
 140  format('C',3(1X,F7.2),1X,F7.2,1X,F8.3)

      end
