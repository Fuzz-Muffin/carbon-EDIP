
      subroutine writexyz

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
 100  call system("/bin/rm -f in.xyz")

 110  continue
      open(unit=8,file='in.xyz',position='append')

      tnow=(t-t0)*timetau

      write(8,*) natom
      write(8,120) tnow,tempk

      do 111 i=1,natom
        if ((xepka).and.(i.eq.natom)) then
          write(8,141) (x(i,ind),ind=1,3)
          goto 111
        end if
        write(8,140) (x(i,ind),ind=1,3)
 111  continue
      close(unit=8)

 115  format(3(F12.5,1X))
 120  format('time=',f8.4,' ps   temp=',f7.1,' K')
 130  format('C',3(1X,F7.2),1X,F7.2,1X,F8.3,1X,F3.1)
 140  format('C',3(1X,F7.2))
 141  format('Xe',3(1X,F7.2))

      end
