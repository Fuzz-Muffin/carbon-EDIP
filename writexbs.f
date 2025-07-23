
      subroutine writexbs

      include "common.f"

      dimension xx(NMAX,3)

      if (avas) then
        call writeavas
        return
      end if

      if (vmd) then
        call writexyz
        return
      end if

      if (ovito) then
        call writeovito
        return
      end if

c +------------------------------------------------+
c | Decide to print in.bs, in.mv or nothing at all |
c +------------------------------------------------+
      if (nsnap.ne.0) then
        do i=1,natom
          do ind=1,3
            if (xbspbc) then
              xx(i,ind)=x(i,ind)-box(ind)*nint(x(i,ind)/box(ind))
            else
              xx(i,ind)=x(i,ind)
            end if
          end do
        end do

        if (numstep.eq.0)                 goto 100
        if (mod(numstep,abs(nsnap)).eq.0) goto 110
      end if
      return

c +-----------------+
c | Write out in.bs |
c +-----------------+
 100  call system("/bin/rm -f in.mv")
      open(unit=8,file='in.bs',status='unknown')
       do i=1,natom
         if (.not.mass) then
           if (i.ne.ipka) write(8,120) (xx(i,ind),ind=1,3)
           if (i.eq.ipka) write(8,125) (xx(i,ind),ind=1,3)
         else
           if (fracmass(i).lt.1.02) then
             write(8,120) (xx(i,ind),ind=1,3)
           else
             write(8,125) (xx(i,ind),ind=1,3)
           end if
         end if
       end do

       write(8,*) 'spec C  0.2 red'
       write(8,*) 'spec C2 0.2 blue'
       if (xepka) write(8,*) 'spec C2 0.4 blue'
       write(8,*) 'bonds C* C* 0 ',bondcutoff,' 0.07 grey'
       write(8,*) 'inc 3'
       write(8,*) '* box=',(real(box(ind)),ind=1,3)
       write(8,*) 'tmat 1 0 0 0 0 1 0 1 0'
      close(unit=8)

      return

c +-----------------+
c | Write out in.mv |
c +-----------------+
 110  if (nsnap.lt.0) call system("/bin/rm -f in.mv")
      open(unit=8,file='in.mv',status='unknown',position='append')
       write(8,130) numstep,t*timetau,tempk
       write(8,140) ((xx(i,ind),ind=1,3),i=1,natom)
      close(unit=8)

 120  format('atom C ',3F9.3)
 125  format('atom C2 ',3F9.3)
 130  format(/,'frame step=',I6,', time=',f7.3,' ps, temp=',f7.1,' K')
 140  format(3F9.3)

      end
