      subroutine writeovito

      ! Write output coordinates for ovito, modified from writeavas.f
      ! -fil

      include "common.f"
      
      real*8:: xx(NMAX,3)

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
 100  call system("/bin/rm -f ovito.xyz")

 110  continue
      open(unit=8,file='ovito.xyz',position='append')
      
      ! Here we use the extended xyz format
      
      tnow=(t-t0)*timetau
      
      write(8,*) natom
      write(8,100) 
      write(8,110) (box(ind),ind=1,3)
      write(8,120) 
      write(8,130) 

      do i=1,natom
        do ind=1,3
          xx(i,ind)= x(i,ind)-box(ind)*nint(x(i,ind)/box(ind))
          xx(i,ind)=xx(i,ind)+box(ind)*0.5
        end do
        write(8,200) atom(i),i,(xx(i,ind),ind=1,3),z(i)
      end do
      close(unit=8)

 100  format("Lattice=""",$)
 110  format(F6.5," 0.0 0.0 0.0 ",F6.5," 0.0 0.0 0.0 ",F6.5,$)
 120  format(" """,$)
 130  format(' Properties=species:I:1:id:I:1:pos:R:3:coord:R:1')

 200  format(I6,1X,I6,3(1X,F6.5),1X,F6.5)

      end
