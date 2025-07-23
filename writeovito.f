      subroutine writeovito

      ! Write output coordinates for ovito, modified from writeavas.f
      ! -fil

      include "common.f"
      
      real*8:: xx(3)

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
      write(8,200) 
      write(8,210) (box(ind),ind=1,3)
      write(8,220) 
      write(8,230) 

      do i=1,natom
        do ind=1,3
          xx(ind)= x(i,ind)-box(ind)*0.5     ! shift centre of box to origin
          xx(ind)=xx(ind)-box(ind)*nint(xx(ind)/box(ind))  ! apply pbc wrapping
          xx(ind)=xx(ind)+box(ind)*0.5       ! shift origin back to centre of box
        end do
        write(8,300) i,(xx(ind),ind=1,3),z(i)
      end do
      close(unit=8)

 200  format("Lattice=""",$)
 210  format(F10.3," 0.0 0.0 0.0 ",F10.3," 0.0 0.0 0.0 ",F10.3,$)
 220  format(" """,$)
 230  format(' Properties=species:S:1:id:I:1:pos:R:3:coord:R:1')

 300  format('C ',I6,3(1X,F10.3),1X,F3.1)

      end
