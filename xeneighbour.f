
      subroutine xeneighbour

      include "common.f"

      parameter (rskin=0.10)
      dimension xold(nmax,3), widthcell(3), ipos(3)
      save      xold

c +-----------------------------+
c | Create file for diagnostics |
c +-----------------------------+
      if (numstep.eq.0) then 
        open(unit=11,file='neighbour.txt',status='unknown')
        write(11,*) 'Using xeneighbour routine'
        write(11,*) 'xebox=',xebox
      end if

c +--------------------------------------------+
c |  Check to see if list needs recalculation  |
c +--------------------------------------------+
      if ((numstep.ne.0).and.(.not.special)) then
        rmovemax=0.0
        do 100 i=1,natom
          ! Atoms not counted in neighbour list don't trigger rebuilds
          do ind=1,3
            if (abs(x(i,ind)).gt.0.5d0*xebox) goto 100
          end do        

          rmove=0.0
          do ind=1,3
            rmove=rmove + (x(i,ind)-xold(i,ind))**2
          end do
          rmovemax=max(rmovemax,rmove)
 100    continue

        if (dsqrt(rmovemax).lt.rskin) return
      end if

c +------------------------+
c |  Clear existing lists  |
c +------------------------+
      rcut=c0+2.0*rskin
      if (special) rcut=c0

      npair=0
      do i=1,natom
        num(i)=0
      end do

      if (cellneighbour) then
c +-------------------------------------------------+
c |  O(N) loop to create the list of indexed pairs  |
c +-------------------------------------------------+
        ! Identify the number of cells for each dimension
        do ind=1,3
          numcells(ind)=int(xebox/(rcut+0.1))
          widthcell(ind)=xebox/real(numcells(ind))

          if (numcells(ind).gt.MAXCELL) stop 'Increase MAXCELL'
        end do

        ! Clear array containing number of atoms in each cell
        do ix=1,numcells(1)
          do iy=1,numcells(2)
            do iz=1,numcells(3)
              ncell(ix,iy,iz)=0
            end do
          end do
        end do

        ! Populate atoms inside xebox into the cells (map coordinate into [0,xebox) first)
        do 110 i=1,natom
          do ind=1,3
            xx= x(i,ind) + 0.5*xebox
            if ((xx.lt.0.0d0).or.(xx.gt.xebox)) goto 110

            ipos(ind)=int(xx/widthcell(ind))+1
            if (ipos(ind).gt.numcells(ind)) ipos(ind)=numcells(ind)
          end do

          numincell= ncell(ipos(1),ipos(2),ipos(3)) + 1
          if (numincell.gt.MAXINCELL) stop 'Increase MAXINCELL'

          ncell(ipos(1),ipos(2),ipos(3))= numincell
          icell(ipos(1),ipos(2),ipos(3),numincell) = i
 110    continue

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ix,iy,iz,iloop,i)
        ! Loop over cells
        do ix=1,numcells(1)
          do iy=1,numcells(2)
            do iz=1,numcells(3)
              ! Loop over atoms in the cell
              do iloop=1,ncell(ix,iy,iz)
                i=icell(ix,iy,iz,iloop)

                call xecheckcell(i,ix,iy,iz,iloop+1) ! other atoms in the same cell

                ! This group of three have no unique inverses
                call xecheckcell(i,ix+1,iy,iz,1)     ! (1,0,0)
                call xecheckcell(i,ix,iy+1,iz,1)     ! (0,1,0)
                call xecheckcell(i,ix,iy,iz+1,1)     ! (0,0,1)

                ! This group of three have a single unique inverse each
                call xecheckcell(i,ix+1,iy+1,iz,1)   ! (1,1,0)
                call xecheckcell(i,ix+1,iy,iz+1,1)   ! (1,0,1)
                call xecheckcell(i,ix,iy+1,iz+1,1)   ! (0,1,1)

                call xecheckcell(i,ix-1,iy+1,iz,1)   ! (-1,1,0)
                call xecheckcell(i,ix+1,iy,iz-1,1)   ! (1,0,-1)
                call xecheckcell(i,ix,iy-1,iz+1,1)   ! (0,-1,1)

                ! This direction has three unique inverses
                call xecheckcell(i,ix+1,iy+1,iz+1,1)   ! (1,1,1)

                call xecheckcell(i,ix-1,iy+1,iz+1,1)   ! (-1,1,1)
                call xecheckcell(i,ix+1,iy-1,iz+1,1)   ! (1,-1,1)
                call xecheckcell(i,ix+1,iy+1,iz-1,1)   ! (1,1,-1)
              end do
            end do
          end do
        end do
!$OMP END PARALLEL DO

c +---------------------------------------------------+
c |  O(N^2) loop to create the list of indexed pairs  |
c +---------------------------------------------------+
      else
!$OMP PARALLEL DO PRIVATE(i,j)
        do i=1,natom-1
          do j=i+1,natom
            call checkrij(i,j)
          end do
        end do
!$OMP END PARALLEL DO
      end if

c +-----------------------+
c | Check no NNN overflow |
c +-----------------------+
      numflag=0
      maxnumi=0
      do i=1,natom
        maxnumi=max(num(i),maxnumi)
        if (num(i).gt.NNN) then
          numflag=1
          write(6,*) 'atom',i,' num(i)=',num(i)
        end if
      end do
      if (numflag.eq.1) stop 'Increase NNN'
      write(11,*) 'numstep/maxnumi=',numstep,maxnumi
      flush(11)

c +----------------------------+
c |  Save current coordinates  |
c +----------------------------+
      do i=1,natom
        do ind=1,3
          xold(i,ind)=x(i,ind)
        end do
      end do

      end

