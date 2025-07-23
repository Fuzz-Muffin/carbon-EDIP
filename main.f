
      program carbon_edip

      include "common.f"

      call printtime(1)

c +-------------------+
c | Read control file |
c +-------------------+
      call defaults
      call readinput
      call checkinput
c disable stress calcs in parallelized version (not correctly coded)
      calcstress = .false.   
      call printinput

c +-------------------------------------+
c | Read coordinate file and initialise |
c +-------------------------------------+
      call constants
      call readcoords
      call readmasses
      call init
      call density
      call printstatus

      if (special) call runspecial

c +-----------+
c | Main loop |
c +-----------+
      do ipass=1,npass
        h=varlist(ipass,1)
        nstep=varlist(ipass,2)
        startt=varlist(ipass,3)
        itherm=varlist(ipass,4)
        igr=varlist(ipass,5)
        imsd=varlist(ipass,6)

        tav=0.0
        vflag=1.0
        nrescale=0

        do nloop=1,nstep
          call varystep
          call verlet
          call distribution
          call reflect
          call therm
          call conductivity

          call fastatom
          call writexbs
          call printstatus
          call writetheta
          call writegr

          if (nsnap.ne.0) then
            if (mod(numstep,abs(nsnap)).eq.0) call writecoords
          end if

          tnow=timetau*(t-t0)
          if ((timemax.gt.0.0).and.(tnow.gt.timemax)) then
            call printav
            call resetmsd
            goto 100
          end if
        end do

        call printav
        call resetmsd
      end do

 100  continue

c +-------------+
c | Cleaning up |
c +-------------+
      call writestress
      call writecoords
      call printrings
      call printfinish
      call printtime(2)

      end
 
