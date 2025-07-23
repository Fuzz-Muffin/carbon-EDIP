
      subroutine printtime(iflag)

      include "common.f"

      character*24 date
      real         tdiff
      integer      tstart,time,tend,rate
      save         tstart

      call fdate(date)

      ncpu=1
!$    ncpu=OMP_GET_MAX_THREADS()

      if (iflag.eq.1) then
        call system_clock(tstart)
        write(6,100) date
        write(6,105) ncpu
      else
        if (static) numstep=1
        call system_clock(tend,rate)
        tdiff=real(tend-tstart)/real(rate)

        write(6,110) date
        write(6,120) tdiff/60.0
        write(6,130) int(tdiff/real(natom*numstep)*1e6)
        write(6,140) (t-t0)*timetau
      end if

 100  format(/,'Simulation begun ',A24)
 105  format('Running on ',I2,' cpus')
 110  format(/,'Simulation completed ',A24)
 120  format('  wall-clock time: ',F7.2,' minutes')
 130  format('  performance: ',i6,' microseconds/atom/timestep')
 140  format('  simulation time: ',F7.2,' ps')

      end

