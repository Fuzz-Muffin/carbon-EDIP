
      subroutine init
c     use IFPORT

      include "common.f"

      logical keep
      dimension xcopy(NMAX,3), vcopy(NMAX,3)

c +-----------------+
c | g(r) statistics |
c +-----------------+
      do ix=1,NGR
        nframe(ix)=0
        do i=1,numbin
          nbin(ix,i)=0.0
        end do
      end do

c +------------------+
c | Angle statistics |
c +------------------+
      do i=1,NGR
        do j=1,3
          do k=1,numbin2
            nbin2(i,j,k)=0
          end do
        end do
      end do

c +---------------+
c | Stress Tensor |
c +---------------+
      do ind1=1,3
	do ind2=1,3
          totstr(ind1,ind2)= 0.0
	end do
      end do

c +-----------------+
c | Kronecker Delta |
c +-----------------+
      do i=1,NNN
        do j=1,NNN
          if (i.ne.j) kron(i,j)=0.0
          if (i.eq.j) kron(i,j)=1.0
        end do
      end do

c +--------------------------+
c | Skip some random numbers |
c +--------------------------+
      do loop=1,nrandom
        call random_number(r1)
!	r1=ran()   ! for use with gfortran
      end do

c +---------------------------------------+
c | Set the temperature conversion factor |
c +---------------------------------------+
      ndim=3
      tfac=2.0d0/ndim / BOLTZMANN

c +-------------------------+
c | Initialise some numbers |
c +-------------------------+
      vflag=0.0d0
      nstress=0
      numstep=0
      iseed=31415
      elosttherm=0.0d0
      t=t0

c +--------------------------------+
c | Delete atoms inside a cylinder |
c +--------------------------------+
      if (numcylinder.gt.0) then
        do i=1,natom
          do ind=1,3
            xcopy(i,ind)= x(i,ind)
            vcopy(i,ind)=vx(i,ind)
          end do
        end do
        ncopy=natom

        natom=0
        do i=1,ncopy
        keep=.true.
          do loop=1,numcylinder
            dxx=xcopy(i,1) - cylinder(loop,1)
            dyy=xcopy(i,2) - cylinder(loop,2)
  
            dxx=dxx - box(1)*nint(dxx/box(1))
            dyy=dyy - box(2)*nint(dyy/box(2))
            dist=sqrt(dxx**2 + dyy**2)
            if (dist.lt.cylinder(loop,3)) keep=.false.
          end do

          if (keep) then
            natom=natom+1
            do ind=1,3
               x(natom,ind)=xcopy(i,ind)
              vx(natom,ind)=vcopy(i,ind)
            end do
          end if
        end do

        write(6,*) 'Created ',numcylinder,' cylinders'
        write(6,*) 'Deleted ',ncopy-natom,' atoms'
      end if

c +-------------------------+
c | Set initial temperature |
c +-------------------------+
      if (temp_start.gt.0.0) then
       write(6,*) ''
       write(6,*) 'Setting velocities from a Maxwell-Boltzmann'
       write(6,*) 'distribution corresponding to a temperature'
       write(6,*) 'of',temp_start,'Kelvin'
       write(6,*) ''
       do i=ifix+1,natom
         do ind=1,3
           ! Generate number from normal distribution (Box-Muller method)
           call random_number(r1)
           call random_number(r2)
!	   r1=ran()   ! for use with gfortran
!	   r2=ran()   ! for use with gfortran
           r_normal=sqrt(-2.0*log(r1))*cos(2*pi*r2)

           ! Scale normal distribution; sigma=sqrt(kT/m)
           vx(i,ind)=sqrt(0.025*temp_start/300.0) * r_normal
         end do
       end do

       ! Subtract off net momentum
       do ind=1,3
         vsum=0.0
         summ=0.0
         do i=ifix+1,natom
           vsum=vsum + fracmass(i)*vx(i,ind)
           summ=summ + fracmass(i)
         end do
         deltav=vsum/summ

         do i=ifix+1,natom
           vx(i,ind)=vx(i,ind) - deltav
         end do
       end do

      end if

c +------------------------+
c | Setup for a Xenon atom |
c +------------------------+
      if (xepka) then
        ! Default Xe mass is Xe-133 (not the natural abundance value of 131.29)
        if (.not.mass) fracmass(natom)=132.90d0/12.01d0

        if (epka.gt.0.1d0) then
          ipka=natom
          zbl=.true.
        end if
      end if

c +------------------------------+
c | Create Primary Knock-On Atom |
c +------------------------------+
      if ((ipka.ge.ifix+1).and.(ipka.le.natom)) then
        variablestep=.true.

        ! Normalise the unit vector
        unit=0.0
        do ind=1,3
          unit= unit + xpka(ind)**2
        end do
        unit=sqrt(unit)

        ! Assign PKA with the appropriate velocity
        velocity= sqrt(2.0*epka/fracmass(ipka))
        do ind=1,3
          vx(ipka,ind)= velocity*xpka(ind)/unit
        end do

        ! Subtract off net momentum
        do ind=1,3
          vsum=0.0
          summ=0.0
          do i=ifix+1,natom
            vsum=vsum + fracmass(i)*vx(i,ind)
            summ=summ + fracmass(i)
          end do
          deltav=vsum/summ

          do i=ifix+1,natom
            vx(i,ind)=vx(i,ind) - deltav
          end do
        end do

        write(6,*)
        write(6,*) '# Creating Primary Knock-on Atom'
        write(6,*) '# Atom Number',ipka
        write(6,*) '# Energy [eV]',epka
        write(6,*) '# Direction',(real(xpka(ind)),ind=1,3)
        write(6,*) '# Setting Momentum to zero'
        write(6,*)
      end if

c +--------------------------------------------+
c | Store initial Potential and Kinetic Energy |
c +--------------------------------------------+
      call volume
      write(6,*)'Initializing...'
      call force
c     call energy
      call properties
      call writexbs
c     call stress(0)

      pestart= pe
      ekestart= eke
      end 
