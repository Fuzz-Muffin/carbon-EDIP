
      subroutine defaults

      include "common.f"

c These defaults are according to 
c   a) carbon/edip/param,
c   b) carbon/edip/run4/all4o, and
c   c) carbon/edip/potential

c Two-body potential parameters
      aa=20.0853862863495d0
      bb=0.827599951322299d0
      beta=0.0490161172713279d0
      sigma=1.25714643580808d0
      a1=1.89225338775144d0
      a2=0.169794491000172d0

c Three-body potential parameters
      gamma=1.2406975366223d0
      xlam=53.7116179513016d0
      qq=3.5d0
      xmu=0.25d0
      pilam=0.0d0

c Counting the coordination
      zlow=1.46d0
      zhigh=2.274108d0
      zalpha=1.544583d0

      flow=zlow
      fhigh=2.0d0
      falpha=zalpha

c Pi-bonding parameters
      zdih=0.30d0
      zrep=0.06d0
      zrep2=0.06d0
      c0=3.2d0

c Coordination cutoff
      bondcutoff=1.85d0

c Primary Knock-on Atom 
      ipka=0
      epka=0.0
      do ind=1,3
        xpka(ind)=1.0
      end do

      c60pka=.false.
      xepka=.false.
      xebox=-1.0d0

c Control File defaults
      npass=0

      varlist(0,1)=0.05d0 ! h
      varlist(0,2)=1000   ! nstep
      varlist(0,3)=300d0  ! temp
      varlist(0,4)=1      ! therm
      varlist(0,5)=1      ! gr
      varlist(0,6)=1      ! msd

      computeshear=.false.
      shear=0.0

      nolinear=.false.
      nodihedral=.false.
      norepulsion=.false.

      usedihedral2=.false.
      usedihedral4=.false.
      newgcutoff=.false.
      steepestdescent=.false.
      velocityverlet=.false.
      xbspbc=.false.
      sphereconstrain=.false.
      sphereradius=-1.0
      cellneighbour=.false.
      zbl=.false.
      avas=.false.
      vmd=.false.

      numcylinder=0

      temp_start=-1.0

      timemax=-1.0

      variablestep=.false.

      mass=.false.

      cyclo=-1.0

      smin=0.0
      smax=-1.0
      calcstress=.false.

      norings=.true.
      special=.false.
      static=.false.

      ifix=0

      nswap=0
      nslab=10

      nrandom=0

      end
