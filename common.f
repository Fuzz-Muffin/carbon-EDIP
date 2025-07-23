
      implicit double precision (a-h,o-z)

      double precision mv

      character*80 namespecial

      logical finiteforce
      logical finiteforce2

      logical static
      logical special
      logical norings
      logical computeshear
      logical usedihedral2
      logical usedihedral4
      logical newgcutoff
      logical calcstress
      logical steepestdescent
      logical velocityverlet
      logical xbspbc
      logical cellneighbour
      logical sphereconstrain
      logical slabconstrain
      logical zbl
      logical avas
      logical vmd
      logical ovito
      logical variablestep
      logical xepka
      logical c60pka

      logical nolinear
      logical nodihedral
      logical norepulsion

      logical mass

C OpenMP parallelization related items
      integer*4 ncpu, icpu
      integer*4 OMP_GET_MAX_THREADS
      integer*4 OMP_GET_THREAD_NUM
      parameter (MAXCPU=24)

C Derived reduced unit conversion factors
C (energies in eV, distances in Angstroms, mass in 12.01 amu)
      parameter (pi=3.141592653589793238462643d0)
      parameter (AMASS=12.01d0*1.66D-27)
      parameter (BOLTZMANN=8.617D-5)      !!! consider 0.025eV/300K !!!
      !parameter (TIMETAU=1e12*1e-10*(AMASS/1.6e-19)**0.5)
      !definition of timetau has moved to constants.f

      parameter (tau1=-6.0d0*2.5d0, tau2=6.0d0)

      parameter (NMAX=120000,NNN=50,NP=0.5*NNN*NMAX)
      parameter (NUMBIN=600,WIDTHBIN=0.01,NGR=2,NVAR=6)
      parameter (NUMBIN2=181,WIDTHBIN2=1.0)
      parameter (MAXCELL=100,MAXINCELL=50)
      parameter (MAXSLAB=1000,MAXCYL=300)

c Parameters for EDIP
      common /PARAM1/ aa,bb,beta,sigma,a1,a2
      common /PARAM2/ qq,xlam,xmu,gamma,pilam
      common /PARAM3/ zlow,zhigh,zalpha
      common /PARAM4/ flow,fhigh,falpha
      common /PARAM5/ zdih,zrep,zrep2,c0,bondcutoff

c Parameters for control file
      common /ENTRY/  varlist(0:100,NVAR),npass,nsnap

      common /ELOST/  elosttherm,nrescale
      common /STATE/  u2,u3,udih,pe,tempk,eke,mv,pestart,ekestart,tav
      common /STEPS/  t0,t,h,numstep,istep,nstep,nloop,nprint,itherm
      common /EXTRA/  timetau,tfac,vflag,startt,ntakof,iseed,imsd
      common /EXTRA2/ timemax,ipass

      common /SHEAR1/ computeshear
      common /SHEAR2/ shear

      common /FLAGS1/ nolinear,nodihedral,norepulsion,usedihedral2
      common /FLAGS2/ norings,static,special
      common /FLAGS3/ newgcutoff,calcstress,usedihedral4
      common /FLAGS4/ steepestdescent,xbspbc,cellneighbour
      common /FLAGS5/ namespecial,velocityverlet,variablestep
      common /FLAGS6/ zbl,avas,vmd,ovito,mass,xepka,c60pka

      common /NATOM/  natom
      common /BOXL/   box(3)
      common /ZZZZ/   z(NMAX), zz(NMAX)
      common /ZDERV/  dzdx(NNN,3), dzdxx(NNN,NNN,3)
      common /ZFORC/  finiteforce(NNN), finiteforce2(NNN,NNN)
      common /COORD/  x(NMAX,3), vx(NMAX,3), x0(NMAX,3), fracmass(NMAX)
      common /DRXYZ/  dr(NMAX,NNN),dx(NMAX,NNN,3),dxdr(NMAX,NNN,3)
      common /STRES/  totstr(3,3),pressure,biaxial,smin,smax,nstress
      common /BOXVOL/ vol
      common /COMM1/  fx(NMAX,3)

      common /CELLS/  ncell(MAXCELL,MAXCELL,MAXCELL), numcells(3),
     *                icell(MAXCELL,MAXCELL,MAXCELL,MAXINCELL)
      common /PAIRS/  ipair(NP),jpair(NP),inum(NP),jnum(NP),npair
      common /INEAR/  near(NMAX,NNN),num(NMAX),kron(NNN,NNN)

      common /ZCYCL/  f(NMAX,NNN), df(NMAX,NNN), dzz(NMAX,NNN)
      common /GATOM/   g(NMAX),  g2(NMAX),  g3(NMAX)
      common /DGATM/  dg(NMAX), dg2(NMAX), dg3(NMAX)
      common /COM01/  cyclo, temp_start

      common /GRDAT/ nbin(NGR,NUMBIN),nframe(NGR),igr
      common /ANGLE/ nbin2(NGR,4,NUMBIN2)
      common /FIXED/ ifix,nrandom
      common /CONDU/ nslab,nswap
      common /SPHER/ sphereradius, sphereconstrain, slabconstrain
      common /PKEVNT/ xpka(3),epka,xebox,ipka

      common /CYLDAT/ cylinder(MAXCYL,3),numcylinder

c Arrays used for reducing energies and forces in parallel sections
      common /EREDC/  u2i(NMAX,MAXCPU), u3i(NMAX,MAXCPU), peatom(NMAX)
      common /FREDC/  fxx(NMAX,3,MAXCPU)

!$OMP THREADPRIVATE (/ZDERV/,/ZFORC/)
