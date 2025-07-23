
      subroutine cutoff

      include "common.f"

c +-------------------------------------------+
c | Compute spherical part of Z for all atoms |
c +-------------------------------------------+
      do i=1,natom
	zz(i)=0
        if (cyclo.gt.0.0) zz(i)=cyclo

        do jj=1,num(i)
	  rij=dr(i,jj)
	  if (rij.gt.zhigh-0.001) then
	    dzz(i,jj)=0.0d0

	  elseif (rij.lt.zlow+0.001) then
	    dzz(i,jj)=0.0d0
	    zz(i)=zz(i)+1.0

	  else
	    frac=(rij-zlow)/(zhigh-zlow)
	    recip= 1.0d0/(1.0d0-1.0d0/frac**3)
	    expz= exp(zalpha*recip)

	    zz(i)=zz(i) + expz
            dzz(i,jj)=-3.0d0/frac**4/(zhigh-zlow) * expz*zalpha*recip**2
	  end if
	end do
      end do

c +----------------------------------------------------------------------+
c | Compute f(r) for the Zpi cutoff (only needed if g(Zi) less than 4.0) |
c +----------------------------------------------------------------------+
      do i=1,natom
 	if (zz(i).lt.4.0) then
        do jj=1,num(i)
	  rij=dr(i,jj)
	  if (rij.gt.fhigh-0.001) then
	    f(i,jj)=0.0
	    df(i,jj)=0.0

	  elseif (rij.lt.flow+0.001) then
	    f(i,jj)=1.0d0
	    df(i,jj)=0.0d0

	  else
	    frac=(rij-flow)/(fhigh-flow)
	    recip= 1.0d0/(1.0d0-1.0d0/frac**3)
	    expf= exp(falpha*recip)

	    f(i,jj)=  expf
	    df(i,jj)= -3.0d0/frac**4/(fhigh-flow) * expf*falpha*recip**2
	  end if
	end do
 	end if
      end do

c +--------------------------------------+
c | Compute g3(Zi) for sp2 contributions |
c +--------------------------------------+
      do i=1,natom
	if ((zz(i).lt.2.0).or.(zz(i).gt.4.0)) then
 	   g3(i)=0.0d0
 	  dg3(i)=0.0d0

	else
	  arg= zz(i)-3.0d0
	   g3(i)=           (arg*arg - 1.0d0)**2
	  dg3(i)= 4.0d0*arg*(arg*arg - 1.0d0)
	end if
      end do

c +-------------------------------------+
c | Compute g2(Zi) for sp contributions |
c +-------------------------------------+
      do i=1,natom
         if ((zz(i).lt.1.0).or.(zz(i).gt.3.0)) then
	   g2(i)=0.0d0
	  dg2(i)=0.0d0

	else
	  arg= zz(i)-2.0d0
	   g2(i)=           (arg*arg - 1.0d0)**2
	  dg2(i)= 4.0d0*arg*(arg*arg - 1.0d0)
	end if
      end do

c +---------------------------------------+
c | Compute g(Zi) for any pi-contribution |
c +---------------------------------------+
      do i=1,natom
        if (zz(i).gt.3.0) then
	   g(i)=  g3(i)
	  dg(i)= dg3(i)

	elseif ((zz(i).gt.2.0).or.(.not.newgcutoff)) then
 	   g(i)= 1.0d0
 	  dg(i)= 0.0d0

	else
	   g(i)=  g2(i)
	  dg(i)= dg2(i)
	end if
      end do

c      do i=1,natom
c 	if (zz(i).lt.3.0) then
c 	  g(i)=1.0
c 	  dg(i)=0.0
c
c 	elseif  (zz(i).gt.4.0) then
c	  g(i)=0.0
c	  dg(i)=0.0
c
c	else
c	  arg= zz(i)-3.0
c	  g(i)=          (arg*arg - 1.0)**2
c	  dg(i)= 4.0*arg*(arg*arg - 1.0)
c	end if
c      end do

      end
