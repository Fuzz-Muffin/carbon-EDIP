
      subroutine reflect

      include "common.f"

      dimension v1(2), v2(2), v3(2), v4(2)

c +---------------------------------+
c | Decision to perform reflections |
c +---------------------------------+
      if (numcylinder.eq.0) return

c +-----------------------------------------------+
c | Test to see if any atoms are inside cylinders |
c +-----------------------------------------------+
      do i=1,natom
        do loop=1,numcylinder
          dxx=x(i,1) - cylinder(loop,1)
          dyy=x(i,2) - cylinder(loop,2)

          dxx=dxx - box(1)*nint(dxx/box(1))
          dyy=dyy - box(2)*nint(dyy/box(2))
          dist=sqrt(dxx**2 + dyy**2)
          if (dist.lt.cylinder(loop,3)) then
            vflag=1.0

            ! Compute theta relative to centre of cylinder
            theta=atan2(dyy,dxx)

            ! Initialise V1
            v1(1)= vx(i,1)
            v1(2)= vx(i,2)

            ! Rotate velocity components of V1 clockwise by theta to get V2
            v2(1)= cos(theta)*v1(1) + sin(theta)*v1(2)
            v2(2)=-sin(theta)*v1(1) + cos(theta)*v1(2)

            ! Reflect x-component of V2 to get V3
            v3(1)= -v2(1)
            v3(2)=  v2(2)

            ! Rotate V3 anti-clockwise by theta to get V4
            v4(1)= cos(theta)*v3(1) - sin(theta)*v3(2)
            v4(2)= sin(theta)*v3(1) + cos(theta)*v3(2)

            ! Assign elements of V4 to atom i
            vx(i,1)=v4(1)
            vx(i,2)=v4(2)
            !write(6,*) 'Step ',numstep,' Reflecting atom ',i
          end if
        end do
      end do

      end
