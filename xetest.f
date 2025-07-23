        program xetest

        include "common.f"

        natom=2
        do ind=1,3
          x(1,ind)=0.0d0
          x(2,ind)=0.0d0
          box(ind)=100.0d0
        end do

        x(2,1)=2.0d0
        do loop=1,10
          do ind=1,3
            fx(1,ind)=0.0d0
            fx(2,ind)=0.0d0
          end do
          
          call xeforce
          x(2,1)=x(2,1)+0.0001d0
        end do

        end
