
      subroutine readmasses

      include "common.f"

      character*80 line

c +--------------------------------------------------------------------------------------+
c | Set default value of fracmass to unity (i.e. 98.9% C-12 and 1.1% C-13) for all atoms |
c +--------------------------------------------------------------------------------------+
      do i=1,natom
        fracmass(i)=1.0
      end do

      if (.not.mass) return

c +--------------------+
c | Parse 'masses.txt' |
c +--------------------+
      open(unit=7,file='masses.txt',status='old')
 100  read(7,105,end=110) line
 105  format(a80)

      read(line,*) iatom,atom_mass
      
      ! Atom number of zero denotes setting default mass for all atoms
      if (iatom.eq.0) then
        do i=1,natom
          fracmass(i)=atom_mass/12.01
        end do

      ! Atom number of -1 denotes setting mass of range of atoms
      elseif (iatom.eq.-1) then
        read(line,*) idummy,iatom,jatom,atom_mass
        do i=iatom,jatom
          fracmass(i)=atom_mass/12.01
        end do

      ! Otherwise, setting mass of individual atom
      elseif ((iatom.gt.0).and.(iatom.le.natom)) then
        fracmass(iatom)=atom_mass/12.01

      ! Error if iatom out-of-range
      else
        stop 'atom out of range in readmasses.f'
      end if
 
      goto 100
 110  continue
      close(unit=7)

      end

