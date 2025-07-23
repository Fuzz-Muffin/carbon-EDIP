      parameter (E_CHAR = 1.602D-19)
      parameter (EPSILON = 8.85D-12)


C ZBL parameters for Xe-C interaction
      z_c1 = 0.02817d0
      z_c2 = 0.28022d0
      z_c3 = 0.50986d0
      z_c4 = 0.18175d0

      z_1 = 6.0d0
      z_2 = 54.0d0
      a_u = 0.46850d0 / (z_1**0.23d0 + z_2**0.23d0)
      z_d1 = 0.20162d0/a_u  
      z_d2 = 0.40290d0/a_u
      z_d3 = 0.94229d0/a_u
      z_d4 = 3.19980d0/a_u

C Lennard-Jones parameters for Xe-C interaction
C From Simonyan et al, J. Chem. Phys. Vol 114, p. 4180 (2001)
      sigmaxe=3.332d0
      epsxe=132.31 * 8.617d-5

C Parameters for switching function
C From gnuplot file: zbl.gp
      zskin=1.0/8.0d0
      zcut=2.7d0
      zshift=0.07d0
