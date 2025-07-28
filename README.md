# carbon-EDIP
### Nigel A. Marks

The carbon environment dependent interaction potential (EDIP), often called the EDIP potential (yes, like ATM machine...).

## Relevant publications
These are the publications that describe the carbo-EDIP potential. If you are using this potential, please consider citing these works.

- N.A. Marks, "Generalizing the environment-dependent interaction potential for carbon", 2000, Phys. Rev. B 63, 035401, DOI [10.1103/PhysRevB.63.035401](https://doi.org/10.1103/PhysRevB.63.035401)

- N.A. Marks, "Modelling diamond-like carbon with the environment-dependent interaction potential", 2002, J. Phys.: Condens. Matter 14 2901, DOI: [10.1088/0953-8984/14/11/308](https://doi.org/10.1088/0953-8984/14/11/308)


### Some details
Size of arrays:
```
      parameter (NMAX=200000,NNN=50,NP=0.5*NNN*NMAX)

      common /ZDERV/  dzdx(NNN,3), dzdxx(NNN,NNN,3)
      common /ZFORC/  finiteforce(NNN), finiteforce2(NNN,NNN)
      common /DRXYZ/  dr(NMAX,NNN),dx(NMAX,NNN,3),dxdr(NMAX,NNN,3)
      common /INEAR/  near(NMAX,NNN),num(NMAX),kron(NNN,NNN)
      common /ZCYCL/  f(NMAX,NNN), df(NMAX,NNN), dzz(NMAX,NNN)

      common /PAIRS/  ipair(NP),jpair(NP),inum(NP),jnum(NP),npair

NMAX*NNN:

/DRXYZ/  5*NMAX*NNN   (doubles)
/INEAR/  1*NMAX*NNN   (integer)
/ZCYCL/  3*NMAX*NNN   (doubles)
/PAIRS/  2*NMAX*NNN   (doubles)
```

assuming 8 bytes per value, we have `88*NMAX*NNN` bytes

```
For NMAX=200000 and NNN=50,    NMAX*NNN=10 million
For NMAX=200000 and NNN=50, 88*NMAX*NNN=880 million
For NMAX=400000 and NNN=50, 88*NMAX*NNN=1.76 billion (~1.8 GB)
```

There are several other large arrays, compilation fails with
`NMAX=400000`, but ok with `NMAX=350000`.  Looks like we have a 
maximum memory allocation of 2GB.  Pretty small ...

