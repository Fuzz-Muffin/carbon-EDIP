
OBJ = checkcell.o\
checkinput.o\
checkrij.o\
conductivity.o\
constants.o\
coordination.o\
cutoff.o\
defaults.o\
density.o\
dihedral.o\
dihedral2.o\
dihedral4.o\
distance.o\
distnab.o\
distribution.o\
edgetherm.o\
energy.o\
fastatom.o\
force.o\
init.o\
linear.o\
main.o\
makespecial.o\
neighbour.o\
pair.o\
pair_zbl.o\
parse.o\
printav.o\
printfinish.o\
printinput.o\
printrings.o\
printstatus.o\
printtime.o\
properties.o\
readcoords.o\
readinput.o\
readmasses.o\
reflect.o\
remspace.o\
repulsion.o\
resetmsd.o\
runspecial.o\
spherecoord.o\
sphereforce.o\
stress.o\
therm.o\
triple.o\
volume.o\
varystep.o\
verlet.o\
writeavas.o\
writecoords.o\
writegr.o\
writeovito.o\
writestress.o\
writetheta.o\
writexbs.o\
writexyz.o\
xecheckcell.o\
xeforce.o\
xeneighbour.o

GOAL = edip
INCL = common.f zblpar.f
FFLAGS = -O2 -fopenmp
#FFLAGS = -O2 -ip -fpe0 -openmp
#FFLAGS = -O3 -h omp
#FC = ftn
FC = gfortran

$(GOAL): $(OBJ) ; \
	$(FC) $(OBJ) $(FFLAGS) -o $(GOAL)

$(OBJ): $(INCL)

.f.o:; $(FC) $(FFLAGS) -c $<

clean:; rm *.o


