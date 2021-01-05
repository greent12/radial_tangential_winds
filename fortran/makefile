#Makefile for vortex recentering routines

#Set these for your fortran90 compiler (F90) and openmp flag (OMP)
F90   =  ifort
OMP   =  fopenmp
FLAGS =  -O2
ECCFLAGS =  -I/Users/tylergreen/eccodes-2.19.1-Source/build/include -L/Users/tylergreen/eccodes-2.19.1-Source/build/lib -leccodes_f90 
EXE   = recenter.x
 
OBJECTS = read_testdata.o ideal_angle.o constants.o kinds.o recenter.o \
          calc_distance_from_point.o uv_to_rt.o maskmod.o

${EXE} : ${OBJECTS}
	${F90} -o recenter.x ${OBJECTS} ${FLAGS} -${OMP}

recenter.o: recenter.f90 read_testdata.o kinds.o constants.o ideal_angle.o calc_distance_from_point.o uv_to_rt.o maskmod.o
	${F90} -c recenter.f90 -${OMP}

read_testdata.o: read_testdata.f90
	${F90} -c read_testdata.f90

ideal_angle.o: ideal_angle.f90 read_testdata.o kinds.o constants.o
	${F90} -c ideal_angle.f90 ${FLAGS} -${OMP}

constants.o: constants.f90
	${F90} -c constants.f90

kinds.o: kinds.f90
	${F90} -c kinds.f90

calc_distance_from_point.o: calc_distance_from_point.f90 kinds.o constants.o
	${F90} -c calc_distance_from_point.f90 

uv_to_rt.o: uv_to_rt.f90 kinds.o constants.o
	${F90} -c uv_to_rt.f90

maskmod.o: maskmod.f90 kinds.o constants.o
	${F90} -c maskmod.f90

gribmod.o: grib_mod.f90 
	${F90} -c grib_mod.f90 ${ECCFLAGS}
clean:
	rm ${EXE} ${OBJECTS} *.mod 

run:
	./${EXE}

