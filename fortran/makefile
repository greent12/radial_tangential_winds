FC   =  gfortran
FLAGS=  -O2
OBJECTS = read_testdata.o recenter.o

recenter.x: ${OBJECTS}
	${FC} -o recenter.x ${OBJECTS}

recenter.o: recenter.f90
	${FC} -c recenter.f90

read_testdata.o: read_testdata.f90
	${FC} -c read_testdata.f90

clean:
	rm *.o *.mod	

run:
	./recenter.x
