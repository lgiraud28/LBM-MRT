#makefile
FF77 = gfortran
CC = g++
FFLAGS =   -w -O3
CCFLAGS = -std=c++98 -pedantic -O3 -Wall
LIB = -lgfortran -lm -lrt
.SUFFIXES: .f .cpp .o
ALL = sim

all : $(ALL)
dep= LBM.o main.o seconds.o 

clean:
	rm -f *.o *~ $(ALL) *.vtk

LBM.o: LBM.cpp LBM.h
main.o: main.cpp LBM.h
seconds.o: seconds.cpp LBM.h


sim : $(dep)
	$(CC)    $(CCFLAGS) $(dep) $(LIB)  -o $@

.f.o:	
	$(FF77)  ${FFLAGS} $< -c
.cpp.o:
	$(CC)    $(CCFLAGS) $< -c