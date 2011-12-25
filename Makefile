.SUFFIXES: .c .f90 .f .for

F90C  = gfortran
mycc  = gcc
mycpp = g++

C_COMPILE_FLAGS = -O3 -g -Wall
F_COMPILE_FLAGS = -O3 -g -Wall

CFLAGS = ${C_COMPILE_FLAGS} 
FFLAGS = ${F_COMPILE_FLAGS}

.c.o:
	${mycc} -o $@ -c ${CFLAGS} $<
.f.o:
	${F90C} -o $@ -c ${FFLAGS} $<
.for.o:
	${F90C} -o $@ -c ${FFLAGS} $<
.f90.o:
	${F90C} -o $@ -c ${FFLAGS} $<
.cpp.o:
	${mycpp} -o $@ -c ${CFLAGS} $<


OBJS2D =  5pt2d.o fsls.o rcm.o

OBJS3D =  7pt3d.o fsls.o rcm.o
	
2d : ${OBJS2D}
	${mycpp} -o 5pt ${OBJS2D} -lm -llapack
	
3d : ${OBJS3D}
	${mycpp} -o 7pt ${OBJS3D} -lm -llapack

clean :
	-rm -f *.o
	-rm -f *~
	-rm -f 5pt 7pt

