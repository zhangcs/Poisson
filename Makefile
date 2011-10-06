.SUFFIXES: .c .f90 .f .for

F90C = gfortran
mycc = gcc

C_COMPILE_FLAGS = -O3 -g -Wall
F_COMPILE_FLAGS = -O3 -g -Wall

CFLAGS = ${C_COMPILE_FLAGS} -lm
FFLAGS = ${F_COMPILE_FLAGS}

.c.o:
	${mycc} -o $@ -c ${CFLAGS} $<
.f.o:
	${F90C} -o $@ -c ${FFLAGS} $<
.for.o:
	${F90C} -o $@ -c ${FFLAGS} $<
.f90.o:
	${F90C} -o $@ -c ${FFLAGS} $<

OBJS2D =  5pt2d.o fsls.o

OBJS3D =  7pt3d.o fsls.o
	
2d : ${OBJS2D}
	${mycc} -o 5pt ${OBJS2D} -lm 
	
3d : ${OBJS3D}
	${mycc} -o 7pt ${OBJS3D} -lm

clean :
	-rm -f *.o
	-rm -f *~
	-rm -f 5pt 7pt

