CC = gcc
#CFLAGS =  -g -O3 -mtune=pentium4 -mmmx -msse -msse2 -mfpmath=sse,387
CFLAGS = -g
INCLUDE = -I$(SUPERLU_INC) -I/usr/people/jbardhan/include -I../fftsvd
FLAGS = -D_MALLOC #-DMVTIME
LDFLAGS = $(CFLAGS)
LIBDIR = -L/usr/people/jbardhan/lib
SUPERLU_INC = /usr/people/jbardhan/src/SuperLU/SuperLU_4.3/SRC
SUPERLU_LIB = /usr/people/jbardhan/src/SuperLU/SuperLU_4.3/lib
LIBS = -L$(SUPERLU_LIB) -lsuperlu_4.3 -llapack -lrefblas -lfftw3 -lfftw3f -lm -lgfortran
