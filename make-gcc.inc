CC = gcc
#CFLAGS =  -g -O3 -mtune=pentium4 -mmmx -msse -msse2 -mfpmath=sse,387
CFLAGS = -g
INCLUDE = -I$(SUPERLU_INC) -I../FFTSVDcp
FLAGS = -DMVTIME
LDFLAGS = $(CFLAGS)
LIBDIR =
SUPERLU_INC = /usr/local/include/SuperLU
SUPERLU_LIB = /usr/local/lib
LIBS = -L$(SUPERLU_LIB) -lsuperlu -llapack -lblas -lfftw3 -lfftw3f -lg2c
