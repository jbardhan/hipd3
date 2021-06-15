CC = gcc
CFLAGS = -g -O3 -D_MALLOC
INCLUDE = -I$(HOME)/include/ -I../fftsvd
FLAGS = -DMVTIME
LDFLAGS = $(CFLAGS)
LIBDIR = -L$(HOME)/lib/
LIBS =  -lsuperlu_4.3 -llapack -lblas -lfftw3 -lfftw3f
