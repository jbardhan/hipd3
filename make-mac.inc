CC = gcc
CFLAGS = -Wall -fast -g -O3 -ftree-vectorize
INCLUDE = -I/opt/local/include -I../fftsvd
FLAGS = 
LDFLAGS = $(CFLAGS)
LIBS = -L/opt/local/lib -lsuperlu_4.3 -lfftw3 -lfftw3f -llapack -lblas -lm
