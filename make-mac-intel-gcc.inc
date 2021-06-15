CC = gcc
CFLAGS = -Wall -g -O3
INCLUDE = -I$(SUPERLU_INC) -I../FFTSVDcp -I../fftw-3.0.1/api
FLAGS = -DMVTIME -DPWD=\"`pwd`\"
LDFLAGS = $(CFLAGS)
LIBDIR =
SUPERLU_INC = ../SuperLU_3.0/SRC
SUPERLU_LIB = ../SuperLU_3.0
LIBS = -L/usr/local/lib -L$(SUPERLU_LIB) -lsuperlu_3.0  -lfftw3 -lfftw3f -framework vecLib
