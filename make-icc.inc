CC = icc
CFLAGS = -O3 -g -xW -unroll -ip
INCLUDE = -I/usr/local/include -I/usr/local/include/SuperLU -I/programs/i386/lib/intel/ipp41/ia32_itanium/include -I../FFTSVDcp
FLAGS = -DINTEL_C #-DMVTIME
LDFLAGS = $(CFLAGS) -static-libcxa -i-static
LIBDIR = -L/usr/local/lib -L/programs/i386/lib/intel/ipp/5.0/ia32/lib -L/programs/i386/lib/intel/mkl/8.0.1/lib/32
LIBS = -lsuperlu /usr/lib/liblapack.a /usr/lib/libblas.a -lg2c -lfftw3 -lfftw3f -lippsemerged -lippsmerged -lippcore

