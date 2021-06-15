CC = icc
CFLAGS = -O3 -axiMKW -tpp7 -unroll -ipo -ipo_obj -openmp
FLAGS = -DMVTIME -DOMP
LDFLAGS = $(CFLAGS)
LIBS =
