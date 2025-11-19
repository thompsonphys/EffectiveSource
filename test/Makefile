PREFIX = $(CONDA_PREFIX)
GCC = $(PREFIX)/bin/clang
CFLAGS  = -g -Wall -std=gnu99 -I../ -I$(PREFIX)/include
LDFLAGS = -L$(PREFIX)/lib -Wl,-rpath,$(PREFIX)/lib
LIBS = -lm -lgsl -lgslcblas

DEPS    = ../effsource.h

all : kerr-equatorial kerr-circular

kerr-equatorial : kerr-equatorial.o kerr-equatorial-coeffs.o kerr-equatorial-dtcoeffs.o kerr-equatorial-dttcoeffs.o effsource.o decompose.o $(DEPS)
	$(GCC) $(CFLAGS) -o kerr-equatorial kerr-equatorial-coeffs.o kerr-equatorial-dtcoeffs.o kerr-equatorial-dttcoeffs.o kerr-equatorial.o effsource.o decompose.o $(LDFLAGS) $(LIBS)

kerr-circular : kerr-circular.o effsource.o decompose.o $(DEPS)
	$(GCC) $(CFLAGS) -o kerr-circular kerr-circular.o effsource.o decompose.o $(LDFLAGS) $(LIBS)

effsource.o : effsource.c
	$(GCC) $(CFLAGS) -O3 -c -o effsource.o effsource.c 

decompose.o : decompose.c
	$(GCC) $(CFLAGS) -O3 -c -o decompose.o decompose.c

kerr-equatorial.o : ../kerr-equatorial.c
	$(GCC) $(CFLAGS) -O3 -c -o kerr-equatorial.o ../kerr-equatorial.c

kerr-equatorial-coeffs.o : ../kerr-equatorial-coeffs.c
	$(GCC) $(CFLAGS) -O0 -c -o kerr-equatorial-coeffs.o ../kerr-equatorial-coeffs.c

kerr-equatorial-dtcoeffs.o : ../kerr-equatorial-dtcoeffs.c
	$(GCC) $(CFLAGS) -O0 -c -o kerr-equatorial-dtcoeffs.o ../kerr-equatorial-dtcoeffs.c

kerr-equatorial-dttcoeffs.o : ../kerr-equatorial-dttcoeffs.c
	$(GCC) $(CFLAGS) -O0 -c -o kerr-equatorial-dttcoeffs.o ../kerr-equatorial-dttcoeffs.c

kerr-circular.o : ../kerr-circular.c
	$(GCC) $(CFLAGS) -O3 -c -o kerr-circular.o ../kerr-circular.c

.PHONY : clean
clean :
	-rm -rf kerr-equatorial kerr-circular
	-rm -rf *.o *.dSYM
