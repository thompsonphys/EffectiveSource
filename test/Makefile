PREFIX = $(CONDA_PREFIX)
GCC = $(PREFIX)/bin/gcc
CFLAGS  = -g -Wall -std=gnu99 -I../ -I$(PREFIX)/include
LDFLAGS = -L$(PREFIX)/lib -Wl,-rpath,$(PREFIX)/lib
LIBS = -lm -lgsl -lgslcblas

DEPS    = ../effsource.h

all : kerr-equatorial kerr-circular point-check mmode-check

PYTHON = python3

POINT_OBJS = point-check.o kerr-equatorial-ctx.o kerr-equatorial-ctx-ld.o kerr-equatorial-point-ctx.o kerr-equatorial-coeffs-ctx.o kerr-equatorial-dtcoeffs-ctx.o kerr-equatorial-dttcoeffs-ctx.o

point-check : $(POINT_OBJS)
	$(GCC) $(CFLAGS) -o point-check $(POINT_OBJS) $(LDFLAGS) $(LIBS)

point-check.o : point-check.c ../effsource_equatorial.h ../effsource_equatorial_point.h
	$(GCC) $(CFLAGS) -O2 -c -o point-check.o point-check.c

MMODE_OBJS = mmode-check.o kerr-equatorial-mmode.o kerr-equatorial-ctx.o kerr-equatorial-ctx-ld.o kerr-equatorial-point-ctx.o kerr-equatorial-coeffs-ctx.o kerr-equatorial-dtcoeffs-ctx.o kerr-equatorial-dttcoeffs-ctx.o

mmode-check : $(MMODE_OBJS)
	$(GCC) $(CFLAGS) -o mmode-check $(MMODE_OBJS) $(LDFLAGS) $(LIBS)

mmode-check.o : mmode-check.c ../effsource_equatorial.h ../effsource_equatorial_point.h ../effsource_equatorial_mmode.h
	$(GCC) $(CFLAGS) -O2 -c -o mmode-check.o mmode-check.c

kerr-equatorial-mmode.o : ../kerr-equatorial-mmode.c ../effsource_equatorial_mmode.h ../effsource_equatorial_point.h ../effsource_dd.h
	$(GCC) $(CFLAGS) -O2 -c -o kerr-equatorial-mmode.o ../kerr-equatorial-mmode.c

kerr-equatorial-ctx.o : ../kerr-equatorial-ctx.c ../effsource_equatorial.h
	$(GCC) $(CFLAGS) -O2 -c -o kerr-equatorial-ctx.o ../kerr-equatorial-ctx.c

kerr-equatorial-ctx-ld.c : ../kerr-equatorial-ctx.c ../codegen/long_double_calc.py
	$(PYTHON) ../codegen/long_double_calc.py ../kerr-equatorial-ctx.c kerr-equatorial-ctx-ld.c

kerr-equatorial-ctx-ld.o : kerr-equatorial-ctx-ld.c ../effsource_equatorial.h
	$(GCC) $(CFLAGS) -O2 -c -o kerr-equatorial-ctx-ld.o kerr-equatorial-ctx-ld.c

kerr-equatorial-point-ctx.o : ../kerr-equatorial-point-ctx.c ../effsource_equatorial_point.h ../effsource_equatorial.h ../effsource_dd.h
	$(GCC) $(CFLAGS) -O2 -c -o kerr-equatorial-point-ctx.o ../kerr-equatorial-point-ctx.c

kerr-equatorial-coeffs-ctx.o : ../kerr-equatorial-coeffs-ctx.c
	$(GCC) $(CFLAGS) -O0 -c -o kerr-equatorial-coeffs-ctx.o ../kerr-equatorial-coeffs-ctx.c

kerr-equatorial-dtcoeffs-ctx.o : ../kerr-equatorial-dtcoeffs-ctx.c
	$(GCC) $(CFLAGS) -O0 -c -o kerr-equatorial-dtcoeffs-ctx.o ../kerr-equatorial-dtcoeffs-ctx.c

kerr-equatorial-dttcoeffs-ctx.o : ../kerr-equatorial-dttcoeffs-ctx.c
	$(GCC) $(CFLAGS) -O0 -c -o kerr-equatorial-dttcoeffs-ctx.o ../kerr-equatorial-dttcoeffs-ctx.c

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

.PHONY : clean all
clean :
	-rm -rf kerr-equatorial kerr-circular point-check mmode-check kerr-equatorial-ctx-ld.c
	-rm -rf *.o *.dSYM
