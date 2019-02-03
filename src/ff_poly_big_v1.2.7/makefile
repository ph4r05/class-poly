CC = gcc
CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -m64 -pedantic -std=gnu99
LDFLAGS = -static
INCLUDES = -I/usr/local/include
LIBS = -lgmp -lm
INSTALL_ROOT = /usr/local

HEADERS = asm.h ff.h ffmontgomery64.h ff2k.h ffext.h ffpolyfromroots.h ffpolybig.h ffpolysmall.h ffpoly.h cstd.h polyparse.h
OBJECTS = ff.o ff2k.o  ffext.o  ffpolyfromroots.o ffpolysmall.o ffpolybig.o ffpoly_big.o polyparse.o

all: libff_poly_big.a

clean:
	rm -f *.o 
	rm -f libff_poly_big.a

install: all
	cp -v ff_poly.h ff_poly_big.h $(INSTALL_ROOT)/include
	mkdir -p $(INSTALL_ROOT)/include/ff_poly
	cp -v $(HEADERS) $(INSTALL_ROOT)/include/ff_poly
	cp -v libff_poly_big.a $(INSTALL_ROOT)/lib

##### ff_poly_big library (requires linking with zn_poly)

libff_poly_big.a: $(OBJECTS)
	ar -r libff_poly_big.a $(OBJECTS)
	ranlib libff_poly_big.a
	
##### executables

##### FF_POLY C modules

ff.o: ff.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ff2k.o : ff2k.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ffext.o: ffext.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ffpoly_big.o: ffpoly.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ffpolyfromroots.o: ffpolyfromroots.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ffpolybig.o: ffpolybig.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ffpolysmall.o: ffpolysmall.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

polyparse.o: polyparse.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ff_poly.h : $(HEADERS)
	touch $@
