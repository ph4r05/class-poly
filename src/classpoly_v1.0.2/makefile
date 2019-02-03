CC = gcc
CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -m64 -pedantic -std=gnu99
LDFLAGS = -static
INCLUDES = -I/usr/local/include
LIBS = -lff_poly_big -lzn_poly -lgmp -lm

HEADERS = bipoly.h bitmap.h class_inv.h classpoly.h classpoly_crt.h classpoly_inv.h phi_gcd.h phi_eval.h phi_fj_strings.h phi_poly.h polycosts.h velu.h \
		   bigX1.h crt.h cstd.h ecurve.h evec.h findcurve.h mpzutil.h iqclass.h pickprimes.h qform.h table.h tecurve.h tecurvecosts.h torcosts.h
OBJECTS = bipoly.o class_inv.o classpoly.o classpoly_crt.o classpoly_inv.o compute_classpoly.o phi_eval.o phi_gcd.o phi_poly.o polycosts.o velu.o \
		  crt.o ecurve.o ecurve_ladic.o evec.o findcurve.o mpzutil.o pickprimes.o iqclass.o qform.o table.o tecurve.o
PROGRAMS = classpoly

all: $(PROGRAMS)

clean:
	rm -f *.o
	rm -f $(PROGRAMS)

##### executables

classpoly: classpoly_program.o $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $< $(OBJECTS) $(LIBS)

classpoly_program.o: classpoly_program.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

##### C modules
	
bipoly.o: bipoly.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

class_inv.o: class_inv.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

classpoly.o: classpoly.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

classpoly_crt.o: classpoly_crt.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

classpoly_inv.o: classpoly_inv.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

compute_classpoly.o: compute_classpoly.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

phi_eval.o: phi_eval.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

phi_gcd.o: phi_gcd.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

phi_poly.o: phi_poly.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

polycosts.o: polycosts.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

velu.o: velu.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

crt.o: crt.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ecurve.o: ecurve.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ecurve_ladic.o: ecurve_ladic.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

evec.o: evec.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

findcurve.o: findcurve.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

mpzutil.o: mpzutil.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

pickprimes.o: pickprimes.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

iqclass.o: iqclass.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

qform.o: qform.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

tecurve.o: tecurve.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

table.o: table.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
