In order to build classpoly, you need to install the following libraries:

1) gmp.lib, available at http://gmplib.org
2) zn_poly.lib, available at http://web.maths.unsw.edu.au/~davidharvey/code/zn_poly/
3) ff_poly_big.lib, available at http://math.mit.edu/~drew/

Note that ff_poly_big and classpoly require a 64-bit environment.

The classpoly makefile expects the library files to be installed to
/usr/local/lib and header files in /usr/local/include.

Before running the classpoly program, you will also need to download
some modular polynomials, which you can find at

    http://math.mit.edu/~drew/SmallModPolys.html
    
You will need to download at least phi_j.tar.  If you are only interested
in computing Hilbert class polynomials (as described in [1]), this is all you
need, but if you want to compute class polynomials for other modular functions
(as described in [2]), you will need to download the corresponding modular
polynomials.  The classpoly program will look for these polynomials
in $(HOME)/phi_files.  If you want to change this (e.g. to a directory in
/usr/local), modify the function phi_dir() in phi_poly.h.

The command line syntax for classpoly is

  classpoly D [inv P filename verbosity]

D: an imaginary quadratic discriminant
   (if you specify a positive D, it will automatically be negated).
    
inv: an integer from the table below that identifies a class invariant

P: an integer modulo which the class polynomial is to be computed.
   (if 0 or unspecified, it will be computed over Z)
   
filename: the name of the output file.  The default is H_%d.txt,
          where %d is replaced by -D.
	  
verbosity: specify of 1 or 2 to see more details of the computation,
	   specify -1 to supress all details
	   
Table of class invariants supported by classpoly:

0 = j
1 = f (Weber function, see section 3 of [2])
2 = f^2
5 = gamma_2 (cube-root of j)
6 = w_{2,3} (double eta-quotient, see section 3 of [2])
9 = w_{3,3}
10 = w_{2,5}
11 = t (related to Ramanujan function, see section 4.4 in [2])
12 = t^2
14 = w_{2,7}
15 = w_{3,5}
21 = w_{3,7}
23 = w_{2,3}^2
24 = w_{2,5}^2
26 = w_{2,13}
27 = w_{2,7}^2
28 = w_{3,3}^2
100+N = A_N (for N=3,5,7,11,13,17,19,23,29,31,41,47,59,71) (Atkin functions)
400+N = w_N^s (for N=3,5,7,13, s = 24/gcd(24,N-1)) (single-eta quotients)
500+p1*p2 = w_{p1,p2}^s, s = 24/gcd(24,(p1-1)*(p2-1)) (double-eta quotients)
   (p1,p2)=(2,3),(2,5),(2,7),(2,13),(3,3),(3,5),(3,7),(3,13),(5,7)


If you use this software in your research, please cite one or both of the
following papers, which describe the algorithms and their implementation
in detail.
in detail.

[1] Andrew V. Sutherland, "Computing Hilbert class polynomials with
the Chinese Remainder Theorem", Math. Comp. 80 (2011), 501-538.

[2] Andreas Enge and Andrew V. Sutherland, "Class invariants by the
CRT method", ANTS IX, LNCS 6197 (2010), 142-156.
