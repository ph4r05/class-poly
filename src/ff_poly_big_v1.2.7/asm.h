#ifndef _ASM_INCLUDE_
#define _ASM_INCLUDE_


/*
    Copyright 2011-2014 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#ifdef __cplusplus
extern "C" {
#endif

#define _asm_div_q_q(q,r,x,y)				asm ("divq %4" :"=a"(q) ,"=d"(r) : "0"(x), "1"(r), "rm"(y))
#define _asm_mult_1_1(z1,z0,x0,y0)		asm ("mulq %3" :"=a"(z0) ,"=d"(z1) : "a"(x0), "rm"(y0))
#define _asm_mult_2_2_1(z1,z0,x1,x0,y0)	asm ("mulq %3" :"=a"(z0) ,"=d"(z1) : "a"(x0), "rm"(y0));(z1)+=(y0)*(x1)
#define _asm_addto_2_2(z1,z0,x1,x0)		asm ("addq %3,%0;adcq %5,%1":"=r"(z0),"=r"(z1): "0"(z0), "rim"(x0),  "1"(z1), "rim"(x1):"cc")
#define _asm_addto_2_1(z1,z0,x0)			asm ("addq %3,%0;adcq $0,%1":"=r"(z0),"=r"(z1): "0"(z0), "rim"(x0),  "1"(z1):"cc")
#define _asm_addto_3_3(z2,z1,z0,x2,x1,x0)	asm ("addq %4,%0;adcq %6,%1;adcq %8,%2":"=r"(z0),"=r"(z1),"=r"(z2): "0"(z0), "rim"(x0),  "1"(z1), "rim"(x1), "2"(z2), "rim"(x2) :"cc")
#define _asm_addto_3_2(z2,z1,z0,x1,x0)		asm ("addq %4,%0;adcq %6,%1;adcq 0,%2":"=r"(z0),"=r"(z1),"=r"(z2): "0"(z0), "rim"(x0),  "1"(z1), "rim"(x1), "2"(z2) :"cc")
#define _asm_subfrom_2_2(z1,z0,x1,x0)		asm ("subq %3,%0;sbbq %5,%1":"=r"(z0),"=r"(z1): "0"(z0), "rim"(x0),  "1"(z1), "rim"(x1):"cc")
// increment needs to propogate the carry - performance not critical here anyway
#define _asm_inc_2(z1,z0)				asm ("addq $1,%0;adcq $0,%1":"=r"(z0),"=r"(z1): "0"(z0), "1"(z1):"cc")
//#define _asm_shiftl_2(z1,z0)				asm ("shlq %0,1;rclq %1,1":"=r"(z0),"=r"(z1): "0"(z0), "1"(z1):"cc")
//#define _asm_shiftr_2(z1,z0)				asm ("shrq %1,1;rcrq %0,1":"=r"(z0),"=r"(z1): "0"(z0), "1"(z1):"cc")


#define _asm_mult_3_2_1(z2,z1,z0,x1,x0,y0)	{ register unsigned long __u; \
									   _asm_mult_1_1 (__u,z0,x0,y0); \
									   _asm_mult_1_1 (z2,z1,x1,y0); \
									   _asm_addto_2_1 (z2,z1,__u); }

// This function assumes that x[1] and y[1] < 2^31
static inline void _asm_mult_3_2_2 (unsigned long z[3], unsigned long x[2], unsigned long y[2])
{
	register unsigned long U, V, R0,R1;
	
	R1 = 0;
	_asm_mult_1_1(R0,z[0],x[0],y[0]);
	_asm_mult_1_1(U,V,x[0],y[1]);
	_asm_addto_2_2(R1,R0,U,V);
	_asm_mult_1_1(U,V,x[1],y[0]);
	_asm_addto_2_2(R1,R0,U,V);
	z[1] = R0;
	z[2] = x[1]*y[1]+R1;
}

// This function assumes that x[1] and y[1] < 2^31
static inline void _asm_mult_3_2_2r (unsigned long z2, unsigned long z1, unsigned long z0, unsigned long x[2], unsigned long y[2])
{
	register unsigned long U, V, R0,R1;
	
	R1 = 0;
	_asm_mult_1_1(R0,z0,x[0],y[0]);
	_asm_mult_1_1(U,V,x[0],y[1]);
	_asm_addto_2_2(R1,R0,U,V);
	_asm_mult_1_1(U,V,x[1],y[0]);
	_asm_addto_2_2(R1,R0,U,V);
	z1 = R0;
	z2 = x[1]*y[1]+R1;
}


// This function assumes that x[1] < 2^31.  For no obvious reason, this is slower than multiplying?!
static inline void _asm_square_3_2 (unsigned long z[3], unsigned long x[2])
{
	register unsigned long U, V, R0,R1;
	
	R1 = 0;
	_asm_mult_1_1(R0,z[0],x[0],x[0]);
	_asm_mult_1_1(U,V,x[0],x[1]);
	_asm_addto_2_2(R1,R0,U,V);
	_asm_addto_2_2(R1,R0,U,V);
	z[1] = R0;
	z[2] = x[1]*x[1]+R1;
}

#ifdef __cplusplus
}
#endif

#endif
