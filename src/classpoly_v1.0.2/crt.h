#ifndef _CRT_INCLUDE_
#define _CRT_INCLUDE_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CRT_DIR					crt_dir()		
#define CRT_DIR_NAME				"temp"
extern char _crt_dir_str[1024];
static inline char *crt_dir (void)
{
	char *s;
	if ( _crt_dir_str[0] ) return _crt_dir_str;
	s = getenv("HOME");
	if ( s ) sprintf (_crt_dir_str, "%s/%s", s, CRT_DIR_NAME);
	else strcpy (_crt_dir_str, CRT_DIR_NAME);
	return _crt_dir_str;
}

#define CRT_MAX_LEVELS		20
#define CRT_BASE_N			8
#define CRT_MAX_MEM			2000000000L		// keep this reasonably small if you plan to run multiple jobs on the same machine
#define CRT_BATCH_SIZE			4096
#define CRT_MAX_JOBS			1024

#define CRT_MAX_COUNT			4096			// # of crt coefficients stored per record --- lower value uses less memory, higher value uses fewer files
#define CRT_BUF_SIZE			((CRT_MAX_COUNT+1)*sizeof(unsigned long))

struct crt_tree {
	unsigned long *m;
	unsigned long *a;
	int *cnts;							// cnts at level 0, contains lens[0] entries
	int lens[CRT_MAX_LEVELS];			// The top level has length 1 and each level below doubles the length of the level above, but it is convenient to have these around
	mpz_t *M[CRT_MAX_LEVELS];			// M[i] contains lens[i] entries
	mpz_t *A[CRT_MAX_LEVELS];			// A[i] contains lens[i] entries (we don't need to allocate space for all levels, we could overwrite, but it is simpler and faster to use the extra space (< 25% increase overall)).
	mpz_t *X[CRT_MAX_LEVELS];			// workspace used by uneval - only allocated if uneval is called
	mpz_t *B;							// always has length n regardless of the number of levels
	mpz_t M0, M1;						// work variables
	mpz_t MM;						// product of all the m[i]
	int levels;							// height of the tree (top level is implicit)
	int n;							// number of moduli
};
typedef struct crt_tree crt_tree_t[1];

// ecrt context is used for explicit crt mod P
struct ecrt_context {
	unsigned long *m;					// n moduli m_1,..., m_n, with M=prod m_i
	unsigned long *a;					// n values a_i = 1/M_i mod m_i where M_i = M/m_i = prod_{j!=i} m_j
//	mpz_t *d;							// n values d_i = a_i M_i mod P where a_i = 1/M_i mod m_i and M_i = prod_{j\ne i} m_i
	mpz_t P;							// output modulus
	mpz_t MP;						// prod m_i mod P
	mpz_t X,Y, Z;						// work variables
	int n;							// number of moduli
	int k;							// number of coefficients
	int j, j0;							// next coefficient to enumerate (after finalizing)
	int delta;							// ceil(lg(n)+2)
	int Climbs;						// number of limbs in C
	int Cbytes;						// Climbs * sizeof(mp_limb_t) for convenience
	int slimbs;						// number of limbs in s (currently always 2)
	int sbytes;						// slimbs * sizeof(mp_limb_t)
	mp_limb_t *Cdata;
	mp_limb_t *sdata;
	char prefix[64];
	int jobs, jobid;						// if jobs=0 then no job processing is performed (i.e. all ecrt precomp, update and postcomp happens in a single process)
	FILE *fp;							// file pointer for job data (if used)
	int batching;						// batching is only supported with job processing
	long fpbase;						// file offset of data stored as blocks of CRT_BATCH_SIZE*Cbytes +CRT_BATCH_SIZE*sbytes
	FILE *infp[CRT_MAX_JOBS];			// used to merge data from multiple jobs
	mp_limb_t *inCdata;
	mp_limb_t *insdata;
};
typedef struct ecrt_context ecrt_context_t[1];



struct crt_file {
	FILE *fp;
	unsigned long id1;
	unsigned long id2;
	int offset;
	int count;
};
typedef struct crt_file crt_file_t[1];

int crt_file_create (crt_file_t fp, unsigned long id1, unsigned long id2, int offset, int count);
void crt_file_write (crt_file_t fp, unsigned long m, unsigned long c[], int n);
void crt_file_close (crt_file_t fp);
int crt_file_open (crt_file_t fp, unsigned long id1, unsigned long id2, int offset);
int crt_file_read (crt_file_t fp, unsigned long *m, unsigned long c[], int n);
void crt_file_delete (crt_file_t fp);

int crt_create_files (crt_file_t **fpp, int id1, int id2, int n);									// creates ceil(n/CRT_MAX_COUNT) crt files sufficient to handle n coefficients
void crt_write_files (crt_file_t *fp, unsigned long m, unsigned long c[], int n);					// writes n coefficient values c[] mod m to crt files	
void crt_close_files (crt_file_t *fp, int n);
int crt_process_files (long *totbits, long *maxbits, int id1, int id2, int n,  unsigned long m[], int k,	// totbits and maxbits are optional, if specified set to total # bits and bits in largest coefficient, respectively
				 int (*callback)(mpz_t C, int i, void *ctx), void *ctx);						// k moduli, n coefficients, calls callback once for each coefficient.  Stops and returns 0 if callback returns false.
																				// Destroys crt files upon successful completion.

int crt_mpz_ui (mpz_t C, mpz_t M, mpz_t C1, mpz_t M1, unsigned long c2, unsigned long m2);	// overlap is not permitted

// assumes C represents an integer of absolute value < M/4
static inline void crt_mpz_mod_to_Z (mpz_t C, mpz_t M)
	{ if ( mpz_sizeinbase(C,2) >= mpz_sizeinbase(M,2)-2 ) { mpz_sub(C,M,C); mpz_neg(C,C); } }
	
// assumes C represents a rational number num/den with |num| < sqrt(M/2) and 0<den<sqrt(M/2) and gcd(num,den)=1.  Returns 0 if this is not true.
int crt_mpz_mod_to_Q (mpz_t num,  mpz_t den, mpz_t C, mpz_t M, mpz_t w[5]);

void crt_tree_init (crt_tree_t t, unsigned long m[], int n);
void crt_tree_clear (crt_tree_t t);
void crt_tree_eval (mpz_t C, unsigned long c[], crt_tree_t t);
void crt_tree_uneval (unsigned long c[], mpz_t C, crt_tree_t t);
static inline void crt_tree_M (mpz_t M, crt_tree_t t) { mpz_set(M,t->MM); }
static inline void crt_tree_eval_to_Z (mpz_t C, unsigned long c[], crt_tree_t t)
	{ crt_tree_eval (C, c, t); crt_mpz_mod_to_Z (C, t->MM); }

void ecrt_init (ecrt_context_t ecrt, unsigned long m[], int k, int n, mpz_t P,					// k moduli, n coefficients, jobs may be zero, jobid ranges from 1 to jobs
		       int jobs, int jobid, char *prefix);											// if jobs is > 0 then jobid and prefix are required, otherwise they are ignored
void ecrt_clear (ecrt_context_t ecrt);
void ecrt_update (ecrt_context_t ecrt, int i, unsigned long c[], int n);							// ith moduli, k coefficient values
void ecrt_finalize (ecrt_context_t ecrt);													// returns k, should only be used for single jobs that do not require batching, otherwise dump and merge
void ecrt_reset (ecrt_context_t ecrt);													// reset for a new round of computations using the same moduli (only supported when batching is not used)
void ecrt_dump (ecrt_context_t ecrt);													// dumps intermediate (nonfinalized) coefficient sums to binary file specified by prefix, jobs, and jobid in ecrt_init
void ecrt_merge (ecrt_context_t ecrt, int jobs, char *prefix);									// merges dumped coefficient sums and sets ecrt for coefficient enumeration
static inline int ecrt_coeff_count (ecrt_context_t ecrt) { return ecrt->k; }
int ecrt_next_coeff (mpz_t C,ecrt_context_t ecrt);											// returns jth coefficient in C, must be called after finalizing

/*
xcrt code is untested and unused at present (too slow)

// xcrt context is used for explicit crt mod P mod m_i
struct xcrt_context {
	unsigned long *m;					// n moduli m_1,...,m_n with M=prod m_i
	double *mi;						// n floating-point approximations to 1/m_i
	unsigned long *a;					// n values a_i = 1/M_i mod m_i  where M_i = M/m_i = prod_{j!=i} m_j
	unsigned long *d;					// n*n values d_ij = (M_i mod P) mod m_j				note this differs from d in ecrt, we don't include a_i in d_ij
	unsigned long *e;					// n values e_i = M mod P mod m_i
	mpz_t P;							// output modulus
	int n;							// number of moduli
	int k;							// number of coefficients
	int delta;							// ceil(lg(n)+2)
	unsigned long *c;					// n*k values c_ij where c_ij is the jth coefficient C_j mod p_i
	unsigned long *t;					// n*k values t_ij = a_i*c_ij 
	unsigned long *s;					// k values s_j = sum_i floor (2^delta*t_ij / p_i)  (accumulated in 64 bits -- assumes delta+log_2(n)+log_2 max_pi < 64 bits!
};
typedef struct xcrt_context xcrt_context_t[1];

void xcrt_init (xcrt_context_t xcrt, unsigned long m[], int n, int k, mpz_t P);
void xcrt_clear (xcrt_context_t ctx);
void xcrt_reduce (xcrt_context_t ctx);


static inline void xcrt_update (xcrt_context_t xcrt, int i, unsigned long c[], int k)
{
	register int j;
	
	for ( j = 0 ; j < k ; j++ ) xcrt->c[i*k+j] = c[j];
}

static inline void xcrt_fetch (xcrt_context_t xcrt, int i, unsigned long c[], int k)
{
	register int j;
	
	for ( j = 0 ; j < k ; j++ ) c[j] = xcrt->c[i*k+j];
}
*/

#ifdef __cplusplus
}
#endif

#endif
