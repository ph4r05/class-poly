#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <gmp.h>
#include "mpzutil.h"
#include "crt.h"
#include "cstd.h"

/*
    Copyright 2010-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

char _crt_dir_str[4096];

#define CRT_FILE_FORMAT_STRING		"%s/CRT_%lu_%lu_%d.crt"

static unsigned long items[CRT_MAX_COUNT+1];		// shared, avoid allocating this on the stack


int crt_mpz_mod_to_Q (mpz_t num,  mpz_t den, mpz_t C, mpz_t M, mpz_t w[5])
{
	__mpz_struct *a0, *a1, *b0, *b1, *B, *q, *t;

	a0 = w[0];  b0 = w[1];  q = w[2];  B = w[3];  t = w[4];
	a1 = num;  b1 = den;
	mpz_div_2exp (B, M, 1);
	mpz_sqrt (B, B);
	mpz_set (a0, M);  mpz_set (a1, C);
	mpz_set_ui (b0, 0); mpz_set_ui (b1, 1);
	for (;;) {
		if ( mpz_cmpabs(b1,B) > 0 ) return 0;
		if ( mpz_cmpabs(a1,B) < 0 ) {
			mpz_gcd (t, a1, b1);
			if ( mpz_cmp_ui(t,1) == 0 ) {
				if ( mpz_sgn(b1) < 0 ) { mpz_neg(a1,a1); mpz_neg(b1,b1); }
				return 1;
			}
		}
		mpz_fdiv_q (q, a0, a1);
		mpz_mul (t, q, a1);  mpz_sub (t, a0, t); mpz_set (a0, a1);  mpz_set (a1, t); 	// a0,a1 = a1,a0-q*a1
		mpz_mul (t, q, b1);  mpz_sub (t, b0, t); mpz_set (b0, b1);  mpz_set (b1, t); 	// b0,b1 = b1,b0-q*b1
	}
}


int crt_file_create (crt_file_t fp, unsigned long id1, unsigned long id2, int offset, int count)
{
	char buf[256];
	
	if ( count < 1 || count > CRT_MAX_COUNT ) { printf ("Invalid count=%d, CRT_MAX_COUNT=%d\n", count, CRT_MAX_COUNT); return 0; }
	sprintf (buf, CRT_FILE_FORMAT_STRING, CRT_DIR, id1, id2, offset);
	fp->fp = fopen (buf, "wb");
	if ( ! fp->fp ) { printf ("Error creating file %s\n", buf); return 0; }
	fp->count = count;
	fp->id1 = id1;
	fp->id2 = id2;
	fp->offset = offset;
	fp->count = count;
	memset(items,0,sizeof(items));
	items[0] = id1;
	items[1] = id2;
	items[2] = offset;
	items[3] = count;
	if ( fwrite (items, sizeof(unsigned long), CRT_MAX_COUNT+1, fp->fp) != CRT_MAX_COUNT+1 ) { printf ("File i/o error in crt_file_write\n"); abort(); }
	return 1;
}

void crt_file_write (crt_file_t fp, unsigned long m, unsigned long c[], int n)
{
	if ( n != fp->count ) { printf ("Count mismatch in crt_file_write, %d != %d\n", n, fp->count); abort(); }
	items[0] = m;
	memcpy (items+1,c,n*sizeof(unsigned long));
	if ( fwrite(items, sizeof(unsigned long), n+1, fp->fp) != n+1 ) { printf ("File i/o error in crt_file_write\n"); abort(); }
	return;
}

void crt_file_close (crt_file_t fp)
{
	fclose (fp->fp);
	fp->fp = 0;
}

int crt_file_open (crt_file_t fp, unsigned long id1, unsigned long id2, int offset)
{
	char buf[256];
	
	sprintf (buf, CRT_FILE_FORMAT_STRING, CRT_DIR, id1, id2, offset);
	fp->fp = fopen (buf, "rb");
	if ( ! fp->fp ) { printf ("Error opening file %s\n", buf); return 0; }
	if ( fread(items,sizeof(unsigned long),CRT_MAX_COUNT+1,fp->fp) != CRT_MAX_COUNT+1 ) { printf ("File i/o error in crt_file_read\n"); return 0; }
	if ( items[0] != id1 || items[1] != id2 || items[2] != offset ) { printf ("Invalid data in CRT file %s\n", buf); return 0; }
	fp->id1 = id1; fp->id2 = id2; fp->offset = offset;
	fp->count = items[3];
	if ( fp->count > CRT_MAX_COUNT ) { printf ("Invalid CRT file %s, count %d > CRT_MAX_COUNT\n", buf, fp->count); return 0; }
	return fp->count;
}

int crt_file_read (crt_file_t fp, unsigned long *m, unsigned long c[], int n)
{
	int cnt;

	if ( n != fp->count ) { printf ("Count mismatch in crt_file_read, %d != %d\n", n, fp->count); exit(0); }
	cnt = fread(items,sizeof(unsigned long),n+1,fp->fp);
	if ( cnt && cnt != n+1 ) {printf ("Partial i/o error in crt_file_write\n"); exit(0); }
	if ( ! cnt ) return 0;
	*m = items[0];
	memcpy (c,items+1,n*sizeof(unsigned long));
	return n;
}

void crt_file_delete (crt_file_t fp)
{
	char buf[256];
	int sts;
	
	if ( fp->fp ) crt_file_close(fp);
	sprintf (buf, CRT_FILE_FORMAT_STRING, CRT_DIR, fp->id1, fp->id2, fp->offset);
	if ( (sts=remove(buf)) ) printf ("error %d occured while deleting file %s\n", sts, buf);
}

// creates an array of crt files sufficient to handle n coefficients
int crt_create_files (crt_file_t **fpp, int id1, int id2, int n)
{
	int i, k;
	
	k = ui_ceil_ratio(n,CRT_MAX_COUNT);
	*fpp = malloc(k*sizeof(crt_file_t));
	for ( i = 0 ; i < k ; i++ )
		if ( ! crt_file_create ((*fpp)[i], id1, id2, i*CRT_MAX_COUNT, ((i+1)*CRT_MAX_COUNT< n ? CRT_MAX_COUNT : n-i*CRT_MAX_COUNT)) ) { err_printf ("Unable to create crt output file in temp directory\n"); abort(); }
	return k;
}

void crt_write_files (crt_file_t *fp, unsigned long m, unsigned long c[], int n)
{
	int i,k;
	
	k = ui_ceil_ratio(n,CRT_MAX_COUNT);
	for ( i = 0 ; i < k ; i++ )
		crt_file_write (fp[i], m, c+i*CRT_MAX_COUNT, ((i+1)*CRT_MAX_COUNT<n ? CRT_MAX_COUNT : n-i*CRT_MAX_COUNT));
}

void crt_close_files (crt_file_t *fp, int n)
{
	int i,k;
	
	k = ui_ceil_ratio(n,CRT_MAX_COUNT);
	for ( i = 0 ; i < k ; i++ ) crt_file_close (fp[i]);
	mem_free (fp);
}

int crt_process_files (long *totbits, long *maxbits, int id1, int id2, int n, unsigned long m[], int k, int (*callback)(mpz_t C, int i, void *ctx), void *ctx)
{
	mpz_t X, Y;
	crt_file_t fp;
	crt_tree_t t;
	unsigned long *r, *s, *a, p;
	long bits;
	register int i, j;
	int sts, cnt, soff, scnt, offset, batch;

	mpz_init (X);  mpz_init (Y);
	if ( (long)k*(long)CRT_MAX_COUNT*sizeof(unsigned long) > CRT_MAX_MEM ) {
		out_printf ("Exceeded CRT_MAX_MEM in crt_process_files, switching to batching\n");
		batch = CRT_MAX_MEM / ((long)k*sizeof(unsigned long));
		if ( batch < 1 ) { err_printf ("Can't handle k=%d moduli, even with batch size of 1\n", k);  mpz_clear(X);  mpz_clear(Y);  return 0; }
		s = malloc (CRT_MAX_COUNT*sizeof(unsigned long));
	} else {
		batch = CRT_MAX_COUNT;
		s = 0;
	}
	r = malloc(k*batch*sizeof(unsigned long));
	a = malloc(k*sizeof(unsigned long));
	crt_tree_init (t, m, k);
	sts = 0;
	if ( totbits ) *totbits =0;  if ( maxbits ) *maxbits = 0;
	for ( offset = 0 ; offset < n ; offset += cnt ) {
		cnt = crt_file_open (fp, id1, id2, offset);
		if ( ! cnt ) { err_printf ("Error opening crt file id1=%d, id2=%d, offset=%d\n", id1, id2, offset); goto done; }
		for ( soff = 0 ; soff < cnt ; soff += scnt ) {
			scnt = cnt-soff;
			if ( scnt > batch ) scnt = batch;
			if ( cnt > scnt ) {
				for ( i = 0 ; i < k && crt_file_read(fp, &p, s, cnt) ; i++ ) {
					if ( p != m[i] ) { err_printf ("Error, unexpected modulus %lu != %lu at record %d\n", p, m[i], i); goto done; }
					memcpy (r+i*scnt, s+soff, scnt*sizeof(unsigned long));
				}
			} else {
				for ( i = 0 ; i < k && crt_file_read(fp, &p, r+i*scnt, scnt) ; i++ )
					if ( p != m[i] ) { err_printf ("Error, unexpected modulus %lu != %lu at record %d\n", p, m[i], i); goto done; }	// read cnt coeffcients for each of k moduli, store in r[]
			}
			if ( i != k ) { err_printf ("Error, found %d != %d records in crt file\n", i, k);  goto done; }		
			for ( j = 0 ; j < scnt ; j++ ) {
				for ( i = 0 ; i < k ; i++ ) a[i] = r[i*scnt+j];
				crt_tree_eval (X, a, t);
				crt_mpz_mod_to_Z (X, t->MM);
				bits = mpz_sizeinbase(X,2);
				if ( ! (*callback) (X, offset+j, ctx) ) goto done;
				if ( totbits ) *totbits += bits;
				if ( maxbits && bits > *maxbits ) *maxbits = bits;
			}
		}
		crt_file_delete(fp);
	}
	if ( offset != n ) { err_printf ("Error, found a total of %d != %d coefficients in crt files\n", offset, n); goto done; }
	sts = 1;
done:
	crt_tree_clear (t);
	free (r);  free (a);
	mpz_clear (X);  mpz_clear(Y);
	return sts;	
}

void crt_tree_init (crt_tree_t t, unsigned long m[], int n)
{
	register int i, j, k;
	register mpz_t *B;
	
	for ( i = 0, j = n ; j > CRT_BASE_N ; i++, j = j-j/2 );
	if ( i >= CRT_MAX_LEVELS ) { printf ("Exceeded CRT_MAX_LEVELS=%d for n=%d in crt_tree_init\n", CRT_MAX_LEVELS, n); exit (0); }
	t->levels = i;
	t->n = n;
	t->m = (unsigned long *)malloc(n*sizeof(*t->m));
	for ( i = 0 ; i < n ; i++ ) t->m[i] = m[i];
	t->a = (unsigned long *)malloc(n*sizeof(*t->a));
	t->cnts = (int *)malloc((1<<t->levels)*sizeof(*t->cnts));
	t->lens[t->levels] = 1;  t->cnts[0] = n;
	// Compute cnts by splitting n down the tree.  We only need the cnts at the base, so we overwrite at each level.
	for ( i = t->levels-1; i >= 0 ; i-- ) {
		t->lens[i] = 2*t->lens[i+1];
		for ( j = t->lens[i+1]-1 ; j >= 0 ; j-- ) {
			k = t->cnts[j];
			t->cnts[2*j] = k/2;
			t->cnts[2*j+1] =  k - k/2;
		}
	}
	if ( t->levels ) {
		t->M[0] = (mpz_t *)malloc((1<<(t->levels+1))*sizeof(mpz_t));
		t->A[0] = (mpz_t *)malloc((1<<(t->levels+1))*sizeof(mpz_t));
		for ( i = 1 ; i < t->levels ; i++ ) { t->M[i] = t->M[i-1]+t->lens[i-1]; t->A[i] = t->A[i-1]+t->lens[i-1]; }
	}
	mpz_init(t->M0);  mpz_init(t->M1);  mpz_init(t->MM);
	t->B = (mpz_t *)malloc(n*sizeof(mpz_t));
	B = t->B;  k = 0;
	for ( i = 0 ; i < t->lens[0] ; i++ ) {
		mpz_set_ui(t->M0,m[k]);
		for ( j = 1 ; j < t->cnts[i] ; j++ ) mpz_mul_ui(t->M0,t->M0,m[k+j]);
		for ( j = 0 ; j < t->cnts[i] ; j++ ) {  mpz_divexact_ui(t->M1,t->M0,m[k+j]); mpz_init_set(*B++,t->M1); }
		if ( t->levels ) { mpz_init_set (t->M[0][i], t->M0); mpz_init2 (t->A[0][i], mpz_sizeinbase(t->M0,2)); }
		k += j;
	}
	for ( i = 1 ; i < t->levels ; i++ )
		for ( j = 0 ; j < t->lens[i] ; j++ ) { mpz_mul(t->M1,t->M[i-1][2*j],t->M[i-1][2*j+1]); mpz_init_set(t->M[i][j], t->M1); mpz_init2(t->A[i][j],mpz_sizeinbase(t->M1,2)); }
	if ( t->levels ) {
		mpz_mul(t->MM, t->M[t->levels-1][0], t->M[t->levels-1][1]);
	} else {
		mpz_mul_ui(t->MM,t->B[0], t->m[0]);
	}
			
	// compute a[i] = (prod_{j!=i}m[j]) mod m[i]
	if ( t->levels ) {
		i = t->levels-1;
		mpz_mod(t->A[i][0],t->M[i][1],t->M[i][0]);  mpz_mod(t->A[i][1],t->M[i][0],t->M[i][1]);
		for ( i-- ; i >= 0 ; i-- ) {
			for ( j = 0 ; j < t->lens[i] ; j+=2 ) {
				mpz_mod(t->M0,t->A[i+1][j/2],t->M[i][j]);					// reduce parent complement mod left child
				mpz_mul(t->M1,t->M0,t->M[i][j+1]);						// multiply by right child
				mpz_mod(t->A[i][j],t->M1,t->M[i][j]);						// reduce mod left child
				mpz_mod(t->M0,t->A[i+1][j/2],t->M[i][j+1]);					// reduce parent complement mod right child
				mpz_mul(t->M1,t->M0,t->M[i][j]);							// multiply by left child
				mpz_mod(t->A[i][j+1],t->M1,t->M[i][j+1]);					// reduce mod right child
			}
		}
		k = 0;
		for ( i = 0 ; i < t->lens[0] ; i++ ) {
			for ( j = 0 ; j < t->cnts[i] ; j++ ) {
				mpz_mul(t->M0,t->A[0][i],t->B[k]);
				t->a[k] = ui_inverse(mpz_fdiv_ui(t->M0,m[k]), m[k]);
				k++;
			}
		}
	} else {
		for ( i = 0 ; i < t->n ; i++ ) t->a[i] = ui_inverse(mpz_fdiv_ui (t->B[i],m[i]), m[i]);
	}
	for ( i = 0 ; i < t->levels ; i++ ) t->X[i] = 0;								// will get allocated in uneval if called
}

void crt_tree_clear (crt_tree_t t)
{
	register int i, j;

	free(t->a);  free(t->m);  free (t->cnts);
	if ( t->levels ) {
		for ( i = 0 ; i < t->levels ; i++ ) for ( j = 0 ; j < t->lens[i] ; j++ ) { mpz_clear(t->A[i][j]); mpz_clear(t->M[i][j]); if ( t->X[i] ) mpz_clear(t->X[i][j]); }
		free (t->A[0]);  free (t->M[0]);   if ( t->X[0] ) free (t->X[0]);
	}
	for ( i = 0 ; i < t->n ; i++ ) mpz_clear(t->B[i]);
	free(t->B);
	mpz_clear(t->M0);  mpz_clear(t->M1);  mpz_clear (t->MM);
}


void crt_tree_eval (mpz_t C, unsigned long c[], crt_tree_t t)
{
	register int i, j, k;
	
	if ( ! t->levels ) {
		mpz_set_ui(t->M1,c[0]);  mpz_mul_ui(t->M0,t->M1,t->a[0]); mpz_mul(C,t->M0,t->B[0]);
		for ( i = 1 ; i < t->cnts[0] ; i++ ) {
			mpz_set_ui(t->M1,c[i]);  mpz_mul_ui(t->M0,t->M1,t->a[i]);  mpz_addmul(C,t->M0,t->B[i]);
		}
		mpz_mod(C,C,t->MM);
		return;
	}
	k = 0;
	for ( i = 0 ; i < t->lens[0] ; i++ ) {
		mpz_set_ui(t->M1,c[k]);  mpz_mul_ui(t->M0,t->M1,t->a[k]); mpz_mul(t->A[0][i],t->M0,t->B[k]);
		for ( j = 1 ; j < t->cnts[i] ; j++ ) {
			mpz_set_ui(t->M1,c[k+j]);  mpz_mul_ui(t->M0,t->M1,t->a[k+j]);  mpz_addmul(t->A[0][i],t->M0,t->B[k+j]);			
		}
		k += j;
	}
	for ( i = 1 ; i < t->levels ; i++ ) {
		for ( j = 0 ; j < t->lens[i] ; j++ ) {
			mpz_mul(t->A[i][j], t->A[i-1][2*j], t->M[i-1][2*j+1]);
			mpz_addmul(t->A[i][j], t->A[i-1][2*j+1], t->M[i-1][2*j]);
		}
	}
	mpz_mul(C,t->A[i-1][0],t->M[i-1][1]);
	mpz_addmul(C,t->A[i-1][1],t->M[i-1][0]);
	mpz_mod(C,C,t->MM);
}

void crt_tree_uneval (unsigned long c[], mpz_t C, crt_tree_t t)
{
	register int i, j, k;

	if ( ! t->levels ) {
		for ( i = 0 ; i < t->n ; i++ ) c[i] = mpz_fdiv_ui(C,t->m[i]);
		return;
	}
	if ( ! t->X[0] ) {
		t->X[0] = (mpz_t *)malloc((1<<(t->levels+1))*sizeof(mpz_t));
		for ( j = 0 ; j < t->lens[0] ; j++ ) mpz_init2 (t->X[0][j], mpz_sizeinbase(t->M[0][j],2));
		for ( i = 1 ; i < t->levels ; i++ ) {
			t->X[i] = t->X[i-1]+t->lens[i-1];
			for ( j = 0 ; j < t->lens[i] ; j++ ) mpz_init2 (t->X[i][j], mpz_sizeinbase(t->M[i][j],2));
		}
	}
	i = t->levels-1;
	mpz_mod(t->X[i][0],C,t->M[i][0]);  mpz_mod(t->X[i][1],C,t->M[i][1]);
	for ( i-- ; i>=0 ; i-- ) for ( j = 0 ; j < t->lens[i] ; j++ ) mpz_mod(t->X[i][j],t->X[i+1][j/2],t->M[i][j]);
	for ( i = j = 0 ; i < t->lens[0] ; i++ ) for ( k = 0 ; k < t->cnts[i] ; j++, k++ ) c[j] = mpz_fdiv_ui(t->X[0][i],t->m[j]);
}


// overlap is not permitted, returns 1 if C is stable (i.e. there is an integer |X| <= M1 such that X = C mod M and X = C1 mod M1)
int crt_mpz_ui (mpz_t C, mpz_t M, mpz_t C1, mpz_t M1, unsigned long c2, unsigned long m2)
{
	static int init;
	static mpz_t X, Y;
	unsigned long r;
	long a, b;
	
	if ( ! init ) { mpz_init(X); mpz_init(Y); init = 1; }
	if ( mpz_cmp_ui(M1,1)==0 ) { mpz_set_ui(M,m2); mpz_set_ui(C,c2); return 0; }
	r = mpz_fdiv_q_ui(C,M1,m2);
	ui_gcd_ext (m2, r, &a, &b);
	if ( b < 0 ) { mpz_mul_ui(M,C,-b);  mpz_neg(M,M); } else { mpz_mul_ui(M,C,b); }
	if ( a < 0 ) { mpz_add_ui(C,M,-a); } else { mpz_sub_ui(C,M,a); }
	mpz_neg(C,C);
	mpz_mul_ui(M,C,m2);
	mpz_mul(C,M,C1);
	if ( b < 0 ) { mpz_mul_ui(M,M1,-b); mpz_neg(M,M); } else { mpz_mul_ui(M,M1,b); }
	mpz_mul_ui(M,M,c2);
	mpz_add(C,C,M);
	mpz_mul_ui(M,M1,m2);
	mpz_mod(C,C,M);
	if ( mpz_cmp(C,C1)==0 ) return 1;
	mpz_sub(X,C,M); mpz_sub(Y,C1,M1);
	if ( mpz_cmp(X,Y)==0 ) return 1;
	return 0;
}

// helper function to handle ecrt datafile creation and initialization
void _ecrt_file_init (ecrt_context_t ecrt)
{
	char filename[256];
	uint64_t x;
	
	if ( ! ecrt->jobs ) { err_printf("_ecrt_file_init called when jobs was not specified\n"); exit (0); }
	sprintf (filename, "%s_%d_%d.ecrt", ecrt->prefix, ecrt->jobs, ecrt->jobid);
	if ( ! (ecrt->fp = fopen(filename, "wb+")) ) { err_printf ("Unable to create file %s!\n", filename); exit(0); }
	// write four 32-bit integers containing n, k, P limbs, and MP limbs
	x = ecrt->n;  fwrite (&x,8,1,ecrt->fp);  x = ecrt->k; fwrite (&x,8,1,ecrt->fp); x = ecrt->P->_mp_size;  fwrite (&x,8,1,ecrt->fp);  x = ecrt->MP->_mp_size;  fwrite (&x,8,1,ecrt->fp);
	if ( fwrite (ecrt->P->_mp_d, sizeof(mp_limb_t), ecrt->P->_mp_size, ecrt->fp) != ecrt->P->_mp_size ) { err_printf ("File write error in ecrt_dump\n"); exit (0); }					// write P limbs
	if ( fwrite (ecrt->MP->_mp_d, sizeof(mp_limb_t), ecrt->MP->_mp_size, ecrt->fp) != ecrt->MP->_mp_size ) { err_printf ("File write error in ecrt_dump\n"); exit (0); }				// write MP limbs
	ecrt->fpbase = ftell(ecrt->fp);
	dbg_printf ("Created ecrt file %s\n", filename);
}

// Performs initialization tasks unrelated to the CRT tree
void _ecrt_init (ecrt_context_t ecrt, int n, int k, mpz_t P, mpz_t MP, int jobs, int jobid, char *prefix)
{
	int nbits;
	long mem;
	int i;

	if ( prefix ) {
		if ( strlen(prefix) >= sizeof(ecrt->prefix) ) { err_printf ("prefix too long in _ecrt_init\n"); exit (0); }
		strcpy(ecrt->prefix,prefix);
	} else {
		 ecrt->prefix[0] = '\0';
	}
	if ( jobs > 1 && ! prefix ) { err_printf ("prefix required in _ecrt_init for jobs=%d\n", jobs); exit (0); }
	if ( jobs > CRT_MAX_JOBS ) { err_printf ("Exceeded CRT_MAX_JOBS=%d\n", CRT_MAX_JOBS); exit (0); }
	if ( jobs <= 1 ) { ecrt->jobs = 0; ecrt->jobid = 0; } else { ecrt->jobs = jobs; ecrt->jobid = jobid; }
	if ( ecrt->jobs && (ecrt->jobid < 1 || ecrt->jobid > ecrt->jobs) ) { err_printf ("invalid jobid=%d for jobs=%d in _ecrt_init\n", jobid, jobs); exit (0); }
	
	ecrt->n = n;
	ecrt->k = k;
	nbits = ui_len(n);
	ecrt->Climbs = ui_ceil_ratio (mpz_sizeinbase(P,2) + 63 + nbits, GMP_NUMB_BITS);  ecrt->slimbs = 2;
	info_printf ("ecrt init %d of %d jobs %d primes %d coefficients, P has %lu bits, Climbs=%d, slimbs=%d\n", ecrt->jobid, ecrt->jobs, ecrt->n, ecrt->k, mpz_sizeinbase(P,2), ecrt->Climbs, ecrt->slimbs);
	mpz_init2(ecrt->X,(ecrt->Climbs+ecrt->slimbs)*GMP_NUMB_BITS); mpz_init2(ecrt->Y,(ecrt->Climbs+ecrt->slimbs)*GMP_NUMB_BITS); mpz_init_set(ecrt->P,P); mpz_init_set(ecrt->MP,MP);
	ecrt->Cbytes = ecrt->Climbs * sizeof(mp_limb_t);
	ecrt->sbytes = ecrt->slimbs * sizeof(mp_limb_t);
	mem =  k * (ecrt->Climbs+ecrt->slimbs) * sizeof(mp_limb_t);
	info_printf ("ecrt mem = %d*(%ld+%ld) + %d*%ld = %ld\n",k, ecrt->Climbs*sizeof(mp_limb_t), ecrt->slimbs*sizeof(mp_limb_t), n, sizeof(*ecrt->a)+sizeof(*ecrt->m),
			 k*(ecrt->Climbs+ecrt->slimbs)*sizeof(mp_limb_t) + n*(sizeof(*ecrt->a)+sizeof(*ecrt->m)));
	if ( mem > CRT_MAX_MEM ) {
		out_printf ("ecrt mem=%ld exceeds CRT_MAX_MEM, using batching\n", mem);
		ecrt->batching = 1;
		ecrt->Cdata = (mp_limb_t *) calloc (CRT_BATCH_SIZE,ecrt->Cbytes);
		ecrt->sdata = (mp_limb_t *) calloc (CRT_BATCH_SIZE,ecrt->sbytes);
	} else {
		ecrt->batching = 0;  ecrt->fp = 0;
		if ( ecrt->jobs ) k = CRT_BATCH_SIZE*ui_ceil_ratio (k, CRT_BATCH_SIZE);		// round up to a multiple of batch size if we know that we will be dumping
		ecrt->Cdata = (mp_limb_t *) calloc (k,ecrt->Climbs*sizeof(mp_limb_t));
		ecrt->sdata = (mp_limb_t *) calloc (k,ecrt->slimbs*sizeof(mp_limb_t));
	}
	if ( ecrt->batching && ! ecrt->jobs ) { out_printf ("ecrt mem=%ld exceeds max CRT mem so batching is required.  Please use job processing\n", mem); exit (0); }
	if ( ecrt->batching ) {
		_ecrt_file_init (ecrt);
		// write out zero values to initialize (this could be avoided but makes life simpler and catches any file allocation issues up front)
		for ( i = 0 ; i < ecrt->k ; i += CRT_BATCH_SIZE ) {
			if ( fwrite(ecrt->Cdata, ecrt->Cbytes, CRT_BATCH_SIZE, ecrt->fp) != CRT_BATCH_SIZE ) { err_printf ("file initialization error in ecrt_post_init\n"); exit (0); }
			if ( fwrite(ecrt->sdata, ecrt->sbytes, CRT_BATCH_SIZE, ecrt->fp) != CRT_BATCH_SIZE ) { err_printf ("file initialization error in ecrt_post_init\n"); exit (0); }
		}
	}
	for ( i = 0 ; i < ecrt->jobs ; i++ ) ecrt->infp[i] = 0;
	ecrt->delta = nbits+1;
	ecrt->j = -1;
}


// Agorithm 2.3 of Hilbert CRT paper.
void ecrt_init (ecrt_context_t ecrt, unsigned long m[], int n, int k, mpz_t P, int jobs, int jobid, char *prefix)
{
	mpz_t *W, X;
	register int i, j;
	
	// allocate temporary space used to compute a_i, delay all other allocations until we free W
	j = n+ui_len(n)+3;
	W = (mpz_t *)malloc(j*sizeof(*W));  for ( i = 0 ; i < j ; i++ ) mpz_init(W[i]);
	ecrt->a = (unsigned long *)malloc (n*sizeof(*ecrt->a));
	ui_crt_coeff (ecrt->a, m, n, W);														// Compute a[i] = M_i mod m[i] where M_i = prod_{j\ne i} m_i
	for ( i = j-1 ; i >= 0 ; i-- ) mpz_clear(W[i]);
	free (W);

	// Compute M mod P iteratively, don't bother with a tree (unless P is huge, it wouldn't make much of a difference)
	mpz_init_set_ui(X,m[0]);  
	for ( i = 1 ; i < n ; i++) { mpz_mul_ui(X,X,m[i]); mpz_mod(X,X,P); }
	_ecrt_init (ecrt, n, k, P, X, jobs, jobid, prefix);
	mpz_clear(X);
	
	ecrt->m = (unsigned long *)malloc (n*sizeof(*ecrt->m));
	for ( i = 0 ; i < n ; i++ ) ecrt->m[i] = m[i];
//	ecrt->d = (mpz_t *)malloc (n*sizeof(*ecrt->d));  for ( i = 0 ; i < n ; i++ ) mpz_init2(ecrt->d[i],mpz_sizeinbase(P,2));

	for ( i = 0 ; i < n ; i++ ) {
		ecrt->a[i] = ui_inverse(ecrt->a[i],m[i]);												// a holds 1/M_i mod m[i]
//		mpz_set_ui(ecrt->X,m[i]); mpz_invert(ecrt->X,ecrt->X,P);
//		mpz_mul (ecrt->X,ecrt->X,ecrt->MP); mpz_mod(ecrt->X,ecrt->X,P);						// compute M_i mod P using division mod P not a CRT tree (faster for log P = o(log^3 maecrt->Xm))
//		mpz_mul_ui(ecrt->X,ecrt->X,ecrt->a[i]); mpz_mod (ecrt->d[i],ecrt->X,P);					// d[i] = a_iM_i mod P
//gmp_printf ("a_%d = %lu, d_%d = %Zd\n", i, ecrt->a[i], i, ecrt->d[i]);
	}
}

void ecrt_clear (ecrt_context_t ecrt)
{
	register int i;
	
	if ( ecrt->fp ) fclose(ecrt->fp);
	for ( i = 0 ; i < ecrt->jobs ; i++ ) if ( ecrt->infp[i] ) fclose (ecrt->infp[i]);
	mpz_clear (ecrt->MP); mpz_clear (ecrt->X); mpz_clear (ecrt->Y); mpz_clear(ecrt->P);
//	for ( i = 0 ; i < ecrt->n ; i++ ) mpz_clear (ecrt->d[i]);  free (ecrt->d);
	free (ecrt->a); free (ecrt->Cdata); free (ecrt->sdata);
	memset(&ecrt,0,sizeof(ecrt));
}

/*
	To avoid disk thrashing when multiple jobs are running on the same machine, acquire an exclusive lock
	on a file in the local directory (note this only works when other jobs are running in the same directory)
*/
FILE *_ecrt_lock (void)
{
	FILE *fp;
	
	fp = fopen("ecrt.lck", "w");
	if ( !fp ) { err_printf ("Error opening ecrt.lck\n"); exit(0); }
	if (lockf(fileno(fp),F_LOCK,0) < 0 ) { err_printf ("flock failed with errno=%d\n",errno); exit(0); }
	return fp;
}

void _ecrt_unlock (FILE *fp)
{
	if ( lockf(fileno(fp),F_ULOCK,0) < 0 ) { err_printf ("flock failed with errno=%d\n",errno); exit(0); }
	fclose(fp);
}

// Algorithm 2.4 of Hilbert CRT paper, augmented to support batching
void ecrt_update (ecrt_context_t ecrt, int i, unsigned long c[], int k)
{
	time_t start, end;
	FILE *lockfp;
	long pos, totbytes;
	register int j0,j;
	register unsigned long a, m;

	if ( k != ecrt->k ) { err_printf ("Error: k=%d does not match ecrt->k=%d in ecrt_update\n", k, ecrt->k); exit (0); }
	if ( ecrt->j >= 0 ) { err_printf ("Error: call to ecrt_update while coefficient enumeration is in progress\n"); exit (0); }
	if ( i > ecrt->n ) { err_printf ("Error: call to ecrt_update with modulus index %d > modulus count %d\n", i, ecrt->n); exit (0); } 
	a = ecrt->a[i];  m = ecrt->m[i];
	// compute d[i] = a_iM_i = a(M/m) mod P
	mpz_set_ui(ecrt->X,m); mpz_invert(ecrt->X,ecrt->X,ecrt->P);
	mpz_mul (ecrt->X,ecrt->X,ecrt->MP); mpz_mod(ecrt->X,ecrt->X,ecrt->P);					// compute M_i mod P using division mod P
	mpz_mul_ui(ecrt->X,ecrt->X,a); mpz_mod (ecrt->Y,ecrt->X,ecrt->P);						// Y = d[i] = a_iM_i mod P
	if ( ecrt->batching ) {
		if ( ! ecrt->fp ) { err_printf ("Error, null file pointer in ecrt_update while using batching\n"); exit (0); }
		// acquire exclusive lock to avoid disk thrashing
		start = time(0);
		lockfp = _ecrt_lock();
		end = time(0);
		if ( end > start+1 ) out_printf ("Waited %ld secs for lock in ecrt_update\n", end-start);
		fseek(ecrt->fp, ecrt->fpbase, SEEK_SET);
		for ( j0=j= 0 ; j< k ; j0++,j++) {
			if ( j0 == CRT_BATCH_SIZE ) {
				if ( fwrite (ecrt->Cdata, ecrt->Cbytes, CRT_BATCH_SIZE, ecrt->fp) != CRT_BATCH_SIZE ) { err_printf ("file write error %d in ecrt_update\n", errno); exit (0); }
				if ( fwrite (ecrt->sdata, ecrt->sbytes, CRT_BATCH_SIZE, ecrt->fp) != CRT_BATCH_SIZE ) { err_printf ("file write error %d in ecrt_update\n", errno); exit (0); }
				j0 = 0; 
			}
			if ( j0 == 0 ) {
				pos = ftell(ecrt->fp);
				if ( fread (ecrt->Cdata, ecrt->Cbytes, CRT_BATCH_SIZE, ecrt->fp) != CRT_BATCH_SIZE ) { err_printf ("file read error %d in ecrt_update\n", errno); exit (0); }
				if ( fread (ecrt->sdata, ecrt->sbytes, CRT_BATCH_SIZE, ecrt->fp) != CRT_BATCH_SIZE ) { err_printf ("file read error %d in ecrt_update\n", errno); exit (0); }
				fseek (ecrt->fp, pos, SEEK_SET);
			}
			mpz_mul_ui(ecrt->X,ecrt->Y,c[j]);		// dont' reduce mod P to save time
			if ( mpn_add (ecrt->Cdata+j0*ecrt->Climbs, ecrt->Cdata+j0*ecrt->Climbs, ecrt->Climbs, ecrt->X->_mp_d, ecrt->X->_mp_size) ) { err_printf ("overflow on C in ecrt"); exit (0); }
			mpz_set_ui(ecrt->X,c[j]);  mpz_mul_ui(ecrt->X,ecrt->X,a);  mpz_mul_2exp(ecrt->X,ecrt->X,ecrt->delta);  mpz_fdiv_q_ui(ecrt->X,ecrt->X,m);	// compute floor((2^delta*c[j]*a[i])/m[i])
			if ( mpn_add (ecrt->sdata+j0*ecrt->slimbs, ecrt->sdata+j0*ecrt->slimbs, ecrt->slimbs, ecrt->X->_mp_d, ecrt->X->_mp_size) ) { err_printf ("overflow on s in ecrt"); exit (0); }
		}
		if ( j0 ) {
			// we write a full batch at the end, even though the trailing entries may be unused (these will be zero)
			if ( fwrite (ecrt->Cdata, ecrt->Cbytes, CRT_BATCH_SIZE, ecrt->fp) != CRT_BATCH_SIZE ) { err_printf ("file write error %d in ecrt_update\n", errno); exit (0); }
			if ( fwrite (ecrt->sdata, ecrt->sbytes, CRT_BATCH_SIZE, ecrt->fp) != CRT_BATCH_SIZE ) { err_printf ("file write error %d in ecrt_update\n", errno); exit (0); }
		}
		_ecrt_unlock(lockfp);
		end = time(0);
		totbytes =  CRT_BATCH_SIZE*  ui_ceil_ratio(k,CRT_BATCH_SIZE) * (ecrt->Cbytes+ecrt->sbytes);
		if ( end > start+1 )
			out_printf ("Wrote %d coefficients (%ld bytes of data) for modulus m[%d]=%ld to disk in %ld secs, throughput %.2f MB/s\n", k, totbytes, i, ecrt->m[i], end-start, (double)totbytes / (1000000.0*(end-start)));
	} else {
		for ( j = 0 ; j < k ; j++ ) {
			mpz_mul_ui(ecrt->X,ecrt->Y,c[j]);		// don't reduce mod P to save time
			if ( mpn_add (ecrt->Cdata+j*ecrt->Climbs, ecrt->Cdata+j*ecrt->Climbs, ecrt->Climbs, ecrt->X->_mp_d, ecrt->X->_mp_size) ) { err_printf ("overflow on C in ecrt"); exit (0); }
			mpz_set_ui(ecrt->X,c[j]);  mpz_mul_ui(ecrt->X,ecrt->X,a); mpz_mul_2exp(ecrt->X,ecrt->X,ecrt->delta);  mpz_fdiv_q_ui(ecrt->X,ecrt->X,m);	// compute floor((2^delta*c[j]*a[i])/m[i])
			if ( mpn_add (ecrt->sdata+j*ecrt->slimbs, ecrt->sdata+j*ecrt->slimbs, ecrt->slimbs, ecrt->X->_mp_d, ecrt->X->_mp_size) ) { err_printf ("overflow on s in ecrt"); exit (0); }
		}
	}
}

void ecrt_reset (ecrt_context_t ecrt)
{
	if ( ecrt->batching ) { err_printf ("Reset not permitted in ecrt when batching is used\n"); exit (0); }
	memset (ecrt->Cdata, 0, ecrt->Cbytes);
	memset (ecrt->sdata, 0, ecrt->sbytes);
	ecrt->j = -1;
}

// assumes X->_mp_alloc >= n and does not bother verifying this!  returned X will always be positive.
static inline void _mpn_limbs_to_mpz (mpz_t X, mp_limb_t *s, int n)
	{ register int i;  for ( i = n-1 ; i >= 0 && ! s[i] ; i-- ); for ( X->_mp_size=i+1 ; i >= 0 ; i-- ) X->_mp_d[i] = s[i]; }

// assumes X->_mp_size <= n and does not bother verifying this!  also assumes X is positive
static inline void _mpn_mpz_to_limbs (mp_limb_t *s, int n, mpz_t X)
	{ register int i;  for ( i = 0 ; i < X->_mp_size ; i++ ) s[i]=X->_mp_d[i];  while ( i < n ) s[i++] = 0; }

void ecrt_merge (ecrt_context_t ecrt, int jobs, char *prefix)
{
	char filename[256];
	uint64_t x, Plimbs, MPlimbs;
	mpz_t P, MP;
	int i, n, k, nbits;

	if ( jobs < 1 || jobs > CRT_MAX_JOBS ) { err_printf ("Invalid jobs=%d in ecrt_merge\n", jobs); exit(0); }
	ecrt->jobs = jobs;
	ecrt->batching = 1;
	for ( i = 0 ; i < jobs ; i++ ) {
		sprintf (filename, "%s_%d_%d.ecrt", prefix, jobs, i+1);
		ecrt->infp[i] = fopen (filename, "rb");
		if ( ! ecrt->infp[i] ) { err_printf ("Error opening ecrt job file %s\n", filename); exit (0); }
		dbg_printf ("Opened %s\n", filename);
		// read four 64-bit integers countaining n, k, and P limbs,
		if ( fread (&x,8,1,ecrt->infp[i]) != 1 ) { err_printf ("File read error in ecrt_merge\n"); exit (0); }
		n = (int)x;
		if ( fread(&x,8,1,ecrt->infp[i]) != 1 ) { err_printf ("File read error in ecrt_merge\n"); exit (0); }
		k = (int)x;
		if ( fread(&Plimbs,8,1,ecrt->infp[i]) != 1 ) { err_printf ("File read error in ecrt_merge\n"); exit (0); }
		if ( fread(&MPlimbs,8,1,ecrt->infp[i]) != 1 ) { err_printf ("File read error in ecrt_merge\n"); exit (0); }
		if ( !i ) {
			mpz_init2(P,Plimbs*GMP_NUMB_BITS); mpz_init2(MP,Plimbs*GMP_NUMB_BITS);   P->_mp_size = Plimbs;  MP->_mp_size = MPlimbs;
			if ( fread (P->_mp_d, sizeof(mp_limb_t), Plimbs, ecrt->infp[i]) != Plimbs ) { err_printf ("File read error in ecrt_merge\n"); exit (0); }
			if ( fread (MP->_mp_d, sizeof(mp_limb_t), MPlimbs, ecrt->infp[i]) != MPlimbs ) { err_printf ("File read error in ecrt_merge\n"); exit (0); }
			ecrt->n = n;  ecrt->k = k;  mpz_init_set(ecrt->P,P); mpz_init_set(ecrt->MP,MP);
		} else {
			if ( n != ecrt->n || k != ecrt->k || Plimbs != P->_mp_size || MPlimbs != MP->_mp_size ) {
				err_printf ("Inconsistent ecrt data in file %d of %d: n=%d vs %d, k=%d vs %d, Psize=%d vs %lu, MPsize=%d vs %lu\n",
						 i, jobs, ecrt->n, n, ecrt->k, k, ecrt->P->_mp_size, Plimbs, ecrt->MP->_mp_size, MPlimbs);
				exit (0);
			}
			if ( fread (P->_mp_d, sizeof(mp_limb_t), Plimbs, ecrt->infp[i]) != Plimbs ) { err_printf ("File read error in ecrt_merge\n"); exit (0); }
			if ( fread (MP->_mp_d, sizeof(mp_limb_t), MPlimbs, ecrt->infp[i]) != MPlimbs ) { err_printf ("File read error in ecrt_merge\n"); exit (0); }
			if ( mpz_cmp(P,ecrt->P) || mpz_cmp(MP,ecrt->MP) ) { err_printf ("Inconsistent P or MP in ecrt data file %d of %d\n", i, jobs); exit (0); }
		}
	}
	ecrt->fpbase = ftell(ecrt->infp[0]);
	nbits = ui_len(ecrt->n);
	ecrt->Climbs = ui_ceil_ratio (mpz_sizeinbase(ecrt->P,2) + 63 + nbits, GMP_NUMB_BITS);  ecrt->slimbs = 2;
	info_printf ("ecrt merge %d jobs %d primes %d coefficients, P has %lu bits, Climbs=%d, slimbs=%d\n", ecrt->jobs, ecrt->n, ecrt->k, mpz_sizeinbase(P,2), ecrt->Climbs, ecrt->slimbs);
	mpz_init2(ecrt->X,(ecrt->Climbs+ecrt->slimbs)*GMP_NUMB_BITS); mpz_init2(ecrt->Y,(ecrt->Climbs+ecrt->slimbs)*GMP_NUMB_BITS);
	ecrt->Cbytes = ecrt->Climbs * sizeof(mp_limb_t);
	ecrt->sbytes = ecrt->slimbs * sizeof(mp_limb_t);
	mpz_clear(P); mpz_clear (MP);
	ecrt->Cdata = (mp_limb_t *)malloc(ecrt->Cbytes*CRT_BATCH_SIZE); ecrt->inCdata = (mp_limb_t *)malloc(ecrt->Cbytes*CRT_BATCH_SIZE);
	ecrt->sdata = (mp_limb_t *)malloc(ecrt->sbytes*CRT_BATCH_SIZE); ecrt->insdata = (mp_limb_t *)malloc(ecrt->sbytes*CRT_BATCH_SIZE);
	ecrt->delta = nbits+1;
	ecrt->j0 = ecrt->j = 0;
}

// dumps current coefficient sums -- required for job processing when batching is not being used
void ecrt_dump (ecrt_context_t ecrt)
{
	int i;

	if ( ! ecrt->jobs ) { err_printf ("ecrt_dump called without job processing enabled\n"); return; }
	if ( ! ecrt->batching ) {
		dbg_printf ("Opening file in ecrt_dump\n");
		_ecrt_file_init(ecrt);
		for ( i = 0 ; i < ecrt->k ; i += CRT_BATCH_SIZE ) {
			if ( fwrite(ecrt->Cdata+i*ecrt->Climbs, ecrt->Cbytes, CRT_BATCH_SIZE, ecrt->fp) != CRT_BATCH_SIZE ) { err_printf ("file write error in ecrt_dump\n"); exit (0); }
			if ( fwrite(ecrt->sdata+i*ecrt->slimbs, ecrt->sbytes, CRT_BATCH_SIZE, ecrt->fp) != CRT_BATCH_SIZE ) { err_printf ("file write error in ecrt_dump\n"); exit (0); }
		}
	}
	fclose (ecrt->fp);  ecrt->fp = 0;
}

static inline void _ecrt_finalize_coeff (ecrt_context_t ecrt, int j)
{
	register unsigned long tf;
	
	 tf = (3UL << (ecrt->delta-2));		// tf = (2^delta*3)/4

	_mpn_limbs_to_mpz (ecrt->Y, ecrt->sdata+j*ecrt->slimbs, ecrt->slimbs);
	mpz_add_ui(ecrt->X,ecrt->Y,tf);  mpz_fdiv_q_2exp(ecrt->X,ecrt->X,ecrt->delta);						// add 2^delta*3/4 to s[j], divide by delta, and take the floor
	_mpn_limbs_to_mpz (ecrt->Y, ecrt->Cdata+j*ecrt->Climbs, ecrt->Climbs);
	mpz_mul(ecrt->X,ecrt->X,ecrt->MP); mpz_sub(ecrt->Y,ecrt->Y,ecrt->X); mpz_mod(ecrt->Y,ecrt->Y,ecrt->P);
	_mpn_mpz_to_limbs (ecrt->Cdata+j*ecrt->Climbs, ecrt->Climbs, ecrt->Y);
}

// Algorithm 2.5 of Hilbert CRT paper
void ecrt_finalize (ecrt_context_t ecrt)
{
	register int j;

	if ( ecrt->batching ) { err_printf ("Invalid call to ecrt_finalize when batching is being used\n"); exit (0); }
	for ( j = 0 ; j < ecrt->k ; j++ ) _ecrt_finalize_coeff(ecrt, j);
	ecrt->j = 0;
}

int ecrt_next_coeff (mpz_t C, ecrt_context_t ecrt)
{
	register int i, j;
	
	if ( ecrt->j < 0 ) { err_printf ("invalid call to ecrt_next_coeff, you must call ecrt_finalize or ecrt_merge first\n"); exit (0); }
	if ( ecrt->j >= ecrt->k ) return 0;
	if ( C->_mp_alloc < ecrt->Climbs ) mpz_realloc2(C,ecrt->Climbs*GMP_NUMB_BITS);
	if ( ecrt->batching ) {
		if ( ! ecrt->j0 || ecrt->j0 == CRT_BATCH_SIZE ) {
			if ( (j=fread(ecrt->Cdata, ecrt->Cbytes, CRT_BATCH_SIZE, ecrt->infp[0])) != CRT_BATCH_SIZE ) { err_printf ("file read error in ecrt_merge, got %d C values\n", j); exit (0); }
			if ( (j=fread(ecrt->sdata, ecrt->sbytes, CRT_BATCH_SIZE, ecrt->infp[0])) != CRT_BATCH_SIZE ) { err_printf ("file read error in ecrt_merge, got %d s values\n", j); exit (0); }		
			for ( i = 1 ; i < ecrt->jobs ; i++ ) {
				if ( fread(ecrt->inCdata, ecrt->Cbytes, CRT_BATCH_SIZE, ecrt->infp[i]) != CRT_BATCH_SIZE ) { err_printf ("file read error in ecrt_merge\n"); exit (0); }
				if ( fread(ecrt->insdata, ecrt->sbytes, CRT_BATCH_SIZE, ecrt->infp[i]) != CRT_BATCH_SIZE ) { err_printf ("file read error in ecrt_merge\n"); exit (0); }		
				for ( j = 0 ; j < CRT_BATCH_SIZE ; j++ ) {				// note that this will not hurt anything if we go past k in the last batch, the values will all be zero
					if ( mpn_add_n (ecrt->Cdata+j*ecrt->Climbs, ecrt->Cdata+j*ecrt->Climbs, ecrt->inCdata+j*ecrt->Climbs, ecrt->Climbs) )  { err_printf ("C overflow in ecrt_next_coeff\n"); exit (0); }
					if ( mpn_add_n (ecrt->sdata+j*ecrt->slimbs, ecrt->sdata+j*ecrt->slimbs, ecrt->insdata+j*ecrt->slimbs, ecrt->slimbs) )  { err_printf ("C overflow in ecrt_next_coeff\n"); exit (0); }
				}
			}
			for ( j = 0 ; j < CRT_BATCH_SIZE ; j++ ) _ecrt_finalize_coeff(ecrt,j);
			ecrt->j0 = 0;
		}
		_mpn_limbs_to_mpz(C, ecrt->Cdata+ecrt->j0*ecrt->Climbs, ecrt->Climbs);
		ecrt->j0++; ecrt->j++;
	} else {
		_mpn_limbs_to_mpz(C, ecrt->Cdata+ecrt->j*ecrt->Climbs, ecrt->Climbs);
		ecrt->j++;
	}
	return 1;
}

/*
	The xcrt code is untested and unused -- in preliminary testing it appears to be too slow to be useful.

void xcrt_init (xcrt_context_t xcrt, unsigned long m[], int n, int k, mpz_t P)
{
	mpz_t *W, M,Y,Z;
	double Mbits, mbits, b;
	register unsigned long a;
	register int i, j;

	// allocate temporary space used to compute a_i, delay all other allocations until we free W
	j = n+ui_len(n)+3;
	xcrt->a = (unsigned long *) malloc (n*sizeof(*xcrt->a));
	W = (mpz_t *)malloc(j*sizeof(*W));  for ( i = 0 ; i < j ; i++ ) mpz_init(W[i]);
	ui_crt_coeff (xcrt->a, m, n, W);														// Compute a[i] = M_i mod m_i  where M_i = M/m_i
	for ( i = j-1 ; i >= 0 ; i-- ) mpz_clear(W[i]);
	free (W);

	xcrt->delta = ui_len(n)+2;
	for ( Mbits = mbits = 0, i = 0 ; i < n ; i++ ) { b = log2(m[i]);  Mbits += b;  if ( b > mbits ) mbits = b; }
	if ( xcrt->delta + mbits + log2(n) > 63 ) { err_printf ("XCRT moduli too large, max log m_i = %f, n=%d\n", mbits, n); exit (0); }
	// Compute M mod P iteratively, don't bother with a tree (unless P is huge, it wouldn't make much of a difference)
	mpz_init2(M,(int)Mbits+1);  mpz_init2(Y,(int)Mbits+1);  mpz_init2(Z,mpz_sizeinbase(P,2));
	mpz_set_ui(M,m[0]);
	for ( i = 1 ; i < n ; i++) { mpz_mul_ui(M,M,m[i]); mpz_mod(M,M,P); }

	xcrt->n = n;
	xcrt->k = k;
	xcrt->m = (unsigned long *) malloc(n*sizeof(*xcrt->m));
	xcrt->mi = (double *)malloc(n*sizeof(*xcrt->mi));
	xcrt->d = (unsigned long *) malloc (n*n*sizeof(*xcrt->d));
	xcrt->e = (unsigned long *) malloc (n*sizeof(*xcrt->e));
	xcrt->c = (unsigned long *) malloc (n*k*sizeof(*xcrt->e));
	xcrt->t = (unsigned long *) malloc (n*k*sizeof(*xcrt->e));
	xcrt->s = (unsigned long *) malloc (k*sizeof(*xcrt->e));
	mpz_mod(Z,M,P);																	// Z = M mod P

	for ( i = 0 ; i < n ; i++ ) {
		xcrt->m[i] = m[i];
		xcrt->mi[i] = 1.0/((double)m[i]);
		xcrt->a[i] = ui_inverse(xcrt->a[i],m[i]);												// a_i = 1/M_i mod m_i
		mpz_divexact_ui (Y,M,m[i]);  mpz_mod(Y,Y,P);
		for ( j = 0 ; j < n ; j++ ) xcrt->d[j*n+i] = mpz_fdiv_ui(Y,m[j]);								// d_ji = (M_i mod P) mod m_j
		xcrt->e[i] = mpz_fdiv_ui(Z,m[i]);													// e_j = M mod P mod m_j
	}
	mpz_init_set(xcrt->P,P);
	mpz_clear(M); mpz_clear (Y); mpz_clear(Z);
}

void xcrt_clear (xcrt_context_t xcrt)
{
	mpz_clear (xcrt->P);
	free (xcrt->m);  free (xcrt->a);  free (xcrt->d);  free (xcrt->e);  free(xcrt->c);  free (xcrt->t);  free (xcrt->s);
}

void xcrt_reduce (xcrt_context_t xcrt)
{
	register long w;
	register int i, i2, j, k, n;
	
	k = xcrt->k; n = xcrt->n;
	for ( i = 0 ; i < n ; i++ ) for ( j = 0 ; j < xcrt->k ; j++ ) xcrt->t[i*k+j] = ui_mod (xcrt->a[i]*xcrt->c[i*k+j], xcrt->m[i]);
	for ( j = 0 ; j < k ; j++ ) {
		xcrt->s[j] = 0;
		for ( i = 0 ; i < xcrt->n ; i++ ) xcrt->s[j] += (xcrt->t[i*k+j]<<xcrt->delta) * xcrt->m[i];
		xcrt->s[j] = 0.75 + (xcrt->s[j] >> xcrt->delta);
	}
	// now do a naive matrix multiplication: c (nxk) = d (nxn) * t (nxk)
	for ( i = 0 ; i < n ; i++ ) {
		for ( j = 0 ; j < k ; j++ ) {
			xcrt->c[i*k+j] = 0;
			for ( i2 = 0 ; i2 < xcrt->n ; i2++ ) xcrt->c[i*k+j] += xcrt->d[i*n+i2] * xcrt->t[i2*k+j];
		}
	}
	for ( i = 0 ; i < xcrt->n ; i++ ) {
		for ( j = 0 ; j < xcrt->k ; j++ ) {
			w = (long)xcrt->c[i*k+j] - (long) xcrt->e[i] * xcrt->s[j];
			xcrt->c[i*k+j] = i_mod(w,xcrt->m[i]);
		}
	}
}
*/

