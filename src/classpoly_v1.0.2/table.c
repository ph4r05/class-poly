#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include "cstd.h"
#include "table.h"

/*
    Copyright 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

/*
	We allocate three sections of memory:
	
	(1) the hash table - 64 bits per entry
	(2) the lookaside table - 128 bits per entry
	(3) overflow lists - 128 bits per entry
	
	In most cases, only (1) gets used, and this is what we want to keep in cache.   We don't really
	care about the size of the rest.
	
	Currently the number of overflow entries is equal to twice the size of the hash table.  If we exceed this we
	are getting a whole lot of collisions and should increase the table size anway.
*/

void _table_alloc (int bits)
{
	if ( bits > 30 ) { err_printf ("table_alloc: bits=%d too large!\n", bits);  exit (0); }
	htab = malloc(sizeof(*htab)*(1UL<<bits));				// use malloc rather than mem_alloc - we don't need the memory initialized
	if ( ! htab ) { printf ("Unable to allocate %ld bytes of memory.\n", sizeof(*htab)*(1<<bits));  abort(); }
	htab_list = malloc(sizeof(*htab_list)*(1UL<<(bits-1)));	// ditto
	if ( ! htab_list ) { printf ("Unable to allocate %ld bytes of memory.\n", sizeof(*htab_list)*(1UL<<(bits-1)));  abort(); }
	htab_next = 1;
	htab_end = 1UL<<(bits-1);
	htab_bits = bits;
//printf("table alloced %ld bytes\n", (1<<bits)*(sizeof(*htab)+sizeof(*htab_list)/2));
}

void _table_list_extend ()
{
	htab_end = (3*htab_end)/2;
	htab_list = realloc (htab_list, htab_end*sizeof(*htab_list));
	if ( ! htab_list ) { printf ("Unable to allocate %ld bytes of memory.\n", htab_end*sizeof(*htab_list));  abort(); }
}

void table_free (void)
{
	free (htab);
	free (htab_list);
	htab_bits = 0;
}


void table_init (int bits)
{
	if ( bits > htab_bits ) {
		if ( htab_bits ) { /*printf ("Table realloced from %d to %d bits\n", bits);*/ table_free(); }
		_table_alloc(bits);
	}
	htab_mask = (1UL<<bits)-1UL;
	memset (htab, 0, (1UL<<bits)*sizeof(*htab));
	htab_next = 1;
}
