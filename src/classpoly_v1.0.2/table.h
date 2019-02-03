#ifndef _TABLE_INCLUDE_
#define _TABLE_INCLUDE_

#include <stdlib.h>
#include <stdint.h>

/*
    Copyright 2007-2014 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#ifdef __cplusplus
extern "C" {
#endif

/*
	Simple table lookup/insert code.  This is designed for small tables
        that fit in cache memory, using 64 bits per entry (128 bits per overflow entry).
	This requires that the caller to either store group entries in a seperate list
	or reconstruct them as required.  The table is also optimized for unsuccessful
	searches, which are expected to be typical (e.g. during a BSGS search).
	
	Each table entry is a pair of 32-bit values, a datum and an index.  The key is
	a 32-bit value (which may be zero), typically coming from a hash function.
	It need not uniquely identity the table entry, but the hope is that key collisions
	are rare (much rarer than table collisions, which may occur fairly often), since
	a key collision will then require the caller to check the corresponding group
	entries.
	
	The datum value can be any NON-ZERO 32-bit value.  The datum values
	need not be unique, but typically they will be.  The lookup function returns a
	list of all datum values in the table whose key matches a specified key.
	
	When used in a BSGS search, datum will be typically be the baby step exponent,
	which will be unique, the the index will be a hash of the group element.
	
	Note that the table size used in a particular search may be any power of 2 smaller
	up to the allocated table size - this is desirable as it reduces the cose of initialization
	and improves locality when this is  small without creating too many
	collisions.  A load factor of around 0.5 works well.
*/

#define TABLE_MAX_MATCHES	512

struct htab_list_item {
	uint32_t key;
	uint32_t datum;
	uint32_t next;
} *htab_list;

uint32_t *htab;

uint32_t htab_next, htab_end;
int htab_bits;								// allocated table size (max)
uint32_t htab_mask;						// mask for initialized table size

void _table_alloc (int bits);
static inline void table_alloc (int bits)
	{ if ( bits <= htab_bits ) return;  _table_alloc(bits); }
void _table_list_extend(void);
void table_init (int bits);						// automatically allocs if necessary
void table_free (void);

void _table_list_insert (struct htab_list_item *pFirst, uint32_t dataum, uint32_t key);
int _table_list_matches (struct htab_list_item *pNext, uint32_t *data, uint32_t key);

// IMPORTANT: dataum must be non-zero and index must be less than 2^n where n is the 
// bits paramater specified to table_init.  These constraints are not verified.
static inline uint32_t table_insert (uint32_t datum, uint32_t key)
{
	register struct htab_list_item *htab_ptr;
	register uint32_t index;

// printf("insert %d, %d\n", key, datum);	
	index = key&htab_mask;
	if ( htab_next == htab_end ) _table_list_extend();
	htab_ptr = htab_list+htab_next;
	htab_ptr->datum = datum;
	htab_ptr->key = key;
	htab_ptr->next = htab[index];
	htab[index] = htab_next++;
	return htab_ptr->next;
}

static inline int _table_get_matches (uint32_t data[TABLE_MAX_MATCHES], uint32_t key, uint32_t next)
{
	register int i;
	
	i = 0;
	do {
		if ( htab_list[next].key == key ) {
			if ( i >= TABLE_MAX_MATCHES ) { printf ("Exceeded MAX_MATCHES=%d for key=%u\n", i, key); abort(); }
			data[i++] = htab_list[next].datum;
		}
		next = htab_list[next].next;
	} while ( next );
	return i;	

}

static inline int table_lookup (uint32_t data[TABLE_MAX_MATCHES], uint32_t key)
{
	register uint32_t index, next;

//printf("lookup %d\n", key);	
	index = key&htab_mask;
	next = htab[index];
	if ( ! next ) return 0;
	return _table_get_matches (data, key, next);
}


// does a combined insert/lookup, returning a list of data for any previously inserted entries with the same key value
static inline int table_insert_matches (uint32_t data[TABLE_MAX_MATCHES], uint32_t datum, uint32_t key)
{
	register uint32_t next;

//printf("insert/lookup %d, %d\n", key, datum);	
	next = table_insert(datum, key);
	if ( ! next ) return 0;
	return _table_get_matches (data, key, next);
}

#ifdef __cplusplus
}
#endif

#endif
