/* Copyright (c) 2011 Lee Nau
 * Released under MIT License (see LICENSE file) */

#ifndef TEIR_H
#define TEIR_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>
#include "seqs.h"
#include "sa.h"
#include "rt.h"
#include "util.h"

// defines
#define MAX_EP_LENGTH 128             // an EP can only contain 64 characters
#define MAX_EP_SET_SIZE (1024*1024)   // allow a megabyte (~1M) EPs
#define MAX_MP_LENGTH 128             // a  MP can only contain 128 characters
#define MAX_MP_SET_SIZE (1024*128)    // allowed total number of MPs
#define MAX_OFFSET_LIST_SIZE 1024     // allow offset lists to have size 1024
#define MAX_STACK_SIZE 1024           // allow stacks to have up to size 1024
#define MAX_OMP_THREADS 32            // max number of allowed openmp threads

// type of data structure to use for indexing sequences and lookups
//#define SA  // use the suffix array data structure
//#define RT  // use the radix trie data structure
#define HYBRID
//#define LINEAR

// which version of the algorithm to use
#define NOVERSION 0
#define SEQVERSION 1
#define OCCVERSION 2

// data structures/types
// elementary/maximal pattern type, including the string, and occurrence counts
typedef struct {
	char s[MAX_EP_LENGTH];
	uint32_t num_occrs;  // total number of occurrences
	uint32_t seq_occrs;  // unique sequence occurrences ("support")
} pattern; 

// global teiresias variables
static pattern EPset_pfless[MAX_EP_SET_SIZE];  // the elementary pattern set
static pattern EPset_sfless[MAX_EP_SET_SIZE];  // EP set, sorted sf-wise less
static uint32_t EPmap[MAX_EP_SET_SIZE];  // map EPset_pfless -> EPset_sfless
static uint32_t EPsetsize = 0;  // (current) size of EP set(s)
static pattern stack[MAX_OMP_THREADS][MAX_STACK_SIZE];  // MP stack(s)
static int64_t stack_top_index[MAX_OMP_THREADS];  // current top of stack
static pattern MPset[MAX_MP_SET_SIZE];  // verified maximal patterns
static uint32_t MPsetsize = 0;  // maximal pattern set size 
static const char alphabet[] = "ACGT";  // DNA alphabet
static const uint32_t ALPHABET_SIZE = 4; 
static uint32_t num_equiv_sets = 1;  // number of unique equivalence sets in EP
static uint32_t equiv_set_first[MAX_EP_SET_SIZE];  // starts of equiv sets in EP
static uint32_t equiv_set_last[MAX_EP_SET_SIZE];   // ends of equiv sets in EP
static uint32_t equiv_set_sizes[MAX_EP_SET_SIZE];  // size of equiv sets in EP

// teiresias core functions
void scan(rt_node *rt, sequence seqs[], uint32_t *SAs[], const uint32_t numseqs,
		const uint32_t L, const uint32_t W, const uint32_t K, uint8_t version); 
void extend(rt_node *rt, sequence seqs[], uint32_t *SAs[], 
		const uint32_t numseqs, const uint32_t L, const uint32_t W, 
		const uint32_t K, pattern p, uint8_t version);
void convolve(rt_node *rt, sequence seqs[], uint32_t *SAs[], 
		const uint32_t numseqs, const uint32_t L, const uint32_t W, 
		const uint32_t K, uint8_t version);

// EP set functions
void EPsetadd(const pattern EP);   // add an element to the EP set
void EPsetprint();    // print all elements in the EP set
void EPsetsizeprint();  // print size of EP set
void EPsortpfless();  // prefix-wise less sort for EP set
void EPsortsfless();  // suffix-wise less sort for EP set
void EPmapsets();     // map prefix-wise less set to suffix-wise less set
void EPmapprint();    // print mapping from prefix-wise less to suffix-wise
void EPequiv();       // find and print information about equivalence sets

// linear search
uint32_t linear_num_occr_exact(const sequence seqs[], const uint32_t numseqs,
		        const char *needle);

// stack functions
void stack_init();  // initialize stacks
void stack_push(int tid, pattern new_top);
void stack_pop(int tid);
pattern stack_top(int tid);
uint32_t stack_size(int tid);

// convolution operations
void conv_op(const char *a, const char *b, const uint32_t L, char *conv);

// MP set functions 
void MPsetadd(sequence seqs[], uint32_t *SAs[], const uint32_t numseqs,
		const pattern newMP, uint8_t version);
bool is_maximal(const pattern P);
void writeMPs(const char *filename);

#endif  // TEIR_H
