/* Copyright (c) 2011 Lee Nau
 * Released under MIT License (see LICENSE file) */

#ifndef RADIXTREE_H
#define RADIXTREE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include "seqs.h"
#include "util.h"

#define RT_NUM_BRANCHES 5  // {A,C,G,T,N}
#define A_IDX 0
#define C_IDX 1
#define G_IDX 2
#define T_IDX 3
#define N_IDX 4
#define DOT_IDX 5
#define RT_MAX_DEPTH 20  // maximum depth the radix trie may have

static uint64_t RT_SIZE_BYTES = 0;  // tally of bytes used for this RT

typedef uint32_t seqvec_t;

// data stored in each trie node
typedef struct {
	uint32_t count;  // total number of occurrences of a word
	seqvec_t seqvec;  // sequence vector for sequence occrs
} rt_data;

// trie node itself
typedef struct rt_node {
    rt_data data; 
	struct rt_node *branch[RT_NUM_BRANCHES];
} rt_node;

rt_node*	 new_rt_node();
void         rt_add(rt_node *node, char *s, rt_data data);
rt_data      rt_get(rt_node *node, char *s, uint32_t seq_num);
void		 rt_inc(rt_node *node, char *s, uint32_t seq_num); 
rt_node*	 rt_build(sequence seqs[], const uint32_t numseqs,
				const uint32_t min_word_length, const uint32_t max_word_length);
rt_data		 rt_search(rt_node *node, char *s);     // normal word search
rt_data	     rt_re_search(rt_node *node, char *s);  // regex search
void		 rt_report_size();  // report stats on memory usage for RT
void		 seqvec_set(seqvec_t *vec, uint32_t seq);  // set a single bit
void         seqvec_print(const seqvec_t vec);  // print hex values for seqvec
uint32_t     rt_seqvec_countseqs(const seqvec_t vec);  // number of bits set 

#endif // RADIXTREE_H
