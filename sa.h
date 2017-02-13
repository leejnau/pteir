/* Copyright (c) 2011 Lee Nau
 * Released under MIT License (see LICENSE file) */

#ifndef SA_H
#define SA_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>    // isalpha
#include <stdbool.h>  // bool type
#include "seqs.h"
#include "divsufsort.h"
#include "util.h"  // min

// data structures/types
typedef struct {
	int64_t left, right;  // use integer because left,right can be negative
} lr_pair;

// functions
uint32_t sa_num_occr(sequence seqs[], uint32_t *SAs[], uint32_t numseqs, 
				  char *target);
uint32_t sa_seq_occr(sequence seqs[], uint32_t *SAs[], uint32_t numseqs, 
				  char *target);
uint32_t sa_re_num_occr(sequence seqs[], uint32_t *SAs[], uint32_t numseqs, 
					 char *target);
uint32_t sa_re_seq_occr(sequence seqs[], uint32_t *SAs[], uint32_t numseqs,
					 char *target);
void print_SAs(uint32_t *SAs[], sequence seqs[], const uint32_t numseqs);
void create_SAs(uint32_t *SAs[], sequence seqs[], const uint32_t numseqs); 
void free_SAs(uint32_t *SAs[], const uint32_t numarrays); 

// needed for building offset lists in teir.c
lr_pair lr_subword_search(char *seq, uint32_t *SA, char *word, uint32_t *i);
bool re_streq(const char *s1, const char *s2, const uint32_t length); 

#endif // SA_H
