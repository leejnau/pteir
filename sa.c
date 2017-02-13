/* Copyright (c) 2011 Lee Nau
 * Released under MIT License (see LICENSE file) */

#include "sa.h"

// compare two strings for equality, faster than strcmp
bool streq(char *w1, char *w2, const uint32_t length) {
	for (uint32_t i = 0; i < length; i++) {
		if (w1[i] != w2[i]) 
			return false;
	}
	return true;
}

// compare two strings for regular expression equality
bool re_streq(const char *s1, const char *s2, const uint32_t length) {
	for (uint32_t i = 0; i < length; i++) {
		if (s1[i] == '\0') {   // s1 is always the original sequence
			//printf("s1 is NULL so returning false\n");
			return false;   // always false if end is reached before s1 over
		}
		if (s2[i] == '.') {    // only s2 can have . regex characters in it
			//printf("s2 is dot so continuing\n");
			continue;       // if so skip because no match is necessary
		}
		else if (s1[i] == s2[i]) {   // if equal continue (obviously)
			//printf("s1 == s2, continuing\n");
			continue;
		}
		else if (s1[i] != s2[i]) {  // if not equal, no match
			//printf("s1 != s2, returning false\n");
			return false;
		}
		else 
			printf("ERROR: UNKNOWN CASE! in re_streq\n");
	}	
	return true;
}

// find the left+right pair that marks the SA boundaries for a word
lr_pair lr_search(char *seq, uint32_t *SA, char *word, uint32_t length,
     			  int64_t low, int64_t high) {
	lr_pair lr;

	if (high < low) {
		lr.left = lr.right = -1;  // signifies improper bounds (left > right)
		return lr;
	}

	int64_t mid = low + ((high - low)/2);
	if (mid </*=*/ 0) {  // if mid is negative then there are no occurrences
		lr.left = lr.right = -1;
		return lr;
	}
	if (strncmp(seq+SA[mid], word, length) > 0) {
		return lr_search(seq, SA, word, length, low, mid-1);
	} else if (strncmp(seq+SA[mid], word, length) < 0) {
		return lr_search(seq, SA, word, length, mid+1, high);
	} else {
		uint32_t left = mid, right = mid;
		while (left > low && streq(seq+SA[left-1], word, length)) {left--;}
		while (right < high && streq(seq+SA[right+1], word, length)) {right++;}
		lr.left = left;
		lr.right = right;
		return lr;
	}
}

// the total number of occurrences of a target word within a SINGLE sequence
uint32_t occr(/*char *seq,*/ sequence *s, uint32_t *SA, char *target) {
	//lr_pair lr = lr_search(seq, SA, target, strlen(target), 0, strlen(seq)-1);
	lr_pair lr = lr_search(s->seq, SA, target, strlen(target), 0, s->length-1);
	if (lr.left == -1 || lr.right == -1) 
		return 0;
	else 
		return (lr.right-lr.left+1);
}

// the total number of occurrences of a target word within ALL sequences
uint32_t sa_num_occr(sequence seqs[], uint32_t *SAs[], uint32_t numseqs, 
				  char *target) {
	uint32_t total_num_occr = 0;
	for (uint32_t i = 0; i < numseqs; i++) {
		total_num_occr += occr(&seqs[i], SAs[i], target);
	}
	return total_num_occr;
}

// the number of unique sequences that a target word appears in ("support")
uint32_t sa_seq_occr(sequence seqs[], uint32_t *SAs[],uint32_t numseqs, 
				  char *target) {
	uint32_t total_seq_occr = 0;
	for (uint32_t i = 0; i < numseqs; i++) {
		if (occr(&seqs[i], SAs[i], target) > 0) {  // TODO: optimize me!
			++total_seq_occr;
		}
	}
	return total_seq_occr;
}

// search for just the first subword of real characters (residues) in a regex
lr_pair lr_subword_search(char *seq, uint32_t *SA, char *word, uint32_t *i) {

	lr_pair lr;   // return value
	uint32_t word_index = *i;
	uint32_t word_length = strlen(word);
	int64_t low = 0;                 // lower search bound
	int64_t high = strlen(seq) - 1;  // upper search bound
	// count contiguous alpha chars
	while (isalpha(word[word_index]) && word_index < word_length) {
		word_index++;  // increment index into subword
	}
	lr = lr_search(seq, SA, word, word_index, low, high);
	*i = word_index;
	return lr;
}

// occurrences of a regular expression in a SINGLE sequence
uint32_t re_occr(sequence *s, uint32_t *SA, char *word) {
	
	uint32_t num_occr = 0;
	uint32_t word_index = 0;  // index into the first '.' char of the word
	lr_pair lr = lr_subword_search(s->seq, SA, word, &word_index);
	
	if (lr.left == -1 || lr.right == -1)
		return 0;  // no occurrences
		
	for (uint32_t x = lr.left; x <= lr.right; x++) {
		if (re_streq(s->seq+SA[x]+word_index, word+word_index, 
					strlen(word+word_index))) {
			++num_occr;
		}
	}
	return num_occr;
}

// count the overall number of occurrences of target in a set of sequences
uint32_t sa_re_num_occr(sequence seqs[], uint32_t *SAs[], uint32_t numseqs, 
					 char *target) {
	uint32_t total_re_num_occr = 0;
	for (uint32_t i = 0; i < numseqs; i++) {
		total_re_num_occr += re_occr(&seqs[i], SAs[i], target); 
	}
	return total_re_num_occr;
}

// count the number of sequences that target occurs in
uint32_t sa_re_seq_occr(sequence seqs[], uint32_t *SAs[], uint32_t numseqs, 
					 char *target) {
	uint32_t total_re_seq_occr = 0;
	uint32_t word_index = 0;
	lr_pair lr;
	for (uint32_t i = 0; i < numseqs; i++) {
		lr = lr_subword_search(seqs[i].seq, SAs[i], target, &word_index);
		for (uint32_t x = lr.left; x <= lr.right; x++) { // for each SA entry
			if (re_streq(seqs[i].seq+SAs[i][x]+word_index, target+word_index,
						strlen(target+word_index))) {
				++total_re_seq_occr;
				break;  // only needs to occur in a single entry of this SA
			}
		}
	}
	return total_re_seq_occr;
}

void print_SA(uint32_t *SA, char *seq, uint32_t SA_size) {
	uint32_t i;
	for (i = 0; i < SA_size; i++) {
		printf("SA[%u]\tsuffix [%u]:\t%s\n", i, SA[i], seq+SA[i]);
	}
}

void print_SAs(uint32_t *SAs[], sequence seqs[], const uint32_t numseqs) {
	for (uint32_t i = 0; i < numseqs; i++) {
		printf("----------------------------------------\n");
		printf("Suffix Array %u\n", i);
		printf("----------------------------------------\n");
		print_SA(SAs[i], seqs[i].seq, seqs[i].length);
	}
}

void create_SAs(uint32_t *SAs[], sequence seqs[], const uint32_t numseqs) {
	printf("Creating %u suffix arrays...", numseqs);
	uint32_t i;

	// construct suffix arrays in parallel
    #pragma omp parallel for shared(seqs,SAs) private(i)
	for (i = 0; i < numseqs; i++) {
		SAs[i] = malloc(sizeof(uint32_t) * (seqs[i].length+1));  // +1 for null
		if (divsufsort(seqs[i].seq, SAs[i], seqs[i].length) != 0) {
			printf("Error: Unable to create suffix array for sequence: %u\n",i);
		} 
	}
	#pragma omp barrier
	printf("done creating suffix arrays\n");
}

void free_SAs(uint32_t *SAs[], const uint32_t numarrays) {
	uint32_t i;
	for (i = 0; i < numarrays; i++) {
		free(SAs[i]);
	}
}

