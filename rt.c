/* Copyright (c) 2011 Lee Nau
 * Released under MIT License (see LICENSE file) */

#include <stdlib.h>
#include <stdio.h>
#include "rt.h"

// set a specific sequence bit in the sequence vector to true
void seqvec_set(seqvec_t *vec, uint32_t seq) {
	(*vec) |= (1 << seq);  // set bit "seq"
}

// print the binary contents of a sequence vector (hex for now)
void seqvec_print(const seqvec_t vec) {
	printf("%x\n", vec);
}

// return an integer that is the total number of bits set for this vector
uint32_t rt_seqvec_countseqs(const seqvec_t vec) {
	uint32_t bits = 0;  // total number of bits (seqs) set, to be returned
	for (uint32_t i = 0; i < sizeof(seqvec_t)*8; i++) {  // for each bit
		if (vec & (1 << i))
			bits++;
	}
	return bits;
}

// return the proper radix trie branch index corresponding to a given character
int rt_char2idx(const char c) {
	switch (c) {
		case 'A':  return A_IDX;
		case 'C':  return C_IDX;
		case 'G':  return G_IDX;
		case 'T':  return T_IDX;
		case 'N':  return N_IDX;    // N means any character (nucleotide)
		case '.':  return DOT_IDX;  // any character, but different from N_IDX
		default:   printf("ERROR: Nonstd char in char2idx, aborting\n");exit(1);
	}
}

// report memory usage of the radix trie
void rt_report_size() {
	printf("Radix Trie Occupies: %lu bytes (%f KB) [%f MB]\n", 
			RT_SIZE_BYTES, RT_SIZE_BYTES/1024.0, RT_SIZE_BYTES/(1024.0*1024.0));
}

rt_node *new_rt_node() 
{
    rt_node *new_node;
    int i;
	RT_SIZE_BYTES += sizeof(rt_node);  // memory counter
    new_node = (rt_node *) malloc (sizeof (rt_node));
    new_node->data.count = 0;
	new_node->data.seqvec = 0;
    for (i = 0; i < RT_NUM_BRANCHES; i++)
		new_node->branch[i] = NULL;
    return new_node;
}

// build a radix trie using multiple sequences, but only increase total count
rt_node* rt_build(sequence seqs[], const uint32_t numseqs, 
			  const uint32_t min_word_length, const uint32_t max_word_length) {
	rt_node *root = new_rt_node();
	char next_word[max_word_length+1];
	uint32_t i, j, s, limit, count = 0, tid;

	printf("Building Radix Trie to depth: %d...\n", RT_MAX_DEPTH);

	///#pragma omp parallel for shared(root,seqs) private(i,tid)
	for (s = 0; s < numseqs; s++) {
		tid = omp_get_thread_num();  // get current thread id
		printf("Adding sequence: %u to trie in thread: %u\n", s, tid);
		for (i = 0; i < seqs[s].length; i++) {
			limit = min(max_word_length, seqs[s].length - i + 1);
			for (j = min_word_length; j </*=*/ limit; j++) {
				memset(next_word, 0, max_word_length+1);  // reset next word
				strncpy(next_word, &(seqs[s].seq[i]), j);  // read in next word
				///#pragma omp critical
				{
					rt_data data = rt_get(root, next_word, s);  // word exists
					//rt_data data = rt_search(root, next_word);
					if (data.count == 0) {  // add new node
						rt_data add_data;  // temporary spot for new node data
						add_data.count = 1;
						add_data.seqvec = 0;
						seqvec_set(&(add_data.seqvec), s);
						rt_add(root, next_word, add_data);
					} /*else {  // data already exists (update stats)
						rt_data tmp = rt_search(root, next_word);
						rt_inc(root, next_word, s);	  // increment stats
					}*/
				} // end pragma omp critical
			}
		}
	}
	///#pragma omp barrier
	
	return root;
}

/* add string s to the trie and attach data to it */
void rt_add(rt_node *node, char *s, rt_data data) {
    if (!*s) {
		node->data.count = data.count;
		node->data.seqvec = data.seqvec;
	} else {
		int idx = rt_char2idx(*s);
		if (!node->branch[idx])
	    	node->branch[idx] = new_rt_node();
		rt_add(node->branch[idx], s+1, data);
    }
}

/* return data associated with a node, and increment counts if appropriate */
rt_data rt_get(rt_node *node, char *s, uint32_t seq_num) {
	rt_data get_data;
    if (node == NULL) {
		get_data.count = 0;  // no occurrences
		get_data.seqvec = 0;  // appears in no sequences	
		return get_data;
	}
	if (*s == '\0') {
		node->data.count++;  // increase the statistics for this sequence
		seqvec_set(&(node->data.seqvec), seq_num);
		return node->data;
    }
    int idx = rt_char2idx(*s);
	return rt_get(node->branch[idx], s+1, seq_num);
}

/* increment data associated with a node */
void rt_inc(rt_node *node, char *s, uint32_t seq_num) {
	if (node == NULL) {
		return;
	}
	if (*s == '\0') {
		node->data.count++;
		seqvec_set(&(node->data.seqvec), seq_num);
	}
	int idx = rt_char2idx(*s);
	return rt_inc(node->branch[idx], s+1, seq_num);
}

/* this does NOT increase node stats - returns number of occurrences */
rt_data rt_search(rt_node *node, char *s) {

	if (strlen(s) > RT_MAX_DEPTH) {
		printf("ERROR: trying to look for pattern beyond RT Max Depth\n");
		exit(1);
	}
	rt_data search_data;
	if (node == NULL) {
		search_data.count = 0;
		search_data.seqvec = 0;
		return search_data;
	}
	if (*s == '\0') {
		return node->data;
	}
	int idx = rt_char2idx(*s);
	return rt_search(node->branch[idx], s+1);
}

/* regular expression search for a pattern s - does NOT increase node stats */
rt_data rt_re_search(rt_node *node, char *s) {
	if (strlen(s) > RT_MAX_DEPTH) {
		printf("ERROR: trying to look for pattern beyond RT Max Depth\n");
		exit(1);
	}
	rt_data search_data, temp_data;
	search_data.count = 0;
	search_data.seqvec = 0;

	if (node == NULL) {
		return search_data;  // return zeros
	}

	if (*s == '\0') {
		return node->data;
	}
	
	// recurse on dots (check all branches)
	int idx = rt_char2idx(*s);
	// DOT_IDX indicates '.' (any character), so search ALL branches
	if (idx == DOT_IDX) {
		for (int i = 0; i < RT_NUM_BRANCHES; i++) {
			temp_data = rt_re_search(node->branch[i], s+1);
			search_data.count  += temp_data.count;   // accumulate brnch cnts
			search_data.seqvec |= temp_data.seqvec;  // OR sequence vectors
		}
		return search_data;

	} else {  // otherwise it's a literal char, so search just that branch
		return rt_re_search(node->branch[idx], s+1);
	}
}

