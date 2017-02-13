/* Copyright (c) 2011 Lee Nau
 * Released under MIT License (see LICENSE file) */

#include "teir.h"

// return the occurrences of a needle n in all sequences up to threshold K
uint32_t linear_num_occr_threshold(const sequence seqs[], 
		const uint32_t numseqs, const char *needle, const uint32_t K) {
	uint32_t num_occrs = 0, n;
	// for each sequence
	for (uint32_t curseq = 0; curseq < numseqs; curseq++) {
		// run a sliding window over the haystack sequence
		for (uint32_t s = 0; s < seqs[curseq].length-strlen(needle) + 1; s++) {
			for (n = 0; n < strlen(needle); n++) {
				if (needle[n] == '.')
					continue;
				else if (needle[n] != seqs[curseq].seq[s+n]) 
					break;
			}
			if (n == strlen(needle)) {  // only inc num_occr if n has been inc
				++num_occrs;
				if (num_occrs >= K) {  // check against threshold K
					return num_occrs;  // return if sufficient support is found
				}
			}
		}
	}
	return num_occrs;  // if insufficient support is there return occrs anyway
}

// return the occurrences of a needle in all sequences (the exact number)
uint32_t linear_num_occr_exact(const sequence seqs[], const uint32_t numseqs,
		const char *needle) {
	uint32_t num_occrs = 0, n;
	// for each sequence
	for (uint32_t curseq = 0; curseq < numseqs; curseq++) {
		// run a sliding window over the haystack sequence
		for (uint32_t s = 0; s < seqs[curseq].length-strlen(needle) + 1; s++) {
			for (n = 0; n < strlen(needle); n++) {
				if (needle[n] == '.')
					continue;
				else if (needle[n] != seqs[curseq].seq[s+n]) 
					break;
			}
			if (n == strlen(needle)) {  // only inc num_occr if n has been inc
				++num_occrs;
			}
		}
	}
	return num_occrs;
}

// count the number of regular characters (not dots) in the string p
uint32_t num_reg_chars(const char *p) {
	uint32_t num_chars = 0, c;
	for (c = 0; c < strlen(p); c++) {
		if (isalpha(p[c])) {
			++num_chars;
		}
	}
	return num_chars;
}

// create and return a string of dots with num_dots dots
char* dots(uint32_t num_dots) {
	char *dotstring = malloc(num_dots+1);  // +1 for null
	for (uint32_t i = 0; i < num_dots; i++) {
		dotstring[i] = '.';
	}
	return dotstring;
}

// run the scan stage of the TEIRESIAS algorithm
void scan(rt_node *rt, sequence seqs[], uint32_t *SAs[], const uint32_t numseqs,
		  const uint32_t L, const uint32_t W, const uint32_t K, uint8_t version)
{

	printf("Starting Teiresias Scan...\n");
	uint32_t i, tid;  // index for alphabet and thread id
	pattern p;  // pattern for consideration (reused)
	#pragma omp parallel for shared(rt,seqs,SAs) private(i,p,tid)
	for (i = 0; i < strlen(alphabet); i++) {  // consider one char at a time
		tid = omp_get_thread_num();
		printf("Scan for character [%c] from thread [%u]\n",
				alphabet[i], tid);
		
		// check for support for a given single alphabet character
		p.s[0] = alphabet[i];  p.s[1] = '\0';
#if defined(SA)
		p.num_occrs = sa_re_num_occr(seqs, SAs, numseqs, p.s);
		if (version == SEQVERSION) {
			p.seq_occrs = sa_re_seq_occr(seqs, SAs, numseqs, p.s);
		}
#elif defined(RT)
		rt_data search_data = rt_re_search(rt, p.s);
		p.num_occrs = search_data.count;
		if (version == SEQVERSION) {
			p.seq_occrs = rt_seqvec_countseqs(search_data.seqvec);
		}
#elif defined(HYBRID)
		if (strlen(p.s) < RT_MAX_DEPTH) {
			rt_data search_data = rt_re_search(rt, p.s);
			p.num_occrs = search_data.count;
			if (version == SEQVERSION) {
				p.seq_occrs = rt_seqvec_countseqs(search_data.seqvec);
			}
		} else {
			p.num_occrs = sa_re_num_occr(seqs, SAs, numseqs, p.s);
			if (version == SEQVERSION) {
				p.seq_occrs = sa_re_seq_occr(seqs, SAs, numseqs, p.s);
			}
		}
#elif defined(LINEAR)
		// for either OCCVERSION or SEQVERSION, num_occrs is needed for conv
		p.num_occrs = linear_num_occr_threshold(seqs, numseqs, p.s, K);
		if (version == SEQVERSION)
			p.seq_occrs = linear_seq_occr_threshold(seqs, numseqs, p.s, K);
#endif // defined(SA)

		if (version == SEQVERSION) {
			if (p.seq_occrs >= K)  {  // proper "seqs" support
				extend(rt, seqs, SAs, numseqs, L, W, K, p, version);
			}
		} else if (version == OCCVERSION) {
			if (p.num_occrs >= K) {   // alternate "occrs" support
				extend(rt, seqs, SAs, numseqs, L, W, K, p, version);
			}
		}
	}
}

// extend portion of the scan phase - extend candidate elementary patterns
void extend(rt_node *rt, sequence seqs[],uint32_t *SAs[],const uint32_t numseqs,
		    const uint32_t L, const uint32_t W, const uint32_t K, pattern p,
			uint8_t version) {

	pattern pprime[ALPHABET_SIZE];  // array of candidate extension patterns
	uint32_t pprime_lastchar = 0;
	char dots[MAX_EP_LENGTH];
	uint32_t d;
	const uint32_t A = num_reg_chars(p.s);  // A from teiresias algorithm

	// check for termination criteria
	if (A == L) {   // reached correct number of residues
		// calculate any necessary lists, exact num_occrs, etc. here 
		// can put off some of these calculations until now
#if defined(LINEAR)
		// up until now only the threshold K was considered - now find exact
		p.num_occrs = linear_num_occr_exact(seqs, numseqs, p.s);
#endif // defined(LINEAR)

		EPsetadd(p);  // found a new elementary pattern to be added to set
		return;
	}

	for (uint32_t i = 0; i <= W - strlen(p.s) - L + A; i++) {
		for (d = 0; d < i; d++) 
			dots[d] = '.';
		dots[d] = '\0';
		
		pprime_lastchar = strlen(p.s) + strlen(dots);  // last char position
		for (uint32_t a = 0; a < strlen(alphabet); a++) {
			strcpy(pprime[a].s, p.s);   // copy in pattern to be extended
			strcat(pprime[a].s, dots);  // concatenate dots to the end of p'
			pprime[a].s[pprime_lastchar] = alphabet[a];  // append final char
			pprime[a].s[pprime_lastchar + 1] = '\0';  // always null terminate
#if defined(SA)
			pprime[a].num_occrs = sa_re_num_occr(seqs,SAs,numseqs,pprime[a].s);
			if (version == SEQVERSION)
				pprime[a].seq_occrs = 
					sa_re_seq_occr(seqs,SAs,numseqs,pprime[a].s);
#elif defined(RT)
			rt_data search_data;
			search_data = rt_re_search(rt, pprime[a].s);
			pprime[a].num_occrs = search_data.count; 
			if (version == SEQVERSION) {
				pprime[a].seq_occrs = rt_seqvec_countseqs(search_data.seqvec);
			}
#elif defined(HYBRID)
			if (strlen(pprime[a].s) < RT_MAX_DEPTH) {
				rt_data search_data;
				search_data = rt_re_search(rt, pprime[a].s);
				pprime[a].num_occrs = search_data.count;
				if (version == SEQVERSION) {
					pprime[a].seq_occrs = \
						rt_seqvec_countseqs(search_data.seqvec);
				}
			} else {
				pprime[a].num_occrs = \
					sa_re_num_occr(seqs, SAs, numseqs, pprime[a].s);
				if (version == SEQVERSION) {
					pprime[a].seq_occrs = \
						sa_re_seq_occr(seqs, SAs, numseqs, pprime[a].s);
				}
			}
#elif defined(LINEAR)
			pprime[a].num_occrs = \
				linear_num_occr_threshold(seqs, numseqs, pprime[a].s, K);
#endif  // defined(SA)

		}
		// extend if support is sufficient
		for (uint32_t a = 0; a < strlen(alphabet); a++) {
			if (version == SEQVERSION) {
				if (pprime[a].seq_occrs >= K) {
					extend(rt, seqs, SAs, numseqs, L, W, K, pprime[a], version);
				}
			} else if (version == OCCVERSION) {
				if (pprime[a].num_occrs >= K) {
					extend(rt, seqs, SAs, numseqs, L, W, K, pprime[a], version);
				}
			}
		}
	}
}

// returns the LENGTH of prefix of a string, given parameter L (L-1 characters)
const uint32_t prefix(const char *s, const uint32_t L) {
    uint32_t num_chars = 0;
    uint32_t length = 0;
    while (num_chars != L-1) {
        if (isalpha(s[length]))
            ++num_chars;
        ++length;
    }
	return length;
}

// returns the INDEX of suffix of a string, given parameter L (L-1 characters)
const uint32_t suffix(const char *s, const uint32_t L) {
    uint32_t num_chars = 0;
    uint32_t length = strlen(s)-1;  // start at end of the string
    while (num_chars != L-1) {
        if (isalpha(s[length]))
            ++num_chars;
        --length;
    }
	return length+1;  // +1 because using 0-base indexing for c strings
}

// initialize stacks (separate stack for each OpenMP thread)
void stack_init() {
	for (int tid = 0; tid < MAX_OMP_THREADS; tid++) {
		stack_top_index[tid] = -1;  // init to -1 (means an empty stack)
	}
}

// push an element onto the stack for thread tid
void stack_push(int tid, pattern new_top) {
	++stack_top_index[tid];
	if (stack_top_index[tid] >= MAX_STACK_SIZE) {
		printf("Error: Stack size exceeded!\n");
		exit(1);
	}
	stack[tid][stack_top_index[tid]] = new_top;
}

// pop an element off of the top of the stack for thread tid
void stack_pop(int tid) {
	if (stack_top_index[tid] >= 0) {  // can make top=-1, meaning stack is empty
		--stack_top_index[tid];
	} else {
		printf("Error: Trying to pop empty stack in tid: %d\n", tid);
		exit(1);
	}
}

// return the top of the stack for thread tid
pattern stack_top(int tid) {
	if (stack_top_index[tid] >= 0) {
		return stack[tid][stack_top_index[tid]];
	} else {
		printf("Error: trying to get top of empty stack in tid: %d\n", tid);
		exit(1);
	}
}

// return stack size for thread tid
uint32_t stack_size(int tid) {
	// stack_top = 0 means there is still one element
	return stack_top_index[tid]+1; 
}

// actual convolution of two strings (a,b) with given parameter L: conv 
void conv_op(const char *a, const char *b, const uint32_t L, char *conv) {
	uint32_t strlen_a = strlen(a);  // want the string length of a
	uint32_t prefix_length_b = prefix(b, L);  // want first part of b

	if (strlen_a + strlen(b+prefix_length_b) >= MAX_MP_LENGTH) {
		printf("Error: convolved string is greater than max possible length\n");
		exit(1);
	}

	strcpy(conv, a);  // first copy in the entirety of a, which is the "prefix"
	conv[strlen_a] = '\0';  // strncpy doesn't necessarily null-terminate
	strcat(conv, b+prefix_length_b);  // only copy in last part of b "suffix"
}

// determine if two strings are convolvable 
bool are_convolvable(const char *a, const char *b, const uint32_t L) {
	uint32_t prefix_length = prefix(b, L);  // want first part of b
	uint32_t suffix_index  = suffix(a, L);  // want last part of a
	
	// if bytes are equal (memcmp == 0), they are convolvable
	if (memcmp(a+suffix_index, b, prefix_length) == 0)
		return true;
	else 
		return false;
}

// calculate occurrences of R and T during left and right extensions
void calc_occrs(rt_node *rt, sequence seqs[], uint32_t *SAs[], 
		const uint32_t numseqs, pattern *R, pattern *T) {
#if defined(SA)
	// calculate occurrences for R and T
	R->num_occrs = sa_re_num_occr(seqs, SAs, numseqs, R->s);
	R->seq_occrs = sa_re_seq_occr(seqs, SAs, numseqs, R->s);
	T->num_occrs = sa_re_num_occr(seqs, SAs, numseqs, T->s);
#elif defined(RT)
	rt_data search_data;
	search_data = rt_re_search(rt, R->s);
	R->num_occrs = search_data.count;
	R->seq_occrs = rt_seqvec_countseqs(search_data.seqvec);
	search_data = rt_re_search(rt, T->s);
	T->num_occrs = search_data.count;
#elif defined(HYBRID)
	if (strlen(R->s) < RT_MAX_DEPTH) {
		rt_data search_data = rt_re_search(rt, R->s);
		R->num_occrs = search_data.count;
		R->seq_occrs = rt_seqvec_countseqs(search_data.seqvec);
	} else {
		R->num_occrs = sa_re_num_occr(seqs, SAs, numseqs, R->s);
		R->seq_occrs = sa_re_seq_occr(seqs, SAs, numseqs, R->s);
	}
	if (strlen(T->s) < RT_MAX_DEPTH) {
		rt_data search_data = rt_re_search(rt, T->s);
		T->num_occrs = search_data.count;
	} else {
		T->num_occrs = sa_re_num_occr(seqs, SAs, numseqs, T->s);
	}
#elif defined(LINEAR)
	R->num_occrs = linear_num_occr_exact(seqs, numseqs, R->s);
	T->num_occrs = linear_num_occr_exact(seqs, numseqs, T->s);
#endif // defined(SA)
}

// the convolution stage of the teiresias algorithm
void convolve(rt_node *rt, sequence seqs[], uint32_t *SAs[], 
		const uint32_t numseqs, 
		const uint32_t L, const uint32_t W, const uint32_t K, uint8_t version) {
	printf("Starting Teiresias Convolution...\n");
	EPsortpfless();  // prefix-wise less sort
	EPsortsfless();  // suffix-wise less sort
	EPmapsets();     // map prefix-wise less set to suffix-wise less set
	EPsetsizeprint();  // print number of EPs
	EPequiv();       // calculate and print equivalence classes for EP
	pattern R;  // convolved pattern that is a candidate for stack push
	pattern T;  // current top of stack
	uint32_t pfless_begin = 0;  // the first element of the pfless set
	char R_temp[MAX_MP_LENGTH];  // temporary string pointer for pattern R
	uint32_t P;  // index for pattern P currently being processed in below loop
	int tid;   // thread id for current OpenMP thread
	stack_init(); // set up the stacks (one per thread)
	for (uint32_t cur_set = 0; cur_set < num_equiv_sets; cur_set++) {
		printf("Processing equivalent set number %u\n", cur_set);
		
		#pragma omp parallel 
		{
			// process each equivalent set simultaneously
			#pragma omp for private(R,T,R_temp,P,tid) schedule(dynamic)
			for (P = equiv_set_first[cur_set]; P <= equiv_set_last[cur_set]; 
					P++) {
				tid = omp_get_thread_num();  // unique tid per thread
				/*printf("Thread: [%d] is processing P[%u]:  %s\n", 
						tid, P, EPset_pfless[P].s);*/

				// push each prefix-wise less EP element onto stack
				stack_push(tid, EPset_pfless[P]);
#ifdef DEBUG
				printf("Candidate EP [%u]: %s\n", P, EPset_pfless[P].s);
				if (P % 10 == 0) {  // report progress every 100th pattern
					printf("Progress: %f%% done\n",
							((float)P/(float)EPsetsize)*100.0);
				}
#endif
	
				while (stack_size(tid) != 0) {  // while the stack is not empty
start:			// goto label
					if (stack_size(tid) == 0) {  // not in pseudocode
						break;  // if stack is already empty, break out of loop
					}
					T = stack_top(tid);  // assign top of stack for candidate
			
					// LEFT EXTEND
					for (uint32_t pfless = pfless_begin; 
					///for (uint32_t pfless = 0; 
							pfless < EPsetsize; pfless++) {
						if (are_convolvable(EPset_pfless[pfless].s, T.s, L)) {
							conv_op(EPset_pfless[pfless].s, T.s, L, R.s);

							calc_occrs(rt, seqs, SAs, numseqs, &R, &T); 

							if (R.num_occrs == T.num_occrs) { // T non-maximal
								stack_pop(tid);
							}
							
							if (version == SEQVERSION && 
									R.seq_occrs >= K && is_maximal(R)) {
								stack_push(tid, R);
								goto start;
							}
							if (version == OCCVERSION && 
								     R.num_occrs >= K && is_maximal(R)) {
								stack_push(tid, R);
								goto start;
							}

							// new - not in pseudocode
							if (R.num_occrs == T.num_occrs) {  
								goto start;
							}
						} // end if are_convolvable
					}  // END LEFT EXTEND

					// RIGHT EXTEND
					for (uint32_t sfless = 0; sfless < EPsetsize; sfless++) {
						// skip those empty sfless elements
						if (EPset_sfless[sfless].s[0] == '\0')
							continue;
						if (are_convolvable(T.s, EPset_sfless[sfless].s, L)) {
							conv_op(T.s, EPset_sfless[sfless].s, L, R.s);
						
							calc_occrs(rt, seqs, SAs, numseqs, &R, &T); 

							if (R.num_occrs == T.num_occrs) { // T non-maximal
								stack_pop(tid);
							}

							if (version == SEQVERSION && 
									R.seq_occrs >= K && is_maximal(R)) {
								stack_push(tid, R);
								goto start;
							}
							if (version == OCCVERSION && 
									R.num_occrs >= K && is_maximal(R)) {
								stack_push(tid, R);
								goto start;
							}

							// new - not in pseudocode
							if (R.num_occrs == T.num_occrs) {  
								goto start;
							}
						} // end if are_convolvable
					}  // END RIGHT EXTEND

					stack_pop(tid);
			
					if (is_maximal(T)) {
#ifdef DEBUG
						// print out maximal patterns
						printf("%s: occrs: %u seqs: %u\n", 
								T.s, T.num_occrs, T.seq_occrs);
#endif
						MPsetadd(seqs, SAs, numseqs, T,version);// add to MP set
					}
				}  // end while stack is not empty
			}  // end for each candidate EP
		} // end omp parallel section
		#pragma omp barrier  // finish each set of EPs before moving on

		// remove elements from sfless set that no longer need to be processed	
		for (pfless_begin =  equiv_set_first[cur_set]; 
			 pfless_begin <= equiv_set_last[cur_set]; pfless_begin++) {
			EPset_sfless[EPmap[pfless_begin]].s[0] = '\0';
		}
	}  // end for each equivalent set 
	printf("Total number of MPs: %u\n", MPsetsize);
}

// add an element to the Elementary Pattern set (includes critical section)
void EPsetadd(const pattern EP) {
	#pragma omp critical  // only allow one EP to added at a time 
	{
		EPset_pfless[EPsetsize] = EP; // add to the prefix-wise less EP set 
		EPset_sfless[EPsetsize] = EP; // add to the suffix-wise less EP set 
		EPsetsize++;
		if (EPsetsize >= MAX_EP_SET_SIZE) {
			printf("Error: exceeded EP set maximum size!\n");
			exit(1);
		}
	}
}

// print the Elementary Pattern (EP) set(s) - first pfless then sfless
void EPsetprint() {
	printf("%-16s\t%16s\n", "pf-less", "sf-less"); // header for ordered EP sets
	for (uint32_t i = 0; i < EPsetsize; i++) {
		printf("%-16s\t", EPset_pfless[i].s);
		printf("%16s\n", EPset_sfless[i].s);
	}
}

// indicate how many Elementary Patterns were found overall
void EPsetsizeprint() {
	printf("Total Elementary Patterns Found: %u\n", EPsetsize);
}

// add an element to the Maximal Pattern set (includes critical section)
void MPsetadd(sequence seqs[], uint32_t *SAs[], const uint32_t numseqs,
		const pattern P, uint8_t version) {
	
	pattern newMP;  // the new maximal pattern to be added
	strcpy(newMP.s, P.s);   // copy in pattern string
	newMP.num_occrs = P.num_occrs;  // overall occrs
	if (version == SEQVERSION)
		newMP.seq_occrs = P.seq_occrs;  // unique seq occrs
	else if (version == OCCVERSION)
		newMP.seq_occrs = 1;  // OCCVERSION only counts 1 sequence
		
	#pragma omp critical  // only update one MP at a time
	{
		MPset[MPsetsize] = newMP;
		MPsetsize++;
		if (MPsetsize >= MAX_MP_SET_SIZE) {
			printf("Error: exceeded MP set maximum size!\n");
			exit(1);
		}
	}
}

// prefix-wise less comparison, as required by the stdlib qsort function
int pf_less_cmp(const void *a, const void *b) {
	const char *x = (*(pattern *)a).s;
	const char *y = (*(pattern *)b).s;
	int i = 0;
	while (i < min(strlen(x),strlen(y))) {
		if (x[i] == y[i]) {  // if characters are identical
			++i;  // increment location in both patterns
		} else if ((x[i] == '.') && (y[i] != '.')) {  // found a dot in x
			return 1;  // true  (y <pf x)  - means x is greater
		} else if ((x[i] != '.') && (y[i]  == '.')) {  // found a dot in y
			return -1;  // false (x <pf y)  - means x is less
		} else if (x[i] != y[i]) {  // no dots found and x != y
			++i;  // increment location
		} else {
			break;  // should not reach here
		}
	}
	//return (strlen(y) < strlen(x));  // if no other condition, return shorter
	return 0;  // equilavent patterns
}

// suffix-wise less comparison, as required by the stdlib qsort function
int sf_less_cmp(const void *a, const void *b) {
	const char *x = (*(pattern *)a).s;
	const char *y = (*(pattern *)b).s;
	int i = strlen(x)-1;  // start end of first string
	int j = strlen(y)-1;  // start end of second string
	while (i >= 0 && j >= 0) {
		if (x[i] == y[j]) {  // if characters are identical
			--i;  // decrement end of x
			--j;  // decrement end of y
		} else if ((x[i] == '.') && (y[j] != '.')) {  // found a dot in x
			return 1;  // true (y <sf x) - means x is greater
		} else if ((x[i] != '.') && (y[j] == '.')) {  // found a dot in y
			return -1;  // false (x <sf y) - means x is less
		} else if (x[i] != y[j]) {  // no dots found and x != y
			--i;  // decrement end of x
			--j;  // decrement end of y
		} else {
			break;  // should not reach here
		}
	}
	//return (strlen(y) < strlen(x));  // if no other condition, return shorter
	return 0;
}

// sort elements in EPset_pfless according to prefix-wise less 
void EPsortpfless() {
	printf("Sorting EP pf-less set according to prefix-wise less\n");
	qsort(EPset_pfless, EPsetsize, sizeof(pattern), pf_less_cmp); //stdlib qsort
}

// sort elements in EPset_sfless according to suffix-wise less
void EPsortsfless() {
	printf("Sorting EP sf-less set according to suffix-wise less\n");
	qsort(EPset_sfless, EPsetsize, sizeof(pattern), sf_less_cmp); //stdlib qsort
}

// map EPset_pfless -> EPset_sfless
void EPmapsets() {
	for (uint32_t pf = 0; pf < EPsetsize; pf++) {
		for (uint32_t sf = 0; sf < EPsetsize; sf++) {
			if (strcmp(EPset_pfless[pf].s, EPset_sfless[sf].s) == 0) {  // match
				EPmap[pf] = sf;
			} 
		}
	}
}

// find information about the number and type of equivalence classes in EP
void EPequiv() {
	for (uint32_t e = 0; e < MAX_EP_SET_SIZE; e++) 
		equiv_set_sizes[e] = 1;  // all sets are size 1 to begin with
	memset(equiv_set_first, 0, sizeof(uint32_t)*MAX_EP_SET_SIZE);  // init to 0
	memset(equiv_set_last,  0, sizeof(uint32_t)*MAX_EP_SET_SIZE);  // init to 0
	equiv_set_sizes[0] = 1;  // first set always has size 1 to begin with
	equiv_set_first[0] = 0;  // first index is always the first EP index: 0
	equiv_set_last[0]  = 0;  // last index begins as 0 as well
	for (uint32_t e = 1; e < EPsetsize; e++) {
		// if the EPs are not set equal, then increase sets and new first index
		if (pf_less_cmp(EPset_pfless[e-1].s, EPset_pfless[e].s) != 0) {
			++num_equiv_sets;  // increase number of equivalent sets
			equiv_set_first[num_equiv_sets-1] = e;  // set new first index
			equiv_set_last[num_equiv_sets-1] = e;   // first could be last also
		} else {  // the current equivalent set grows since EP are set equal
			equiv_set_sizes[num_equiv_sets-1]++;  // increase size of this set
			equiv_set_last[num_equiv_sets-1] = e;  // set new last index
		}
	}
	printf("Number of equivalent sets in EP: %u\n", num_equiv_sets);
}

// print the integer mapping between EPset_pfless -> EPset_sfless
void EPmapprint() {
	printf("%-16s\t%-s\t%16s\n", "pf-less", "map", "sf-less"); 
	for (uint32_t i = 0; i < EPsetsize; i++) {
		printf("%-16s\t", EPset_pfless[i].s);
		printf("%-u\t", EPmap[i]);
		printf("%16s\n", EPset_sfless[i].s);
	}
}

// determine if a candidate pattern is in fact maximal by looking at existing MP
bool is_maximal(const pattern P) {
	
	for (uint32_t m = 0; m < MPsetsize; m++) {  // for each existing MP
		if (strcmp(P.s, MPset[m].s) == 0) {  // strings are identical
			return false;  // cannot be maximal if identical to another MP
		}

		if (strlen(P.s) > strlen(MPset[m].s)) {  // candidate string longer
			continue;
		}

		// only compare substrings which overlap
		for (uint32_t mstart = 0; mstart < strlen(MPset[m].s)-strlen(P.s)+1; 
				mstart++) {  // maximal start loc
				
			// check if these substrings are regex equal
			if (re_streq(MPset[m].s+mstart, P.s, strlen(P.s))) {
				// pattern is subsumed, so check num occrs for each
				if (P.num_occrs <= MPset[m].num_occrs) {
					return false;  // offset list size of P is smaller or equal
				}
			}
		}
	}
	return true;  // if all else fails, it is maximal
}

// write Maximal Pattern set to a file (disk)
void writeMPs(const char *filename) {
	printf("Writing Maximal Patterns to output file: %s\n", filename);
	FILE *MPoutfile = fopen(filename, "w");
	if (MPoutfile == NULL) {
		printf("Error: Unable to open output file for writing: %s\n",
				filename);
		exit(1);
	}
	for (uint32_t m = 0; m < MPsetsize; m++) {
		fprintf(MPoutfile, "%u\t%u\t%s\n", 
				MPset[m].num_occrs, MPset[m].seq_occrs, MPset[m].s);
	}
	fclose(MPoutfile);
}
