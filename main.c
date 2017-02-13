/* Copyright (c) 2011 Lee Nau
 * Released under MIT License (see LICENSE file) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>  // gettimeofday
#include "rt.h"
#include "sa.h"
#include "divsufsort.h"
#include "seqs.h"
#include "teir.h"

void print_cmdline_error() {
	printf("Error: improper command line arguments\n");
	printf("Usage:\n");
	printf("-f input_file_name\n");
	printf("-L parameter_value_for_L\n");
	printf("-W parameter_value_for_K\n");
	printf("-K parameter_value_for_W\n");
	printf("-v [seq|occ] sequence support or occurrence support\n");
	exit(1);
}

int main(int argc, char *argv[])
{
	char infilename[256];  // input file name
	char outfilename[256];   // output file name
	uint32_t L,W,K;  // teiresias parameters
	uint8_t version = NOVERSION;  // start as non-seq and non-occ version
	char strver[256]; // string of the version type

	if (argc != 11) {  // program name + 10 args
		print_cmdline_error();
	}

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-f") == 0) {  // input file
			++i;  strcpy(infilename, argv[i]);
		} else if (strcmp(argv[i], "-L") == 0) {
			++i;  L = atoi(argv[i]);
		} else if (strcmp(argv[i], "-W") == 0) {
			++i;  W = atoi(argv[i]);
		} else if (strcmp(argv[i], "-K") == 0) {
			++i;  K = atoi(argv[i]);
		} else if (strcmp(argv[i], "-v") == 0) {
			++i;  
			if (strcmp(argv[i], "seq") == 0) {
				version = SEQVERSION;
			} else {  // default to "occ" version
				version = OCCVERSION;		
			}
			++i;
		} else {
			print_cmdline_error();
		}
	}
	
	if (version == OCCVERSION)
		strcpy(strver, "occ");
	else 
		strcpy(strver, "seq");
	printf("Running TEIRESIAS with the following arguments\n");
	printf("file: %s, L: %u, W: %u, K: %u, version: %s\n",
			infilename, L, W, K, strver);

	sequence seqs[MAX_SEQS];
	uint32_t num_seqs = 0;
	
	struct timeval begin_t, end_t;

	// read in sequence information
	gettimeofday(&begin_t, NULL);
	num_seqs = readseqs(infilename, seqs);
	gettimeofday(&end_t, NULL);
	printf("readseqs took %lu s and %lu us\n",
			end_t.tv_sec - begin_t.tv_sec, end_t.tv_usec - begin_t.tv_usec);

	// data structures
	uint32_t *SAs[num_seqs];  // suffix arrays (one per sequence)
	rt_node *rt;  // radix trie (just one for all sequences)
#if defined(SA)
	create_SAs(SAs, seqs, num_seqs); // create the suffix arrays
#elif defined(RT)
	rt = rt_build(seqs, num_seqs, 1, RT_MAX_DEPTH);
	rt_report_size();
#elif defined(HYBRID)
	gettimeofday(&begin_t, NULL);
	rt = rt_build(seqs, num_seqs, 1, RT_MAX_DEPTH);
	gettimeofday(&end_t, NULL);
	printf("RT construct took: %lu s and %lu us\n", 
			end_t.tv_sec - begin_t.tv_sec, end_t.tv_usec - begin_t.tv_usec);
	rt_report_size();
	gettimeofday(&begin_t, NULL);
	create_SAs(SAs, seqs, num_seqs);
	gettimeofday(&end_t, NULL);
	printf("SA construct took: %lu s and %lu us\n",
			end_t.tv_sec - begin_t.tv_sec, end_t.tv_usec - begin_t.tv_usec);
#elif defined(LINEAR)
	printf("No data structure construction, conducting linear searches\n");
#endif // defined(SA)
	
	// run teiresias
	gettimeofday(&begin_t, NULL);
	scan(rt, seqs, SAs, num_seqs, L, W, K, version);
	gettimeofday(&end_t, NULL);
	printf("Scan phase took: %lu s and %lu us\n", 
			end_t.tv_sec - begin_t.tv_sec, end_t.tv_usec - begin_t.tv_usec);

	gettimeofday(&begin_t, NULL);
	convolve(rt, seqs, SAs, num_seqs, L, W, K, version);
	gettimeofday(&end_t, NULL);
	printf("Convolution phase took %lu s and %lu us\n",
		end_t.tv_sec - begin_t.tv_sec, end_t.tv_usec - begin_t.tv_usec);

	gettimeofday(&begin_t, NULL);
	sprintf(outfilename, "%s_L%u_W%u_K%u_%s.out", infilename, L, W, K, strver);
	writeMPs(outfilename);
	gettimeofday(&end_t, NULL);
	printf("Writing MPs to file took %lu s and %lu us\n",
		end_t.tv_sec - begin_t.tv_sec, end_t.tv_usec - begin_t.tv_usec);

	return 0;
}

