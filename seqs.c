/* Copyright (c) 2011 Lee Nau
 * Released under MIT License (see LICENSE file) */

#include "seqs.h"

void initseqs(sequence seqs[]) {
	for (uint32_t i = 0; i < MAX_SEQS; i++) {
		seqs[i].header = NULL;
		seqs[i].seq    = NULL;
		seqs[i].length = 0;
	}
}

uint32_t readseqs(const char *filename, sequence seqs[]) {
	printf("Reading input sequence file: %s\n", filename);
	int fd = open(filename, O_RDONLY);  // open file readonly
	if (fd == -1) {
		printf("Error: Unable to open input fasta file: %s\n", filename);
		exit(1);
	}
	struct stat fdinfo;  // file descriptor information
	stat(filename, &fdinfo);  // gather information about the file
	char *filebuf;   // buffer into which the file data is stored from mmap
	uint32_t cur_seq_bytes = 0;  // track how much memory current seq using

	int s = -1;  // total num sequences / current sequence
	const uint32_t LINE_SIZE = 256;
	const uint32_t HEADER_SIZE = 256;
	char line[LINE_SIZE+2];  // +2 for newline and null

	initseqs(seqs);  // initialize all of the sequences

	// memory map the fasta file for reading only
	filebuf = mmap(0, fdinfo.st_size, PROT_READ, MAP_PRIVATE, fd, 0);

	// read sequence data into seqs
	for (int i = 0; i < fdinfo.st_size; i++) {
		uint32_t seqidx = 0;
		if (filebuf[i] == '\n') {  // found solitary newline
			continue;  // go to beginning of for-loop	
		}
		if (filebuf[i] == '>') {  // found header
			++s;  // increment current sequence (to 0 initially, from -1)
			seqs[s].header = malloc(HEADER_SIZE);
			cur_seq_bytes = (1024*1024);  // start at 1MB allocated
			seqs[s].seq = malloc(cur_seq_bytes);
			while (filebuf[i] != '\n') {
				seqs[s].header[seqidx] = filebuf[i];  ++i, ++seqidx;
			}
			seqs[s].header[seqidx] = '\0';  // append null
			printf("Read header: [%s]\n", seqs[s].header);
		} else {  // found sequence line
			seqidx = seqs[s].length;  // reset temporary pointer into sequence
			while (filebuf[i] != '\n') {
				seqs[s].seq[seqidx] = filebuf[i];  ++i, ++seqidx; // copy in seq
				seqs[s].length++;
				if (seqs[s].length == cur_seq_bytes) {
					cur_seq_bytes += (1024*1024);  // add another 1MB chunk
					seqs[s].seq = realloc(seqs[s].seq, cur_seq_bytes);
				}
			}
			seqs[s].seq[seqs[s].length] = '\0';  // append null
			//printf("Sequence [%u] is now: [%s]\n", s, seqs[s].seq);
		}
	}

	close(fd);  // close the file
	munmap(filebuf, fdinfo.st_size);  // unmap the buffer
	return (s+1);  // total number of sequences is s+1 (since 0-based index)
}

void freeseqs(const char *filename, sequence seqs[]) {
	printf("Freeing sequences data for file: %s\n", filename);
	for (uint32_t i = 0; i < MAX_SEQS; i++) {
		if (seqs[i].header != NULL) 
			free(seqs[i].header);
		if (seqs[i].seq != NULL)
			free(seqs[i].seq);
		seqs[i].length = 0;
	}
}

