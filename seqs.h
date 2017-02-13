/* Copyright (c) 2011 Lee Nau
 * Released under MIT License (see LICENSE file) */

#ifndef SEQS_H
#define SEQS_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <sys/stat.h> 
#include <sys/mman.h>  // mmap
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

// constants
#define MAX_SEQS 4096  // maximum number of possible sequences allowed

// data structures/types
typedef struct {
	uint32_t length;  // length of this sequence
	char *seq;  // actual raw sequence data, no newlines
	char *header;  // header for this sequence
} sequence;

// functions
uint32_t readseqs(const char *filename, sequence seqs[]);
void freeseqs(const char *filename, sequence seqs[]);

#endif // SEQS_H
