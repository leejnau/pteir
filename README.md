# pteir
### A Scalable Multicore Implementation of the Teiresias Algorithm
### 2/12/2017

Introduction
--------------------------------------------------------------------------------
The Teiresias motif discovery algorithm was published by IBM Research in 1998 (see references below). This is an open-source implementation of that algorithm, written in C and using the OpenMP parallel programming model. It supports multiple foundational data structures and both stages of the algorithm, scanning and convolution, are parallelized.

This code was originally published on google code (now in archive):  https://code.google.com/archive/p/pteir/
It has now been added to github as a public repository.

This code was published at the Bioinformatics Open Source Conference 2011:  https://www.open-bio.org/wiki/BOSC_2011

The code was utilized in my MS thesis:  https://etd.ohiolink.edu/pg_10?1630639973388::NO:10:P10_ETD_SUBID:61769

Building
--------------------------------------------------------------------------------
The included Makefile should be all that is needed to build this software.
A working OpenMP implementation is required, as well as the gcc compiler.
Simply type "make" in the code directory.

Running
--------------------------------------------------------------------------------
After building the software, the user may run it, using these parameters:
  * -f filename (this is the input file name in FASTA format)
  * -L parameterval (this is the L parameter as described in the papers)
  * -W parameterval (this is the W parameter as described in the papers)
  * -K parameterval (this is the K parameter as described in the papers)
  * -v version (this can be either "seq" or "occ" depending on preference)

For instance, here is a sample configuration:
```bash
./teir -f ecoli.fa -L 3 -W 10 -K 50000 -v occ
```

References
--------------------------------------------------------------------------------
The following publications describe the algorithm in detail:
  * Rigoutsos, I. and A. Floratos, Combinatorial Pattern Discovery in Biological Sequences: the TEIRESIAS Algorithm. Bioinformatics, 14(1), January 1998.  https://academic.oup.com/bioinformatics/article/14/1/55/267266/Combinatorial-pattern-discovery-in-biological

  * Rigoutsos, I. and A. Floratos. Motif Discovery Without Alignment Or Enumeration. Proceedings 2nd Annual ACM International Conference on Computational Molecular Biology (RECOMB '98), New York, NY. March 1998.  http://dl.acm.org/citation.cfm?id=279118
