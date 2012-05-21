
# This affects the memory usage of the program
# we use 1 byte for every 4 bp in kmers. Ideally
# this parameter should be a multiple of 4.
# Actual maximum kmer size is 1 less.
MAX_KMER_SIZE = 32

CC = g++
CXX = g++
INCLUDES = -I.
CXXFLAGS = -c -Wall -Wno-reorder $(INCLUDES) -DMAX_KMER_SIZE=$(MAX_KMER_SIZE) -fPIC
LDFLAGS =
LDLIBS  = -lm -lz

all: CXXFLAGS += -O3
all: target


debug: CXXFLAGS += -g -O0
debug: LDFLAGS += -g
debug: target

profile: CXXFLAGS += -p -g -O2
profile: LDFLAGS += -p -g
profile: clean
profile: target

target: BFGraph

OBJECTS =   Kmer.o KmerIterator.o KmerIntPair.o hash.o fastq.o FilterReads.o BuildContigs.o SimplifyGraph.o KmerMapper.o CompressedSequence.o Contig.o

testread: testread.o $(OBJECTS)
	$(CC) $(INCLUDES) $(LDFLAGS) $(LDLIBS) $(OBJECTS) testread.o -o testread

BFGraph: BFGraph.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) BFGraph.o $(LDFLAGS) $(LDLIBS) -o BFGraph


BFGraph.o: BFGraph.cpp
FilterReads.o: FilterReads.cpp  BloomFilter.hpp
BuildContigs.o: BuildContigs.cpp  BloomFilter.hpp
SimplifyGraph.o: SimplifyGraph.cpp
KmerIntPair.o: KmerIntPair.cpp
fastq.o: fastq.hpp fastq.cpp 
kmer.o: kmer.hpp kmer.cpp
KmerIterator.o: KmerIterator.hpp KmerIterator.cpp
KmerMapper.o: KmerMapper.cpp KmerMapper.hpp
CompressedSequence.o: CompressedSequence.cpp CompressedSequence.hpp
Contig.o: Contig.cpp Contig.hpp

#BloomFilter.o: BloomFilter.hpp
hash.o: hash.hpp hash.cpp	

clean:
	rm -rf *.o
	rm -rf BFGraph
