CXX      = g++
CXXFLAGS = -O3 -std=c++17

# Adjust these to where the libraries are installed/built.
# gbwtgraph: https://github.com/jltsiren/gbwtgraph
# gbwt:      https://github.com/jltsiren/gbwt
# sdsl-lite: https://github.com/simongog/sdsl-lite
GBWTGRAPH_DIR ?= $(HOME)/lib/gbwtgraph
GBWT_DIR      ?= $(HOME)/lib/gbwt
SDSL_DIR      ?= $(HOME)/lib/sdsl-lite

INCLUDES = \
    -I$(GBWTGRAPH_DIR)/include \
    -I$(GBWT_DIR)/include \
    -I$(SDSL_DIR)/include

LIBDIRS = \
    -L$(GBWTGRAPH_DIR)/lib \
    -L$(GBWT_DIR)/lib \
    -L$(SDSL_DIR)/lib

LIBS = -lgbwtgraph -lgbwt -lsdsl -lz -ldivsufsort -ldivsufsort64

all: bin/filter_paths bin/trace_haplotypes

bin/filter_paths: filter_paths.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

bin/trace_haplotypes: trace_haplotypes.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LIBDIRS) $(LIBS)

clean:
	rm -f bin/filter_paths bin/trace_haplotypes
