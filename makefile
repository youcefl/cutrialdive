# Youcef Lemsafer, 2026.04.21

CXX ?= g++
CC ?= gcc

BIN ?= ./bin

SIEVER_FLAG ?= -DCUTRIALDIVE_HAS_PRIMESIEVE

GMP_INCLUDE ?=
GWNUM_DIR ?= ../p95v3019b21.source/gwnum

CXXFLAGS_TO_ADD ?=

CFLAGS ?= -march=native -mtune=native -O3 -I$(GWNUM_DIR) $(GMP_INCLUDE)
CXXFLAGS ?= -std=c++23 -march=native -mtune=native -O3 -fopenmp $(CXXFLAGS_TO_ADD) $(SIEVER_FLAG) $(GMP_INCLUDE)
PRIMESIEVE_LIB_FLAGS ?= -L/usr/lib/x86_64-linux-gnu/ -lprimesieve
GMP_LIB_FLAGS ?= -L/usr/lib/x86_64-linux-gnu/ -lgmp
GWNUM_LIB_FLAGS ?= -L../p95v3019b21.source/gwnum -l:gwnum.a
LIB_FLAGS_TO_ADD ?=
LIB_FLAGS ?= $(PRIMESIEVE_LIB_FLAGS) $(GMP_LIB_FLAGS) $(LIB_FLAGS_TO_ADD) $(GWNUM_LIB_FLAGS)


OBJS = $(BIN)/smarandache.o $(BIN)/gw_utility.o \
       $(BIN)/prp.o \
	   $(BIN)/timer.o

all: $(BIN)/smarandachization

$(BIN)/gw_utility.o: gw_utility.c gw_utility.h
	$(CC) $(CFLAGS) -c -o $@ $<

$(BIN)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BIN)/smarandachization: $(OBJS) main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIB_FLAGS)

clean:
	rm -f $(OBJS)
	rm -f $(BIN)/smarandachization
