# Youcef Lemsafer, 2026.04.21

CXX ?= g++
CC ?= gcc

OUT_DIR ?= ./bin

GMP_DIR ?= /usr/local
GMP_INCLUDE_DIR ?= $(GMP_DIR)/include
GMP_LIB_DIR ?= $(GMP_DIR)/lib
GWNUM_DIR ?= ../p95v3019b21.source/gwnum

PRIMESIEVE_INCLUDE_DIR ?= /usr/local/include/primesieve
SIEVER_FLAGS ?= -DCUTRIALDIVE_HAS_PRIMESIEVE -I$(PRIMESIEVE_INCLUDE_DIR)
PRIMESIEVE_LIB_FLAGS ?= -L/usr/lib/x86_64-linux-gnu/ -lprimesieve

CXXFLAGS_TO_ADD ?=

CFLAGS ?= -march=native -mtune=native -O3 -I$(GWNUM_DIR) -I$(GMP_INCLUDE_DIR) $(SIEVER_FLAGS)
CXXFLAGS ?= -std=c++23 -march=native -mtune=native -O3 -fopenmp $(CXXFLAGS_TO_ADD) $(SIEVER_FLAGS) -I$(GMP_INCLUDE_DIR)
GMP_LIB_FLAGS ?= -L$(GMP_LIB_DIR) -lgmp
GWNUM_LIB_FLAGS ?= -L../p95v3019b21.source/gwnum -l:gwnum.a
LIB_FLAGS_TO_ADD ?=
LIB_FLAGS ?= $(PRIMESIEVE_LIB_FLAGS) $(GMP_LIB_FLAGS) $(LIB_FLAGS_TO_ADD) $(GWNUM_LIB_FLAGS)


OBJS = $(OUT_DIR)/smarandache.o $(OUT_DIR)/gw_utility.o \
       $(OUT_DIR)/prp.o \
	   $(OUT_DIR)/timer.o

all: $(OUT_DIR)/cutrialdive

$(OUT_DIR)/gw_utility.o: gw_utility.c gw_utility.h
	$(CC) $(CFLAGS) -c -o $@ $<

$(OUT_DIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OUT_DIR)/cutrialdive: $(OBJS) main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIB_FLAGS)

clean:
	rm -f $(OBJS)
	rm -f $(OUT_DIR)/cutrialdive
