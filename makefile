# Youcef Lemsafer, 2026.04.21

CXX ?= g++
CC ?= gcc

OUT_DIR ?= ./bin
SRC_DIR ?= ./src

CUTRIALDIVE_VERSION = \"0.0.1\"

GMP_DIR ?= /usr/local
GMP_INCLUDE_DIR ?= $(GMP_DIR)/include
GMP_LIB_DIR ?= $(GMP_DIR)/lib
GWNUM_DIR ?= ../p95v3019b21.source/gwnum

PRIMESIEVE_INCLUDE_DIR ?= /usr/local/include/primesieve
SIEVER_FLAGS ?= -DCUTRIALDIVE_HAS_PRIMESIEVE -I$(PRIMESIEVE_INCLUDE_DIR)
PRIMESIEVE_LIB_FLAGS ?= -L/usr/lib/x86_64-linux-gnu/ -lprimesieve

CXXFLAGS_TO_ADD ?=

CC_OPTIM_FLAGS ?= -O3 -march=native -mtune=native

COMMON_FLAGS = $(CC_OPTIM_FLAGS) -DCUTRIALDIVE_VERSION="$(CUTRIALDIVE_VERSION)" -I$(GMP_INCLUDE_DIR) $(SIEVER_FLAGS)

CFLAGS ?= $(COMMON_FLAGS) -I$(GWNUM_DIR)
CXXFLAGS ?= -std=c++23 -fopenmp $(COMMON_FLAGS) $(CXXFLAGS_TO_ADD)

GMP_LIB_FLAGS ?= -L$(GMP_LIB_DIR) -lgmp
GWNUM_LIB_FLAGS ?= -L../p95v3019b21.source/gwnum -l:gwnum.a
LIB_FLAGS_TO_ADD ?=
LIB_FLAGS ?= $(PRIMESIEVE_LIB_FLAGS) $(GMP_LIB_FLAGS) $(LIB_FLAGS_TO_ADD) $(GWNUM_LIB_FLAGS)


OBJS = $(OUT_DIR)/autotests.o \
	   $(OUT_DIR)/command_line.o \
       $(OUT_DIR)/gw_utility.o \
       $(OUT_DIR)/mode.o \
       $(OUT_DIR)/modular_arithmetic_detail.o \
       $(OUT_DIR)/prp.o \
	   $(OUT_DIR)/smarandache.o \
       $(OUT_DIR)/timer.o \
	   $(OUT_DIR)/trial_factoring.o \
	   $(OUT_DIR)/trial_factoring_options.o

all: $(OUT_DIR)/cutrialdive

$(OUT_DIR)/gw_utility.o: $(SRC_DIR)/gw_utility.c $(SRC_DIR)/gw_utility.h
	$(CC) $(CFLAGS) -c -o $@ $<

$(OUT_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OUT_DIR)/cutrialdive: $(OBJS) $(SRC_DIR)/main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIB_FLAGS)

clean:
	rm -f $(OBJS)
	rm -f $(OUT_DIR)/cutrialdive
