# Youcef Lemsafer
# 2021.12.23

CC = gcc
CXX = g++
GMP_DIR = /usr/local
GWNUM_DIR = /usr/local/gwnum
CFLAGS = -fPIC -O3 -I$(GMP_DIR)/include -I$(GWNUM_DIR)
CXXFLAGS = -fPIC -O3 -I$(GMP_DIR)/include -I$(GWNUM_DIR)
NVCC = nvcc
NVCCFLAGS = -std=c++14 -O3 -I$(GMP_DIR)/include -I$(GWNUM_DIR)

OBJS = gw_utility.o cutrialdiv.o prp.o

all: cutrialdiv

gw_utility.o: gw_utility.h gw_utility.c
	$(CC) $(CFLAGS) -c -o gw_utility.o gw_utility.c

prp.o: prp.cpp
	$(CXX) $(CXXFLAGS) -c -o prp.o prp.cpp

cutrialdiv.o: cutrialdiv.cu hgint.hpp
	$(NVCC) $(NVCCFLAGS) -o cutrialdiv.o -c cutrialdiv.cu

cutrialdiv: $(OBJS)
	$(NVCC) -o cutrialdiv $(OBJS) -L$(GMP_DIR)/lib -L$(GWNUM_DIR) -lgwnum -lgmp -ldl

clean:
	rm cutrialdiv $(OBJS)
