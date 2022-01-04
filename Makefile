# Youcef Lemsafer
# 2021.12.23

CXX=nvcc
CXFLAGS=-std=c++14 -O3 -DYAMPP_HAS_UINT128 -I/media/E/dev/yampp

all: cutrialdiv

cutrialdiv: cutrialdiv.cu
	$(CXX) $(CXFLAGS) -o cutrialdiv cutrialdiv.cu

clean:
	rm cutrialdiv
