# Youcef Lemsafer
# 2021.12.23

GMP_DIR = /usr/local
NVCC = nvcc
NVCCFLAGS = -std=c++14 -O3 -I$(GMP_DIR)/include

all: cutrialdiv

cutrialdiv: cutrialdiv.cu
	$(NVCC) $(NVCCFLAGS) -o cutrialdiv cutrialdiv.cu -L$(GMP_DIR)/lib -lgmp

clean:
	rm cutrialdiv
