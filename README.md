[![CI](https://github.com/youcefl/cutrialdive/actions/workflows/build-and-test-cpu.yml/badge.svg)](https://github.com/youcefl/cutrialdive/actions/workflows/build-and-test-cpu.yml) [![CI](https://github.com/youcefl/cutrialdive/actions/workflows/build-gpu.yml/badge.svg)](https://github.com/youcefl/cutrialdive/actions/workflows/build-gpu.yml) [![Sponsor](https://img.shields.io/badge/Sponsor-%E2%9D%A4-red.svg)](https://github.com/sponsors/youcefl)

<h1 align="center">cutrialdive</h1>
<p align="center"><em>GPU-accelerated trial factoring for some number sequences</em></p>

## What it is

`cutrialdive` is a command-line trial-factoring tool for integer sequences whose terms can be evaluated modulo a prime through a recurrence relation.

It is primarily designed for GPU-accelerated trial factoring on recent NVIDIA GPUs using CUDA, but can also be built as a CPU-only binary when no GPU is available.

The current implementation supports Mersenne numbers, Proth numbers, and Smarandache numbers in several bases, and is designed to make it easy to add new sequence families.


## Supported integer sequences

The integer sequences already supported are

| Sequence | Definition | Value of command line switch --num-seq |
|----------|------------|----------------------------------------|
| $M_n$ | ${2}^{n}-1$ | mersenne |
| $Sm_b(n)$ | Concatenation of the first n integers written in base b | smarandache,b |
| $Pr_k(n)$ | $k\cdot{2}^{n}+1$ | proth,k |
| $R_k(n)$ | $k\cdot{2}^{n}-1$ | riesel,k |

## Getting started

Prebuilt binaries are available on the Releases page.  
Alternatively, the project can be built from source with CMake.

Below are a few examples of use to get you started.

$Sm(n) = 12345...\overline{n}$ is the concatenation of the base 10 representation of the first n integers $\ge 1$.

To trial-factor the Smarandache numbers $Sm(10)=12345678910,\dots,Sm(999)=12345\ldots998999$
using all primes in the range $\left[0,2^{37}\right[$:

```bash
$ cutrialdive --num-seq smarandache,10 -s 10 -e 1000 --tf-bits 37 --output sm10_10_1000_37b_f.txt
Trial factoring starts, results will be written to file `sm10_10_1000_37b_f.txt'
TF range is [2, 137438953472[.
Checkpoint file will be written to `sm10_10_1000_37b_f.txt.chkpnt' (period is 360s).
Use `--resume "sm10_10_1000_37b_f.txt.chkpnt"' to resume execution.
Using GPU device 0: NVIDIA GeForce RTX 4060 Ti (15.6 GiB)
Sieving for primes using up to 4 threads.
[100.00% done][137438953447][Elapsed 00:02:08]                                                                                  
Copying prime data to device took 16.16s (cumulated time)
[Factoring took 128.21s (sieve: 72.06s)]
Trial factoring took 128.75s

```

The output is a file where each line is of the form `<index> - [<factors>]`:

```bash
$ head -3 sm10_10_1000_37b_f.txt && echo "..." && tail -3 sm10_10_1000_37b_f.txt
10 - 2, 5, 1234567891
11 - 3, 7, 13, 67, 107, 630803
12 - 2, 3, 2437, 2110805449
...
997 - 13, 284201
998 - 2, 3, 3881, 7927, 1386139
999 - 3, 23, 227, 31139, 632629
```

One can specify a device id using option `--device` and tweak various other parameters, e.g.

```bash
$ cutrialdive --device 0 --threads 8 --num-seq smarandache,10 -s 10 -e 1000 --tf-start $((2**37)) --tf-end $((2**38)) --output sm10_10_1000_37b_38b_f.txt
Trial factoring starts, results will be written to file `sm10_10_1000_37b_38b_f.txt'
TF range is [137438953472, 274877906944[.
Checkpoint file will be written to `sm10_10_1000_37b_38b_f.txt.chkpnt' (period is 360s).
Use `--resume "sm10_10_1000_37b_38b_f.txt.chkpnt"' to resume execution.
Using GPU device 0: NVIDIA GeForce RTX 4060 Ti (15.6 GiB)
Sieving for primes using up to 8 threads.
[100.00% done][274877906899][Elapsed 00:02:01]                                                                                  
Copying prime data to device took 16.72s (cumulated time)
[Factoring took 121.88s (sieve: 82.06s)]
Trial factoring took 122.45s

```

trial factors the same sequence on device 0 but using the prime numbers in range $\left[{2}^{37}, {2}^{38}\right[$

## Some results

The tool has already been used to obtain the following results.

### Large probable primes found

Here $Sm_b(n)$ denotes the base-b Smarandache number obtained by concatenating the base-b representations of $1,2,\dots,n$. For base 10, I simply write $Sm(n)$.

| Number | Size in base 10 | Date |
|--------|-----------------|------|
| $\dfrac{Sm(58628)}{{2}^2\cdot3\cdot1597\cdot2473879}$ | 282023 | Jun 2026 |
| $\dfrac{Sm(52936)}{{2}^3\cdot13\cdot83\cdot3165054961}$ | 253560 | Jun 2026 |
| $\dfrac{Sm(43213)}{7\cdot83}$ | 204956 | Jun 2026 |
| $\dfrac{Sm(42558)}{2\cdot3\cdot7\cdot42709}$ | 201677 | Jun 2026 |
| $\dfrac{Sm(42193)}{106739}$ | 199854 | Jun 2026 |
| $\dfrac{Sm(26896)}{{2}^4\cdot67\cdot283\cdot31546981}$ | 123361 | Jun 2026 |
| $\dfrac{Sm_2(18653)}{13}$ | 74367 | May 2026 |
| $\dfrac{Sm_2(17564)}{{2}^2\cdot3\cdot5\cdot397752637}$ | 69440 | May 2026 |
| $\dfrac{Sm_2(17466)}{2\cdot2459}$ | 69004 | May 2026 |
| $\dfrac{Sm(13365)}{{3}^3\cdot5\cdot1597}$ | 55713 | Jan 2023 |

### Systematic trial factoring results

I systematically trial-factored:

- $Sm_b(n)$ for all $b \in \lbrace2,3,4,5,6,7,8,9,11,12\rbrace$ and all $n \le 1000$, up to ${2}^{42}$
- $Sm_{16}(n)$ for all $n \le 2000$, up to ${2}^{43}$.

This search led to **5862 new factor submissions** to FactorDB and **more than 200 new primality certificates**.


## Performances

The table below gives an idea of the performance one can expect from the tool.

All benchmarks were run on the Proth sequence

$$
    2718281828459 \times 2^n + 1,\qquad n \in [50000000, 50000000+N[
$$

with `--threads 8`, N being the "Candidate count".


| TF bounds                | Number of primes | Candidate count | GPU | Time (s) | $(candidate \times prime)/s$ |
|--------------------------|------------------|-----------------|-----|----------|------------|
| $\left[0, 2^{35}\right[$ | $1'480'206'279$ | $5e5$ | GeForce RTX 5090 | 531 | $1.39 \times {10}^{12}$ |
| $\left[0, 2^{36}\right[$ | $2'874'398'515$ | $5e5$ | GeForce RTX 5090 | 826 | $1.74 \times {10}^{12}$ |
| $\left[0, 2^{37}\right[$ | $5'586'502'348$ | $5e5$ | GeForce RTX 5090 | 1409 | $1.98 \times {10}^{12}$ |
| $\left[0, 2^{40}\right[$ | $41'203'088'796$ | $5e5$ | GeForce RTX 5090 | 9021 | $2.28 \times {10}^{12}$ |
| $\left[0, 2^{35}\right[$ | $1'480'206'279$ | $5e5$ | GeForce RTX 5060 Ti | 2464 | $300.37 \times {10}^{9}$ |
| $\left[0, 2^{36}\right[$ | $2'874'398'515$ | $5e5$ | GeForce RTX 5060 Ti | 3870 | $371.37 \times {10}^{9}$ |
| $\left[0, 2^{37}\right[$ | $5'586'502'348$ | $5e5$ | GeForce RTX 5060 Ti | 6550 | $426.45 \times {10}^{9}$ |
| $\left[0, 2^{40}\right[$ | $41'203'088'796$ | $5e5$ | GeForce RTX 5060 Ti | 42000 | $490.51 \times {10}^{9}$ |

In practical terms, factoring 500,000 Proth candidates to 40 bits completes in 2.5 hours
on an RTX 5090 and 11.7 hours on an RTX 5060 Ti.

On these benchmarks, throughput improves with larger trial-factoring bounds, indicating
that fixed overheads are amortized efficiently as the number of sieved primes grows.
On the RTX 5090, the implementation reaches up to $2.28 \times {10}^{12}$
 $(candidate \times prime)/s$, while the RTX 5060 Ti reaches up to $4.91 \times {10}^{11}$
 $(candidate \times prime)/s$.


## Building from source

### Build requirements

Version numbers below indicate either a minimum required version or a known working version used for development/testing.

#### Required for all builds
- [CMake](https://cmake.org/) (>= 3.25.1)
- A C/C++ compiler (tested with gcc/g++ 12.2)
- [GMP](https://gmplib.org/) (tested with 6.3.0)
- [primesieve](https://github.com/kimwalisch/primesieve) (tested with 12.13)

#### Required for GPU builds
CUDA compiler and runtime (tested with 13.2.78)
Older CUDA versions may work, but can produce slower code.

#### Required if PRP testing is enabled
[GWNUM](https://www.mersenne.org/download/) (tested with 30.19b21)

#### Optional, for testing only
- [PARI/GP](https://pari.math.u-bordeaux.fr/download.html) (tested 2.15.2)

#### Notes on tests

The test suite depends on Catch2, but you should not need to install it manually: CMake fetches Catch2 automatically when tests are enabled.

Tests that depend on PARI/GP are not run by default. They can be run separately using the [pari] tag.

### Typical build commands


The examples below assume that the project sources are located in a directory named `cutrialdive_sources`.

#### GPU build without PRP support
```bash
cmake -S cutrialdive_sources -B build-gpu \
  -DCMAKE_BUILD_TYPE=Release \
  -DCUTRIALDIVE_ENABLE_GPU=ON \
  -DCMAKE_CUDA_ARCHITECTURES=89 \
  -DGMP_ROOT=/path/to/gmp/install \
  -Dprimesieve_DIR=/path/to/primesieve/install/lib/cmake/primesieve

cmake --build build-gpu -j
```

Run the test suite:
```bash
build-gpu/tests/cutrialdive_tests "[pari]"
```

Run PARI/GP based tests (optional):
```bash
build-gpu/tests/cutrialdive_tests "[pari]"
```

#### CPU-only build without PRP support
```bash
cmake -S cutrialdive_sources -B build-cpu \
  -DCMAKE_BUILD_TYPE=Release \
  -DGMP_ROOT=/path/to/gmp/install \
  -Dprimesieve_DIR=/path/to/primesieve/install/lib/cmake/primesieve

cmake --build build-cpu -j
```

When `CUTRIALDIVE_ENABLE_GPU=OFF` (which is the default), the project builds a CPU-only binary
with the same command-line interface (except for the GPU-only switches) but without CUDA support.

#### GPU build with PRP support enabled
```bash
cmake -S cutrialdive_sources -B build-gpu-prp \
  -DCMAKE_BUILD_TYPE=Release \
  -DCUTRIALDIVE_ENABLE_GPU=ON \
  -DCMAKE_CUDA_ARCHITECTURES=89 \
  -DCUTRIALDIVE_ENABLE_PRP=ON \
  -DGMP_ROOT=/path/to/gmp/install \
  -Dprimesieve_DIR=/path/to/primesieve/install/lib/cmake/primesieve \
  -DGWNUM_INCLUDE=/path/to/gwnum \
  -DGWNUM_LIB=/path/to/gwnum/gwnum.a

cmake --build build-gpu-prp -j
```

#### Choosing CMAKE_CUDA_ARCHITECTURES

The value of CMAKE_CUDA_ARCHITECTURES must match the target GPU architecture(s).
For example, 89 targets Ada Lovelace GPUs. If you want one binary to support multiple GPU generations, you can pass a semicolon-separated list such as:

```bash
-DCMAKE_CUDA_ARCHITECTURES="89;103;120"
```
Adjust this list to the GPUs you actually want to support.


## What motivated the project

In late 2021 I became interested in Smarandache numbers. Since finding prime Smarandache numbers is difficult
(a search coordinated on mersenneforum found none for \(n \le 10^6\); see [here](https://www.mersenneforum.org/node/15124?p=661563#post661563)),
I shifted focus to a more tractable problem: finding prime cofactors of Smarandache numbers.

I first wrote a sieve and a naive CPU trial division program and gradually improved them, but the project then remained dormant for about three years.
When I came back to it, I decided to push the idea much further.

Because large probable-prime cofactors of base-10 Smarandache numbers are themselves quite rare, I extended the search to Smarandache numbers in other bases.
That in turn motivated a more general framework: rather than writing a tool for one specific sequence, I wanted a trial-factoring engine that could be reused across several families of recurrence-defined integer sequences.


## Roadmap

Here is a list of points I would like to work on  for the upcoming versions

1. **Closed-form modular evaluation for Smarandache numbers**  
   Use the rational closed form of $Sm_b(n)$ to compute $Sm_b(n) \bmod p$ in $O(\log n)$ big integer operations via modular exponentiation and modular inverses, replacing the current big-integer limb reduction. This is expected to speed up Smarandache searches significantly at large $n$.
2. **Reverse Smarandache numbers support**: support for sequence $Smr(n) = \overline{n}...4321$
3. **Support for stateful sequences**: support for stateful sequences would make
it possible to support sequences for which residue propagation needs more than
the residue of previous term e.g. Fibonacci, Lucas, Smarandache-Wellin, etc.



## Support

If this tool helps your research or saves you compute time, consider 
[sponsoring](https://github.com/sponsors/youcefl) its development. 
Sponsorship helps fund GPU compute time, hardware for testing, and 
continued feature development.

