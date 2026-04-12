# pas-core-math

A faithful port of the [CORE-MATH](https://core-math.gitlabpages.inria.fr/) binary32 (single-precision floating-point) C library to Free Pascal.

## Overview

CORE-MATH is a research project from INRIA providing correctly-rounded implementations of standard math functions. This project ports the entire binary32 function set to Free Pascal, targeting **bit-exact results** matching the C reference implementation across all 2³² possible single-precision float inputs.

**Key properties:**
- Correctly-rounded results (round-to-nearest-even) for all inputs
- Bit-exact agreement with the C reference library
- Exhaustively verified via full 2³² input space testing
- x86_64 Linux target with inline assembly for performance-critical operations

## Implemented Functions

42 functions implemented across 5 phases:

| Phase | Functions |
|-------|-----------|
| 1 (14) | `rsqrtf`, `tanhf`, `atanpif`, `cospif`, `acosf`, `cbrtf`, `sinpif`, `atanf`, `asinf`, `acospif`, `log2f`, `asinpif`, `tanpif`, `coshf` |
| 2 (17) | `logf`, `exp2f`, `log1pf`, `exp2m1f`, `expm1f`, `exp10f`, `log10f`, `erfcf`, `log2p1f`, `erff`, `sinhf`, `expf`, `atanhf`, `exp10m1f`, `log10p1f`, `asinhf`, `acoshf` |
| 3 (2)  | `tgammaf`, `lgammaf` |
| 4 (5)  | `hypotf`, `atan2f`, `atan2pif`, `powf`, `compoundf` |
| 5 (4)  | `sinf`, `cosf`, `sincosf`, `tanf` |

All functions use the `pcr_` prefix (Pascal Correctly Rounded). The C reference declarations use the `cr_` prefix.

## Repository Layout

<pre>
/
├── <a href="src/">src/</a>
│   ├── <a href="src/pascoremath.pas">pascoremath.pas</a>            # Main library — 42 pcr_* functions
│   ├── <a href="src/pascoremathtypes.pas">pascoremathtypes.pas</a>       # Shared types, bit-cast helpers, MulWide (x86_64 ASM)
│   ├── <a href="src/pascoremathhelperfuncs.pas">pascoremathhelperfuncs.pas</a> # Primitives: fmaf, fabsf, sqrtf, etc.
│   ├── <a href="src/hexfloat.pas">hexfloat.pas</a>               # Utility to parse C99 hex float literals
│   ├── <a href="src/ccoremath.pas">ccoremath.pas</a>              # External declarations for C reference (cr_* functions)
│   ├── <a href="src/pascoremath.inc">pascoremath.inc</a>            # Shared FPC compiler directives
│   └── <a href="src/tests/">tests/</a>
│       ├── <a href="src/tests/TestHarness.pas">TestHarness.pas</a>        # Exhaustive 2³² correctness tester
│       ├── <a href="src/tests/Benchmark.pas">Benchmark.pas</a>          # Throughput benchmark (C CORE-MATH vs PCM)
│       ├── <a href="src/tests/BenchmarkFPC.pas">BenchmarkFPC.pas</a>       # Throughput benchmark (FPC builtins vs PCM)
│       └── <a href="src/tests/build.sh">build.sh</a>               # Build script
├── <a href="install_dependencies.sh">install_dependencies.sh</a>        # Install FPC, GCC, and other dependencies
├── <a href="LICENSE">LICENSE</a>
└── <a href="README.md">README.md</a>
</pre>

## Requirements

- **Free Pascal Compiler** (FPC) 3.2.2 or later
- **GCC** (to compile the C reference library for testing/benchmarking)
- **x86_64 Linux** (inline assembly in `pascoremathtypes.pas` uses x86_64 instructions including `MUL`, `BSF`, and `BSR`)

## Building

```bash
cd src/tests
bash build.sh
```

This compiles both the Pascal library and the C reference library, then links the test and benchmark binaries into `bin/`.

## Running

### Benchmark (PCM vs C CORE-MATH)

Measures throughput (Mops/s) for each function, comparing the C reference and PCM implementations:

```bash
LD_LIBRARY_PATH=src/ bin/Benchmark
```

### Benchmark (PCM vs FPC builtins)

Compares PCM against Free Pascal's built-in math functions (no C dependency required):

```bash
bin/BenchmarkFPC
```

### Exhaustive Correctness Test

Tests every possible 32-bit float input against the C reference. This is a full 2³² = 4,294,967,296 input sweep and takes several hours to complete:

```bash
LD_LIBRARY_PATH=src/ bin/TestHarness
```

Any mismatch is reported as a failure with the input value and both outputs.

## Benchmark Results

### PCM vs FPC builtins (50 million calls per function, x86_64 Linux)

FPC's built-in math functions operate on `Double` internally even for `Single` inputs. PCM targets `Single` precision throughout, which is why it is substantially faster:

| Function   | FPC (Mops/s) | PCM (Mops/s) | Speedup |
|------------|-------------|--------------|---------|
| `sinf`     | 45.9        | 133.7        | 2.9×    |
| `cosf`     | 44.9        | 130.9        | 2.9×    |
| `tanf`     | 17.4        | 129.2        | 7.4×    |
| `asinf`    | 8.5         | 403.2        | 47.4×   |
| `acosf`    | 8.3         | 190.8        | 23.0×   |
| `atanf`    | 28.2        | 310.6        | 11.0×   |
| `sinhf`    | 5.4         | 359.7        | 66.6×   |
| `coshf`    | 5.5         | 400.0        | 72.7×   |
| `tanhf`    | 31.5        | 287.4        | 9.1×    |
| `asinhf`   | 27.8        | 234.7        | 8.4×    |
| `acoshf`   | 6.3         | 186.6        | 29.6×   |
| `atanhf`   | 7.1         | 357.1        | 50.3×   |
| `expf`     | 17.2        | 182.5        | 10.6×   |
| `logf`     | 47.6        | 347.2        | 7.3×    |
| `log2f`    | 44.1        | 365.0        | 8.3×    |
| `log10f`   | 14.5        | 378.8        | 26.1×   |
| `atan2f`   | 17.7        | 94.2         | 5.3×    |
| `hypotf`   | 65.4        | 196.9        | 3.0×    |
| `powf`     | 1.8         | 106.8        | 59.3×   |
| `sincosf`  | 18.2        | 108.0        | 5.9×    |

**On average, PCM is ~22.9× faster than FPC builtins** (arithmetic mean over 20 functions).

## Correctness Guarantee

The test harness performs an exhaustive bit-level comparison between `pcr_*` and `cr_*` for all 2³² single-precision float values (including NaN, infinities, subnormals, and signed zeros). A function is considered correct only when zero mismatches are found across the full input space.

## References

- CORE-MATH project: https://core-math.gitlabpages.inria.fr/
- CORE-MATH paper: Sibidanov, Zimmermann et al. — *Towards a correctly-rounded and fast libm for binary32*
- IEEE 754-2019 standard for floating-point arithmetic
