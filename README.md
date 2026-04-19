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

## Functions Beyond FPC's Math Unit

FPC's `Math` unit does not provide single-function equivalents for 18 of the 42 implemented functions. The table below shows the closest FPC compound expression for each, together with why the pas-core-math version is preferable.

| Function | Computes | Nearest FPC expression | Why pas-core-math is better |
|---|---|---|---|
| `pcr_acospif` | acos(x) / π | `ArcCos(x) / Pi` | Correctly rounded in one step; no rounding loss from the division |
| `pcr_asinpif` | asin(x) / π | `ArcSin(x) / Pi` | Same as above |
| `pcr_atanpif` | atan(x) / π | `ArcTan(x) / Pi` | Same as above |
| `pcr_atan2pif` | atan2(y,x) / π | `ArcTan2(y,x) / Pi` | Same as above |
| `pcr_cospif` | cos(x · π) | `Cos(x * Pi)` | Avoids cancellation error near x = 0.5; returns exact 0.0 there |
| `pcr_sinpif` | sin(x · π) | `Sin(x * Pi)` | Returns exact 1.0 at x = 0.5, exact 0.0 at integers |
| `pcr_tanpif` | tan(x · π) | `Tan(x * Pi)` | Returns +Inf at x = 0.5; FPC returns a large finite value |
| `pcr_cbrtf` | x^(1/3) | `Power(x, 1/3)` | Handles negative inputs correctly; correctly rounded |
| `pcr_rsqrtf` | 1 / √x | `1 / Sqrt(x)` | Single correctly-rounded operation instead of two |
| `pcr_exp2f` | 2^x | `Power(2, x)` | Faster and correctly rounded |
| `pcr_exp2m1f` | 2^x − 1 | `Power(2, x) - 1` | Accurate near x = 0 where subtraction cancels |
| `pcr_exp10f` | 10^x | `Power(10, x)` | Faster and correctly rounded |
| `pcr_exp10m1f` | 10^x − 1 | `Power(10, x) - 1` | Accurate near x = 0 |
| `pcr_expm1f` | e^x − 1 | `Exp(x) - 1` | Accurate near x = 0 where subtraction cancels |
| `pcr_log1pf` | ln(1 + x) | `Ln(1 + x)` | Accurate near x = 0 where addition cancels |
| `pcr_log2p1f` | log₂(1 + x) | `Log2(1 + x)` | Accurate near x = 0 |
| `pcr_log10p1f` | log₁₀(1 + x) | `Log10(1 + x)` | Accurate near x = 0 |
| `pcr_compoundf` | (1 + x)^n | `Power(1 + x, n)` | Accurate for small x |

The following four functions have **no FPC equivalent** at all:

| Function | Computes |
|---|---|
| `pcr_erff` | Error function erf(x) |
| `pcr_erfcf` | Complementary error function erfc(x) |
| `pcr_lgammaf` | Natural logarithm of the Gamma function |
| `pcr_tgammaf` | Gamma function Γ(x) |

## Repository Layout

<pre>
/
├── <a href="src/">src/</a>
│   ├── <a href="src/pascoremath32.pas">pascoremath32.pas</a>          # Main library — 42 pcr_* functions (binary32)
│   ├── <a href="src/pascoremathtypes.pas">pascoremathtypes.pas</a>       # Shared types, bit-cast helpers, Mulu64u64 (x86_64 ASM)
│   ├── <a href="src/pascoremathhelperfuncs.pas">pascoremathhelperfuncs.pas</a> # Primitives: fmaf, fabsf, sqrtf, etc.
│   ├── <a href="src/hexfloat.pas">hexfloat.pas</a>               # Utility to parse C99 hex float literals
│   ├── <a href="src/ccoremath32.pas">ccoremath32.pas</a>            # External declarations for C reference (cr_* functions, binary32)
│   ├── <a href="src/pascoremath.inc">pascoremath.inc</a>            # Shared FPC compiler directives
│   └── <a href="src/tests/">tests/</a>
│       ├── <a href="src/tests/TestHarness32.pas">TestHarness32.pas</a>      # Exhaustive 2³² correctness tester
│       ├── <a href="src/tests/Benchmark32.pas">Benchmark32.pas</a>        # Throughput benchmark (C CORE-MATH vs pas-core-math)
│       ├── <a href="src/tests/BenchmarkFPC32.pas">BenchmarkFPC32.pas</a>     # Throughput benchmark (FPC builtins vs pas-core-math)
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

To enable AVX2 and tune for a modern Intel/AMD core:

```bash
bash build.sh -dAVX2 -CfAVX2 -CpCOREI -OpCOREI
```

This compiles both the Pascal library and the C reference library, then links the test and benchmark binaries into `bin/`.

## Running

### Benchmark (pas-core-math vs C CORE-MATH)

Measures throughput (Mops/s) for each function, comparing the C reference and pas-core-math implementations:

```bash
LD_LIBRARY_PATH=src/ bin/Benchmark32
LD_LIBRARY_PATH=src/ bin/Benchmark32 sinf   # run a single function (case-insensitive exact match)
```

### Benchmark (pas-core-math vs FPC builtins)

Compares pas-core-math against Free Pascal's built-in math functions (no C dependency required):

```bash
bin/BenchmarkFPC32
bin/BenchmarkFPC32 sinf                     # run a single function (case-insensitive exact match)
```

### Exhaustive Correctness Test

Tests every possible 32-bit float input against the C reference. This is a full 2³² = 4,294,967,296 input sweep and takes several hours to complete:

```bash
LD_LIBRARY_PATH=src/ bin/TestHarness32
LD_LIBRARY_PATH=src/ bin/TestHarness32 --func sinf            # run a single function
LD_LIBRARY_PATH=src/ bin/TestHarness32 --func sinf --pct 1    # single function, 1% sampling for quick checks
```

Any mismatch is reported as a failure with the input value and both outputs.

## Benchmark Results

### pas-core-math vs FPC builtins (50 million calls per function, x86_64 Linux)

FPC's built-in math functions operate on `Double` internally even for `Single` inputs. pas-core-math targets `Single` precision throughout, which, in part, explains why it is substantially faster.

**On average, pas-core-math is ~2.4× faster than FPC builtins in Windows. In Linux, pas-core-math is 6-20x times faster than FPC builtins**.

## Correctness Guarantee

The test harness performs an exhaustive bit-level comparison between `pcr_*` and `cr_*` for all 2³² single-precision float values (including NaN, infinities, subnormals, and signed zeros). A function is considered correct only when zero mismatches are found across the full input space.

## References

- CORE-MATH project: https://core-math.gitlabpages.inria.fr/
- CORE-MATH paper: Sibidanov, Zimmermann et al. — *Towards a correctly-rounded and fast libm for binary32*
- IEEE 754-2019 standard for floating-point arithmetic
