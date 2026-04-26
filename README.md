# pas-core-math

A faithful port of the [CORE-MATH](https://core-math.gitlabpages.inria.fr/) C library to Free Pascal, covering both **binary32** (single precision) and **binary64** (double precision).

> **Status:** binary32 is complete and exhaustively verified across all 2³² inputs. binary64 is undergoing correctness and performance work.

License: [MIT](LICENSE).

## Overview

CORE-MATH is a research project from INRIA providing correctly-rounded implementations of standard math functions. This project ports the CORE-MATH function set to Free Pascal, targeting **bit-exact results** matching the C reference implementation.

**Key properties:**
- Correctly-rounded results (round-to-nearest-even) for all inputs
- Bit-exact agreement with the C reference library (binary32: verified across all 2³² inputs; binary64: see [Correctness Guarantee](#correctness-guarantee))
- x86_64 Linux target with inline assembly for performance-critical operations

### Naming convention

| Prefix | Precision | Reference |
|---|---|---|
| `pcr_*f` | binary32 (`Single`) | C: `cr_*f` |
| `pcr_*`  | binary64 (`Double`) | C: `cr_*` |

`pcr_` stands for *Pascal Correctly Rounded*.

## Implemented Functions

### binary32 — 42 functions (complete)

| Phase | Functions |
|-------|-----------|
| 1 (14) | `rsqrtf`, `tanhf`, `atanpif`, `cospif`, `acosf`, `cbrtf`, `sinpif`, `atanf`, `asinf`, `acospif`, `log2f`, `asinpif`, `tanpif`, `coshf` |
| 2 (17) | `logf`, `exp2f`, `log1pf`, `exp2m1f`, `expm1f`, `exp10f`, `log10f`, `erfcf`, `log2p1f`, `erff`, `sinhf`, `expf`, `atanhf`, `exp10m1f`, `log10p1f`, `asinhf`, `acoshf` |
| 3 (2)  | `tgammaf`, `lgammaf` |
| 4 (5)  | `hypotf`, `atan2f`, `atan2pif`, `powf`, `compoundf` |
| 5 (4)  | `sinf`, `cosf`, `sincosf`, `tanf` |

### binary64 — 41 functions ported (in progress)

| Phase | Functions |
|-------|-----------|
| 1 (13) | `rsqrt`, `tanh`, `atanpi`, `cospi`, `acos`, `cbrt`, `sinpi`, `atan`, `asin`, `acospi`, `log2`, `asinpi`, `tanpi` |
| 2 (18) | `log`, `exp2`, `log1p`, `exp2m1`, `expm1`, `exp10`, `log10`, `erfc`, `log2p1`, `erf`, `sinh`, `exp`, `atanh`, `exp10m1`, `log10p1`, `asinh`, `acosh`, `cosh` |
| 3 (2)  | `tgamma`, `lgamma` |
| 4 (4)  | `hypot`, `atan2`, `atan2pi`, `pow` |
| 5 (4)  | `sin`, `cos`, `sincos`, `tan` |

A small number of binary64 functions still have edge-case mismatches (NaN sign and signed-zero corners; `max_ulp = 0` — i.e., they are correctly rounded for the numeric value).

## Functions Beyond FPC's Math Unit

FPC's `Math` unit does not provide single-function equivalents for 18 of the implemented function families. The table below lists each family by its base name (read with the `f` suffix for binary32, without for binary64 — e.g. `acospi` → `pcr_acospif` and `pcr_acospi`).

| Family | Computes | Nearest FPC expression | Why pas-core-math is better |
|---|---|---|---|
| `acospi` | acos(x) / π | `ArcCos(x) / Pi` | Correctly rounded in one step; no rounding loss from the division |
| `asinpi` | asin(x) / π | `ArcSin(x) / Pi` | Same as above |
| `atanpi` | atan(x) / π | `ArcTan(x) / Pi` | Same as above |
| `atan2pi` | atan2(y,x) / π | `ArcTan2(y,x) / Pi` | Same as above |
| `cospi` | cos(x · π) | `Cos(x * Pi)` | Avoids cancellation error near x = 0.5; returns exact 0.0 there |
| `sinpi` | sin(x · π) | `Sin(x * Pi)` | Returns exact 1.0 at x = 0.5, exact 0.0 at integers |
| `tanpi` | tan(x · π) | `Tan(x * Pi)` | Returns +Inf at x = 0.5; FPC returns a large finite value |
| `cbrt` | x^(1/3) | `Power(x, 1/3)` | Handles negative inputs correctly; correctly rounded |
| `rsqrt` | 1 / √x | `1 / Sqrt(x)` | Single correctly-rounded operation instead of two |
| `exp2` | 2^x | `Power(2, x)` | Faster and correctly rounded |
| `exp2m1` | 2^x − 1 | `Power(2, x) - 1` | Accurate near x = 0 where subtraction cancels |
| `exp10` | 10^x | `Power(10, x)` | Faster and correctly rounded |
| `exp10m1` | 10^x − 1 | `Power(10, x) - 1` | Accurate near x = 0 |
| `expm1` | e^x − 1 | `Exp(x) - 1` | Accurate near x = 0 where subtraction cancels |
| `log1p` | ln(1 + x) | `Ln(1 + x)` | Accurate near x = 0 where addition cancels |
| `log2p1` | log₂(1 + x) | `Log2(1 + x)` | Accurate near x = 0 |
| `log10p1` | log₁₀(1 + x) | `Log10(1 + x)` | Accurate near x = 0 |
| `compound` | (1 + x)^n | `Power(1 + x, n)` | Accurate for small x (binary32 only at present) |

The following four families have **no FPC equivalent** at all:

| Family | Computes |
|---|---|
| `erf` | Error function erf(x) |
| `erfc` | Complementary error function erfc(x) |
| `lgamma` | Natural logarithm of the Gamma function |
| `tgamma` | Gamma function Γ(x) |

## Repository Layout

<pre>
/
├── <a href="src/">src/</a>
│   ├── <a href="src/pascoremath32.pas">pascoremath32.pas</a>          # Main library — binary32 (42 pcr_*f functions)
│   ├── <a href="src/pascoremath64.pas">pascoremath64.pas</a>          # Main library — binary64 (41 pcr_* functions)
│   ├── <a href="src/pascoremathtypes.pas">pascoremathtypes.pas</a>       # Shared types, bit-cast helpers, Mulu64u64 (x86_64 ASM)
│   ├── <a href="src/pascoremathhelperfuncs.pas">pascoremathhelperfuncs.pas</a> # Primitives: fmaf, fabsf, sqrtf, etc.
│   ├── <a href="src/hexfloat.pas">hexfloat.pas</a>               # Utility to parse C99 hex float literals
│   ├── <a href="src/ccoremath32.pas">ccoremath32.pas</a>            # External declarations for C reference (cr_*f)
│   ├── <a href="src/ccoremath64.pas">ccoremath64.pas</a>            # External declarations for C reference (cr_*)
│   ├── <a href="src/pascoremath.inc">pascoremath.inc</a>            # Shared FPC compiler directives
│   ├── <a href="src/inc_64/">inc_64/</a>                    # 47 binary64 port/const includes (one per family)
│   └── <a href="src/tests/">tests/</a>
│       ├── <a href="src/tests/TestHarness32.pas">TestHarness32.pas</a>      # Exhaustive 2³² binary32 correctness tester
│       ├── <a href="src/tests/Benchmark32.pas">Benchmark32.pas</a>        # binary32 throughput (C CORE-MATH vs pas-core-math)
│       ├── <a href="src/tests/BenchmarkFPC32.pas">BenchmarkFPC32.pas</a>     # binary32 throughput (FPC builtins vs pas-core-math)
│       ├── <a href="src/tests/TestHarness64.pas">TestHarness64.pas</a>      # binary64 sampling correctness tester
│       ├── <a href="src/tests/Benchmark64.pas">Benchmark64.pas</a>        # binary64 throughput (C CORE-MATH vs pas-core-math)
│       ├── <a href="src/tests/BenchmarkFPC64.pas">BenchmarkFPC64.pas</a>     # binary64 throughput (FPC builtins vs pas-core-math)
│       └── <a href="src/tests/build.sh">build.sh</a>               # Build script (builds both binary32 and binary64)
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
bash build.sh -dAVX2 -CfAVX2 -CpCOREAVX -OpCOREAVX
```

This compiles both the Pascal library (binary32 and binary64) and the C reference libraries, then links the test and benchmark binaries into `bin/`.

## Running

### Benchmark (pas-core-math vs C CORE-MATH)

Measures throughput (Mops/s) for each function, comparing the C reference and pas-core-math implementations:

```bash
LD_LIBRARY_PATH=src/ bin/Benchmark32
LD_LIBRARY_PATH=src/ bin/Benchmark32 sinf   # run a single function (case-insensitive exact match)

LD_LIBRARY_PATH=src/ bin/Benchmark64
LD_LIBRARY_PATH=src/ bin/Benchmark64 sin    # binary64 single-function run
```

### Benchmark (pas-core-math vs FPC builtins)

Compares pas-core-math against Free Pascal's built-in math functions (no C dependency required):

```bash
bin/BenchmarkFPC32
bin/BenchmarkFPC32 sinf

bin/BenchmarkFPC64
bin/BenchmarkFPC64 sin
```

#### Reducing benchmark variance

Run-to-run numbers — especially for single-function runs — vary with CPU frequency scaling, core migration, and background noise. For more reproducible results, pin the benchmark to one core by prepending `taskset -c 1 env`:

```bash
taskset -c 1 env LD_LIBRARY_PATH=src/ bin/Benchmark32 sinf
taskset -c 1 bin/BenchmarkFPC32 sinf
```

In our measurements, pinning alone cuts variance roughly in half (e.g., ±14% → ±5% spread over 5 runs). On bare metal, also set the CPU governor to `performance` for the most stable numbers:

```bash
sudo cpupower frequency-set -g performance
```

### Correctness Tests

**binary32 — exhaustive 2³² sweep.** Tests every possible 32-bit float input against the C reference. Full sweep takes several hours:

```bash
LD_LIBRARY_PATH=src/ bin/TestHarness32
LD_LIBRARY_PATH=src/ bin/TestHarness32 --func sinf            # single function
LD_LIBRARY_PATH=src/ bin/TestHarness32 --func sinf --pct 1    # 1% sampling for quick checks
```

**binary64 — sampling-based.** A 2⁶⁴ exhaustive sweep is infeasible; the harness uses a mix of structured patterns and random sampling. Use `--pct N` to control coverage:

```bash
LD_LIBRARY_PATH=src/ bin/TestHarness64 --pct 1     # quick regression check (recommended)
LD_LIBRARY_PATH=src/ bin/TestHarness64 --pct 100   # full configured sample set (long)
LD_LIBRARY_PATH=src/ bin/TestHarness64 --func sin --pct 1
```

Any mismatch is reported as a failure with the input value and both outputs.

## Benchmark Results

> Benchmark numbers below are for the **binary32** port. binary64 benchmarking is in progress; run `bin/Benchmark64` and `bin/BenchmarkFPC64` for current numbers on your hardware.

### binary32 — pas-core-math vs FPC builtins (50 million calls per function, x86_64 Linux)

FPC's built-in math functions operate on `Double` internally even for `Single` inputs. pas-core-math targets `Single` precision throughout, which, in part, explains why it is substantially faster.

**On average, pas-core-math is ~2.4× faster than FPC builtins on Windows. On Linux, pas-core-math is 6-20× faster than FPC builtins.**

This Single-vs-Double-internal asymmetry does **not** apply to binary64 — FPC's `Math` unit is natively `Double`, so the binary64 speedup story will be different (and is expected to be flatter). Numbers will be added once the binary64 port stabilises.

## Correctness Guarantee

**binary32.** The test harness performs an exhaustive bit-level comparison between `pcr_*f` and `cr_*f` for all 2³² single-precision float values (including NaN, infinities, subnormals, and signed zeros). A function is considered correct only when zero mismatches are found across the full input space.

**binary64.** A 2⁶⁴ exhaustive sweep is not feasible. Verification combines: (a) random sampling across the input space, (b) structured patterns covering NaN, ±∞, subnormals, signed zeros, and exact-result inputs (integers, powers of two, transcendental landmarks), and (c) regression corners accumulated from previously-found mismatches. Most binary64 functions show zero mismatches at `--pct 1`; a small number still report mismatches with `max_ulp = 0` (i.e., the numeric value is correctly rounded but a NaN sign or signed-zero edge case differs from the C reference). These are tracked and being resolved.

## References

- CORE-MATH project: https://core-math.gitlabpages.inria.fr/
- CORE-MATH paper: Sibidanov, Zimmermann et al. — *Towards a correctly-rounded and fast libm for binary32*
- IEEE 754-2019 standard for floating-point arithmetic
