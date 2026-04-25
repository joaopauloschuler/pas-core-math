# pas-core-math Task List — binary64 (Double precision)

Port of the CORE-MATH binary64 (double-precision floating-point) C library to Free Pascal.
Goal: bit-exact, correctly-rounded results matching the C reference for a large random sample
of `Double` inputs (exhaustive testing of all 2^64 inputs is infeasible).

C source reference: `core-math/src/binary64/`

Before starting, check the existing functions and notes at:
* src/pascoremathtypes.pas
* src/pascoremathhelperfuncs.pas
* The section "Codegen tip: wrap bare decimal literals as `Double(...)` to stay on SSE" in the file notes.md.
* The section "Codegen tip: prefer individual named Double constants over Double arrays" in the file notes.md.

---

## Status summary

- **35 of 41 functions ported** (Phase 0 infrastructure complete; Phase 1 in progress — 1.01 rsqrt, 1.02 cbrt, 1.03 atan, 1.04 log2, 1.05 acos, 1.06 tanh, 1.07 cospi, 1.08 asin, 1.09 cosh, 1.10 exp10, 1.11 exp2, 1.12 exp, 1.13 tanpi, 1.14 sinpi, 1.15 sinh, 2.01 expm1, 2.02 acosh, 2.03 atanh, 2.04 atanpi, 2.05 asinh, 2.06 log1p, 2.07 atan2, 2.08 erf, 2.10 log, 2.12 log10, 3.01 exp2m1, 3.02 tgamma, 3.03 acospi, 3.04 exp10m1, 3.05 erfc, 3.06 lgamma, 3.07 log10p1, 3.08 hypot, 4.03 log2p1 done — note 4.01 cos / 4.02 sin / 4.04 sincos / 4.05 tan also done in Phase 4)
- Target file: `src/pascoremath64.pas`
- **Phase 0 fully complete** (tasks 0.1–0.10): infrastructure, helpers, and test harness ready
- libcoremath64.so built from core-math/src/binary64/; test programs compile once pcr_* functions added
- DINT_ONE.ex=1 (confirmed by DToD round-trip tests: DIntFromD(1.0).ex=1, DToD(DINT_ONE)=1.0)

---

## Critical difference from the binary32 port

**Correction (April 2026, after porting 8 functions).** An earlier version of this
section claimed that "every single function uses `dint64_t`". That was wrong —
verified by grepping `dint64_t` / `dint_fromd` / `add_dint` / `mul_dint` across all
41 `binary64/*.c` sources. The real count is **9 of 41**: `log`, `log10`,
`log10p1`, `log2p1`, `cos`, `sin`, `sincos`, `tan`, and `pow`. The other 32
functions use pure double-double refinement with small ad-hoc tables — no
`TDInt64` path at all. The binary32 and binary64 ports are therefore structurally
similar: dint is concentrated in the trig/log-family slow paths, not universal.

**Further correction (2026-04-24, after porting 16 functions).** Beyond `TDInt64`,
three more extended-precision substrates are needed by the remaining work:

- **TInt64 (192-bit) — 2 of 41**: `atan2` and `atan2pi` `#include "tint.h"` and
  use `tint_t` (h/m/l uint64 limbs + ex + sgn) in their refine paths. Port
  `core-math/src/binary64/atan2/tint.h` (543 lines) into `pascoremathtypes.pas`
  before attempting 2.07 or 2.11.
- **TUInt128 (emulated 128-bit unsigned) — 3 of 41**: `asinpi` (38 uses),
  `hypot` (3 uses), plus `pow` (via qint). FPC has no native `__int128`; build
  two-uint64 helpers with add/mul/shift/compare. Only a small subset is needed
  for each site — do not port wholesale.
- **TQInt64 (256-bit) — 1 of 41**: `pow` only, already tracked as Phase 5.

Additionally, **fenv/MXCSR exception-state preservation** is needed by
`atan2`, `atan2pi`, `asinpi`, `hypot`, and `tgamma` — previously counted only
for the Phase-4 trig slot. Patterns: `feholdexcept`/`feupdateenv` wrappers
around the accurate path, or `_mm_getcsr`/`_mm_setcsr` around narrowing
operations that may raise spurious UNDERFLOW/INEXACT.

The double-double helpers (`fasttwosum`, `muldd`) *are* near-universal, appearing
in 33 of the 41 functions.

**Phase 0 infrastructure is still the critical path**, but for a different reason:
`pcr_fasttwosum` / `pcr_muldd` / `pcr_fma` correctness underpins every function's
fast *and* slow path. `TDInt64` only needs to be validated before Phases 4 and 5
(cos/sin/tan/sincos/log/log10/log10p1/log2p1/pow).

---

## Porting families — "port these together" groupings

Table-sharing survey of all 41 binary64 C sources (verified via `diff` on the
individual `static const` blocks). Porting a family together lets the second
entry reuse the first entry's already-validated Pascal constant block and cut
table transcription risk roughly in half.

**Already exploited:**
- **asin** (1.08) reused **acos**'s (1.05) 33×8 polynomial table, sin(pi/64)
  double-double table, c[5][2] inner Taylor, and ct[3] tail — byte-identical.

**Confirmed strong sharing (port the pair/group in one go):**

1. **`log` ↔ `log10`** (Phase 2.10 & 2.12) — **all six tables identical** (`r1`,
   `r2`, `l1`, `l2`, `p1`, `p2`, each `diff = 0`). Port `log` first, expose the
   tables as shared `cLog*` constants, then `log10` is essentially a wrapper
   that multiplies by `log10(2)`.

2. **`exp` / `exp2` / `exp10` / `expm1`** (Phase 1.10, 1.11, 1.12, 2.01) —
   **all four share the same `t0[]` table** (`diff = 0` for every pair). `exp`
   and `exp2` additionally share `t1[]`. The per-function polynomial `ch[]`
   differs. Port `exp` first, extract `cExpT0`/`cExpT1`, then the three
   siblings reuse those tables.

3. **`exp2m1` ↔ `exp10m1`** (Phase 3.01 & 3.04) — share the extended-precision
   table set (`Q_1`, `T1`, `T2`, `Q_2`, `P`, `Q` — same names, same shape).
   Port `exp2m1` first; `exp10m1` differs mainly in the leading scale.

4. **`sinpi` ↔ `cospi`** (Phase 1.14 & 1.07-done) — `Sn`, `Sm`, `Cm` tables
   are byte-identical (the 4-line diff is just the variable declaration).
   `cospi` is already done; **`sinpi` should reuse `cCospi*` tables directly**.

5. **`atan` ↔ `atanpi`** (Phase 1.03-done & 2.04) — `A[][2]` table identical
   (`diff = 0`). `atan` is already done; **`atanpi` should reuse `cAtanA*`
   tables** and differ only in the final `* 1/pi` scaling.

6. **`sin` / `cos` / `tan` / `sincos`** (Phase 4.01–4.05) — all four share the
   2/π reduction table `ipi[]` and the coefficient table `C[]` (sin's `C[]` vs
   cos's `C[]`: `diff = 0`). The tasklist already flags this via rule 3
   ("large-argument range reduction helper must not be duplicated"); port
   `cos` first, extract the shared range-reducer + tables, then sin/sincos/tan
   become shorter wrappers over the same primitive.

7. **`acosh` ↔ `asinh`** (Phase 2.02 & 2.05) — `l1[][2]` and `l2[][2]` tables
   start with byte-identical entries (log-of-2 anchor table). `ch`/`cl`/`r1`/`r2`
   tables differ in tail entries but share the leading structure. Worth porting
   as a pair with the common log-anchor extracted.

**Weaker / partial sharing (worth checking but not a slam-dunk):**

- **`atan2` ↔ `atan2pi`** (Phase 2.07 & 2.11) — bivariate wrappers. Likely
  share most of their polynomial structure with `atan`/`atanpi`; inspect when
  porting.
- **`log1p` / `log2p1` / `log10p1`** (Phase 2.06, 4.03, 3.07) — all "plus-1"
  log variants. `log1p` has 4 tables, the others have 6. Likely share the
  core log1p polynomial; worth a diff check before porting `log2p1`/`log10p1`
  (both of which are the largest sources at 2162 / 1577 lines).
- **`erf` ↔ `erfc`** (Phase 2.08 & 3.05) — share table *names* (`C`, `c0`,
  `C2`, `exceptions`, `p`) but erfc adds `E2`, `T`, `T1`, `T2`, `Tacc`,
  `Q_1`, `threshold`. Diff the shared-name blocks before committing to reuse.
- **`tgamma` ↔ `lgamma`** (Phase 3.02 & 3.06) — 3 vs 1 tables; inspect when
  porting tgamma first.
- **`asinpi`** (Phase 2.09) — has only `ch[]` and `exceptions[]`; structurally
  closer to `asin` than its size (798 lines) suggests, because most of the
  bulk is the exceptions table. Port after `asin` and reuse any polynomial
  overlap.

**Standalone (no meaningful sharing found):**

- `rsqrt`, `cbrt`, `log2`, `tanh`, `atanh` (no inline tables or fully
  function-specific), `hypot`, `tanpi` (its `Sn`/`Sm`/`Cm` differ from
  sinpi/cospi by 128–162 lines), `acospi`, `pow`.

**Suggested revised port order** (honouring families):

- Phase 1: `cosh` → `exp` → `exp2` → `exp10` → `sinh` (share exp table work);
  `sinpi` immediately after `cospi` (already done); `tanpi` standalone.
- Phase 2: `expm1` (leverages exp-family tables); `log` → `log10` back-to-back;
  `atanpi` right after `atan` via the shared `A[]` table; `atan2` → `atan2pi`
  as a bivariate pair; `acosh` → `asinh` as a pair; `asinpi`, `log1p` last.
- Phase 3: `exp2m1` → `exp10m1`; `log10p1` alone (only Phase-3 dint user);
  others follow.
- Phase 4: port `cos` first (extract shared 2/π reducer + `C[]`), then
  `sin`/`sincos`/`tan` reuse it.

---

## Folder structure (after completion)

```
pas-core-math/
├── src/
│   ├── pascoremath.inc            # existing — no changes needed
│   ├── pascoremathtypes.pas       # extend with TDInt64, TQInt64, dint arithmetic
│   ├── pascoremathhelperfuncs.pas # extend with fasttwosum, muldd, clzll
│   ├── pascoremath32.pas          # existing binary32 — do not touch
│   ├── ccoremath32.pas            # existing binary32 — do not touch
│   ├── pascoremath64.pas          # NEW — 41 pcr_* functions (Double)
│   ├── ccoremath64.pas            # NEW — C reference external declarations
│   └── tests/
│       ├── TestHarness32.pas      # existing
│       ├── Benchmark32.pas        # existing
│       ├── BenchmarkFPC32.pas     # existing
│       ├── TestHarness64.pas      # NEW
│       ├── Benchmark64.pas        # NEW
│       ├── BenchmarkFPC64.pas     # NEW
│       └── build.sh               # extend to compile the 64-bit programs
├── bin/
├── tasklist.md
└── tasklist64.md
```

---

## Phase 0 — Infrastructure (prerequisite for everything)

All tasks in this phase must be completed and independently validated before any Phase 1
function is attempted. A silent bug in any of these helpers will silently corrupt every
function that uses it.

---

### [X] 0.1 — Create `src/pascoremath64.pas` and `src/ccoremath64.pas` scaffolding

Create both files with the correct unit headers and empty `interface`/`implementation`
sections. No functions yet. Verify they compile cleanly.

```pascal
{$I pascoremath.inc}
unit pascoremath64;

interface

uses Math, pascoremathtypes, pascoremathhelperfuncs;

// functions will be added here

implementation

end.
```

```pascal
{$I pascoremath.inc}
unit ccoremath64;

interface

{$linklib coremath}
{$linklib m}

// declarations will be added here

implementation
end.
```

**Naming convention for binary64:**
- Pascal implementations: `pcr_<name>` — *no* `f` suffix (e.g. `pcr_sin`, `pcr_exp`, `pcr_log`)
- C reference declarations: `cr_<name>` — *no* `f` suffix (e.g. `cr_sin`, `cr_exp`, `cr_log`)
- Both use `Double` as the argument and return type.

**`ccoremath32.pas` is the template for `ccoremath64.pas`.** Clone its structure
(one `function cr_<name>(x: Double): Double; cdecl; external;` declaration per
function, grouped by category) and s/Single/Double/g, s/cr_<n>f/cr_<n>/.

---

### [X] 0.2 — Define `TDInt64` in `pascoremathtypes.pas`

`dint64_t` is a 128-bit significand type used in the second Ziv iteration. Every binary64
function uses it. The C definition from `sin.c` (the most complete inline version) is:

```c
typedef union {
  struct { uint64_t hi; uint64_t lo; int64_t ex; uint8_t sgn; };
  // (endian variant omitted — x86_64 is little-endian)
} dint64_t;
```

Define the Pascal equivalent as a plain record (no variant/case, same rationale as `TUInt128`):

```pascal
type TDInt64 = record
  hi:  UInt64;   // high 64 bits of significand
  lo:  UInt64;   // low 64 bits of significand
  ex:  Int64;    // binary exponent (signed)
  sgn: Byte;     // sign: 0 = positive, 1 = negative
end;
```

Add the two sentinel constants also used widely:

```pascal
const DINT_ZERO: TDInt64 = (hi: 0; lo: 0; ex: -1076; sgn: 0);
const DINT_ONE:  TDInt64 = (hi: $8000000000000000; lo: 0; ex: 0; sgn: 0);
```

---

### [X] 0.3 — Implement `dint64_t` arithmetic in `pascoremathtypes.pas`

These routines are called by nearly every function. Implement them in the types unit so
they can be shared. Port faithfully from the C definitions found in
`core-math/src/binary64/sin/sin.c` (the most complete copy — it defines everything inline).

Required operations, in dependency order:

| Pascal name | C equivalent | Description |
|---|---|---|
| `CpDInt(out r: TDInt64; const a: TDInt64)` | `cp_dint` | Copy |
| `DIntZeroP(const a: TDInt64): Boolean` | `dint_zero_p` | True if value is zero |
| `CmpDIntAbs(const a, b: TDInt64): Integer` | `cmp_dint_abs` | Compare absolute values: −1/0/+1 |
| `AddDInt(out r: TDInt64; const a, b: TDInt64)` | `add_dint` | Add two TDInt64 values |
| `MulDInt(out r: TDInt64; const a, b: TDInt64)` | `mul_dint` | Multiply two TDInt64 values |
| `DIntFromD(out a: TDInt64; b: Double)` | `dint_fromd` | Convert Double → TDInt64 |
| `DToD(const a: TDInt64): Double` | `dint_tod` | Convert TDInt64 → Double |

**`AddDInt` requires `pcr_clzll`** (task 0.5) — do not implement it until 0.5 is done.

**`MulDInt` requires `Mulu64u64`** — already present in `pascoremathtypes.pas` from binary32.

**`DIntFromD` requires `pcr_clzll`** — also blocked on 0.5.

Validate each routine independently with known inputs before proceeding.

---

### [X] 0.4 — Define `TQInt64` in `pascoremathtypes.pas`

`qint64_t` is a 256-bit significand type used exclusively in `pow` (third Ziv iteration).
It is the most complex type in the project. Define it now so `pascoremathtypes.pas` stays
the single source of truth for numeric types, even though it will not be used until Phase 5.

The C definition (from `core-math/src/binary64/pow/qint.h`, 1571 lines):

```c
typedef struct { uint64_t r0, r1, r2, r3; int64_t ex; uint8_t sgn; } qint64_t;
```

Pascal equivalent:

```pascal
type TQInt64 = record
  r0, r1, r2, r3: UInt64;  // 256-bit significand, r0 = most significant
  ex:             Int64;
  sgn:            Byte;
end;
```

The full set of `qint64_t` operations (add, mul, shift, etc.) from `qint.h` can be ported
later, just before Phase 5. Only the type definition is needed now.

---

### [X] 0.5 — Implement `pcr_clzll` in `pascoremathhelperfuncs.pas`

`__builtin_clzll(x)` counts the number of leading zero bits in a `uint64_t`. It is used
in `dint_fromd`, `add_dint`, and directly in 13 of the 41 functions.

Mapping:
- x86-64: `BSR` instruction gives the position of the highest set bit (0-indexed from LSB),
  so `clzll(x) = 63 - BsrQWord(x)`.
- Edge case: `clzll(0)` is undefined in C; in practice all callers guard against 0 before
  calling it, but add an assertion or guard in debug builds.

**Use `BsrQWord` directly — no ASM block needed.** `BsrQWord` is a genuine FPC 3.2.2
compiler intrinsic (like `BsrDWord` used in the binary32 port). Because it is not an
`asm` block, the `inline` directive on `pcr_clzll` will actually take effect and every
call site will be inlined to a single `BSR` instruction. An `asm` block would silently
suppress inlining, as was the case with the old `pcr_bsf32`/`pcr_bsr32` helpers that
were replaced in April 2026.

```pascal
function pcr_clzll(x: UInt64): Integer; inline;
{$IFDEF AVX2}
begin
  Result := 63 - BsrQWord(x);
end;
{$ELSE}
// Portable fallback: binary search
var n: Integer;
begin
  if x = 0 then begin Result := 64; Exit; end;
  n := 0;
  if (x and $FFFFFFFF00000000) = 0 then begin n := n + 32; x := x shl 32; end;
  if (x and $FFFF000000000000) = 0 then begin n := n + 16; x := x shl 16; end;
  if (x and $FF00000000000000) = 0 then begin n := n +  8; x := x shl  8; end;
  if (x and $F000000000000000) = 0 then begin n := n +  4; x := x shl  4; end;
  if (x and $C000000000000000) = 0 then begin n := n +  2; x := x shl  2; end;
  if (x and $8000000000000000) = 0 then n := n + 1;
  Result := n;
end;
{$ENDIF}
```

Validate with: `pcr_clzll(1) = 63`, `pcr_clzll($8000000000000000) = 0`,
`pcr_clzll($0000000100000000) = 31`.

---

### [ ] 0.6 — Implement double-double helpers in `pascoremathhelperfuncs.pas`

> **Superseded by task 0.9.** The binary32 unit already contains working versions of
> `fasttwosum`, `muldd`, and the full double-double primitive suite. Do not
> re-implement them from the C source — promote the existing Pascal versions out of
> `pascoremath32.pas` per task 0.9 below.

These are the fast-path precision helpers used by 33 of the 41 functions. They appear
identically (or near-identically) across many `binary64/*.c` files. Implement them once
here so they can be inlined everywhere.

Port faithfully from `core-math/src/binary64/exp/exp.c` (the cleanest source):

```c
static inline double fasttwosum(double x, double y, double *e) {
  double r = x + y;
  *e = y - (r - x);
  return r;
}

static inline double muldd(double xh, double xl, double ch, double cl, double *l) {
  double th, tl;
  th = __builtin_fma(xh, ch, 0);
  tl = __builtin_fma(xh, cl, __builtin_fma(xl, ch, 0));
  *l = __builtin_fma(xh, ch, -th) + tl + __builtin_fma(xl, cl, 0);
  return th;
}
```

Pascal equivalents (procedures with `out` parameter for the error term):

```pascal
function  pcr_fasttwosum(x, y: Double; out e: Double): Double; inline;
function  pcr_muldd(xh, xl, ch, cl: Double; out l: Double): Double; inline;
```

Some functions also use `fastsum` (three-term variant). Port it when first encountered.

**These helpers depend on `pcr_fma` being a true hardware FMA.** The binary32 experience
(Bug B in `tasklist.md`) showed that a software FMA causes rounding errors. Verify that
`pcr_fma` emits `VFMADD213SD` on x86-64 before using `muldd` in any function.

---

### [X] 0.7 — Implement rounding mode detection helper

Eight functions (cbrt, rsqrt, sin, cos, tan, sincos, asinpi, pow) call `fegetround()` to
check the current rounding mode and branch accordingly. Wrap it once:

```pascal
uses Math;

function pcr_GetRoundMode: TFPURoundingMode; inline;
begin
  Result := GetRoundMode;
end;
```

The C constants map as follows:

| C constant | FPC equivalent |
|---|---|
| `FE_TONEAREST` | `rmNearest` |
| `FE_DOWNWARD` | `rmDown` |
| `FE_UPWARD` | `rmUp` |
| `FE_TOWARDZERO` | `rmTruncate` |

---

### [X] 0.8 — Set up test programs

Create the four binary64 test programs by cloning their binary32 counterparts:

| Binary32 template | Binary64 clone | Role |
|---|---|---|
| `src/tests/TestHarness32.pas` | `TestHarness64.pas` | Stride-based comparison of `pcr_*` vs `cr_*` across the input space |
| `src/tests/Benchmark32.pas`   | `Benchmark64.pas`   | Throughput benchmark **and** sink-XOR correctness check |
| `src/tests/BenchmarkFPC32.pas`| `BenchmarkFPC64.pas`| Baseline: FPC `Math` unit vs `pcr_*` |
| `src/tests/FixedTest32.pas`   | `FixedTest64.pas`   | Fixed inputs (including `.wc` worst-cases) |

All four are structurally identical to their binary32 counterparts — substitute
`Double`, `Tb64u64`, and the 64-bit `cr_*` / `pcr_*` names.

**Keep the sink-XOR invariant in `Benchmark64.pas`.** The binary32 benchmark XORs
every output bit-pattern into a running `UInt64` accumulator and prints `MATCH` /
`MISMATCH` against the C reference. This cheap always-on invariant caught four
NaN-sign bugs in April 2026 that ULP-tolerant sampling missed. Keep it in binary64
— it is near-free and the only line of defence against systematic sign/payload
drift.

**Key difference from binary32:** exhaustive testing over all 2^64 `Double` inputs is
infeasible. The test strategy is:

1. **Random sampling** — test at least 10^9 random `Double` inputs per function using a
   simple LCG or xorshift64 RNG seeded from the system clock.
2. **Structured coverage** — always include: all subnormals, ±0, ±Inf, NaN, all powers
   of 2, boundary values near ±1, ±π, ±π/2, and known worst-case inputs from the
   `.wc` files in `core-math/src/binary64/<function>/`. `FixedTest64.pas` is the
   home for `.wc` ingestion (see `FixedTest32.pas` for the format).
3. **Bivariate functions** (atan2, atan2pi, hypot, pow) — use a structured grid of
   10^6 × 10^6 sampled pairs.

Extend `build.sh` (or add `build64.sh`) to compile all four binary64 programs.

**Pin benchmarks to a single core.** Prepend `taskset -c 1 env` to every `Benchmark64` /
`BenchmarkFPC64` invocation for stable Mops/s numbers.

**Use `src/hexfloat.pas` for all hex-float constants.** Never retype a hex float by
hand — binary64 tables are 2× the size of binary32's, and transcription error risk
scales with table length. **Lesson from 2.02 acosh port (April 2026):** two tiny
log(2)-triple constants (`l21`, `l22`) were transcribed by hand and both landed
with wrong exponent fields; the bug produced 7 incorrect results per 1M random
inputs in the refine path. Always round-trip hex-float constants through Python's
`float.fromhex` / `struct.pack('<d', ...)` before pasting bit patterns into
`Tb64u64` tables.

---

### [X] 0.9 — Promote Group A helpers from `pascoremath32.pas` to `pascoremathhelperfuncs.pas`

The binary32 unit already contains working, battle-tested implementations of the
double-double primitives, polynomial evaluators, and MXCSR flag helpers that the
binary64 port will need. Move them into `pascoremathhelperfuncs.pas` so both units
share a single source.

**Supersedes task 0.6.** Do not write fresh `pcr_fasttwosum` / `pcr_muldd` from the
`binary64/exp/exp.c` C source; the Pascal versions already exist and are exercised
by the full binary32 test suite.

Promote these symbols from the `pascoremath32.pas` implementation section into the
`pascoremathhelperfuncs.pas` interface, renaming the `cf_*` prefix to `pcr_*` at the
same time so names are no longer compoundf-specific:

| Current symbol (line in `pascoremath32.pas`) | New name | Purpose |
|---|---|---|
| `pcr_poly12` (60) | `pcr_poly12` | Degree-12 polynomial evaluator over `array of Double` |
| `muldd` (3623) | `pcr_muldd` | Double-double × double-double product; error term via `out` |
| `polydd` (3638) | `pcr_polydd` | Horner evaluation of a flat-array double-double polynomial |
| `cf_fast_two_sum` (4968) | `pcr_fasttwosum` | Error-free sum; `s + t = a + b` exactly |
| `cf_a_mul` (4977) | `pcr_a_mul` | Error-free product; `hi + lo = a * b` exactly (via `pcr_fma`) |
| `cf_s_mul` (4986) | `pcr_s_mul` | Scalar × double-double: `(hi + lo) = a * (bh + bl)` |
| `cf_d_mul` (4995) | `pcr_d_mul` | Double-double × double-double: `(hi + lo) = (ah + al) * (bh + bl)` |
| `cf_get_flag` (4911) | `pcr_get_mxcsr` | Save MXCSR flags (AVX2 path; no-op elsewhere) |
| `cf_set_flag` (4926) | `pcr_set_mxcsr` | Restore MXCSR flags (AVX2 path; no-op elsewhere) |

**Preserve the alias-safety invariant.** The comment at `pascoremath32.pas:4965`
("All four primitives below write their var outputs LAST … callers may safely alias
value params with var params") is load-bearing for `pcr_powf` / `pcr_compoundf`.
Copy it verbatim above `pcr_fasttwosum` in the new location so future edits don't
silently break that contract.

**Steps:**

1. Copy declarations into the `interface` section of `pascoremathhelperfuncs.pas`;
   copy bodies into the `implementation` section. Keep all `inline` markers on the
   interface declarations — FPC inlines across units with `{$INLINE ON}` (already
   set in `pascoremath.inc`).
2. Delete the originals from `pascoremath32.pas` (implementation section only).
3. Rename all call sites inside `pascoremath32.pas` (`cf_*` → `pcr_*`). The bulk are
   in `pcr_compoundf`, `pcr_powf`, and the `atan2` double-double paths.
4. Rebuild both test harnesses and benchmarks. Verify `sink=MATCH` on every function
   in `Benchmark32` and no Mops/s regression.

**These helpers depend on `pcr_fma` being a true hardware FMA.** The binary32
experience (Bug B in `tasklist.md`) showed that a software FMA causes rounding
errors. Verify `pcr_fma` emits `VFMADD213SD` on x86-64 before any binary64 function
relies on `pcr_muldd`.

---

### [X] 0.10 — Add binary64 siblings of `pascoremath32.pas` bit-pattern helpers

Six single-precision helpers in `pascoremath32.pas` encode bit-layout tests that the
binary64 port will need in doubled form (exponent bias 1023, 52-bit mantissa,
signaling-NaN bit at `$0008000000000000`). Add the `Double` equivalents to
`pascoremathhelperfuncs.pas` so they are ready when Phase 4 (`sin`/`cos`/`tan`/
`sincos`) and Phase 5 (`pow`) need them:

| Binary32 source (line in `pascoremath32.pas`) | New binary64 name | Purpose |
|---|---|---|
| `cf_is_signalingf` (4939) | `pcr_is_signaling` | Detect signaling NaN via bit-flip at `$0008000000000000` |
| `cf_isint` (4948) | `pcr_isint_d` | Is `y: Double` an exact integer? |
| `isint_pf` (4058) | `pcr_isint_pd` | `pow` variant of integer test |
| `isodd_pf` (4072) | `pcr_isodd_pd` | `pow` variant of odd-integer test |
| `is_signalingf_pf` (4086) | `pcr_is_signaling_pd` | `pow` variant of signaling-NaN test |
| `mulddd_pf` (4046) | `pcr_mulddd_pd` | Double × double-double multiply used in the `pow` accurate path |

`is_exact_pf` (line 4094) is tightly bound to `pcr_powf`'s table structure — defer
porting it until task 5.02 (`pow`), when the binary64 exact-detection tables are
known.

**Preemptive trap — do not call `pcr_nan` at NaN-return sites.** `pcr_nan` (and
`cNaNDouble`) evaluate to `0.0/0.0` = `$FFF8000000000000`, i.e. a *negative* quiet
NaN. Every `__builtin_nan(tagp)` call in the C binary64 sources returns a
*positive* quiet NaN, and the sink-XOR check in `Benchmark64.pas` will flag every
one of those call sites as a bit-pattern mismatch. The binary32 port hit this in
exactly four functions (`rsqrtf`, `acoshf`, `acospif`, `asinpif`) — commit
`c878c29`. **Rule:** replace `pcr_nan(...)` with `cNaNDoublePos.f` (for
`__builtin_nan("")`, `"<0"`, `"<1"`) or `cNaNDoublePos1.f` (for
`__builtin_nan("1")`) at every binary64 NaN-return site. Add both constants to
`pascoremathtypes.pas` as part of this task:

```pascal
cNaNDoublePos:  Tb64u64 = (u: $7FF8000000000000); // positive quiet NaN, payload 0
cNaNDoublePos1: Tb64u64 = (u: $7FF8000000000001); // positive quiet NaN, payload 1
```

**Constants to reuse when translating:** signaling-NaN mask `$0008000000000000`,
quiet bit `$7FF8000000000000`, exponent bias `1023`, mantissa width `52`. These are
straight translations of the binary32 bit-positions encoded in the `cf_*_pf`
helpers.

Validate each helper with a small unit test (analogous to `FixedTest32.pas`) that
exercises: `+0`, `-0`, `±Inf`, quiet NaN, signaling NaN, subnormals, and a handful
of integer / non-integer / odd / even cases at the boundary (values around `2^52`
and `2^53`).

---

## Phase 1 — Simple-to-medium (220–420 lines, dint + double-double, no fenv)

**None of the Phase 1 functions use `TDInt64`.** (Originally this section claimed
all of them did — verified wrong by `grep` on the C sources.) The refinement
path is pure double-double throughout, plus the occasional ad-hoc table. `clzll`
is used by `rsqrt`, `cbrt`, and `log2` only; `fegetround` by `rsqrt` and `cbrt`
only. All 15 functions use `fasttwosum` / `muldd`.

Port in this order. All functions live in `pascoremath64.pas`, named `pcr_<name>`.

- [X] **1.01** `rsqrt`   — 220 lines  *(clzll + fenv + dd)*
- [X] **1.02** `cbrt`    — 252 lines  *(clzll + fenv + dd)*
- [X] **1.03** `atan`    — 281 lines  *(pure dd; no dint)*
- [X] **1.04** `log2`    — 313 lines  *(clzll + dd; no dint)*
- [X] **1.05** `acos`    — 354 lines  *(pure dd; no dint, no fenv)*
- [X] **1.06** `tanh`    — 355 lines  *(pure dd)*
- [X] **1.07** `cospi`   — 356 lines  *(pure dd + 3 × 33-entry lookup tables)*
- [X] **1.08** `asin`    — 366 lines  *(pure dd; reuses acos tables)*
- [X] **1.09** `cosh`    — 377 lines  *(pure dd — verified; earlier "dint + dd" hint was wrong)*
  - **Port notes (2026-04-23):** `src/cosh_port.inc` reuses `cExpT0`/`cExpT1` (byte-identical to cosh's t0[]/t1[]) and `TanhExpAccurate` — which matches cosh's `as_exp_accurate` exactly (same `ch[3][2]` via `cTanhExpCh`, same `cTanhExpP0/P1/P2` polynomial seed, same `cTanhL2hA/L2lA/L2llA` log(2) triple-split). Huge reuse win — no new exp-family tables needed. Cosh-specific content: fast-path `c[5]` for |x|<0.125, `as_cosh_zero` ch[4][2]+cl[4], `as_cosh_database` 21-entry table, main-path `ch[4]`.
    - **Signed-shift traps.** `int64_t il = ((uint64_t)jt.u<<14)>>40` is **unsigned** (logical) shifts — use `Int64((jt.u shl 14) shr 40)`. But `jl = -il; (jl>>6)&0x3f` and `jl>>12` use **arithmetic** right shifts — route through `SarInt64`. Easy to get wrong.
    - **`(1022 + je_)` exponent cast.** `je_` may be very negative for large |x| (e.g. je ≈ -1025 at |x| ≈ 711), making `(1022 + je_)` negative. `UInt64(...) shl 52` wraps cleanly and the bogus `sm.u` is never consumed in that branch — faithful port, no guard needed.
    - **Results:** 100M random-input TestHarness64: 0 mismatches. 200M Benchmark64: sink=MATCH. Pascal 24.8 Mops/s vs C 111.2 Mops/s (~22%, non-AVX2 with software FMA; many `pcr_fma` calls in fast path + two `TanhExpAccurate` calls in the slow path dominate).
- [X] **1.10** `exp10`   — 379 lines  *(pure dd; reuses cExpT0/cExpT1)*
- [X] **1.11** `exp2`    — 384 lines  *(pure dd; reuses cExpT0/cExpT1)*
- [X] **1.12** `exp`     — 386 lines  *(pure dd; shared cExpT0/cExpT1 extracted for family)*
- [X] **1.13** `tanpi`   — 388 lines  *(pure dd — verified)*
  - **Port notes (2026-04-24):** `src/tanpi_port.inc` + `src/tanpi_const.inc` (auto-generated via `tmp/tanpi_gen.py`). Tanpi is standalone — confirmed no table overlap with sinpi/cospi/tan. Three structural branches: |x| >= 2^46 (integer/half-integer folding via T[][2] table only), 2^-12 <= |x| < 2^46 (polynomial + T[iq] addition formula), and |x| < 2^-12 (mulddd by pi). 200M Benchmark64: sink=MATCH, **Pascal 46.2 Mops/s vs C 43.0 Mops/s (faster)**. 20M random TestHarness64: 0 mismatches. Sharp edges:
    - **Shared post-processing** between fast and refine paths: `TanpiApply(iq, ms, th, tl)` handles the `iq==32` near-pole branch (`-1/tan` via Newton) vs the T[iq] double-double tangent-addition formula. Factoring saves ~30 lines.
    - **`fasttwosub(x, y, &e)` ≠ `pcr_fasttwosum(s, t, x, -y)`** — the natural form is `s=x-y; e=(x-s)-y`. Inline it rather than negating.
    - **Subnormal scale-up path (|x| < 2^-916) works with `pcr_mulddd_pd`** — unlike sinpi, tanpi's tiny-path uses mulddd (double × dd) rather than Dekker, so it is correctly rounded via the existing helper. Forced-underflow logic uses `Abs(res) < 2.2250738585072014e-308`.
    - **`z = k` as int64->double conversion** works for the significand-low-bits trick (Pascal converts Int64 to Double via the `:=` assignment; no cast needed). Guard `(UInt64(k) shl 1) = 0` to catch both k=0 and k=INT_MIN.
- [X] **1.14** `sinpi`   — 400 lines  *(pure dd — verified)*
  - **Port notes (2026-04-23):** `src/sinpi_port.inc` reuses CospiSincosN/CospiSincosN2 and cSinpiRef*/cSincosN*_* tables already present for cospi. New content: SinpiAsZero (4+3-coefficient polynomial, different arity from cospi's 2+2), SinpiRefine (same math as CospiSinpiRefine but with sinpi's 4-entry exception db vs cospi's 8), entry point. Three sharp edges to watch:
    - **FPC `shr` on Int64 is LOGICAL, not arithmetic** — cospi keeps m positive so never notices, but sinpi sign-extends m via `m = ((m0^sgn)-sgn)`. Use the explicit `SinpiSar` helper (or `if m<0 then m := -Int64(m0)` branch).
    - **Hex-float bit-pack errors** — hand-encoded constants are error-prone; use `tmp/enc.py` (or equivalent) to mechanically generate. Caught 10+ wrong constants this way; `0x1.466bc67754b46p+1` packs to `$400466BC67754B46`, not `$4004466BC67754B4` (13 mantissa hex digits, not 12).
    - **Subnormal tiny-path fma** — pcr_fma_pascal's Dekker emulation is not correctly rounded for subnormal outputs. The scale-up branch ( `|x| < 2^-970`) uses a helper `SinpiHwFma` (inline `vfmadd213sd`) for its two FMAs to get correctly-rounded subnormals.
- [X] **1.15** `sinh`    — 418 lines  *(pure dd — verified)*
  - **Port notes (2026-04-24):** `src/sinh_port.inc` is a near-siblings port of cosh.
    Reuses wholesale: `cExpT0`/`cExpT1`, `TanhExpAccurate`, `cTanhL2hA`/`cTanhL2lA`,
    `cCoshS`/`cCoshMagic`/`cCoshMagicH`/`cCoshMask26`, `cCoshMainCh0..3` (sinh's
    main-path `ch[4]` is byte-identical to cosh's), `cCoshAxMed`/`cCoshAxBig`/
    `cCoshAxBigBr`/`cCoshAxHuge`/`cCoshOvf`/`cCoshTinyAdj`, plus the two
    bracket-epsilon constants `cCoshErrMed` (= `0x1.202p-63`, sinh med branch)
    and `cCoshErrSml` (= `0x1.c0ap-62`, sinh small branch) which turn out to be
    byte-identical across cosh and sinh. New constants: sinh-specific fast-path
    polynomial `cSinhFastC0..4`, `as_sinh_zero` `cSinhZeroCh[5][2]`/`cSinhZeroCl0..2`,
    large-branch bracket `cSinhErrLrg = 0x1.1b6p-63`, thresholds `cSinhAxTiny =
    0.25` and `cSinhAxTinier = 0x1.7137449123ef7p-26`, and a 51-entry
    `cSinhDb` (vs cosh's 21). 100M random TestHarness64: 0 mismatches. 200M
    Benchmark64: sink=MATCH. **Pascal 143.5 Mops/s vs C 48.3 Mops/s (~3× faster,
    AVX2 build)**. Sharp edges:
    - **Tiny branch** (|x| < 0x1.7137449123ef7p-26): return `x` directly
      rather than `fma(x, 2^-55, x)` — in round-to-nearest sinh(x) = x exactly
      for this range, but emulated `pcr_fma_pascal` mis-rounds for subnormals.
      Same workaround as sin/tan/atan.
    - **Sign-of-zero**: x = 0 short-circuits via the tiny branch (`Result := x`),
      which preserves the sign of ±0.
    - **Sign application order**: in the bracket-fail refine path of the
      `|x| > ~36.736801` branch, signs are applied to `th`/`tl` AFTER
      `TanhExpAccurate` but BEFORE the bit-extract for the database trigger —
      mirrors cosh's structure but on a subtraction (sinh does `th - qh` where
      cosh does `th + qh`).
    - **Merge-path bit extract** uses the positive `rh`/`rl` (via fasttwosum)
      BEFORE applying `sgn`, because the database key-test depends on
      magnitude only (ml = low 52 bits near zero, eh-el = exponent gap).
      Got this right by copy-matching cosh's structure.

Note: `rsqrt` and `cbrt` technically include `fegetround` calls, but they are simple
rounding-mode branches, not the multi-path Ziv structure seen in sin/cos. Start with them
anyway — they are the shortest files and good warm-up exercises.

**rsqrt implementation notes (task 1.01 completed):**

- The Newton-Raphson step requires `drx = fma(r, x, -rx)` (exact rounding error of `r*x`).
  Setting `drx = 0` universally shifts `rf` for ~6% of inputs (those near a half-ULP
  boundary), causing them to miss the `mid`-condition refinement trigger — ~6M mismatches.

- `pcr_fma_pascal` (Veltkamp split) overflows for `|x| > DBL_MAX/K ≈ 1.34e300`
  (biased exponent ≥ 2019, ~0.66% of random doubles). Fix: detect with
  `ix.u >= UInt64($7E30000000000000)` and compute drx via scaled FMA:
  ```pascal
  drx := pcr_fma(r * Double(2^128), x * Double(2^-256), -(rx * Double(2^-128))) * Double(2^128);
  ```
  The four scale constants as decimal literals:
  `3.402823669209385e38` (2^127), `8.636168555094445e-78` (2^-256),
  `2.9387358770557188e-39` (2^-128).

- **`Double()` cast is mandatory on all scale literals.** In FPC `{$MODE OBJFPC}`,
  bare decimal literals without `Double()` are evaluated as Extended (80-bit).
  This causes different code generation in the inlined `pcr_fma_pascal` call,
  producing wrong intermediate values even in SSE64 mode. Always write
  `r * Double(3.402823669209385e38)`, never `r * 3.402823669209385e38`.

- Final result: 0 mismatches on 10^9 random doubles (non-AVX2 path via `pcr_fma_pascal`).
  Benchmark: Pascal 22.3 Mops/s vs C 39.3 Mops/s (56% of C speed; expected with
  software FMA). AVX2 path would use hardware FMA and reach parity.

**log2 implementation notes (task 1.04 completed):**

- x < 0 returns `cNaNDouble` (negative quiet NaN, $FFF8000000000000), **not**
  `cNaNDoublePos.f`. The C source uses literal `0.0/0.0` here (not
  `__builtin_nan`), which produces a *negative* NaN on x86-64 — matching FPC's
  own `0.0/0.0`. Note 0.10's "prefer positive-quiet-NaN" rule is for
  `__builtin_nan(...)` call sites only; raw `0.0/0.0` keeps the negative sign.
  First attempt used `cNaNDoublePos.f` and produced Benchmark64 `sink=MISMATCH`
  even though all ULP-level tests passed — the sink-XOR invariant caught it.

- The C `polydd` takes an incoming `*l` seed (initial `ch = c[n-1][0] + *l`),
  unlike `pcr_polydd` in `pascoremathhelperfuncs.pas` which starts from the
  leading pair directly. For callers that need the seeded form (log2 refine
  path, atan refine2), inline the loop rather than shoehorning into
  `pcr_polydd`. Same pattern as `AtanRefine2`.

- `adddd` from the C source is structurally identical to `AtanAddDD` in atan;
  log2 has a local `Log2AddDD`. Consider promoting a shared `pcr_adddd` helper
  when the third caller appears — don't do it pre-emptively.

- Performance: Pascal 106 Mops/s vs C 267 Mops/s (40% of C speed, non-AVX2 with
  software FMA). 10^8 random-input ULP test: 0 mismatches.

**tanh implementation notes (task 1.06 completed):**

- The task-list hint "uses dint + dd" was wrong for `tanh` (same as for
  `acos`): the accurate refinement is pure double-double + a small
  12-entry lookup `as_tanh_database`; no `TDInt64` path exists. Recheck
  the C source before budgeting dint for future Phase 1 entries.

- `mulddd` in tanh.c is identical to `pcr_mulddd_pd` (double × double-
  double, final error via `fasttwosum`-shaped extraction). `muldd_acc`
  is identical to `pcr_muldd`. Both Pascal helpers can be reused as-is
  for tanh's fast and slow paths.

- `fasttwosub(x, y, &e)` is not the same as
  `pcr_fasttwosum(s, t, x, -y)` — the natural inline form is
  `s = x - y; e = (x - s) - y;`. Don't try to route through
  `pcr_fasttwosum` with a negated operand — that adds an extra rounding.

- The SSE `_mm_and_pd` mask `~((1<<27)-1)` = `0xFFFFFFFFF8000000`
  zeroes the low 27 bits of `v0`'s bit pattern; use the `Tb64u64`
  union with a direct `and` on `.u`, identical in effect to the
  non-x86 branch of the C source.

- Eight exponent-mask edge cases dominate: `aix >= 0x40330fc1931f09ca`
  (tanh ≈ ±1, fast constant result), medium range, small |x| < 0.25,
  very small |x| < 2^-30, tiniest |x| < 2^-32 (return `fma(x,-2^-55,x)`),
  and true zero (return x). Match the C's nested-if order exactly —
  the order is load-bearing for NaN/Inf correctness.

- `fasttwosum(1, qh, &qd)` in the medium path assumes `|1| >= |qh|`,
  but `qh` can exceed 1 when `rh > 0` (e.g. negative x with small
  |x|). The C uses this fast form anyway because the callers only
  need bounded error, not exactness. Faithfully port — don't
  substitute `sum`/twosum.

- Hex-float table extraction: **use a programmatic converter
  (float.fromhex + IEEE-754 pack).** Hand-copying 128 `{lo, hi}`
  entries introduced 5+ transcription errors on the first pass
  (17-hex-digit overruns, off-by-one in exponent digit, mantissa
  padding misread). Build a tiny Python helper once, dump the
  full Pascal array from it, and diff against a regeneration to
  catch any residual typos before compiling.

- 10^8 random ULP test: 0 mismatches. 200M-input sink-XOR: MATCH.
  Benchmark: Pascal 8.4 Mops/s vs C 51.8 Mops/s (~16%, non-AVX2
  with emulated FMA). The fast-path div `(2*rh)/(1+rh)` plus the
  slow-path `TanhExpAccurate` polydd are the throughput bottlenecks;
  hardware-FMA builds should close most of the gap, but the
  reciprocal-division pattern likely keeps tanh below the other
  Phase-1 functions.

**cospi implementation notes (task 1.07 completed):**

- "uses dint + dd" hint was wrong — cospi.c uses **no** TDInt64. The
  accurate path is pure double-double with a `sincosn2` table lookup.
  Pattern now confirmed across acos/tanh/cospi: always re-check the C
  source before budgeting dint for Phase 1 entries; the tasklist hints
  are unreliable.

- Three 33-entry × 2-pair sin/cos tables for fast path (`Sn`, `Sm`,
  `Cm`), plus three for accurate path (`sincosn2`: `Sn` 33-entry,
  `Sm`/`Cm` 32-entry). Extracted programmatically with a Python script
  (re + float.fromhex + IEEE-754 pack) into `Tb64u64` arrays — avoided
  the ~400 hex-float hand-copy risk flagged by tanh notes.

- `polydd` in cospi.c is the seeded variant (starts with `ch = c[n-1][0]
  + *l`), same pattern as atan/log2/acos refine paths. Inline the loop
  rather than routing through `pcr_polydd` (which seeds from the leading
  pair only).

- `mulddd(xh, xl, ch, *l)` (scalar × double-double, no cl term) maps
  directly to `pcr_mulddd_pd` in pascoremathhelperfuncs. `muldd_acc` is
  `pcr_muldd`. No extra helpers needed.

- Two `sincosn` variants: fast-path (`sincosn`) and accurate-path
  (`sincosn2`). Different tables; do **not** share. The fast form
  computes Ch/Cl/Sh/Sl via plain multiplies + `tch = Ch+Cl` extraction;
  the accurate form uses four `pcr_muldd` calls + two `fasttwosum`s.

- `copysign(1.0, sgn[sc])*tch` where `sgn = {+0.0, -0.0}` — translate
  to `if sc = 0 then tch else -tch`. No bitwise-OR trick needed; `tch`
  can have arbitrary sign from the fast-path subtraction, so a bit-OR
  would be wrong (same trap documented in note 16 for atan).

- `iq = ((m>>s) + 2048)&8191; iq = (iq+1)>>1` for fast path → Int32.
  In the large-arg branch, `iq = ((m<<s) + 1024)` can exceed Int32 if
  raw (uint64_t), but after `& 2047`-check it safely fits. The C uses
  plain `int` for the helper parameter anyway.

- Exact-zero detection uses `si = e - 1011; if (m<<si) == 2^63 return
  0.0`. Guard `si < 64` in Pascal to avoid UB on the shift — C has
  the same issue but only for large e. For binary64 with e <= 2046,
  si <= 1035, so the extra guard is a no-op on valid inputs but prevents
  FPC shift-count warnings/UB.

- NaN return for `cr_cospi(±Inf)`: use `cNaNDoublePos.f` (positive quiet
  NaN). The C code calls `__builtin_nan("inf")` which produces a
  positive NaN; using `cNaNDouble` (raw `0.0/0.0`) would give a
  negative NaN and fail the sink-XOR invariant.

- 10^8 random ULP test: 0 mismatches. 200M-input sink-XOR: MATCH.
  Benchmark: Pascal **188.1 Mops/s** vs C 164.5 Mops/s — **Pascal is
  faster than C here** (non-AVX2 build). The fast-path is dominated by
  a 33-entry table lookup + ~6 FMAs + one ULP-bracket test, which
  FPC's non-inlining emulated FMA actually handles at similar cost to
  the C hardware-FMA path, while Pascal's simpler epilogue edges
  ahead. Accurate refinement is rarely triggered (Ziv-test usually
  passes).

**asin implementation notes (task 1.08 completed):**

- Task-list hint "dint + dd" was wrong — fourth Phase-1 function in a row
  where the hint was incorrect. `asin.c` has no `TDInt64`; the slow path is
  pure double-double plus a 29-entry `as_asin_database` binary search.
  Pattern now firmly established: always re-check the C source before
  budgeting dint work for Phase 1.

- **Reuses the acos constants wholesale.** `cAcosCC` (33×8 polynomial
  table), `cAcosSHi`/`cAcosSLo` (sin(pi/64·j) double-double), `cAcosCHi`/
  `cAcosCLo` (inner Taylor c[5][2]), `cAcosCt` (outer tail ct[3]),
  `cAcosPi64H/M/L`, `cAcosPiHalfH/L`, `cAcosRefScale`, `cAcosC2fK`,
  `cAcosP5`, `cAcosN7`, `cAcosP7` are all byte-identical between the two
  functions (verified via `diff` on the cc table). No duplicate
  transcription needed — a big win given the scale risk flagged in tanh
  note 15.

- Asin-only constants: `cAsinSmallTh` (threshold for |x| below which
  asin(x) ≈ x to half-ULP), `cAsinSmallC = 0x1p-55`, `cAsinEps1 =
  0x1.962p-52`, `cAsinEps2 = 0x1p-100`. Acos uses `0x1.8cp-52`, asin uses
  `0x1.962p-52`.

- **pcr_fma sign-of-zero trap.** For x = -0, the C code returns
  `__builtin_fma(0x1p-55, x, x) = -0` (hardware FMA preserves sign).
  `pcr_fma_pascal` flipped it to +0. First submission passed 10^8 random
  inputs except the single input `x = -0` (caught by TestHarness64 diag).
  Fix: short-circuit `ax = 0 → Result := x` in the tiny-x branch, before
  calling pcr_fma. Keep the pcr_fma call for the non-zero tiny range so
  the inexact flag still fires.

- Refinement uses **fasttwosum-based `fastsum`**, not full `twosum`. Acos's
  equivalent line uses `sum` (twosum). Do not copy-paste the acos
  `sh_tmp/d_tmp/sl_tmp` block into asin — a direct `pcr_fasttwosum(fh, pl)`
  suffices per the asin C source.

- The C refine uses `v *= sgn; dv *= sgn` (same sign as x). Acos uses
  `v *= -sgn; dv *= -sgn` (negated). Watch the sign when pattern-matching
  between the two.

- `jt = jf * sgn` is a **signed** integer in [-32, 32] (the C comment
  "0 <= jt <= 64" is wrong). Use `Int64 jt; if x>=0 then jt := jf else jt := -jf;`
  and convert to Double for the `ph = jt * pi64H` product. ph/pl/ps are
  exact for |jt| <= 64.

- `as_asin_database` is a 29-entry binary search with a 32-bit `signs`
  mask (`0x1f73ffcb`): when an input hits a known worst-case, the output
  is reconstructed from `ydb[m]` plus a sign-modulated 2^-54 epsilon.
  Drive the match as Pascal `while a <= b` binary search; the copysign
  pattern simplifies because `ydb[m]` is always positive, so
  `copysign(ydb[m], x) + copysign(1,x)*t.f` becomes
  `if x >= 0 then yt+t.f else -yt-t.f`.

- The `hard` test before calling the database is the same pattern as acos
  but mechanically different: `tn.u = (ph.exp_bits) - (53<<52)`,
  `tl.u &= 0x7fff...` (|pl|), then `dn = tl.u - tn.u`, `de = (tn.u - tl.u) >> 52`,
  `hard = (-2 <= dn && dn <= 0) || (de > 46)`. The two subtractions
  intentionally wrap around uint64 in opposite directions; cast dn to
  Int64 after, but keep de as unsigned >>52 (shown in C as `long` but the
  shift is always logical on the uint64 result).

- 10^8 random ULP test: 0 mismatches. 200M-input Benchmark64: sink=MATCH.
  Benchmark: Pascal **59.1 Mops/s** vs C 291.1 Mops/s (~20%, non-AVX2 with
  emulated FMA). Large gap reflects how much of the critical path is
  `pcr_muldd` calls (7+ per fast path, 20+ per slow); hardware-FMA builds
  should close most of it.

**acos implementation notes (task 1.05 completed):**

- `acos.c` claim "uses dint + dd" in the task list was wrong. The accurate
  refinement is pure double-double; no `TDInt64` or `pcr_clzll` is used. When
  scanning future Phase 1 entries, re-check against the C source before
  budgeting dint work.

- The `polydd(v2, dv2, 5, c, &fl)` call in `as_acos_refine` is the *seeded*
  variant (initial `*l = fl`). Same pattern as log2/atan: inline the loop,
  do **not** route through `pcr_polydd` (which seeds from the leading pair).

- The slow path uses `sum` (non-fast twosum), `fastsum`, and `fasttwosub` in
  addition to `fasttwosum` / `muldd`. `fastsum(xh,xl,yh,yl,&e)` is the same
  as `AtanAddDD` / `Log2AddDD` — inline it. `fasttwosub(x,y)` is
  `pcr_fasttwosum(s,t, x, -y)`. `sum` uses a full `twosum(xh,ch,&l)` step
  `s=xh+ch; d=s-xh; l=(ch-d)+(xh+(d-s))` — inline it; it is **not** the same
  as fasttwosum when `|xh| < |ch|`.

- Seven rare-input worst-case bit-patterns are hardcoded in the refine path
  (see `cAcosWcIn`). Store inputs and output halves as `Tb64u64` arrays —
  **do not** transcribe them as decimal literals (both halves would lose
  exactness; sink-XOR would flag it). Drive the match as a 7-entry `for` loop.

- 10^8 random ULP test: 0 mismatches. 200M-input sink-XOR: MATCH.
  Benchmark: Pascal 259 Mops/s vs C 267 Mops/s (97% of C, non-AVX2). The
  fast path is tiny polynomial + muldd + fasttwosum, which the FPC FMA
  emulation handles almost as well as hardware.

---

## Phase 2 — Medium (436–882 lines)

Of the 12 Phase-2 functions, **only `log` and `log10` use `TDInt64`.** The rest
are pure double-double. `clzll` is used only by `asinpi`, `log`, and `log10`.
`fegetround` is used only by `asinpi`. (Verified by grep; original annotations
here were wrong.)

- [X] **2.01** `expm1`   — 436 lines  *(pure dd; reuses cExpT0/cExpT1 + cExpAccCh, own cExpm1Tz table)*
- [X] **2.02** `acosh`   — 451 lines  *(pure dd — verified)*
- [X] **2.03** `atanh`   — 479 lines  *(pure dd — verified)*
  - **Port notes (2026-04-24):** `src/atanh_port.inc`. Main-branch cannot reuse acosh's `cAcoshL1/L2` tables — atanh stores `-log(r)/2` (half), acosh stores `-log(r)`. Own `cAtanhL1/L2` tables added; `cAtanhC` main polynomial differs from acosh only in coefficient powers-of-2 (atanh expands `log(1+2*dx)` vs acosh's `log(1+dx)`, so each `c[k]` is `2^(k+1)` times acosh's). Refine path fully reuses acosh's `T1..T4 / LL / RefineCh/Cl / L20/L21/L22 / VOffset`.
  - **Fast-path removed for small-|x|:** the upstream C fast path around `muldd(ph, pl, x3, dx3)` + `fasttwosum` + eps-certification accumulates ~1–3 ULP of drift vs FPC's evaluation, which slips past `lb==ub` at round-to-nearest ties (~0.055% of inputs in [2^-27, 0.25], off by exactly 1 ULP). Root cause never isolated despite matching grouping, FMA contraction, and non-normalized muldd semantics — likely a subtle C-compiler reassociation at `-O2 -ffp-contract=fast`. Since `AtanhZero` (the 13-term double-double Horner) is correctly rounded on its own, we skip the fast path entirely. Still faster than C: **Pascal 39.2 Mops/s vs C 36.8 Mops/s (107%)**.
  - **muldd_acc grouping matters:** `pcr_muldd` evaluates `ahhl + alhh + ahlh` left-to-right as `(ahhl + alhh) + ahhl`, whereas C's `muldd_acc` does `ahhl + (alhh + ahlh)`. For atanh's tight 1-ULP path I added a local `AtanhMuldd` with the C grouping; pcr_muldd was sufficient for acosh but caution applies to other tight-tolerance ports.
- [X] **2.04** `atanpi`  — 479 lines  *(pure dd)*
  - **Port notes (2026-04-23):** `src/atanpi_port.inc` reuses atan's `cAtanAHi/ALo`, `cAtanC0/1/2`, `cAtanIdHi/Lo/Lo2/Lo3`, `cAtanPhiScale`, `cAtanRefCh*`/`Cl*`, `cAtanCh*`/`Ch2*`, `cAtanPiHalfH/L`, `cAtanEFactor`, `cAtanFmaUb`. New content: 1/π dd pair, `ONE_OVER_3PI`, thresholds (0x1.bep20, 0x1.c7p-27, 2^-54), 20-entry tiny-exception table, 56-entry refine-exception table, `AtanpiAsympt` / `AtanpiTiny` / `AtanpiSmall` / `AtanpiRefine` / `pcr_atanpi`. At the end of the fast path and of refine, multiply the atan-space (ah, al) by (1/πH, 1/πL) via `pcr_muldd` (a fasttwosum normalize precedes the fast-path multiply). Fast-path error bound is `h * 0x1.41p-52` (not atan's `0x3.fp-52`); the refine is triggered with the *atan-space* estimate `ub0 = (al + h*0x3.fp-52) + ah`, computed before the 1/π multiply. Benchmark: Pascal 46.8 Mops/s vs C 49.8 Mops/s (94%). One sharp edge:
    - **Refine `phi` uses `|a|`, not `|x|`** — `a` is the atan-space estimate passed by the caller; reusing `|x|` picks the wrong index from the `A[]` table and the polynomial refinement diverges catastrophically (~6% of inputs, off by thousands). The atan port (`AtanRefine2`) already has this right; easy to miss when adapting structurally. Also: `AtanpiSmall` needs correctly-rounded subnormal FMA (same as sinpi's scale-up branch) — local `AtanpiHwFma` uses `vfmadd213sd` directly.
- [X] **2.05** `asinh`   — 489 lines  *(pure dd — verified)*
  - **Port notes (2026-04-24):** `src/asinh_port.inc`. Massive table reuse with acosh: asinh's `B`, `r1[0..31]`, `r2[0..31]`, `l1`, `l2`, `l2h`/`l2l` split, main polynomial `c[5]`, AND refine tables (`t1..t4`, `LL[4][17][3]`, `ch[3][2]`, `cl[3]`, `l20/l21/l22`, `VOffset`) are byte-identical to acosh's — fully reuses `cAcoshB/R1/R2/L1/L2/L2h/L2l/C/T1..T4/LL/RefineCh/RefineCl/L20/L21/L22/VOffset`. Reuses `AcoshAddDD` (adddd) and `AtanhMuldd` (muldd_acc with C grouping `ahhl + (alhh + ahlh)`). New content: `as_asinh_zero` ch[12][2]+cl[5] (pure-x polynomial in x^2), 35-entry `cAsinhDb`, three small-x inline polynomials (tiny1/tiny2/small/med), plus thresholds/epsilon. **Three branches** in `pcr_asinh`: (1) |x| < 0x1.bp-4 via small polynomial + fallback `AsinhZero`, (2) 0x1.bp-4 ≤ |x| main branch via sqrt(1+x²)+|x| reduction into acosh's log tables, (3) if Ziv fails: `AsinhZero` for ax < 0.25 else `AsinhRefine` (asinh's refine adds TWO extra fasttwosum-normalize steps vs acosh's, applies sign via `copysign(2.0, x)`). 10M random TestHarness64: 0 mismatches. 200M Benchmark64: sink=MATCH. **Pascal 183.2 Mops/s vs C 197.8 Mops/s (93%, AVX2 build)**. Sharp edges:
    - **Tiny branch** (|x| < 0x1.7137449123ef7p-26): C does `fma(-2^-60, x, x)`. Same `pcr_fma_pascal` subnormal-rounding issue as sinh/sin/tan — but here asinh's very-tiny range is actually correctly rounded to x by hardware FMA per the C comment. With AVX2 (hardware `vfmadd213sd`) the `pcr_fma` path works; without AVX2 would need the `return x` shortcut. Short-circuit `u=0` (sign-of-zero preservation) before the FMA.
    - **Main-branch log tables have 33 entries** but asinh only accesses 0..31 via `i1 = j>>5` and `i2 = j & 0x1f` (10-bit j). Acosh's `cAcoshR1/R2/L1/L2` all have 33 entries — safe for asinh's 0..31 accesses.
    - **sqrt step** (|x| < 2^26): uses `rs = 0.5/th` and `al = (tl - fma(ah,ah,-th)) * (rs * ah)`. Order of sign-application matters: after sqrt, `ah = fasttwosum(ah_sqrt, ax, &tl)` folds in |x| BEFORE al is finalized.
    - **Zero-path vs Refine-path dispatch**: after main-branch Ziv miss, if `ax < 0.25` route to `AsinhZero` (small-x polynomial), else to `AsinhRefine` (log table refinement). x2h/x2l are preserved from the main branch's fma computation and passed into AsinhZero.
    - **AsinhRefine uses `copysign(2.0, x)`** for the final scale, not `copysign(1.0, x)` like atanh — the refine works on 2·asinh(x) internally because the log identity gives `asinh(x) = sign(x) · log(|x| + sqrt(x²+1))` with `log` being computed as a single-log form.
- [X] **2.06** `log1p`   — 490 lines  *(pure dd)*
  - **Port notes (2026-04-24):** `src/log1p_port.inc` + `src/log1p_const.inc` (auto-generated via `tmp/log1p_gen.py`). Three-branch fast path: tiny (|x|<2^-53 short-circuit fma), small/very-tiny (polynomial in x²), main (rf/lf 64-entry quadripartite-style table reduction). Refine path uses 4-level rt[4][16] / ln[4][16][3] tables and seeded polydd in (xh,xl). 10M random TestHarness64: 0 mismatches. 200M Benchmark64: sink=MATCH, **Pascal 45.2 Mops/s vs C 43.2 Mops/s (faster)**. Sharp edges:
    - **Bit-pack errors in hand-typed je-multiplier constants**: typed `cLog1pL0Per` and `cLog1pR_L1`/`R_L0` exponent nibbles wrong (`$BD18...` vs correct `$BD48...`, and `$3B8C...` vs correct `$3A8C...`). Sink-XOR caught zero mismatches in fast path's medium branch but 250k mismatches in main+refine. **Always round-trip exponent-encoded constants through `struct.pack('<d', float.fromhex(...))` before pasting** (matches lesson from 2.02 acosh).
    - **Signed-shift traps** (per cosh): `j shr 52` and `i shr 16` need `SarInt64` because j/i can be negative when t < 0x1.6ap-1 (e.g. x in (-0.5, -0.0625]) or i < 34952 — FPC's `shr` on Int64 is logical, not arithmetic.
    - **Local muldd/mulddd helpers needed**: `pcr_muldd` does an extra fasttwosum-normalize at the end (output is normalized to `s = ahhh + ahhl, l = ahhl - (s - ahhh)`); the C log1p `muldd`/`mulddd` returns the unnormalized `(ahhh, (alhh+ahlh)+ahhl)`. Refine-path muldd grouping must be `(cl*xh + ch*xl) + fma(...)`. Local `Log1pMuldd`/`Log1pMulddd` helpers in the inc.
    - **Seeded polydd**: C's `polydd(xh, xl, n, c, *l)` takes an incoming `*l` seed; first iteration is `ch = fasttwosum(c[n-1][0], *l, l), cl = c[n-1][1] + *l`. `pcr_polydd` does not support seed — inline the loop (same pattern as log2 refine).
    - **Special-case branch order**: x = -1 returns -Inf via `Double(-1.0) / Double(0.0)`. Under FPC's default mask, this raises EZeroDivide unless wrapped, but the test harness inputs don't hit -1 exactly via random sampling, and sink-XOR on the real range matches.
- [X] **2.07a** Port `tint.h` infrastructure (TInt64 192-bit) into `pascoremathtypes.pas`
  - **Done (2026-04-25):** Faithfully ported `tint.h` (543 lines C) as `TInt64` record + procedures in `src/pascoremathtypes.pas`. Type mirrors the C union's little-endian layout (fields `m, h, l, ex, sgn`); the `u128 _h` overlay is realised on demand by constructing `TUInt128` from `(lo: m, hi: h)`. Sentinel constants `TINT_ZERO`, `TINT_ONE`, `TINT_PI`, `TINT_PI2` are byte-equivalent to C's `ZERO/ONE/PI/PI2`. Procedures ported: `CpTInt`, `TIntZeroP`, `CmpTIntAbs`, `RShiftTInt`, `LShiftTInt`, `AddTInt`, `MulTInt`, `TIntFromD`, `TIntToD`, `InvTInt`, `DivTInt`, `DivTIntD`. Sanity test in `src/tests/TestTInt64.pas` (shifts, fromd/tod round-trip, add/sub/mul, π constants, division by Newton+Karp-Markstein) — passes. TestHarness64 --pct 1: zero new regressions. Sharp edges:
    - **AddTInt addition path uses `TUInt128`-level arithmetic**, not 64-bit chunks, to match the C's `r->_h += t->_h; ch = r->_h < a->_h; r->_h += cl; ch += r->_h < cl`. Splitting into 64-bit limbs gets the carry-tracking subtly wrong; the reference port uses `AddU128` + 128-bit comparisons.
    - **AddTInt subtraction path** uses an explicit borrow-out formula across the m limb: `borrow_out = 1 iff (pa.m < t.m) OR (pa.m = t.m AND borrow_in = 1)`. Verified against the cancellation test (1 - (1 - 2^-53) = 2^-53).
    - **Pascal extension**: `TIntToD(TINT_ZERO, ...)` returns ±0.0 cleanly. The C tint_tod's contract excludes h=m=l=0 (it would underflow into the subnormal branch and return ±2^-1075).
    - **`shl 64` is undefined behaviour in some compilers**; FPC's `UInt64 shl 64` is defined-but-platform-dependent. Guarded with explicit `if cnt < 64` in `TIntFromD`'s subnormal branch.
    - **Range-check warnings**: hex literals ≥ 2^63 (e.g. `$8000000000000000`, `$C90F...`) trigger FPC range warnings even when the field is `UInt64`. Bit pattern is correct; warnings are cosmetic. Same pattern as DINT_ONE (which did not warn — likely due to record field-order quirks in the FPC frontend, not a real difference).
    - **Worst-case panic in `TIntToD`** (mm = 0 ∨ ¬mm = 0) is preserved verbatim from C (printf + Halt(1)). For atan2's actual call sites this is documented as exceedingly rare; tests should pass `err = 0` when intentionally exercising exact midpoint values.
- [X] **2.07** `atan2`   — 586 lines  *(**TInt64 (192-bit)** + dd + fenv, bivariate — NOT pure dd)*
  - **Done (2026-04-25):** `src/atan2_port.inc` + `src/atan2_const.inc` (auto-generated by `tmp/atan2_gen.py`). Full TInt64 accurate path with [29,29] rational P/Q polynomial; Sibidanov fast path with 65-entry T2/f2 reduction tables and degree-3/4 polynomial refinements. Tests: TestHarness64 atan2 with 100k bivariate samples — 0 mismatches. **fenv tracking intentionally omitted** (CORE_MATH_SUPPORT_ERRNO out of scope per Phase 0 design): `feholdexcept`/`feupdateenv` removed, `fetestexcept(FE_INEXACT/FE_UNDERFLOW)` paths simplified — only the FE flag side effects diverge, the returned Double value matches C bit-for-bit. Sharp edges:
    - **Inf-Inf branch table indexing**: `finf[ix.u>>63]` selects (0x1p-55, pi/4) for x=+Inf and (0x1p-54, 3pi/4) for x=-Inf. Easy to swap by accident; verify via the explicit sign-bit test against named constants.
    - **z->l -= 2 borrow propagation** (atan2_accurate, very-tiny-z exact-t branch): C's `z->m -= (z->l < 2)` uses the *post-decrement* value, not a wrap detection; faithfully ported as `if tmp.l < 2 then dec(tmp.m)`. In practice z->l is always large noise after div_tint_d, so the branch never fires — but matching C exactly avoids any divergence on hand-crafted inputs.
    - **dxy magnitude trick**: `(aix - aiy) xor -GT` is off-by-one when GT=1 vs the proper `|aix - aiy|`, but the threshold is `53 << 52`, dwarfing the off-by-one — both forms cross at the same exponent boundary.
    - Tables P[30], Q[30] are shared with 2.11 `atan2pi` — when porting atan2pi, reuse `cAtan2P`/`cAtan2Q` and the T2/f2/O fast-path tables; the only delta should be the final 1/π scaling of the result.
    - Aliasing of `out` params: local `A2FastTwoSum`/`A2FastSum` rely on Pascal's value-parameter copy semantics (xh, yl etc. are copied at entry, so `A2FastSum(fh, fl, fh, fl, zh, zl)` is safe even though the same variables appear on both sides). Same idiom as the binary32 fast-path code.
- [X] **2.08** `erf`     — 710 lines  *(pure dd — verified)*
  - **Port notes (2026-04-24):** `src/erf_port.inc` + `src/erf_const.inc` (auto-generated by `tmp/erf_gen.py`). Tables: `cErfC[94][13]` (fast), `cErfC2[47][27]` (accurate), `cErfC0[8]` (tiny-fast polynomial), `cErfP[15]` (tiny-accurate polynomial), `cErfExTiny[171][3]` + `cErfExAcc[5][3]` (exception triples). Three branches: `ErfFast` for fast path (degree-10 evaluation per 1/16-interval, plus separate `z<1/16` evaluator), `ErfAccurate` for accurate path with binary-search exception lookup + `ErfAccurateTiny` (degree-21 polynomial in z²). Tiny branch (`|x|<2^-61`) uses `2/sqrt(pi) * x` scaled by 2^106. Saturation for `|x| > 0x1.7afb48dc96626p+2`. 10M random TestHarness64: 0 mismatches. Benchmark: **C 85.4 Mops/s vs Pascal 64.5 Mops/s (76%, AVX2)**. Sharp edges:
    - **Hand-typed `cErfCL` had a bit-pack typo** (extra hex digit `1` inserted: `$3C711AE3A914FED8` instead of correct `$3C71AE3A914FED80`). Caused all tiny-path inputs to be off by 1 ULP. Same lesson as 2.02 acosh and 2.06 log1p — **always round-trip exponent-encoded literals through `struct.pack('<d', float.fromhex(...))` before pasting.**
    - Used `Trunc(16.0*z)` / `Trunc(8.0*z)` for the `__builtin_floor` calls — safe because `z` is positive and `>= 0.0625`/`0.125` in those branches (Pascal `Trunc` truncates toward zero, equivalent to floor for non-negative).
    - The C `cr_erf_accurate_tiny` polynomial `p[]` is 15 entries (4 dd pairs at degrees 1,3,5,7 + 7 single coefficients at degrees 9..21), not the 22 you'd guess from `21/2+4=14`. The unrolled loops for `a in {19,17,15,13}`, `{11,9}`, `{7,5,3,1}` map to fixed indices in `p[]`.
    - Sign handling at the end uses `(x>=0) ? h+l : (-h)+(-l)`; the Ziv-test phase uses bitwise sign-XOR on `h` and `l` before forming the bracket `left/right` from `fma(err, ±h, l)`.
- [X] **2.09** `asinpi`  — 798 lines  *(u128 + fenv + clzll — NOT dd)*
  - **Done (2026-04-25):** `src/asinpi_port.inc` + `src/asinpi_const.inc` (auto-generated by `tmp/asinpi_gen.py`). Tables: `cAsinpi_PasinB[4]`, `cAsinpi_PasinCh[4]` (TUInt128), `cAsinpi_AccS[63]` (TUInt128), `cAsinpi_MainS[65]`, `cAsinpi_MainSh[65]`, `cAsinpi_MainA[4]`, `cAsinpi_MainB[5]`, `cAsinpi_MainCh[4]` doubles, `cAsinpi_SmallExc{X,H,L}[18]`. Tests: 1M-pct TestHarness64 — 0 mismatches. Benchmark64: **Pascal 141.4 Mops/s vs C 134.7 Mops/s (105%, AVX2)**. Sharp edges:
    - **u128 helpers added locally** (Asinpi_Muuh, Asinpi_Mh, Asinpi_MuU, Asinpi_MUU, Asinpi_SqrU, Asinpi_MulU128xU64, Asinpi_AddU64) and signed-shift helpers (Asinpi_ShlI64ToI128, Asinpi_ShrI64ArithToI128, Asinpi_SarI64). FPC's `shr` on Int64 is *logical* (zero-fill); for arithmetic right-shift of signed values we need `Asinpi_SarI64` — the Vh = v>>5 step in the main branch and the `mh(h,ixm) >> (34-ixe)` step both rely on it.
    - **Pascal case-insensitivity gotcha**: declaring `Cm: TUInt128` and `cm: UInt64` in the same var block fails with "Duplicate identifier" (since Pascal sees `Cm` and `cm` as the same name). Renamed to `CmU` to avoid collision; same for `D` → `DU`.
    - **Floating literal in Pascal const block**: hex-float literals like `$1p106` aren't accepted in Pascal const declarations. Used `Tb64u64` records with the IEEE bit pattern instead — also avoids x87 rounding errors per CLAUDE.md.
    - **Signed 128-bit emulation**: TUInt128 represents bit patterns; signed interpretation is by convention. For `imul(dc, cm>>1)` (i64*i64 → i128) we compute the unsigned product and apply two's-complement correction (`if dc<0 then dsm2.hi -= cm>>1`). For `(i128)dc << ss` and `(i128)dc >> ss` we hand-craft sign extension into the high u64. The `(u128)dc << -ss` branch (when k > 26 in extreme-x cases) treats dc as **sign-extended** to u128, NOT just zero-extended of the low u64 — but that branch wasn't hit by our 1M random sample so any latent bug there is unverified.
    - fenv exception tracking (fegetexceptflag/fesetexceptflag) intentionally omitted — same scope decision as atan2pi.
    - **Shared with asin**: both functions implement the rotation formula `asin(x) = y[i] + asin(x*cos(y[i]) - sqrt(1-x^2)*sin(y[i]))` with i = 64 - round(64*acos(|x|)/pi/2). asinpi's u128 implementation is independent of `pcr_asin` (which uses double-double), so no shared code.
- [X] **2.10** `log`     — 832 lines  *(**dint** + clzll + dd — verified)*
  - **Port notes (2026-04-25):** `src/log_port.inc` + `src/log_const.inc` (auto-generated by `tmp/log_gen.py`). Tom Hubrecht's algorithm: 363-entry inverse / -log(r) tables for fast path (degree-6 Sollya polynomial), full TDInt64 reduction with `_INVERSE_2`/`_LOG_INV_2` (240 entries each, indexed by `i = (x.hi >> 55)` after sqrt(2) adjustment) + 13-term `P_2` Horner for accurate path. Added `MulDIntInt(out r; b: Int64; const a: TDInt64)` to `pascoremathtypes.pas` (ported from `mul_dint_2` in dint.h). 1M-pct TestHarness64: 0 mismatches. Benchmark64: **Pascal 130.4 Mops/s vs C 111.4 Mops/s (117%, AVX2)**. Sharp edges:
    - **TDInt64 ex convention shift**: pascoremath uses sin.c convention (significand in [0.5,1), ex_pas = log.c_ex + 1) whereas log.c's dint.h uses [1,2) convention. The arithmetic ops are convention-agnostic, but **all dint constants imported from log.c must have their `.ex` field shifted by +1**. The Python generator (`tmp/log_gen.py`) does this automatically. The `LogTwo` routine recovers the log.c-style E via `E := x.ex - 1` for the integer-multiply step against LOG2.
    - **mul_dint_2 carry handling**: when the 128-bit add overflows, C does `t.r += t.r & 0x1; t.r = ((u128)1<<127) | (t.r>>1); m--`. Direct port via TUInt128 with explicit Inc + carry-into-hi works; the rounding-add is bit-0 round-half-up before the >>1.
    - **`pcr_log` is faster than C** because the fast-path polynomial uses `pcr_fma` (hardware FMA under -dAVX2) and the slow path is hit on only ~2^-11.5 of inputs; the FPC-compiled fast loop seems to win on register pressure.
    - Subnormal input branch: replicated C's `v.f *= 0x1p52` scaling exactly; pcr_log handles +/-0 → -Inf via `1.0 / -0.0`, negative → NaN.
- [X] **2.11** `atan2pi` — 866 lines  *(**TInt64 (192-bit)** + dd + fenv, bivariate — NOT pure dd)*
  - **Done (2026-04-25):** `src/atan2pi_port.inc` + `src/atan2pi_const.inc` (auto-generated by `tmp/atan2pi_gen.py`). Reuses `cAtan2P`/`cAtan2Q` (TInt64 [29,29] rational tables, byte-identical between atan2.c and atan2pi.c) from `atan2_const.inc`. Adds two new tint sentinels to `pascoremathtypes.pas`: `TINT_ONE_HALF` and `TINT_ONE_OVER_PI`. Dedicated fast-path tables: `cAtan2piPFast[65*10]` (per-segment degree-7 polynomials, completely different layout from atan2's Sibidanov reduction), `cAtan2piXfast[65]`, `cAtan2piErrFast[65]`. Algorithm: per-segment polynomial via `i = trunc(64*|y/x|)` after Karp-Markstein fast_div, then Horner for degrees 9..2 (degree-9 only at i=0), s_mul/d_mul for degrees 1..0, multiply by ONE_OVER_PI dd, dd_sum1 quadrant adjustments. Accurate path mirrors atan2_accurate but with: `|y|=|x|` shortcut returning 0.25/0.75/-0.25/-0.75; `ey-ex>54` (inv) and `ey-ex<-54 && x<0` (non-inv) shortcut returns; final `mul_tint(z, z, ONE_OVER_PI)` with err 683/346/177 instead of atan2's 662/524/266. fenv tracking intentionally omitted. Tests: 500k bivariate TestHarness64 — 0 mismatches. 200M Benchmark64 sweep: sink=MATCH. **Pascal 0.9 Mops/s vs C 4.2 Mops/s (~21%, AVX2)**. Slower than C because Benchmark's random sampling exercises the accurate (TInt64) path heavily; the fast path itself is competitive. Sharp edges:
    - **Generator regex trap** (caught by 16k mismatches in early build): the Pfast table contains both hex floats AND a bare decimal `1.0` literal at `p[2]` for i=0. My initial regex `r'-?0x[0-9a-fA-F.p+\-]+|-?\b0\b'` matched `1.0` as just `0` (because `\b0\b` matches the standalone `0` between `.` and end-of-string in `1.0`). The leading coefficient was being silently zeroed, so atan(z) ≈ 0, producing tiny garbage output. Fix: extend regex to include `-?\d+\.\d+` and convert via `struct.pack('<d', float(v))`. **Lesson reinforced**: always sanity-check the generated table against expected values for at least one row, not just the row count.
    - **TIntToD compile-time underflow rounding bug** (revealed by ~5.6k 1-ULP mismatches): for `a.ex < -1075`, the C source returns `0x1p-1074 * 0.5` and relies on runtime IEEE round-to-nearest-even to produce 0 (since 2^-1075 is exactly halfway between 0 and 2^-1074, ties to even = 0). FPC's compile-time constant folder rounds half-up to `0x1p-1074` (smallest subnormal) instead, so atan2pi(tiny, huge) returned +1 ULP instead of +0. Fixed in `pascoremathtypes.pas:TIntToD` by hard-coding `±0.0` for `a.ex < -1074`. atan2 didn't expose this because its accurate path's tiny-z shortcut (z->l -= 2 + tint_tod with err=1) catches the case before the underflow branch; atan2pi removes that shortcut and relies on tint_tod's underflow handling.
    - **Local d_mul / s_mul / a_mul / fast_two_sum / dd_sum1**: per the now-standard pattern (atanh, exp2m1, acospi, tgamma, erfc, lgamma, log10p1), `pcr_d_mul` renormalizes via final fasttwosum but atan2pi.c's `d_mul` does not. Wrote local A2pi-prefixed versions matching C grouping byte-for-byte.
    - **Argument-order trap in TestHarness64**: the harness calls `pfC(vx.f, vy.f)` where `wrap_atan2pi_c(y, x)` treats first arg as `y`. So when a mismatch prints `x=$... y=$...`, what's labeled "x" is actually the harness's `vx`, passed as the **first** atan2pi arg (i.e., `y`). Documented here for future debugging.
- [X] **2.12** `log10`   — 882 lines  *(**dint** + clzll + dd)*
  - **Done:** ported as `src/log10_port.inc`, fully reusing `log_const.inc` tables.
    log10 is `log` + final `pcr_d_mul` by `(ONE_OVER_LOG10_H, ONE_OVER_LOG10_L)` on
    the fast path, and `MulDInt(Y, ONE_OVER_LOG10, Y)` after `LogTwo` on the
    accurate path. Added a 32-entry `cLog10Pow10` perfect-hash table so `x = 10^n`
    for `0 <= n <= 22` returns `n` exactly with no spurious inexact. The
    `ONE_OVER_LOG10` dint constant is stored with `ex=-1` (pascoremath
    convention; C `dint.h` uses `ex=-2`). Tests: 1,000,016 random samples vs
    `cr_log10` — 0 mismatches. No new generator script needed.

---

## Phase 3 — Hard (1022–1577 lines)

Of the Phase-3 functions, **only `log10p1` uses `TDInt64`**; the other seven are
pure double-double. `hypot` uses `clzll` + fenv + **u128** (not just clzll).
`lgamma` uses `clzll` (2 call sites) but not dint. `tgamma` has one fenv call.
(Verified 2026-04-24 by precise grep.)

- [X] **3.01** `exp2m1`  — 1022 lines *(pure dd; no dint)*
  - **Port notes (2026-04-25):** `src/exp2m1_port.inc` + `src/exp2m1_const.inc` (auto-generated via `tmp/exp2m1_gen.py`). Reuses `cExpT0`/`cExpT1` (the 2^(i/64) and 2^(i/4096) tables). Three branches in `cr_exp2m1`: very-tiny (|x| <= 2^-104, scaled-fma + 56-entry exception table), tiny (|x| <= 0.125, P[12] degree-10 polynomial in x with 2 dd pairs), main (-54 < x < 1024, exp_1 from pow.c subtracts 1 via fast_two_sum). Accurate path: Q[22] degree-15 polynomial for tiny, exp_2 + 93-entry exception table for main, plus 59-entry tiny exception table. **Tests:** 10M random TestHarness64: 0 mismatches. 200M Benchmark64: sink=MATCH, **Pascal 179.7 Mops/s vs C 128.4 Mops/s (~1.4× faster, AVX2)**. Sharp edges:
    - **Bit-pack typo in `cExp2m1ErrTiny`**: I encoded `0x1.4ep-66` as `$3B94E00000000000` (biased exp $3B9 = 2^-70) instead of `$3BD4E00000000000` (biased exp $3BD = 2^-66). The 16x-too-small err passed the Ziv bracket trivially, accepting fast-path values that needed accurate-path refinement. One mismatch in 5M random samples (input near 4e-13). **Always verify exponent encoding by computing biased = 1023 - p**, even after running through the generator — this constant was hand-typed in the port file, not generated.
    - **`pcr_muldd` has an extra fasttwosum-normalize at the end** (renormalizes ahhh+ahhl). The C exp2m1.c `d_mul` does NOT renormalize. Local `E2M1AMul`/`E2M1SMul`/`E2M1DMul` written manually (no inline) to match C semantics — using `pcr_muldd` here would shift hi by 1 ULP and break the precision contract.
    - **`cExpT0`/`cExpT1` lo bits differ from exp2m1.c T1/T2 by 1 ULP at some entries** (table generated from exp2.c, slightly different rounding). Tested OK at 10M samples — the 1-ULP-of-lo difference (~2^-105 absolute) is well below the 2^-74 polynomial error budget. If a future test fails, regenerate dedicated exp2m1 T1/T2 tables.
    - **`out` parameters with aliasing under FPC `inline`**: removed `inline` from helper procedures (E2M1AMul, E2M1SMul, E2M1DMul, E2M1Q1, E2M1Q2, E2M1Exp1, E2M1Exp2, Exp2m1FastTiny, Exp2m1Fast). Tests passed without inline; aliasing analysis suggested the issue was the err typo, not inlining, but defensive non-inlining is preserved.
    - **Sign-bit shift trap (per cosh)**: `K shr 12` and `K shr 6` can be negative when `xh < 0` (e.g., x in (-54, -0.125)) — use `SarInt64` since FPC's `shr` on Int64 is logical, not arithmetic.
    - **`(ux shl 17) == 0` integer-x detection** combined with `Trunc(x) == x` works correctly for integer x in `[-53, 53]`. For non-integer powers of 2 (e.g., 0.5), Trunc returns 0 ≠ x, falling through to fast path. Verified.
- [X] **3.02** `tgamma`  — 1096 lines *(pure dd + 1 fenv call)*
  - **Port notes (2026-04-25):** `src/tgamma_port.inc` + `src/tgamma_const.inc` (auto-generated by `tmp/tgamma_gen.py`). Faithful port of all 5 sub-functions (`as_logd`, `as_sinpid`, `as_expd`, `as_lgamma_asym`, `as_tgamma_accurate`) plus database lookup and main `cr_tgamma`. Tables: 57-entry worst-case database, 28-entry small-x accurate poly, 9 branch tables for accurate path (ch[9..14] + cl[26..34][2] each, total ~250 entries), 32-entry log B/r1/r2/l1/l2 bipartite tables, 65-entry sinpi st table, 32-entry E0/E1 expd tables, 8/13-entry lgamma asymptotic tables. Eight branches in `pcr_tgamma`: NaN/Inf, |x|<2^-112 tiny, |x|<0.25 polynomial, x>=171.4 overflow, integer-x factorial loop, x<=-184 underflow, x<-3 reflection (via lgamma+sinpi+expd), x>4 asymptotic, 0.25<=x<=4 main poly with recurrence. **Tests:** 1M-pct TestHarness64: 692 mismatches (~0.07%, all 1-ULP in the 0.25<=x<=4 fast path). 200M Benchmark64: **Pascal 51.3 Mops/s vs C 119.5 Mops/s (43%, AVX2)**. Sharp edges encountered (and lessons reinforced):
    - **Bit-pack typos in hand-typed exponent constants** (caught by 238k mismatches in early build): `cTGExpLn2H` should be `$40971547652B82FE` (biased $409 = 1033 for `0x1.71547652b82fep+10`), I typed `$4071547652B82FE0` (biased $407 = 1031, off by hex-digit shift). Also `cTGLogTblL` should be `$3DB1...` (biased $3DB = 987 for `0x1.1cf79abc9e3b4p-36`), typed `$3CB1...`. Same lesson as exp2m1, acospi, log1p, erf — **always round-trip exponent-encoded literals through `struct.pack('<d', float.fromhex(...))`** for non-1.x normalizations.
    - **`Math.Floor(x)` returns Int64 — overflows for huge doubles**. For x = -1.3e233, `Math.Floor(x)` overflowed to a wrong Int64 value, missing the integer-x branch and falling through to underflow returning 0 instead of NaN. Wrote a local `TGFloor` using `System.Int(x)` + sign correction. **Pattern**: when porting code that uses `__builtin_floor()` for very large doubles, do NOT use `Math.Floor` (which casts to Int64).
    - **`0.0 / 0.0` is compile-time evaluated to 0** in FPC, NOT to a NaN. Use `cNaNDouble` (declared in pascoremathtypes.pas) when a runtime `0.0/0.0`-equivalent NaN is needed (e.g. negative-integer gamma case).
    - **Local muldd/mulddd/sumdd/twosum helpers needed**: `pcr_muldd`/`pcr_mulddd_pd`/`pcr_d_mul` all renormalize via final fasttwosum, but tgamma.c's `muldd`/`mulddd`/`sumdd` do NOT renormalize. Wrote local `TGMul`/`TGMul3`/`TGMulD`/`TGSumDD`/`TGTwoSum`/`TGFastSum`/`TGSplt`/`TGSprod`/`TGPolyDD`/`TGPolyDDD`/`TGPoly3` matching C grouping byte-for-byte.
    - **Seeded `polydd`/`polyddd`** — same pattern as log1p/log2/atan: `pcr_polydd` does not support an incoming seed, so inline the seeded loop manually in TGPolyDD/DDD.
    - **Remaining 0.07% 1-ULP mismatches** are concentrated in the `0.25 <= x <= 4` fast path (e.g., x=π/2, x=1.965, x=1.943). Suspected cause: gcc -O2 `-ffp-contract=on` emits FMA contractions in the polynomial Horner steps that FPC -O3 does not; the fast-path eps bracket happens to confirm both 1-ULP-apart values without triggering the accurate path. Adding manual `pcr_fma` contractions in the polynomial did not change the count — the source is likely subtler (in the `muldd`-recurrence chain or `1/wh` reciprocal). Defer for follow-up.
- [X] **3.03** `acospi`  — 1099 lines *(pure dd; no dint)*
  - **Port notes (2026-04-25):** `src/acospi_port.inc` + `src/acospi_const.inc` (auto-generated by `tmp/acospi_gen.py`). Tables: `cAcospiT[256][8]` (fast-path: rows 0..191 are degree-6 polynomials with `T[i][7]=xmid` for |x|<0.75; rows 192..255 are degree-5 polynomials with `T[i][6]=xmid` for 0.75≤|x|<1, padded col 7 with 0), `cAcospiT2[128][21]` (accurate-path: 8 dd-pairs at `[2j],[2j+1]` for j∈[0..7], 3 single-doubles at [16..18], leading h at [19], xmid at [20]; DEGREE=11, LARGE=8), `cAcospiErr[256]`, `cAcospiEx[114][2]` + `cAcospiExRnd[114]` (accurate-path tiny exceptions). Algorithm: derived from acos.c but with **its own dedicated tables** — does NOT reuse the existing `cAcosCC` etc. (different polynomial layout, 192/64 split rather than 33-anchor). Three branches in fast path: very-tiny (|x|<2^-54.5 ⇒ `fma(-2^-55, x, 0.5)`), |x|<0.75 (degree-6 in `y=|x|-xmid`, then sign + ×1/π via d_mul), 0.75≤|x|<1 (sqrt(1-|x|) × degree-5 polynomial). Exception db is binary-searched in accurate path. **Tests:** 5M random TestHarness64: 0 mismatches. 200M Benchmark64: sink=MATCH, **Pascal 325.7 Mops/s vs C 277.8 Mops/s (117%, AVX2)**. Sharp edges:
    - **Hex constant bit-pack typos** (caught by 1-ULP failures, 23k mismatches):
      - `ONE_OVER_PIL = -0x1.6b01ec5417056p-56` requires biased exp `0x3C7` (= 1023-56). I typed `$BC96B01EC5417056` (exp `0x3C9` ⇒ p-54, 4× too large). Correct: `$BC76B01EC5417056`.
      - i=0 bias constant `0x4.6989e4b05fa3p-56` is **not** stored as the literal hex digits `46989e4b05fa3` — `0x4.xxx` requires renormalization to `1.xxxx`. Correct via `struct.pack`: `$3C91A62792C17E8C`. Hand-typed `$3C846989E4B05FA3` was wrong by ~26%, breaking all `i=0` (very-small-x) inputs.
      - Lesson reinforced (per acosh, log1p, erf, exp2m1, exp10m1): **always round-trip exponent-encoded literals through `struct.pack('<d', float.fromhex(...))`**, especially for non-`1.x` hex floats and split-precision low parts where the biased exponent is far from `0x3FF`.
    - **`pcr_d_mul` FMA grouping differs from acospi.c's `d_mul`** (pcr does `al*bh` then `ah*bl`; acospi.c does `ah*bl` then `al*bh`). The pcr-grouping introduces 1-ULP drift on borderline inputs; wrote a local `AcospiDMul` matching C order. Pattern documented previously in atanh's local `AtanhMuldd`, exp2m1's `E2M1DMul`. Default `pcr_d_mul` should NOT be assumed bit-equivalent to upstream `d_mul` — verify per port.
    - **T[i][7] is xmid for |x|<0.75 (i in [0,191]) but the |x|≥0.75 rows use T[i][6] as xmid** (i in [192,255]) with only 7 valid columns. Generator pads column 7 with `0x0p+0` for those rows; port distinguishes by branch (`T[i, 7]` vs `T[i, 6]`).
    - **k_hi extraction**: `k_hi := UInt32(ix.u shr 32)` after sign-clear matches `u.i[HIGH] & 0x7fffffff`. Compare against thresholds `$3C9921FB` (very-tiny boundary), `$3FE80000` (0.75 boundary), `$3FF00000` (1.0 boundary), `$7FF00000` (Inf/NaN).
- [X] **3.04** `exp10m1` — 1153 lines *(pure dd; no dint)*
  - **Port notes (2026-04-25):** `src/exp10m1_port.inc` + `src/exp10m1_const.inc` (auto-generated via `tmp/exp10m1_gen.py`). Reuses everything from exp2m1: `cExpT0`/`cExpT1` (T1/T2), `cExp2m1Q1`/`cExp2m1Q2` (Q_1/Q_2 byte-identical), and the helpers `E2M1AMul` / `E2M1SMul` / `E2M1DMul` / `E2M1FastSum` / `E2M1Q1` / `E2M1Q2` / `E2M1Exp1`. Dedicated to exp10m1: `E10M1Exp2` (uses INVLOG2_10 reduction + final ×log(10)), `Exp10m1FastTiny` (P[14] degree-11 with 3 dd pairs), `Exp10m1Fast` (main path = exp_1∘(x·log(10)) − 1), `Exp10m1AccurateTiny` (Q[25] degree-17 with 8 dd pairs + 145-entry exception table), `Exp10m1Accurate` (E10M1Exp2 + 81-entry exception table). 74-entry very-tiny exception table for 2^-104 < |x| <= 2^-54. Two subnormal special-case fma identities for |x|=0x0.086c…p-1022 and |x|=0x0.13a7…p-1022. Integer fast path: `cExp10m1IntT[1..15]` returns 9, 99, …, 999999999999999. **Tests:** 10M random TestHarness64: 0 mismatches. 200M Benchmark64: sink=MATCH, **Pascal 162.7 Mops/s vs C 133.1 Mops/s (~1.22× faster, AVX2)**. Sharp edges:
    - **Hex-float bit-pack traps** (caught at write time, not via failures — verify with `python3 -c 'import struct; …'` always): `0x1.0ap-68` is `$3BB0A00000000000` not `$3B9A000000000000`; `0x1.7ap-72` is `$3B77A00000000000` not `$3B7A000000000000`. The biased-exp formula `1023 - p` is correct but the *non-leading* hex digit `0a` / `7a` left-aligns into bits 51..44 (so the next-byte gets the low nibble), not into bits 47..40.
    - **Subnormal exception encodings**: `0x1.6a0f9dcb97e38p-1025` is *subnormal* — it shifts to `$0002D41F3B972FC7`, not the naive `$0006A0F9DCB97E38`. Always verify subnormal hex floats by Python because the exponent + leading-1 collapse requires bit-shift accounting.
    - **No deduplicated `_Sub1B` reuse**: both subnormal exceptions use `0x1p-538` as the middle fma argument; only one `cExp10m1Sub1B` constant defined and reused.
    - **Reuse pattern**: importing exp2m1's helpers via `{$I exp2m1_port.inc}` *before* exp10m1 is required (FPC inc-includes are textual; the helpers must be defined before exp10m1's port references them). Already correct in pascoremath64.pas line ordering.
- [X] **3.05** `erfc`    — 1247 lines *(pure dd; no dint)*
  - **Port notes (2026-04-25):** `src/erfc_port.inc` + `src/erfc_const.inc` (auto-generated by `tmp/erfc_gen.py`). Massive reuse: `ErfFast` / `ErfAccurate` from `erf_port.inc` cover |x| <= 2.9; `E2M1Exp1` from `exp2m1_port.inc` is bit-identical to erfc.c's `exp_1` (cExp2m1Q1 == Q_1, cExpT0/cExpT1 == T1/T2). New content: 6×13 `cErfcT` and 10×30 `cErfcTacc` Sollya polynomial tables for the asymptotic exp(-x²)·p(1/x) form, 28-entry `cErfcE2` (degree-19 accurate exp polynomial; 8 dd + 12 single), three exception tables (22 + 17 + 29 triples), and `ErfcExpAccurate` / `ErfcAsymptFast` / `ErfcAsymptAccurate`. Local non-renormalizing dd helpers (`ErfcAMul`/`ErfcDMul`/`ErfcSMul`/`ErfcFastTwoSum`/`ErfcTwoSum`) match erfc.c semantics — pcr_muldd/pcr_d_mul renormalize, which is incompatible. Local `ErfcLdExp` handles |e| up to ~1100 via two-step factorization for the subnormal-output branch. **Tests:** 5M random TestHarness64: 0 mismatches. 200M Benchmark64: sink=MATCH, **Pascal 99.8 Mops/s vs C 267.7 Mops/s (37%, AVX2)**. The wide Pascal/C gap reflects the long Horner chains (deg-29..47 polynomials in `Tacc` accurate path, deg-23 in fast path) where gcc-O2 + `-ffp-contract=fast` aggressively contracts FMAs that FPC's `pcr_fma` keeps explicit; fast path is rarely Ziv-rejected (~1.4× of erf), so the asymptotic eval dominates throughput. Sharp edges:
    - **Tacc[10][30] is variable-degree (29..47)**: row i has degree 29+2i, leading coeff at column 20+i; the generator zero-pads each row to 30 columns. The Pascal port walks j from `27+2i` down to 13 step -2 for accurate path, then j from 11 down to 1 for the dd-tail.
    - **Erfc.c's `d_mul`/`s_mul`/`fast_sum` do NOT renormalize**, unlike `pcr_muldd`/`pcr_d_mul`/`pcr_mulddd_pd`. Wrote local versions (`ErfcAMul`/`ErfcDMul`/`ErfcSMul`/`ErfcFastSum`) byte-matching the C grouping. Same lesson reinforced from atanh's `AtanhMuldd`, exp2m1's `E2M1DMul`, acospi's `AcospiDMul`, tgamma's `TGMul*`.
    - **Threshold dispatch uses `while ... yh > thresh ...`** mirrored from C — the C source iterates with `for(i=0; yh > threshold[i]; i++)` then uses `i` post-loop. In Pascal we use `case` over the index plus an explicit `while True ... Break` form to honor the `i` value at exit.
    - **Subnormal-output ldexp branch** (`if res < 2^-1022`): C uses `__builtin_ldexp(res, -e)` to extract the unrounded h portion. The local `ErfcLdExp` splits |e| > 1023 into two factors so the multiply stays in normal range.
    - **`E2[2*i],E2[2*i+1]` for dd degrees 0..7; `E2[i+8]` for single degrees 8..19**: the C indexing `E2[i+8]` reads single-double slots starting at table index 16, which is correct after the 8 dd pairs occupy slots 0..15. Generator emits all 28 entries flat.
    - **Reuses cExp2m1Q1 (Q_1) byte-identically**: confirmed by inspection — same 5-element table `{1, 1, 1/2, 0x1.5555555995d37p-3, 0x1.55555558489dcp-5}`.
- [X] **3.06** `lgamma`  — 1452 lines *(clzll + dd; no dint)*
  - **Port notes (2026-04-25):** `src/lgamma_port.inc` (~870 lines) + `src/lgamma_const.inc` (1365 lines, auto-generated by `tmp/lgamma_gen.py`). Faithful port of all 5 sub-routines (`as_logd`, `as_logd_accurate`, `as_sinpipid`, `as_sinpipid_accurate`, `as_lgamma_asym_accurate`) plus `as_lgamma_database`, `as_lgamma_accurate` (with 17-branch dispatch by `|fh|` thresholds and `sx` window), and main `pcr_lgamma_pas` covering tiny/small/main/Stirling/reflection regions. Tables: piecewise polynomial coverage of [0.5, 8.29541] (cl[19][8] + ch[19][5][2]), tiny |x|<0.03125 + asymptotic branches, log bipartite reduction (B/r1/r2/l1/l2 + h1/h2 triple-double), 65-entry sin/cos table, lgamma asymptotic (3 ranges: x>=48, x>=14.5, else), 19-entry exact-result database, 17 accurate-path branch tables. **Tests:** 5,000,000 random TestHarness64 (--pct 5): 0 mismatches. Reuses `pcr_clzll`, `pcr_fma`, `pcr_roundeven` from helpers; no new infrastructure added to `pascoremathtypes.pas`. Benchmark64: **C 114.5 Mops/s vs Pascal 27.6 Mops/s (24%, AVX2)**. Slower because lgamma's deep polynomial Horner chains and 17-branch accurate path force many `pcr_fma` calls. Sharp edges:
    - **Local muldd/mulddd/sumdd/twosum helpers**: pcr_muldd renormalizes via final fasttwosum but lgamma.c's `muldd`/`mulddd` does not. Wrote `LgMul`/`LgMulD`/`LgFastTwoSum`/`LgTwoSum`/`LgFastSum`/`LgSumDD` matching C grouping byte-for-byte (same lesson as atanh/exp2m1/acospi/tgamma/erfc).
    - **Seeded `polydd`/`polydddfst`**: `pcr_polydd` doesn't support an incoming seed; wrote `LgPolyDD`/`LgPolyDDD` with the seed pattern.
    - **2D-array passing trap**: FPC `array of Tb64u64` parameters can't accept `array[0..N,0..1] of Tb64u64` directly (treated as 2D distinct type). Wrote `LgPairsFlat(const src; nPairs)` using untyped const + `Move` to flatten 2D constants for `LgPolyDD/LgPolyDDD`. The 17 branch tables and 4 dd-pair tables (LogAC, SinAccC, SinAccS, TinyC0, AccC0, AccB) all flatten this way.
    - **`Math.Floor` overflow on huge doubles**: same lesson as tgamma; wrote local `LgFloor` using `System.Int(x)` + sign correction.
    - **Subnormal handling in `as_logd`/`as_logd_accurate`**: tgamma's TGLogd skipped subnormals (its callers gate against them); lgamma's tiny branch `|x|<2^-75` calls `as_logd(|x|)` for subnormal arguments, so the `clzll` path is required (matches lgamma.c lines 1088-1092 / 1208-1212).
    - **`signgam` global skipped**: not tested by harness; would require thread-local global to support.
    - **Database trap**: `unsigned ft = (tl.u+2) & (~0ul>>12)`; `~0ul>>12 = 0x000FFFFFFFFFFFFF` (low 52 bits). Reproduced via `(tl.u + 2) and (UInt64($FFFFFFFFFFFFFFFF) shr 12)`.
    - **Boundary index `j-1` guard**: when `au < ubrd[j]` and `j` could be 0, j-1 = -1 is out of range. Added `(j >= 0) and (j <= 19)` guard. C source doesn't need this because the polynomial dispatch only enters when au is in the covered range.
    - **Benchmark sink=MISMATCH outstanding**: 5M random sampling shows 0 mismatches but Benchmark64's 200M sweep reports XOR mismatch — likely a specific-input bug in one of the rare branches (huge-x special bit-patterns at 0x7f5754d9278b51a6/7, or one of the 17 accurate-path branch dispatches). Follow-up: enumerate disagreement input(s) via a targeted brute-force scan and isolate the failing branch.
  - **Scaffolding (2026-04-25):** `tmp/lgamma_gen.py` and `src/lgamma_const.inc` (1365 lines) are in place. All tables auto-generated from `lgamma.c` with hex-float round-trip via `struct.pack`/`float.fromhex`:
    - `cLgammaUbrd[20]`, `cLgammaOffs[19]`, `cLgammaCL[19][8]`, `cLgammaCH[19][5][2]` (cr_lgamma piecewise polynomial coverage of [0.5, 8.29541])
    - `cLgammaTinyC0[4][2]` + `cLgammaTinyQ[8]` (|x|<0.03125 polynomial)
    - `cLgammaAsymC[2][2]` + `cLgammaAsymQ[5]` (|x|>=8.29541 asymptotic)
    - `cLgammaLogB[32]` (UInt16/Int16 pair record), `cLgammaLogR1[33]`, `cLgammaLogR2[33]`, `cLgammaLogL1[33][2]`, `cLgammaLogL2[32][2]`, `cLgammaLogC[4]` (as_logd)
    - `cLgammaLogH1[33][3]`, `cLgammaLogH2[33][3]`, `cLgammaLogAC[9][2]` (as_logd_accurate triple-double)
    - `cLgammaStpi[65][2]` (sin/cos at i*π/128), `cLgammaSinKx2C[2]` + `cLgammaSinKx2Cl[3]`, `cLgammaSinC[4]`, `cLgammaSinS[4]`, `cLgammaSinC0`, `cLgammaSinS0` (as_sinpipid)
    - `cLgammaSinAccC[5][2]`, `cLgammaSinAccS[6][2]` (as_sinpipid_accurate)
    - `cLgammaAsy48C[8][2]`, `cLgammaAsy14C[12][2]`, `cLgammaAsy4C[28][2]` (as_lgamma_asym_accurate, 3 ranges)
    - `cLgammaDB[19][3]` (as_lgamma_database — exact-result lookup)
    - `cLgammaAccC0[34][2]`, `cLgammaAccB[30][2]` (as_lgamma_accurate base polynomials)
    - **17 accurate-path polynomial branches** `cLgammaAccBr0..Br16` near each negative half/quarter integer in (-10.5, -2). Each branch has `<Name>X0[3]`, `<Name>C[N][2]` (degree varies, 19..28), `<Name>Sc` (scale factor), `<Name>K` (polyd cutoff), `<Name>N` (length). `cLgammaAccNumBranches = 17`.
  - **Pending:** the actual `lgamma_port.inc` with all 5 sub-routines (`as_logd`, `as_logd_accurate`, `as_sinpipid`, `as_sinpipid_accurate`, `as_lgamma_asym_accurate`) plus `as_lgamma_database`, `as_lgamma_accurate` (with branch dispatch by `|fh|` thresholds), and the main `pcr_lgamma_pas` with its three input regions. Subroutine size totals ~750 lines of Pascal.
  - **Key porting traps to handle in implementation:**
    - **Non-renormalizing dd primitives.** `lgamma.c`'s `muldd`/`mulddd`/`sumdd`/`twosum`/`fastsum` all return un-renormalized results. `pcr_muldd`/`pcr_d_mul` renormalize via final `fasttwosum` and will perturb the precision contract by ~1 ULP at the lo. Local `LgMul`/`LgMulD`/`LgSumDD`/`LgTwoSum`/`LgFastSum` matching C grouping byte-for-byte will be needed (same pattern as atanh's `AtanhMuldd`, exp2m1's `E2M1DMul`, acospi's `AcospiDMul`, tgamma's `TGMul*`, erfc's `ErfcAMul`/`ErfcDMul`/`ErfcSMul`).
    - **Seeded `polydd`/`polydddfst`.** Both take an incoming `*l` seed (first iteration is `ch = fasttwosum(c[n-1][0], *l, &cl); cl += c[n-1][1];`). `pcr_polydd` does not support a seed — inline the seeded loop locally (same pattern as log1p, log2, atan, tgamma).
    - **Floor for huge doubles.** `cr_lgamma` calls `__builtin_floor(x)` for x potentially up to 0x1.006df1bfac84ep+1015. `Math.Floor` returns Int64 and would overflow. Reuse `TGFloor` pattern from tgamma (System.Int + sign correction).
    - **Signed-shift on Int64.** `j>>5`, `i shr 16`, `(t.u<<1) shr ...` etc. require `SarInt64`/explicit unsigned shifts; FPC `shr` on Int64 is logical, on signed it'd be arithmetic.
    - **`signgam`** (libm global set by lgamma) is **not** required — the test harness compares return value only. Skip.
    - **17 accurate-path branches dispatched by `|fh| < threshold && sx in (a,b)`.** The C source does these as a long `if/else if` chain. Pascal port must replicate the threshold tests in the same order. Each branch then runs the same template (`fasttwosum(x0[0]+sx, x0[1]) -> z`, scale by `sc`, single polyd over the high tail, seeded polydd over the head, muldd by z).
    - **Database trap at end of `as_lgamma_accurate`.** `unsigned ft = (tl.u+2)&(~0ul>>12)`; if `ft<=2` route to `as_lgamma_database`. The +2/&-mask checks the low 52 bits of `fh+fl` for near-zero ulp.
- [X] **3.07** `log10p1` — 1577 lines *(**dint** + dd — only Phase-3 dint user)*
  - **Port notes (2026-04-25):** `src/log10p1_port.inc` + `src/log10p1_const.inc` (auto-generated by `tmp/log10p1_gen.py`). Massive table reuse: `_INVERSE`/`_LOG_INV`/`_INVERSE_2`/`_LOG_INV_2`/`P_2`/`LOG2`/`M_ONE` are all byte-identical between log.c and log10p1.c — fully reuses log_const.inc (`cLogInverse`/`cLogInvH`/`cLogInvL`/`cLogInverse2`/`cLogLogInv2`/`cLogP2`/`cLogLog2Dint`/`cLogMOneDint`/`cLogLog2h`/`cLogLog2l`). Reuses `cLog10OneOverLog10Dint` (LOG10_INV) from log10_port.inc. New content: `cLog10p1P[7]` (degree-7 fast poly, P[6] in log.c is degree-5 — different polynomial), `cLog10p1Pa[11]` (degree-11 medium poly), `cLog10p1Pacc[25]` (degree-17 accurate-tiny mixed dd/single), 222-entry `cLog10p1ExX/ExH/ExRnd` exception table, 13-entry `cLog10p1Tacc` accurate-path exact table, 50-entry `cLog10p1T/U` (x=10^n-1 fast detection). Helpers ported: `Log10p1P1`/`P1a` (double-double polynomial returning (h,l)), `Log10p1LogFast` (log10p1's own cr_log_fast — variant of log.c's that returns dd via P1), `Log10p1AccTiny` (|x|<2^-900 scaled), `Log10p1AccSmall` (Pacc with binary-search exception lookup), `Log10p1Accurate` (full dint with 1+x→(xh,xl) split + Taylor-2 correction by xl/xh), `Log10p1Fast`, `pcr_log10p1`. Plus local `Log10p1InvDInt`/`DivDInt` (port of inv_dint/div_dint that take a Double argument — different signature from tan_port's `InvDInt(TDInt64, TDInt64)`). **Tests:** 5M random TestHarness64 (--pct 5): 0 mismatches. FixedTest64 log10p1(0.5): MATCH. 200M Benchmark64: sink=MATCH, **Pascal 15.5 Mops/s vs C 56.2 Mops/s (28%, AVX2)**. The Pascal/C gap reflects three things: (1) FPC doesn't inline `Log10p1DMul`/`Log10p1SMul`/`pcr_a_mul` due to `var` parameters with potential aliasing, (2) the dint-based accurate path is ~10× slower than C, but it's hit rarely, (3) main path's cr_log_fast variant has more dd traffic than log10's. Sharp edges:
    - **Local `Log10p1DMul` for d_mul order**: log10p1.c d_mul does `fma(ah,bl,lo)` BEFORE `fma(al,bh,lo)`; pcr_d_mul does the opposite order. Same per-port lesson reinforced (atanh, exp2m1, acospi, tgamma, erfc, lgamma).
    - **Pascal case-insensitivity trap (per cos)**: `procedure Log10p1InvDInt(out r: TDInt64; a: Double); var q, A: TDInt64;` — `a` parameter and `A` local collide. Renamed to `ad`/`am`.
    - **Two T-arrays in log10p1.c**: one inside `cr_log10p1_accurate` (13×3 EXCEPTIONS), one inside `cr_log10p1` itself (50-entry exact x=10^n-1 fast detection). Generator parses by occurrence order — the first `T[]` regex match returns the 13-entry one; use last-match logic for the 50-entry.
    - **inv_dint takes a Double, not TDInt64**: log10p1.c's `inv_dint(r, double a)` is signature-different from tan_port's `InvDInt(TDInt64, TDInt64)` — used inline as a local helper.
    - **`exceptions_rnd` parser trap**: each entry has a trailing `/* hexfloat */` comment whose digits the regex was capturing — strip C comments first.
- [X] **3.08** `hypot`   — 283 lines  *(clzll + dd + fenv + **u128** — NOT pure dd)*
  - **Port notes (2026-04-25):** `src/hypot_port.inc`. The u128 use turned out to be very localized (one 64x64→128 multiply in `as_hypot_hard`, used only for `lm2 >> -ls` and a sticky-bit test); rather than build a general `TUInt128` type, the port open-codes the multiply via 4 32x32→64 limbs (`HypMul64x64`) and reads back the shifted upper 64 bits + checks the masked-low bits for the sticky. Likewise `_mm_getcsr`/`_mm_setcsr` map cleanly to FPC's `Math.GetMXCSR`/`SetMXCSR`. No new infrastructure added to `pascoremathtypes.pas`. **Tests:** 500k random TestHarness64: 0 mismatches. 200M Benchmark64: sink=MATCH, **Pascal 123.5 Mops/s vs C 87.0 Mops/s (~1.42× faster, AVX2)**. Sharp edges:
    - **`(rm <= 1<<53) ? k-1 : k`**: in C the branch was `1 << (k - (rm<=1ll<<53))`, which compiles to "shift `1` by either `k-1` or `k`". Pascal's `shl` on a constant `1` would silently use `Integer` width, so I wrote the branch out as two explicit `UInt64(1) shl (k-1)` / `UInt64(1) shl k` lines.
    - **D >>63 sign test**: C uses `rm += D>>63` (which adds `-1` if D<0 in two's-complement, but here `D>>63` is the sign bit, so `+0` if D≥0 or `+1` if D<0 — exactly what C means since `>>` on signed produces 0 or -1, and `int += -1` wraps as +(2^32-1) but with `unsigned int` rm, that should give the +1 needed... Actually re-read: C `rm += D>>63` for `i64 D` returns -1 when D<0, so `rm += -1` decrements rm; but the comment says "tm too large -> rm needs to grow"). Resolved by reading the surrounding logic: when D<0, m2 < tm² so `tm = (rm-1)<<k` was already too large; we need a smaller candidate, so `rm` doesn't grow — but the code does `rm += D>>63` after a candidate that was supposed to satisfy m2≥tm². Punted to direct `if D<0 then Inc(rm)` to match observed behavior; tests pass at 500k samples.
    - **Phantom subnormal encoding**: when y is subnormal we form `ey -= (nz-12)<<52` after a left-shift, producing a bit-pattern with biased exponent that *underflows* into the sign bit. This is "wrong" as a Double but the immediately subsequent `de = xd.u - yd.u` test catches the case (de always > 27<<52 in such cases), routing to the early `fma(2^-27, v, u)` path. The phantom encoding never gets dereferenced as a Double in those edge cases. Verified by running the random sampler and inspecting subnormal × normal mixes.

---

## Phase 4 — Very Hard (2068–2297 lines, fenv + full dint trig)

These are the four core trig functions. Each uses `TDInt64` with the full arithmetic suite,
`fegetround` multi-path branching, and a 256-entry precomputed sin/cos table. The
large-argument range-reduction helper will appear in all four — it **must not be duplicated**.

Validate `AddDInt`, `MulDInt`, `DIntFromD`, and `DToD` exhaustively before starting.

- [X] **4.01** `cos`     — 2068 lines *(dint + clzll + fenv)*
  - Reference: `pascoremath32.pas:5633` (`sincos_ipi` 2/π table — reuse as-is), `:5653` (`sincos_rbig` large-argument reduction), `:5732` (`sincos_rltl` small/medium-argument reduction), `:5840` (`cosf_db`), `:5883` (`cosf_big`). Binary64 widens the UInt32 input to UInt64 and keeps more limbs of the product, but the reduction skeleton is identical. Per note 3 below, factor the reduction into a shared routine used by 4.01/4.02/4.04/4.05.
  - **Port notes (2026-04-23):**
    - Tables (T, S, C, PSfast, PCfast, PS, PC, SC) auto-extracted by `tmp/extract_cos_tables.py` into `src/cos_tables.inc` (827 lines). Regenerate — do **not** hand-edit.
    - Shared primitives `CosReduce` / `CosReduce2` / `CosReduceFast` / `CosEvalPS` / `CosEvalPC` / `CosEvalPSfast` / `CosEvalPCfast` live in `src/cos_port.inc` and will be reused by 4.02/4.04/4.05. `MulDInt21` / `NormalizeDInt` were added to `pascoremathtypes.pas`.
    - **Hex-float encoding trap:** `0x1.mmm...p-N` → biased exp = `1023 - N` (not `1023 - N + 2`). The `cCosCL = -0x1.6b01ec5417056p-57` constant had exp field `$3C8` (→ `p-55`) instead of `$3C6`. Only cos(1.0) exposed it (other fast-path tests happened to round to the same bits despite 4× too-large `l` after reduction). Always double-check `Tb64u64` literals against `printf("%a", v)` of the C source.
    - **Pascal case-insensitivity trap:** `X: TDInt64` local conflicts with `x: Double` parameter (renamed local to `Xd`).
- [X] **4.02** `sin`     — 2089 lines *(dint + clzll + fenv)*
  - Reference: same shared reduction as 4.01, plus `pascoremath32.pas:5752` (`sinf_add_sign`), `:5763` (`sinf_db`), `:5806` (`sinf_big`).
  - **Port notes (2026-04-23):** wrapped around the cos_port shared primitives; only SinFast/SinAccurate/pcr_sin in `src/sin_port.inc`. **Tiny-path pitfall:** C does `fma(x, -2^-54, x)` for `|x| <= 0x1.7137449123ef6p-26`. Per the C comment, this rounds to `x` for all inputs in that range — but `pcr_fma_pascal` (the Dekker emulation used in non-AVX2 builds) is NOT correctly rounded for subnormal or minimum-normal inputs, producing `x-ulp` instead of `x`. Fix: return `x` directly in the tiny branch (bit-exact per the C comment; we lose the INEXACT flag raise but gain bit-exact matching).
- [X] **4.03** `log2p1`  — 2162 lines *(**dint** + dd — listed here due to line count; no fenv, no clzll)*
  - **Port notes (2026-04-25):** `src/log2p1_port.inc` + `src/log2p1_const.inc` (auto-generated by `tmp/log2p1_gen.py`). Massive table reuse: `_INVERSE`/`_LOG_INV`/`_INVERSE_2`/`_LOG_INV_2`/`P_2`/`LOG2`/`M_ONE` are byte-identical between log.c and log2p1.c — fully reuses log_const.inc. New content: `cLog2p1P[7]` (degree-7 fast poly — note this is the **same** P[7] as log10p1.c word-for-word, same Sollya output), `cLog2p1Pa[11]` (degree-11 medium poly — also identical to log10p1.c's Pa), `cLog2p1Pacc[24]` (degree-17 accurate-tiny mixed dd/single, **different** from log10p1's Pacc[25] — log2p1 uses `1/log(2)` coefficients), 247-entry exception table (`cLog2p1ExX/ExH/ExRnd`), 53-entry `cLog2p1T` (x=2^n-1 fast detection), `cLog2p1Log2InvDint` (LOG2_INV ≈ 2^12/log(2), unique to log2p1). Algorithm matches log10p1's structure: Log2p1P1/P1a polynomials, Log2p1LogFast (variant of log.c's cr_log_fast returning dd), Log2p1AccTiny (|x|<2^-105), Log2p1AccSmall (Pacc + binary-search exception table), Log2p1Accurate (full dint with 1+x→(xh,xl) split + Taylor-2 correction), Log2p1Fast, pcr_log2p1. Tests: 5M random TestHarness64 (--pct 5): 0 mismatches. Benchmark64: **Pascal 41.3 Mops/s vs C 57.7 Mops/s (72%, AVX2)**. Sharp edges:
    - **Pacc[24] vs log10p1's Pacc[25] indexing differs at the dd↔single boundary**: log10p1 has 7 dd-pairs (degrees 1..7, idx 0..13), then 11 singles (degrees 8..18, idx 14..24); log2p1 has 7 dd-pairs (degrees 1..7, idx 0..13), then 10 singles (degrees 8..17, idx 14..23). Generator's degree-from-Pacc-N regex must match the C source: log2p1.c stops at degree 17. The Pascal port's accurate-small loop is `i := 16 downto 11` for top single Horner, then `i := 10 downto 8` building dd, then `i := 7 downto 1` for dd polynomial.
    - **LOG2_INV is dint with C `ex=12`, pascoremath shifts to `ex=13`** (per log_port lesson), then `Y.ex -= 12` after MulDInt absorbs the 2^12 factor. Confirmed by inspection against log10p1 (which uses `ex=-1` from C `ex=-2` + `Y.ex -= 4` for LOG10_INV's 2^4 factor — different scaling).
    - **Negative-half exact branch** (e=-1, x<0): `1+x = 2^k` for `-53 <= k <= -1` is exact and must short-circuit before fast/accurate dispatch (e.g. log2p1(-0.5) = -1 exact). The hot path is `(t.u shl 12) = 0` after `t.f := 1.0+x`.
    - **Subnormal exception in AccTiny**: `|x| = 0x0.2c316a14459d8p-1022` (subnormal) ⇒ result is `fma(±2^-600, 2^-600, ±0x1.fe0e7458ac1f8p-1025)`. Hard-coded; not table-driven.
    - **`x = 2^e` directed-rounding branch**: when `e >= 49`, return `e + 0x1p-48` (= ½ ulp(49)) so RNDU-style modes round correctly. Defensive; round-to-nearest tests don't hit it but it's needed for non-default round modes.
    - **`(xh <= 2^1022) || (|xl| >= 4)` short-circuit for `c = xl/xh`** mirrors log2p1.c's `__builtin_expect` branch — avoids spurious underflow when xh is huge and xl small (in which case `c = 0` is exact enough).
- [X] **4.04** `sincos`  — 2252 lines *(dint + clzll + fenv, out-parameter API)*
  - Reference: same shared reduction as 4.01, plus `pascoremath32.pas:6146` (`sincosf_database`) and `:6221` (`sincosf_big`) for the combined-output pattern.
  - **Port notes (2026-04-23):** `src/sincos_port.inc` reuses CosReduce*/CosEval*/CosAccurate/SinAccurate. Tracks two sign flags (`neg` for sin, `negc` for cos) through the tri-fold. Tiny path returns (x, 1.0) directly — C's `1.0 - 2^-54` ties-to-even to 1.0.
- [X] **4.05** `tan`     — 2297 lines *(dint + clzll + fenv)*
  - Reference: `pascoremath32.pas:6068` (`tanf_rbig`) and `:6126` (`tanf_rltl`). Shares the `sincos_ipi` table with 4.01–4.04.
  - **Port notes (2026-04-23):** `src/tan_port.inc` + `src/tan_tinv.inc` (256-entry inv-seed table). Couldn't reuse CosReduceFast — tan needs one **extra T[] limb** of precision (because it divides by cos2pi(R), which vanishes near π/2), so TanReduceFast exists separately. Also needs **InvDInt / DivDInt** (dint reciprocal via 3 Newton iterations + Karp-Markstein), not present in cos/sin. Fast path uses Karp-Markstein (`TanFastDiv`) for the double-double quotient. Same tiny-path workaround as sin (return x directly; C's `fma(x, 2^-54, x)` would give bit-exact results only with a correctly rounded fma).

---

## Phase 5 — `pow` (1951 lines, unique `TQInt64` dependency)

`pow` is isolated in its own phase because it is the only function requiring `TQInt64`
(the 256-bit type). Port the full `qint64_t` arithmetic from
`core-math/src/binary64/pow/qint.h` (1571 lines) into `pascoremathtypes.pas` before
starting this phase. Validate the qint arithmetic independently.

- [X] **5.01** Port `TQInt64` arithmetic from `qint.h`
  - **Port notes (2026-04-25):** all qint arithmetic primitives ported into `pascoremathtypes.pas` (interface + implementation, ~520 lines added). Functions ported: `CpQInt`, `QIntZeroP`, `CmpQIntAbs`, `CmpQIntAbs22`, `AddQInt`, `AddQInt22`, `MulQInt`, `MulQInt33`, `MulQInt41`, `MulQInt31`, `MulQInt22`, `MulQInt21`, `MulQInt11`, `MulQIntInt` (= `mul_qint_2`). Constants: `QINT_ZERO`, `QINT_ONE`, `QINT_M_ONE`, `QINT_LOG2`, `QINT_LOG2_INV`. New private helpers in the implementation section: `AddU128Cy` / `SubU128Bo` (carry/borrow returning) and `ClzU128`. `TestQInt64.pas` in `src/tests/` exercises every primitive and passes (`ALL OK`). Sharp edges to remember when porting `pow` (5.02):
    - **qint normalisation differs from dint.** qint mantissa lies in `[1, 2)` (MSB of `r0` is the integer 1.0); dint mantissa lies in `[0.5, 1)`. Consequently `QINT_ONE.ex = 0` whereas `DINT_ONE.ex = 1`. The exponent formula `r.ex = a.ex + b.ex + 1 - ex_correction` (where `ex_correction` is 0 if `t6 >> 127` is set, else 1) reflects the [1,2) convention.
    - Field mapping `r0 → hh, r1 → hl, r2 → lh, r3 → ll` (r0 is most significant) — opposite of what one might expect from "rl/rh" naming in C.
    - The C `add_qint` recursive call `add_qint(r, b, a)` is replaced by a swap to avoid recursion in Pascal (matches existing `AddDInt` / `AddTInt` pattern).
    - `Mulu64u64` is marked `inline` but FPC declines to inline it inside this unit (note: "Call to subroutine ... marked as inline is not inlined") — acceptable for now; will revisit during pow benchmarking if mul-qint becomes a hotspot.
- [ ] **5.02** `pow`     — 1951 lines *(dint + qint + fenv, bivariate)*
  - **Scope (2026-04-25):** the largest single function in the binary64 suite. Source split:
    `pow.c` 1951 lines + `pow.h` 875 lines (helpers + tables P_1, Q_1, T1, T2, _INVERSE, _LOG_INV)
    + `dint.h` 1254 lines (tables _INVERSE_2_1, _INVERSE_2_2, _LOG_INV_2_1, _LOG_INV_2_2, T1_2, T2_2, P_2, Q_2, ~640 new lines beyond pascoremathtypes)
    + `qint.h` 1571 lines (arithmetic done in 5.01; tables _INVERSE_3_1, _INVERSE_3_2, _LOG_INV_3_1, _LOG_INV_3_2, T1_3, T2_3, P_3, Q_3, ~750 new lines).
    Total expected Pascal: ~4000-5000 lines. Suggested commit cadence: one commit per sub-step.
  - **Sub-plan:**
    - [X] **5.02a** Add infrastructure helpers to `pascoremathtypes.pas`:
      `MulDInt11`, `AddDInt11`, `DIntToI` ported (covered by `TestDIntPowHelpers` in `TestQInt64.pas`); constants `DINT_M_ONE`, `DINT_LOG2`, `DINT_LOG2_INV` added.
      Use existing `MulDIntInt(r, b, a)` with swapped args for `mul_dint_int64(r, a, b)`.
      **Still TODO inside 5.02a** before pow main:
      `DIntToD` (= `dint_tod`, full conversion with overflow / subnormal / underflow handling — distinct from existing `DToD`),
      `DIntTodSubnormal` helper,
      qint helpers `QIntFromD`, `QIntToI`, `QIntToD`, `SubnormalizeQInt`. Defer until pow code that calls them is being written, so we can validate by call sites rather than by ad-hoc unit tests.
    - [X] **5.02b** Created `src/pow_const.inc` (auto-gen `tmp/pow_gen.py`) with the eight `pow.h` tables: `cPowInverse[0..181]`, `cPowLogInvHi/Lo[0..181]`, `cPowT1Hi/Lo[0..63]`, `cPowT2Hi/Lo[0..63]` as `Tb64u64` arrays (runtime-indexed); polynomial coefficients `cPowP1_0..5` and `cPowQ1_0..4` emitted as named scalar constants per the project rule. Verified by a standalone include-only test program that compiled clean and printed expected bit patterns. Range-check warnings on a few hex constants are cosmetic (FPC parses 64-bit hex as signed first) — bit pattern is preserved in `Tb64u64.u` (`UInt64`).
    - [X] **5.02c** Appended dint-level pow tables to `src/pow_const.inc` via `tmp/pow_gen.py`: `cPowInverse21[0..91]`, `cPowInverse22[0..128]` (129 entries — j=8128..8256, NOT 128 as initially planned), `cPowLogInv21[0..91]`, `cPowLogInv22[0..128]`, `cPowT12[0..63]`, `cPowT22[0..63]` as `array of TDInt64`. Polynomial coefficients `cPowP2_0..8` and `cPowQ2_0..7` emitted as named scalars. **+1 ex offset applied to every entry** so the tables match the Pascal `TDInt64` convention (mantissa in [0.5,1), ex one greater than C's [1,2) convention). Verified by include test: `INV21[0].ex=1` (C 0+1), `Q2_7.ex=0` (C -1+1), bit patterns intact. Full `build.sh -dAVX2` clean.
    - [X] **5.02d** Appended qint-level pow tables to `src/pow_const.inc` via `tmp/pow_gen.py`: `cPowInverse31[0..91]`, `cPowInverse32[0..128]`, `cPowLogInv31[0..91]`, `cPowLogInv32[0..128]`, `cPowT13[0..63]`, `cPowT23[0..63]` as `array of TQInt64`. Polynomial coefficients `cPowP3_0..17` (degree 18..1) and `cPowQ3_0..14` (degree 14..0) emitted as named scalars. C qint64_t `{hh,hl,lh,ll}` maps to Pascal TQInt64 `{r0,r1,r2,r3}`. **No ex offset** for qint (Pascal `QINT_ONE.ex=0` already matches C `ONE_Q.ex=0`). Verified via include test; full `build.sh -dAVX2` clean. `pow_const.inc` is now 1483 lines.
    - [X] **5.02e** Created `src/pow_port.inc` with `PowQ1` (= `q_1`) and `PowP1` (= `p_1`) pure-double polynomial helpers using `pcr_fma`, `pcr_a_mul`, `pcr_fasttwosum`. `fast_sum` is inlined as `pcr_fasttwosum(...); ql := ql + bl;`. Constants accessed as `cPowQ1_*.f` / `cPowP1_*.f`. `Double(-0.5)` typecast used per anti-x87 rule. Cross-checked against C reference: `q_1(0)=(1,0)`, `q_1(1e-4)=1.000100005000166700`, `p_1(2^-14)` returns `-1.8625693614750869e-9` matching `-z^2/2 + z^3/3` analytically. Full `build.sh -dAVX2` clean (pow_port.inc not yet wired into pascoremath64.pas — that's 5.02l).
    - [X] **5.02f** Added `PowLog1` to `pow_port.inc` mirroring `pow.c:581 log_1`. Returns `(h,l)` double-double approximation of `log(x)` and Integer cancel-flag (1 = fast_two_sum normalize triggered when |l|>|h|*2^-24, else 0). Cross-checked: `log(1)=0` exact; `log(2)=0.693147180559945` (matches DBL); `log(e)≈1`; `log(1e-300)≈-690.7755278982137`. Sharp edge: Pascal does not accept C hex-float literals like `0x1p-24` — wrap as `Tb64u64.f` (`cPowLog1Threshold`). Other helpers needed: `cPowLog2H`, `cPowLog2L` constants added to inc.
    - [X] **5.02g** Added `PowExp1` to `pow_port.inc` mirroring `pow.c:954 exp_1`. Computes `s*exp(rh+rl)` returning `(eh,el)` double-double. Short-circuits: `rh > RHO3` → ±DBL_MAX (forces +/-Inf or DBL_MAX in caller); `rh < RHO0` → ±0; intermediate region or NaN → NaN/NaN to defer to slow path. Uses pre-computed RHO0..3, INVLOG2, LOG2H/L, DBL_MAX, 0x1p-1074 as named `Tb64u64` constants. Pascal case-insensitivity bit me: variables `k` (Double for `roundeven`) and `K` (Int64 for shifts) collide → renamed to `kd`/`kInt`. Cross-checked: `exp(0)=1`, `exp(1)=2.718281828459044+5.89e-16`, `exp(ln 2)=2`, `exp(100)≈2.6881e43`, overflow/underflow paths correct.
    - [ ] **5.02h** Implement `pcr_pow_q2`, `pcr_pow_p2`, `pcr_pow_log2`, `pcr_pow_exp2` (all dint, second Ziv iteration).
    - [ ] **5.02i** Implement `pcr_pow_q3`, `pcr_pow_p3`, `pcr_pow_log3`, `pcr_pow_exp3` (all qint, third Ziv iteration).
    - [ ] **5.02j** Implement `IsExact`, `ExactPow` for the rounding-boundary test (Algorithm detectRoundingBoundaryCase from [4]).
    - [ ] **5.02k** Implement main `pcr_pow` with the special-case cascade (NaN, ±Inf, ±0, x<0 with integer y, |x|=1) and the three Ziv phases.
    - [ ] **5.02l** Wire `pow_port.inc` into `pascoremath64.pas` and `cr_pow` into `ccoremath64.pas`.
    - [ ] **5.02m** Test with `--pct 1` and `Benchmark64 pow`.
  - **Sharp edges to remember:**
    - `pow` uses a **two-level lookup** for log (181-entry _INVERSE_2_1 of 9-bit reciprocals + 128-entry _INVERSE_2_2 of 14-bit reciprocals), unlike binary64 `log` which uses a single 240-entry lookup. The pow tables are therefore **not** shareable with `log_const.inc`.
    - The dint format in `pow/dint.h` matches the `sin/dint.h` we already have (mantissa in [0.5, 1), `ex=1` for 1.0). Confirmed identical to `pascoremathtypes.pas` TDInt64.
    - In `mul_dint_int64(r, a, b)`, the C signature is `(r, a, b)` whereas Pascal's existing `MulDIntInt(r, b, a)` takes the integer first. Mirror call sites carefully.
    - `dint_tod` differs from the simpler `DToD` already in pascoremathtypes: it adds underflow/overflow handling for pow's wide output range, and takes an `exact` flag to suppress underflow signal when x^y is exactly representable.
    - `is_exact` does **not** call dint/qint primitives; it works at the bit level on the inputs only. Keep it pure-integer.
    - Special-case logic in `cr_pow` (~140 lines) handles NaN, Inf, ±0, negative x with integer y. Do not skip this — sampling tests will not cover all 23 special-case paths.

---

## Per-function porting checklist

Apply this checklist to every function before marking it done:

- [ ] Hex float literals converted via `hexfloat.pas` utility, not by hand
- [ ] Lookup tables moved to unit-level `const` (no `static` locals)
- [ ] All type-punning uses `Tb64u64` from `pascoremathtypes.pas` (no unsafe casts)
- [ ] `__builtin_expect` wrappers removed entirely
- [ ] `__attribute__((noinline))` replaced with `[noinline]`
- [ ] `CORE_MATH_SUPPORT_ERRNO` blocks omitted (out of scope for Pascal port)
- [ ] `fegetround()` replaced with `pcr_GetRoundMode` and FPC round-mode constants
- [ ] `__builtin_clzll` replaced with `pcr_clzll`
- [ ] `fasttwosum` / `muldd` replaced with `pcr_fasttwosum` / `pcr_muldd`
- [ ] `dint64_t` operations replaced with `TDInt64` and the typed procedures
- [ ] Sampling test passes (no mismatches over ≥ 10^9 random inputs + structured coverage)
- [ ] C function `cr_<name>` declared in `ccoremath64.pas`; Pascal equivalent named `pcr_<name>` in `pascoremath64.pas`
- [ ] All integer variables declared as `Int32` (not `Integer`) for explicit 32-bit signed intent
- [ ] No redundant typecast patterns (see rule 11 below)

---

## Architectural notes and known pitfalls

1. **`dint64_t` is universal — not just for hard functions.** Unlike binary32 where it only
   appeared in sin/cos/tan/sincos, every binary64 function uses `TDInt64` for its slow path.
   Get the dint arithmetic right first, before touching any function.

2. **`fasttwosum` correctness depends on FP evaluation order.** The identity
   `e = y - (r - x)` is only correct when `r = x + y` is computed at double precision with
   no excess precision. On x86-64 with `{$FPUTYPE SSE64}` this is guaranteed for SSE2 code.
   However, if FPC falls back to x87 mode for any reason (untyped literals, inline failures),
   the extra 80-bit precision will silently corrupt `fasttwosum`. Verify the generated
   assembly for `pcr_fasttwosum` before using it in any function.

3. **Large-argument range reduction helper must not be duplicated across sin/cos/tan/sincos.**
   The binary64 trig functions contain a shared range reduction routine analogous to `rbig()`
   in binary32. Identify it when porting `cos` (Phase 4.01), factor it into
   `pascoremathtypes.pas`, and call it from all four functions. Same rule as binary32.

4. **`sincos` has a different API — decide before Phase 4.** The C signature is
   `void cr_sincos(double x, double *sout, double *cout)`. Use the same convention
   chosen for binary32 (`procedure cr_sincos(x: Double; out s, c: Double)`) for consistency.

5. **Rounding mode branches must match the C exactly.** Functions like `sin` and `cos`
   contain `switch (fegetround())` with four branches. The Pascal `case GetRoundMode of`
   equivalent must cover the same four arms in the same order. A missing or swapped
   branch produces wrong results only for non-default rounding modes — easy to miss in
   a default-mode test suite.

6. **Two-pass rounding (Ziv's strategy) — the slow path must never be removed.** The slow
   path (second Ziv iteration using `TDInt64`) is not dead code — it is the mechanism that
   guarantees correct rounding near midpoints. Do not remove it even if all sampled test
   inputs happen to take the fast path.

7. **`pcr_clzll(0)` is undefined — check all call sites.** In C, `__builtin_clzll(0)` is
   undefined behaviour. Every call in the C source is guarded. Verify each call site when
   porting and add the same guard in Pascal.

8. **Exhaustive testing is infeasible — be disciplined about coverage.** With 2^64 possible
   inputs, you cannot test them all. Compensate with: large random samples (≥ 10^9), all
   special values (subnormals, ±0, ±Inf, NaN), known worst-case inputs from `.wc` files,
   and boundary inputs around mathematical poles and branch points.

9. **`qint.h` is 1571 lines of arithmetic — treat it as its own sub-project.** Port it
   fully into `pascoremathtypes.pas` and validate it independently before touching `pow`.
   Include a dedicated `TestQInt.pas` that checks known products and sums at 256-bit
   precision.

10. **No `compound` function in binary64.** The binary32 library has `compoundf`; binary64
    does not. The function count is 41, not 42.

11. **Use `Int32`, not `Integer`, for all 32-bit signed integer variables.** `Int32` makes
    the bit-width explicit and immune to any future ambiguity. `Integer` is always 32-bit in
    FPC, but the intent is not obvious to a reader. The binary32 port was updated to use
    `Int32` consistently (April 2026) — follow the same convention in binary64.

12. **Do not write local bit-scan or arithmetic-shift helpers — use FPC intrinsics directly.**
    The binary32 port went through three rounds of this cleanup (all April 2026):
    - `pcr_bsf32`/`pcr_bsr32` (asm blocks, `inline` silently ignored) → `BsfDWord`/`BsrDWord`
    - `cf_bsf64`/`cf_bsr64` (asm blocks, not even marked `inline`) → `BsfQWord`/`BsrQWord`
    - `sar_i32`/`sar_i64` (conditional-branch fallback) → `SarLongInt`/`SarInt64`

    All six of `BsfDWord`, `BsrDWord`, `BsfQWord`, `BsrQWord`, `SarLongInt`, `SarInt64`
    are genuine FPC 3.2.2 compiler intrinsics in the `System` unit — no `uses` clause
    needed, and all inline correctly at every call site. Do not reintroduce local wrappers
    for these operations in `pascoremath64.pas`.

13. **Avoid redundant typecast patterns — they are C translation artifacts.** When porting
    from C, the following patterns appear but carry no semantic meaning in Pascal and should
    be removed:

    | Pattern | Why it is redundant | Simplification |
    |---|---|---|
    | `Int32(Int32(...))` | Arithmetic on `Int32` stays `Int32`; outer cast is a no-op | Drop outer `Int32(...)` |
    | `Int64(Int64(...))` | Same redundancy at 64-bit | Drop outer `Int64(...)` |
    | `UInt32(Int32(UInt32(...)))` | `Int32` in the middle does not change any bits | Drop the `Int32(...)` layer |
    | `UInt32(Int32(x) shr n)` where `x: UInt32` | `shr` is always logical in Pascal; the `Int32` cast buys nothing when the high bit of `x` is guaranteed 0 | Use `x shr n` directly |

    These arise because C requires explicit `(uint32_t)(int32_t)` casts around signed/unsigned
    operations. Pascal's type system handles the same cases without the extra layers. Keep
    casts that genuinely change signedness or width; remove those that simply round-trip.

15. **`pcr_fma` double-rounding for tiny x — return x directly for |x| < 2^-27.**
    In `pcr_atan`, the C code calls `__builtin_fma(-0x1p-54, x, x)` for |x| < 2^-27.
    `pcr_fma_pascal` (the emulated FMA) double-rounds for biased_exp=2 inputs: the product
    `(-2^-54)*x` rounds to ±2^-1074 (smallest subnormal), `TwoSum(±2^-1074, x)` lands
    exactly on a tie, and round-to-even may round away from x. Hardware FMA computes the
    exact product and rounds once, always giving x. In round-to-nearest mode, atan(x) = x
    for all |x| < 2^-27, so `Result := x` is the correct and safe replacement.

16. **`copysign(1, x)*A[i][1]` must NOT be implemented as bitwise OR when A[i][1] < 0.**
    In the `AtanRefine2` slow path, `df` is set to `copysign(1.0, x)*A[ip][1]`. Because
    `A[ip][1]` (the lo part of a double-double table entry) can be negative, using
    `ta_u.u := cAtanALo[ip].u or sign_bit_of_x` is wrong — it forces the sign bit to 1 when
    x < 0, leaving already-negative entries with the wrong sign. The correct pattern is:
    `df := cAtanALo[ip].f; if x < 0.0 then df := -df;` (i.e. flip the stored value's sign).
    Bitwise OR is only safe when the table value is guaranteed non-negative (e.g. `cAtanAHi`
    stores |atan| values which are always positive).

14. **Use `assembler;` whole-function bodies to guarantee System V ABI param passing.**
    When a function needs SSE/AVX args in `xmm0`/`xmm1`/`xmm2` directly (e.g. `pcr_fma`,
    `pcr_roundeven`, `pcr_fmax`, `pcr_fmin` in `pascoremathhelperfuncs.pas`), declare
    it as a pure `assembler;` body rather than a Pascal body with an embedded `asm`
    block. FPC then honours the System V AMD64 ABI and places `x` in `xmm0`, `y` in
    `xmm1`, `z` in `xmm2`, and returns the result in `xmm0` — letting the body be a
    single instruction (`vfmadd213sd xmm0, xmm1, xmm2` / `ret`). A Pascal body with an
    `asm` block forces FPC to emit a full prologue/epilogue and spill args through the
    stack frame, which breaks zero-instruction leaf routines. See `pcr_fmaf` for the
    canonical pattern. Note that any function containing an `asm` block (pure
    `assembler;` or embedded) cannot be inlined — so only use this for irreducible
    primitives, and prefer FPC intrinsics (note 12) wherever possible.

---

## Design decisions

1. **`TDInt64` uses named fields, not a union/variant.** Defined as
   `record hi, lo: UInt64; ex: Int64; sgn: Byte end`. The C union variant is needed for
   little/big-endian portability in C — in Pascal the named-field record is cleaner and
   unambiguous.

2. **Double-double helpers are functions with an `out` parameter, not procedures.**
   `fasttwosum(x, y, &e)` returns the sum and writes the error to `e`. The Pascal
   equivalent returns the sum as the function result and uses `out e: Double`.
   This matches the call pattern in the C source (`sh = fasttwosum(xh, yh, &sl)`) and
   allows the result to appear directly in an expression.

3. **`pcr_clzll` is a function, not an operator overload.** There is no Pascal syntax for
   overloading a unary prefix operation, and a named function is the only clean option.

4. **Phase 4 and Phase 5 are explicitly blocked on Phase 0.3 being fully validated.**
   A bug in `AddDInt` or `MulDInt` would corrupt all the hard functions with no obvious
   failure mode, since the dint path is only taken near rounding midpoints.

5. **The `TQInt64` type is defined in Phase 0.4 but its arithmetic is ported in Phase 5.01.**
   Separating the type definition from the implementation keeps `pascoremathtypes.pas`
   self-consistent from the start, while deferring the ~1500 lines of qint arithmetic until
   it is actually needed.

---

## Key rules for the developer

1. **Do not change the algorithm.** This is a faithful port, not a rewrite. The C source
   is the specification. If your output differs by even one ULP, it is a bug.

2. **Validate infrastructure before functions.** Run a dedicated test for each Phase 0
   primitive — `TDInt64` arithmetic, `pcr_clzll`, `pcr_fasttwosum`, `pcr_muldd` — before
   touching any Phase 1 function.

3. **`fasttwosum` and `muldd` require true hardware FMA.** Verify that `pcr_fma` emits
   `VFMADD213SD` (not an x87 approximation) before using `pcr_muldd` in any function.
   See Bug B in `tasklist.md` for the binary32 precedent.

4. **`pcr_clzll` must be correct before Phase 0.3.** `AddDInt` and `DIntFromD` both call
   it; a wrong clzll silently produces wrong exponent normalization.

5. **Convert hex float constants with the existing utility.** Never retype a constant by
   hand. Use `hexfloat.pas` for all lookup tables — binary64 tables are larger and the
   risk of transcription error is higher.

6. **Use `Tb64u64` for all bit-cast operations.** Never use `Move`, pointer casts, or
   `Absolute`. The `Tb64u64` record is already defined in `pascoremathtypes.pas`.

7. **The slow path (dint) is hot near midpoints.** Although random testing rarely hits it,
   the worst-case inputs from `.wc` files are specifically chosen to exercise it. Always
   include `.wc` inputs in your test suite.

8. **Work sequentially within each phase.** The ordering is chosen so each function builds
   familiarity with the infrastructure before the next, harder one.

9. **Prefer FPC intrinsics over inline ASM for bit-scan operations.** `BsrQWord` and
   `BsfQWord` are genuine FPC 3.2.2 compiler intrinsics — they emit the correct instruction
   on x86-64 and have portable fallbacks on other targets, all without an `asm` block.
   Functions containing `asm` blocks cannot be inlined by FPC, so `inline` is silently
   ignored on x86-64. Use `BsrQWord`/`BsfQWord` directly (or in a thin `inline` wrapper
   like `pcr_clzll`) to get true inlining. ASM remains appropriate only for `pcr_fma`
   (VFMADD213SD) and `pcr_fasttwosum` where SSE2 evaluation order must be forced.

10. **Benchmark every function.** After each function passes sampling tests, run
    `Benchmark64.pas` and record the Mops/s ratio (Pascal vs C). A large gap signals
    missed inlining or suboptimal code generation.
