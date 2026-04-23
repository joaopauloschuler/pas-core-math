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

- **8 of 41 functions ported** (Phase 0 infrastructure complete; Phase 1 in progress — 1.01 rsqrt, 1.02 cbrt, 1.03 atan, 1.04 log2, 1.05 acos, 1.06 tanh, 1.07 cospi, 1.08 asin done)
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

The double-double helpers (`fasttwosum`, `muldd`) *are* near-universal, appearing
in 33 of the 41 functions.

**Phase 0 infrastructure is still the critical path**, but for a different reason:
`pcr_fasttwosum` / `pcr_muldd` / `pcr_fma` correctness underpins every function's
fast *and* slow path. `TDInt64` only needs to be validated before Phases 4 and 5
(cos/sin/tan/sincos/log/log10/log10p1/log2p1/pow).

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
scales with table length.

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
- [ ] **1.09** `cosh`    — 377 lines  *(pure dd — verified; earlier "dint + dd" hint was wrong)*
- [ ] **1.10** `exp10`   — 379 lines  *(pure dd — verified)*
- [ ] **1.11** `exp2`    — 384 lines  *(pure dd — verified)*
- [ ] **1.12** `exp`     — 386 lines  *(pure dd — verified)*
- [ ] **1.13** `tanpi`   — 388 lines  *(pure dd — verified)*
- [ ] **1.14** `sinpi`   — 400 lines  *(pure dd — verified)*
- [ ] **1.15** `sinh`    — 418 lines  *(pure dd — verified)*

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

- [ ] **2.01** `expm1`   — 436 lines  *(pure dd)*
- [ ] **2.02** `acosh`   — 451 lines  *(pure dd)*
- [ ] **2.03** `atanh`   — 479 lines  *(pure dd)*
- [ ] **2.04** `atanpi`  — 479 lines  *(pure dd)*
- [ ] **2.05** `asinh`   — 489 lines  *(pure dd)*
- [ ] **2.06** `log1p`   — 490 lines  *(pure dd)*
- [ ] **2.07** `atan2`   — 586 lines  *(pure dd, bivariate)*
- [ ] **2.08** `erf`     — 710 lines  *(pure dd; no dint — earlier "dint only" annotation was wrong)*
- [ ] **2.09** `asinpi`  — 798 lines  *(clzll + fenv + dd; no dint)*
- [ ] **2.10** `log`     — 832 lines  *(**dint** + clzll + dd — first real dint user in the sequence)*
- [ ] **2.11** `atan2pi` — 866 lines  *(pure dd, bivariate; no dint)*
- [ ] **2.12** `log10`   — 882 lines  *(**dint** + clzll + dd)*

---

## Phase 3 — Hard (1022–1577 lines)

Of the Phase-3 functions, **only `log10p1` uses `TDInt64`**; the other seven are
pure double-double. `hypot` uses `clzll` but not dint or fenv. `lgamma` uses
`clzll` (2 call sites) but not dint. (Verified by grep.)

- [ ] **3.01** `exp2m1`  — 1022 lines *(pure dd; no dint)*
- [ ] **3.02** `tgamma`  — 1096 lines *(pure dd; no dint)*
- [ ] **3.03** `acospi`  — 1099 lines *(pure dd; no dint)*
- [ ] **3.04** `exp10m1` — 1153 lines *(pure dd; no dint)*
- [ ] **3.05** `erfc`    — 1247 lines *(pure dd; no dint)*
- [ ] **3.06** `lgamma`  — 1452 lines *(clzll + dd; no dint)*
  - Reference: `pascoremath32.pas:3197` (`lgamma_as_sinpi`) and `:3211` (`lgamma_as_ln`) — existing `Double`-precision auxiliary functions used by `pcr_lgammaf`. Mirror their shape when porting; the coefficient tables (`c_nz1`/`c_nz2`/`c_nz3`, `rn_md`/`rd_md`) will be re-derived from `binary64/lgamma/lgamma.c`.
- [ ] **3.07** `log10p1` — 1577 lines *(**dint** + dd — only Phase-3 dint user)*
- [ ] **3.08** `hypot`   — 283 lines  *(clzll + dd; no dint, no fenv — listed here due to historical grouping, not complexity)*

---

## Phase 4 — Very Hard (2068–2297 lines, fenv + full dint trig)

These are the four core trig functions. Each uses `TDInt64` with the full arithmetic suite,
`fegetround` multi-path branching, and a 256-entry precomputed sin/cos table. The
large-argument range-reduction helper will appear in all four — it **must not be duplicated**.

Validate `AddDInt`, `MulDInt`, `DIntFromD`, and `DToD` exhaustively before starting.

- [ ] **4.01** `cos`     — 2068 lines *(dint + clzll + fenv)*
  - Reference: `pascoremath32.pas:5633` (`sincos_ipi` 2/π table — reuse as-is), `:5653` (`sincos_rbig` large-argument reduction), `:5732` (`sincos_rltl` small/medium-argument reduction), `:5840` (`cosf_db`), `:5883` (`cosf_big`). Binary64 widens the UInt32 input to UInt64 and keeps more limbs of the product, but the reduction skeleton is identical. Per note 3 below, factor the reduction into a shared routine used by 4.01/4.02/4.04/4.05.
- [ ] **4.02** `sin`     — 2089 lines *(dint + clzll + fenv)*
  - Reference: same shared reduction as 4.01, plus `pascoremath32.pas:5752` (`sinf_add_sign`), `:5763` (`sinf_db`), `:5806` (`sinf_big`).
- [ ] **4.03** `log2p1`  — 2162 lines *(**dint** + dd — listed here due to line count; no fenv, no clzll)*
- [ ] **4.04** `sincos`  — 2252 lines *(dint + clzll + fenv, out-parameter API)*
  - Reference: same shared reduction as 4.01, plus `pascoremath32.pas:6146` (`sincosf_database`) and `:6221` (`sincosf_big`) for the combined-output pattern.
- [ ] **4.05** `tan`     — 2297 lines *(dint + clzll + fenv)*
  - Reference: `pascoremath32.pas:6068` (`tanf_rbig`) and `:6126` (`tanf_rltl`). Shares the `sincos_ipi` table with 4.01–4.04.

---

## Phase 5 — `pow` (1951 lines, unique `TQInt64` dependency)

`pow` is isolated in its own phase because it is the only function requiring `TQInt64`
(the 256-bit type). Port the full `qint64_t` arithmetic from
`core-math/src/binary64/pow/qint.h` (1571 lines) into `pascoremathtypes.pas` before
starting this phase. Validate the qint arithmetic independently.

- [ ] **5.01** Port `TQInt64` arithmetic from `qint.h`
- [ ] **5.02** `pow`     — 1951 lines *(dint + qint + fenv, bivariate)*

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
