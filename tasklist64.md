# pas-core-math Task List — binary64 (Double precision)

Port of the CORE-MATH binary64 (double-precision floating-point) C library to Free Pascal.
Goal: bit-exact, correctly-rounded results matching the C reference for a large random sample
of `Double` inputs (exhaustive testing of all 2^64 inputs is infeasible).

C source reference: `core-math/src/binary64/`

---

## Status summary

- **0 of 41 functions ported** (work not yet started)
- Target file: `src/pascoremath64.pas`

---

## Critical difference from the binary32 port

In binary32, `dint64_t` was only needed by the four hardest trig functions (sin/cos/tan/sincos).
In binary64, **every single function** uses `dint64_t` — it is the standard second-pass data
type for Ziv's rounding strategy across the entire library. The double-double helpers
(`fasttwosum`, `muldd`) are similarly universal, appearing in 33 of the 41 functions.

**Phase 0 infrastructure is therefore the critical path.** No function can be correctly ported
until `TDInt64` and its arithmetic are implemented and independently validated.

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

### [ ] 0.1 — Create `src/pascoremath64.pas` and `src/ccoremath64.pas` scaffolding

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

---

### [ ] 0.2 — Define `TDInt64` in `pascoremathtypes.pas`

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

### [ ] 0.3 — Implement `dint64_t` arithmetic in `pascoremathtypes.pas`

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

### [ ] 0.4 — Define `TQInt64` in `pascoremathtypes.pas`

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

### [ ] 0.5 — Implement `pcr_clzll` in `pascoremathhelperfuncs.pas`

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

### [ ] 0.7 — Implement rounding mode detection helper

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

### [ ] 0.8 — Set up test programs

Create `src/tests/TestHarness64.pas`, `src/tests/Benchmark64.pas`, and
`src/tests/BenchmarkFPC64.pas`. These are structurally identical to their binary32
counterparts but use `Double`, `Tb64u64`, and the 64-bit `cr_*` / `pcr_*` functions.

**Key difference from binary32:** exhaustive testing over all 2^64 `Double` inputs is
infeasible. The test strategy is:

1. **Random sampling** — test at least 10^9 random `Double` inputs per function using a
   simple LCG or xorshift64 RNG seeded from the system clock.
2. **Structured coverage** — always include: all subnormals, ±0, ±Inf, NaN, all powers
   of 2, boundary values near ±1, ±π, ±π/2, and known worst-case inputs from the
   `.wc` files in `core-math/src/binary64/<function>/`.
3. **Bivariate functions** (atan2, atan2pi, hypot, pow) — use a structured grid of
   10^6 × 10^6 sampled pairs.

Extend `build.sh` (or add `build64.sh`) to compile all three binary64 programs.

---

## Phase 1 — Simple-to-medium (220–420 lines, dint + double-double, no fenv)

All functions in this phase use `TDInt64` and most use `fasttwosum`/`muldd`, but have
a straightforward one- or two-pass Ziv structure with no `fegetround` branches.
Port in this order. All functions live in `pascoremath64.pas`, named `pcr_<name>`.

- [ ] **1.01** `rsqrt`   — 220 lines  *(uses clzll, fenv — see note below)*
- [ ] **1.02** `cbrt`    — 252 lines  *(uses clzll, fenv)*
- [ ] **1.03** `atan`    — 281 lines  *(uses dint + dd)*
- [ ] **1.04** `log2`    — 313 lines  *(uses clzll + dd)*
- [ ] **1.05** `acos`    — 354 lines  *(uses dint + dd)*
- [ ] **1.06** `tanh`    — 355 lines  *(uses dint + dd)*
- [ ] **1.07** `cospi`   — 356 lines  *(uses dint + dd)*
- [ ] **1.08** `asin`    — 366 lines  *(uses dint + dd)*
- [ ] **1.09** `cosh`    — 377 lines  *(uses dint + dd)*
- [ ] **1.10** `exp10`   — 379 lines  *(uses dint + dd)*
- [ ] **1.11** `exp2`    — 384 lines  *(uses dint + dd)*
- [ ] **1.12** `exp`     — 386 lines  *(uses dint + dd)*
- [ ] **1.13** `tanpi`   — 388 lines  *(uses dint + dd)*
- [ ] **1.14** `sinpi`   — 400 lines  *(uses dint + dd)*
- [ ] **1.15** `sinh`    — 418 lines  *(uses dint + dd)*

Note: `rsqrt` and `cbrt` technically include `fegetround` calls, but they are simple
rounding-mode branches, not the multi-path Ziv structure seen in sin/cos. Start with them
anyway — they are the shortest files and good warm-up exercises.

---

## Phase 2 — Medium (436–882 lines)

- [ ] **2.01** `expm1`   — 436 lines  *(uses dint + dd)*
- [ ] **2.02** `acosh`   — 451 lines  *(uses dint + dd)*
- [ ] **2.03** `atanh`   — 479 lines  *(uses dint + dd)*
- [ ] **2.04** `atanpi`  — 479 lines  *(uses dint + dd)*
- [ ] **2.05** `asinh`   — 489 lines  *(uses dint + dd)*
- [ ] **2.06** `log1p`   — 490 lines  *(uses dint + dd)*
- [ ] **2.07** `atan2`   — 586 lines  *(uses dint + dd, bivariate)*
- [ ] **2.08** `erf`     — 710 lines  *(uses dint only)*
- [ ] **2.09** `asinpi`  — 798 lines  *(uses dint + clzll + fenv)*
- [ ] **2.10** `log`     — 832 lines  *(uses dint + clzll)*
- [ ] **2.11** `atan2pi` — 866 lines  *(uses dint, bivariate)*
- [ ] **2.12** `log10`   — 882 lines  *(uses dint + clzll)*

---

## Phase 3 — Hard (1022–1577 lines)

- [ ] **3.01** `exp2m1`  — 1022 lines *(uses dint)*
- [ ] **3.02** `tgamma`  — 1096 lines *(uses dint + dd)*
- [ ] **3.03** `acospi`  — 1099 lines *(uses dint)*
- [ ] **3.04** `exp10m1` — 1153 lines *(uses dint)*
- [ ] **3.05** `erfc`    — 1247 lines *(uses dint)*
- [ ] **3.06** `lgamma`  — 1452 lines *(uses dint + dd)*
- [ ] **3.07** `log10p1` — 1577 lines *(uses dint)*
- [ ] **3.08** `hypot`   — 283 lines  *(uses dint + dd + fenv — listed here due to fenv complexity)*

---

## Phase 4 — Very Hard (2068–2297 lines, fenv + full dint trig)

These are the four core trig functions. Each uses `TDInt64` with the full arithmetic suite,
`fegetround` multi-path branching, and a 256-entry precomputed sin/cos table. The
large-argument range-reduction helper will appear in all four — it **must not be duplicated**.

Validate `AddDInt`, `MulDInt`, `DIntFromD`, and `DToD` exhaustively before starting.

- [ ] **4.01** `cos`     — 2068 lines *(dint + clzll + fenv)*
- [ ] **4.02** `sin`     — 2089 lines *(dint + clzll + fenv)*
- [ ] **4.03** `log2p1`  — 2162 lines *(dint — listed here due to line count)*
- [ ] **4.04** `sincos`  — 2252 lines *(dint + clzll + fenv, out-parameter API)*
- [ ] **4.05** `tan`     — 2297 lines *(dint + clzll + fenv)*

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
    | `LongWord(Int32(LongWord(...)))` | `Int32` in the middle does not change any bits | Drop the `Int32(...)` layer |
    | `LongWord(Int32(x) shr n)` where `x: LongWord` | `shr` is always logical in Pascal; the `Int32` cast buys nothing when the high bit of `x` is guaranteed 0 | Use `x shr n` directly |

    These arise because C requires explicit `(uint32_t)(int32_t)` casts around signed/unsigned
    operations. Pascal's type system handles the same cases without the extra layers. Keep
    casts that genuinely change signedness or width; remove those that simply round-trip.

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
