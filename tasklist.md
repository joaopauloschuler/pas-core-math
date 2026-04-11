# pas-core-math Task List

Port of the CORE-MATH binary32 (single-precision) library to Free Pascal.
Goal: bit-exact, correctly-rounded results matching the C reference for all 2^32 float inputs.

---

## Status summary

- **All 42 functions ported and committed** (branch: `a1`)
- **Bug-fix round in progress** — detected with `bin/TestHarness --pct 1`
- **10 functions fixed** (8 bugs + sinf/cosf/sincosf/tanf x87 precision fix), committed in two commits (`25d897a`, `24acd3c`)
- **2 functions still failing** at `--pct 1`: `powf` (458 657 mismatches) and `compoundf` (1 903 717 mismatches) — see "Open bugs" section below
- Benchmark sample: `acosf  C=322.6 Mops/s  Pascal=222.2 Mops/s`

---

## Folder structure

```
pas-core-math/
├── src/
│   ├── pascoremath.inc         # compiler directives + CPU/AVX capability flags ({$I pascoremath.inc})
│   ├── pascoremathtypes.pas    # TUInt128, builtins, trig helpers (rbig, etc.)
│   ├── pascoremath.pas         # Pascal implementations (pcr_* functions)
│   ├── ccoremath.pas           # C reference external declarations (cr_* functions)
│   ├── laz-project/
│   │   ├── pas-core-math.lpi
│   │   └── pas-core-math.lps
│   ├── tools/
│   │   └── HexFloatConvert.pas
│   └── tests/
│       ├── TestHarness.pas
│       ├── TestMulWide.pas
│       └── Benchmark.pas       # single-threaded Mops/s comparison: cr_* (C) vs pcr_* (Pascal)
├── bin/
└── tasklist.md
```

---

## Phase 0 — Infrastructure (prerequisite for everything)

- [x] **0.1** Create `src/pascoremath.inc` containing the compiler directives and CPU/AVX
  capability flags. Include it at the top of every unit with `{$I pascoremath.inc}` —
  placed before the `unit` keyword so the mode directives take effect in time:
  ```pascal
  {$I pascoremath.inc}
  unit pascoremathtypes;
  ```
  Content of `pascoremath.inc`:
  ```pascal
  {$IFDEF FPC}
    {$COPERATORS ON}
    {$MODE OBJFPC}
    {$INLINE ON}
    {$MACRO ON}
    {$LONGSTRINGS ON}
    {$CODEPAGE UTF8}
    {$IFDEF CPU32BITS}
      {$DEFINE CPU32}
    {$ENDIF}
    {$IFDEF CPU64BITS}
      {$DEFINE CPU64}
    {$ENDIF}
    {$IFDEF CPUARM}
      {$DEFINE NOTAVX}
    {$ENDIF}
    {$IFDEF CPUAARCH64}
      {$DEFINE NOTAVX}
    {$ENDIF}
    {$IFDEF CPUPOWERPC}
      {$DEFINE NOTAVX}
    {$ENDIF}
    {$IFDEF CPUM68K}
      {$DEFINE NOTAVX}
    {$ENDIF}
    {$IFDEF NOTAVX}
      {$UNDEF AVX}
      {$UNDEF AVX2}
      {$UNDEF AVX512}
      {$UNDEF AVXANY}
    {$ENDIF}
  {$ELSE}
    // AVX code is supported only under FPC
    {$UNDEF AVX}
    {$UNDEF AVX2}
    {$UNDEF AVX512}
    {$UNDEF AVXANY}
    {$IFDEF CPU32BITS}
      {$DEFINE CPU32}
    {$ENDIF}
    {$IFDEF CPU64BITS}
      {$DEFINE CPU64}
    {$ENDIF}
  {$ENDIF}

  {$IFDEF CPU32}
    {$IFDEF AVX}
      {$DEFINE AVX32}
    {$ENDIF}
    {$IFDEF AVX2}
      {$DEFINE AVX32}
    {$ENDIF}
  {$ENDIF}

  {$IFDEF CPU64}
    {$IFDEF AVX}
      {$DEFINE AVX64}
    {$ENDIF}
    {$IFDEF AVX2}
      {$DEFINE AVX64}
    {$ENDIF}
  {$ENDIF}

  {$IFDEF AVX512}
    {$DEFINE AVX64}
    {$UNDEF AVX32}
  {$ENDIF}

  {$IFDEF AVX32}
    {$DEFINE AVXANY}
  {$ENDIF}

  {$IFDEF AVX64}
    {$DEFINE AVXANY}
  {$ENDIF}
  ```
  AVX/AVX2/AVX512 are never defined by the code itself — they must be passed externally
  via the compiler command line (e.g. `-dAVX2`). The block above only derives the
  secondary flags (`AVX32`, `AVX64`, `AVXANY`) and disables AVX on non-x86 targets.

- [x] **0.1b** Define `TUInt128` as a pure record in `pascoremathtypes.pas`:
  ```pascal
  type TUInt128 = record
    lo, hi: UInt64;
  end;
  ```

- [x] **0.2** Implement `MulWide(a, b: UInt64): TUInt128` in `pascoremathtypes.pas`.
  - Primary path: x86-64 inline assembly using the `MUL` instruction (`rdx:rax = rax * src`).
  - Portable fallback: four 32-bit partial products for non-x86-64 targets (ARM, etc.).

- [x] **0.3** Overload `+` for `TUInt128 + UInt64 → TUInt128` (implements the `p1 += p0>>64` pattern as `p1 := p1 + p0.hi`).
  This is a genuine addition with carry propagation — implement as:
  ```pascal
  operator+(const a: TUInt128; b: UInt64): TUInt128; inline;
  begin
    Result.lo := a.lo + b;
    Result.hi := a.hi + UInt64(Result.lo < b);  // carry
  end;
  ```
  Mark `inline` — it is a hot-path operation called inside `MulWide`. Also try marking `MulWide` itself as `inline`; trust the compiler to handle it.

- [x] **0.4** Define type-punning records in `pascoremathtypes.pas`:
  ```pascal
  type Tb32u32 = record case Boolean of
    False: (f: Single);  True: (u: LongWord); end;
  type Tb64u64 = record case Boolean of
    False: (f: Double);  True: (u: UInt64);   end;
  ```

- [x] **0.5** Implement builtin equivalents in `pascoremathtypes.pas`:
  | C | Pascal |
  |---|---|
  | `__builtin_roundeven(x)` | `RoundEven(x: Double): Double` |
  | `__builtin_copysign(x,y)` | `Math.CopySign` or manual |
  | `__builtin_ctz(n)` | `BsfDWord(n)` (BSF instruction) |
  | `__builtin_fmax(x,y)` | NaN-aware `FMax(x,y: Double)` |
  | `__builtin_fmin(x,y)` | NaN-aware `FMin(x,y: Double)` |
  | `__builtin_fma(a,b,c)` | `Math.FMA` (FPC 3.2+); verify rounding |
  | `__builtin_inff()` | `Infinity` (from `Math`) |
  | `__builtin_expect(e,v)` | Delete entirely (branch-hint only) |
  | `__attribute__((noinline))` | `[noinline]` or `{$INLINE OFF}` |

- [x] **0.6** Write a hex-float conversion utility (script or small tool) to convert C99
  hex float literals (e.g. `0x1.62e42fefa39efp-1`) to Pascal `Double` / `UInt64`
  constants. All lookup tables must be converted using this tool, not by hand.

- [x] **0.7** Set up a test harness that compiles and runs both the C reference and the
  Pascal implementation, then compares results bit-for-bit for all 2^32 `Single`
  inputs (exhaustive). For bivariate functions, agree on a sampling strategy.

- [x] **0.8** Set up `Benchmark.pas`: single-threaded, 10 million calls per function, two
  independent timed loops. Reports Mops/s for `cr_*` (C) and `pcr_*` (Pascal) and
  validates agreement via XOR sink. Pattern:
  ```pascal
  sink_c   := 0;
  sink_pas := 0;

  t1 := Now;
  for i := 0 to 10_000_000 - 1 do
  begin
    v.u    := i * stride;
    sink_c := sink_c xor Tb32u32(cr_rsqrtf(v.f)).u;
  end;
  t2 := Now;
  for i := 0 to 10_000_000 - 1 do
  begin
    v.u      := i * stride;
    sink_pas := sink_pas xor Tb32u32(pcr_rsqrtf(v.f)).u;
  end;
  t3 := Now;

  if sink_c = sink_pas then WriteLn('OK') else WriteLn('MISMATCH');
  WriteLn('C:      ', 10_000_000 / MillisecondsBetween(t2,t1) / 1000 :0:2, ' Mops/s');
  WriteLn('Pascal: ', 10_000_000 / MillisecondsBetween(t3,t2) / 1000 :0:2, ' Mops/s');
  ```
  The sink doubles as a quick sanity check — `MISMATCH` means the two implementations
  disagree on at least one input in the sample. Full correctness is validated by
  `TestHarness.pas` (exhaustive 2^32).

---

## Phase 1 — Simple univariate (no u128, no FMA, ≤ 130 lines)

Port in this order. All functions live in `pascoremath.pas`, named `pcr_<name>f`.

- [x] **1.01** `rsqrt`   — 89 lines
- [x] **1.02** `tanh`    — 89 lines
- [x] **1.03** `atanpi`  — 106 lines
- [x] **1.04** `cospi`   — 114 lines
- [x] **1.05** `acos`    — 115 lines
- [x] **1.06** `cbrt`    — 117 lines
- [x] **1.07** `sinpi`   — 117 lines
- [x] **1.08** `atan`    — 118 lines
- [x] **1.09** `asin`    — 120 lines
- [x] **1.10** `acospi`  — 126 lines  *(uses FMA — verify `Math.FMA` correctness first)*
- [x] **1.11** `log2`    — 126 lines
- [x] **1.12** `asinpi`  — 128 lines  *(uses FMA)*
- [x] **1.13** `tanpi`   — 130 lines
- [x] **1.14** `cosh`    — 132 lines

---

## Phase 2 — Medium univariate

- [x] **2.01** `log`       — 133 lines
- [x] **2.02** `exp2`      — 138 lines
- [x] **2.03** `log1p`     — 140 lines
- [x] **2.04** `exp2m1`    — 146 lines
- [x] **2.05** `expm1`     — 146 lines
- [x] **2.06** `exp10`     — 150 lines
- [x] **2.07** `log10`     — 150 lines
- [x] **2.08** `erfc`      — 151 lines
- [x] **2.09** `log2p1`    — 151 lines
- [x] **2.10** `erf`       — 152 lines
- [x] **2.11** `sinh`      — 156 lines
- [x] **2.12** `exp`       — 115 lines *(listed here due to two-pass rounding logic)*
- [x] **2.13** `atanh`     — 158 lines
- [x] **2.14** `exp10m1`   — 158 lines
- [x] **2.15** `log10p1`   — 162 lines
- [x] **2.16** `asinh`     — 168 lines
- [x] **2.17** `acosh`     — 173 lines

---

## Phase 3 — Longer special univariate

- [x] **3.01** `tgamma`  — 205 lines
- [x] **3.02** `lgamma`  — 259 lines

---

## Phase 4 — Bivariate and compound (FMA-heavy)

- [x] **4.01** `hypot`    — 282 lines  *(uses FMA, fenv rounding modes)*
- [x] **4.02** `atan2`    — 231 lines  *(uses FMA)*
- [x] **4.03** `atan2pi`  — 190 lines  *(uses FMA)*
- [x] **4.04** `compound` — 611 lines  *(uses FMA, fenv)*
- [x] **4.05** `pow`      — 325 lines  *(uses FMA, fenv)*

---

## Phase 5 — u128 functions (hardest)

Depends on Phase 0.2–0.3 being fully correct. Validate `MulWide` independently before
starting this phase.

- [x] **5.01** `sin`    — 222 lines
- [x] **5.02** `cos`    — 206 lines
- [x] **5.03** `sincos` — 245 lines
- [x] **5.04** `tan`    — 199 lines

---

## Per-function porting checklist

Apply this checklist to every function before marking it done:

- [ ] Hex float literals converted via the conversion utility (0.6), not by hand
- [ ] Lookup tables moved to unit-level `const` (no `static` locals)
- [ ] All type-punning uses `Tb32u32` / `Tb64u64` records from `pascoremathtypes.pas` (no unsafe casts)
- [ ] `__builtin_expect` wrappers removed entirely
- [ ] `__attribute__((noinline))` replaced with `[noinline]`
- [ ] `CORE_MATH_SUPPORT_ERRNO` blocks omitted (out of scope for Pascal port)
- [ ] Exhaustive test passes (bit-exact match against C reference for all inputs)
- [ ] C function `cr_<name>f` declared in `ccoremath.pas`; Pascal equivalent named `pcr_<name>f` in `pascoremath.pas`

---

## Architectural notes and known pitfalls

1. **`rbig()` must not be duplicated.** The large-argument range-reduction helper `rbig()` is
   byte-for-byte identical in `sin`, `cos`, `tan`, and `sincos`. It must live once in
   `pascoremathtypes.pas` and be called by all four functions in `pascoremath.pas`, not duplicated.
   Duplication would make any future bug fix require four parallel edits.

2. **`sincos` has a different API — an explicit decision is required.** The C signature is
   `void cr_sincosf(float x, float *sout, float *cout)` — two output pointers, no return value.
   In Pascal this must be declared as one of:
   - `procedure cr_sincosf(x: Single; out s, c: Single)` — closest to the C original, or
   - a function returning a small record `TSinCos = record s, c: Single end`.
   Pick one convention and apply it consistently before porting Phase 5.

3. **Rounding mode support needs a Phase 0 infrastructure task.** Phase 4 functions (`hypot`,
   `compound`, `pow`, `atan2`, `atan2pi`) call `fesetround`/`FE_TONEAREST` from `<fenv.h>`.
   The Pascal equivalent is `SetRoundMode` from the `Math` unit. This must be evaluated and
   wrapped in `pascoremathtypes.pas` before Phase 4 starts — it is not safe to defer until
   then.

4. **`BsfDWord` is 32-bit only — use `BsfQWord` for 64-bit values.** The builtins table maps
   `__builtin_ctz(n)` to `BsfDWord`, which is correct when `n` is a `LongWord`. If `ctz` is
   ever called on a `UInt64` value, `BsfQWord` must be used instead. Check the argument type
   at each call site; using the wrong variant silently operates on only the low 32 bits.

5. **`roundeven_finite` is architecturally complex and deserves its own sub-task.** The C
   source has four distinct implementations selected at compile time: AVX, SSE4.1, ARMv8, and
   a portable software fallback using bit manipulation. Task 0.5 lists it as a single line but
   it warrants the same treatment as `MulWide` (task 0.2): implement the software fallback
   first, then add an x86-64 SSE4.1/AVX path, and validate both against the C reference.

6. **Two-pass rounding (Ziv's strategy) — the slow path must never be removed.** Several
   functions (`exp` and others) compute a fast approximation and then check whether the result
   is close enough to a rounding boundary. If it is not, a slower, higher-precision second pass
   is used to resolve the ambiguity. This is not dead code — it is the mechanism that guarantees
   correct rounding. A developer unfamiliar with the pattern may remove the slow path believing
   it is unreachable. It must be kept exactly as ported from the C source.

---

## Design decisions

1. **`TUInt128` is a plain record — no variant/case.** Defined as `record lo, hi: UInt64 end`.
   The variant form was considered but rejected in favour of simplicity and explicit field access.

2. **`MulWide` is a named function, not an overloaded `*` operator.** FPC already owns the
   signature `UInt64 * UInt64 → UInt64` and will not allow a second overload with a different
   return type. A named function is the only option.

3. **Only `+` needs to be overloaded for `TUInt128`, not a full shift operator.** The only
   shift used in the C code is `p0>>64`, which always shifts by exactly 64 — meaning it simply
   reads the high word. The pattern `p1 += p0>>64` therefore becomes `p1 := p1 + p0.hi`, so a
   general 128-bit shift operator is not needed and should not be added.

4. **Phase 5 (sin/cos/tan/sincos) is explicitly blocked on Phase 0.2–0.3 being fully validated.**
   A silent bug in `MulWide` or the `+` overload would corrupt all four functions with no
   obvious failure mode, since the error only manifests for large-argument inputs that trigger
   `rbig()`.

5. **The per-function checklist is intended to double as a PR template.** Copy it into each
   pull request description so reviewers can confirm each step was completed before merging.

---

## Key rules for the developer

1. **Do not change the algorithm.** This is a faithful port, not a rewrite or optimization.
   The C source is the specification. If your output differs by even one ULP, it is a bug.

2. **Convert hex float constants systematically.** Never retype a constant by hand.
   Use the conversion utility from task 0.6 for every lookup table.

3. **`MulWide` must be correct before Phase 5.** Validate it with known pairs covering
   overflow, zero, and maximum-value cases before porting any u128 function.

4. **FMA correctness is critical.** Verify that FPC's `Math.FMA` produces correctly-
   rounded results. If not, implement a software FMA before touching Phase 4.

5. **`__builtin_expect` is a no-op hint.** Simply delete it; do not replace it with
   anything.

6. **Test exhaustively.** For single-argument functions, all 2^32 float inputs is
   feasible. Run exhaustive tests before moving to the next function.

7. **Work sequentially within each phase.** The ordering within each phase is chosen
   to build familiarity gradually. Do not skip ahead.

8. **ASM is allowed and encouraged where beneficial.** Use inline assembly for
   performance-critical infrastructure (`MulWide`, `BsfDWord`/`BsfQWord`, etc.).
   Always provide a portable Pascal fallback for non-x86-64 targets.

9. **Inline everything possible.** Mark all small helpers, wrappers, type-punning
   accessors, and operator overloads as `inline`. The compiler ignores the hint when
   inlining is not beneficial — the cost of marking `inline` unnecessarily is zero.

10. **Benchmark every function.** After each function passes exhaustive testing, run
    `Benchmark.pas` and record the Mops/s ratio (Pascal vs C). A large gap is a signal
    to investigate missed inlining or suboptimal code generation.

---

## Open bugs (detected 2026-04-11, `--pct 1` sampling)

### Bug A — `pcr_compoundf`: returns +Inf for tiny subnormal x with large y

**Symptom:** Large error — C returns `0x3F800000` (1.0), Pascal returns `0x7F800000` (+Inf).
1 903 717 mismatches out of 10 M sampled pairs.

**Example inputs:**
```
x=$000001AD (~6.0e-43, tiny subnormal)   y=$55555703 (~1.47e13, large positive)
x=$0000035A                               y=$555558B1
```

**Expected result:** `compound(x, y) = (1+x)^y ≈ 1.0`
because `x*y ≈ 6e-43 * 1.47e13 ≈ 8.8e-30`, so `exp(y * log(1+x)) ≈ exp(8.8e-30) ≈ 1.0`.

**Root cause (suspected):** The Pascal `pcr_compoundf` implementation takes a wrong code path
for very small subnormal `x` values. Instead of recognising that `(1+x)^y ≈ 1`, it overflows
somewhere in the intermediate computation and returns +Inf.
The C source (`core-math/src/binary32/compound/compoundf.c`, ~1110 lines) has a dedicated
subnormal handling path (search for `// subnormal numbers`); it is likely that the Pascal port
is missing a guard or has an incorrect comparison for this range.

**How to debug:**
- Add `--diag 10` to the TestHarness bivariate output (already supported as of commit `25d897a`)
  and focus on tiny subnormal x (bits 0x0000_0001 to roughly 0x007F_FFFF) with large positive y.
- Compare `pcr_compoundf` step-by-step with C for input `x=$000001AD, y=$55555703`.
- Grep for `subnormal` in `compoundf.c` and verify the corresponding Pascal guard conditions.

---

### Bug B — `pcr_powf`: 1-ULP rounding error for large x, tiny negative y

**Symptom:** Small error — C returns `0x3F7FFFFF` (1.0 − 1 ULP), Pascal returns `0x3F800000` (1.0).
458 657 mismatches out of 10 M sampled pairs.

**Example inputs:**
```
x=$5ACCA329 (~2.88e16, large positive)   y=$B058276B (~−7.86e-10, tiny negative)
x=$5ACCA4D6                              y=$B0582919
(pattern: x large ~2.88e16, y small negative ~−7.86e-10, with slow drift)
```

**Expected result:** `x^y = exp(y * log(x)) ≈ 1 − 2.98e-8`, which is just below the
float32 midpoint between `0x3F7FFFFF` and `0x3F800000`, so the correct rounding is `0x3F7FFFFF`.

**Investigation so far:**
- All lookup tables (`ix_pf`, `lix_pf`, `c_pf`, `ce_pf`, `tb_pf`) were verified bit-exact
  against the C source.
- The fast path (no `accurate2` call) is taken for these inputs — the borderline check
  `((rr.u + 468) & 0xFFFFFFF) <= 936` is False (check value ≈ 268 M >> 936).
- Python simulation of the full fast path (with both true IEEE FMA and 80-bit-extended
  approximation) produces `0x3F7FFFFF` (= correct C answer) for the example input.
- Therefore the Pascal computation must diverge from the simulation somewhere during
  actual FPC code generation — likely due to x87 excess precision in an intermediate
  variable that is nominally `Double` but kept in an 80-bit x87 register.

**Root cause (suspected):** `pcr_fma` is a double-rounding approximation
(`Double(Extended(x)*Extended(y)+Extended(z))`), not a true IEEE FMA. When FPC inlines it,
the `Double(...)` cast may not force a spill to memory, leaving intermediate values at
80-bit extended precision. This excess precision in `z_pf` (the reduced argument) or `h_pf`
propagates through the exponential polynomial and shifts `rr_pf` to the wrong side of the
float32 midpoint.

**How to fix:**
- Option 1 (preferred): Replace `pcr_fma` with a true hardware FMA using FPC inline assembly:
  ```pascal
  function pcr_fma(x, y, z: Double): Double; inline;
  begin
    asm
      vmovsd xmm0, x
      vfmadd213sd xmm0, y, z   // xmm0 = x*y + z
      vmovsd Result, xmm0
    end;
  end;
  ```
  Requires AVX/FMA3 (available on all modern x86_64). Add a `{$IFDEF CPUX86_64}` guard
  with the 80-bit fallback for other architectures.
- Option 2: Force intermediate variables to be stored to memory (defeating x87 register
  caching) by adding `volatile`-style stores — but FPC has no standard mechanism for this.
- Option 3: Compile with `-CfAVX2` or `-CfSSE4` to enable FMA instructions globally (risky
  for portability, and currently `{$FPUTYPE SSE64}` already fails to prevent x87 usage for
  untyped literals).

**Note:** The same `pcr_fma` issue may affect `pcr_compoundf` and potentially other functions
in edge cases not caught at 1% sampling. Fix `pcr_fma` first, then re-run `--pct 1` to check.
