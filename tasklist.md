# pas-core-math Task List

Port of the CORE-MATH binary32 (single-precision) library to Free Pascal.
Goal: bit-exact, correctly-rounded results matching the C reference for all 2^32 float inputs.

---

## Status summary

- **All 42 functions ported and committed** (branch: `a1`)
- **Bug-fix round in progress** — detected with `bin/TestHarness32 --pct 1`
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
│   ├── pascoremath32.pas       # Pascal implementations (pcr_* functions)
│   ├── ccoremath32.pas         # C reference external declarations (cr_* functions)
│   ├── laz-project/
│   │   ├── pas-core-math.lpi
│   │   └── pas-core-math.lps
│   ├── tools/
│   │   └── HexFloatConvert.pas
│   └── tests/
│       ├── TestHarness32.pas
│       ├── TestMulu64u64.pas
│       └── Benchmark32.pas     # single-threaded Mops/s comparison: cr_* (C) vs pcr_* (Pascal)
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

- [x] **0.2** Implement `Mulu64u64(a, b: UInt64): TUInt128` in `pascoremathtypes.pas`.
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
  Mark `inline` — it is a hot-path operation called inside `Mulu64u64`. Also try marking `Mulu64u64` itself as `inline`; trust the compiler to handle it.

- [x] **0.4** Define type-punning records in `pascoremathtypes.pas`:
  ```pascal
  type Tb32u32 = record case Boolean of
    False: (f: Single);  True: (u: UInt32); end;
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

- [x] **0.8** Set up `Benchmark32.pas`: single-threaded, 10 million calls per function, two
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
  `TestHarness32.pas` (exhaustive 2^32).

---

## Phase 1 — Simple univariate (no u128, no FMA, ≤ 130 lines)

Port in this order. All functions live in `pascoremath32.pas`, named `pcr_<name>f`.

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

Depends on Phase 0.2–0.3 being fully correct. Validate `Mulu64u64` independently before
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
- [ ] C function `cr_<name>f` declared in `ccoremath32.pas`; Pascal equivalent named `pcr_<name>f` in `pascoremath32.pas`
- [ ] x87-avoidance pass applied (Phase 7) — `tools/x87_audit.py` reports zero hits for the function. **Cast-selection rule**: pick the cast that matches the assignment target, not the literal in isolation. Accurate-path expressions assigning to `Double` accumulators want `Double(1.0)` even inside a binary32 function; fast-path `Single` helpers (e.g. `ir: array[0..1] of Single`) want `Single(1.0)`. Blanket `Single(...)` would silently downgrade dd/Horner refinement precision.

---

## Architectural notes and known pitfalls

1. **`rbig()` must not be duplicated.** The large-argument range-reduction helper `rbig()` is
   byte-for-byte identical in `sin`, `cos`, `tan`, and `sincos`. It must live once in
   `pascoremathtypes.pas` and be called by all four functions in `pascoremath32.pas`, not duplicated.
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
   `__builtin_ctz(n)` to `BsfDWord`, which is correct when `n` is a `UInt32`. If `ctz` is
   ever called on a `UInt64` value, `BsfQWord` must be used instead. Check the argument type
   at each call site; using the wrong variant silently operates on only the low 32 bits.

5. **`roundeven_finite` is architecturally complex and deserves its own sub-task.** The C
   source has four distinct implementations selected at compile time: AVX, SSE4.1, ARMv8, and
   a portable software fallback using bit manipulation. Task 0.5 lists it as a single line but
   it warrants the same treatment as `Mulu64u64` (task 0.2): implement the software fallback
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

2. **`Mulu64u64` is a named function, not an overloaded `*` operator.** FPC already owns the
   signature `UInt64 * UInt64 → UInt64` and will not allow a second overload with a different
   return type. A named function is the only option.

3. **Only `+` needs to be overloaded for `TUInt128`, not a full shift operator.** The only
   shift used in the C code is `p0>>64`, which always shifts by exactly 64 — meaning it simply
   reads the high word. The pattern `p1 += p0>>64` therefore becomes `p1 := p1 + p0.hi`, so a
   general 128-bit shift operator is not needed and should not be added.

4. **Phase 5 (sin/cos/tan/sincos) is explicitly blocked on Phase 0.2–0.3 being fully validated.**
   A silent bug in `Mulu64u64` or the `+` overload would corrupt all four functions with no
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

3. **`Mulu64u64` must be correct before Phase 5.** Validate it with known pairs covering
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
   performance-critical infrastructure (`Mulu64u64`, `BsfDWord`/`BsfQWord`, etc.).
   Always provide a portable Pascal fallback for non-x86-64 targets.

9. **Inline everything possible.** Mark all small helpers, wrappers, type-punning
   accessors, and operator overloads as `inline`. The compiler ignores the hint when
   inlining is not beneficial — the cost of marking `inline` unnecessarily is zero.

10. **Benchmark every function.** After each function passes exhaustive testing, run
    `Benchmark32.pas` and record the Mops/s ratio (Pascal vs C). A large gap is a signal
    to investigate missed inlining or suboptimal code generation.

---

## Phase 6 — Performance optimization

Benchmark baseline (2026-04-11, FPC 3.2.2 -O2, x86_64 Linux):

| Function   | C (Mops/s) | Pascal (Mops/s) | Ratio |
|------------|-----------|-----------------|-------|
| lgammaf    | 151.5     | 35.3            | 4.3×  |
| atan2f     | 277.8     | 88.5            | 3.1×  |
| tgammaf    | 333.3     | 133.3           | 2.5×  |
| expf       | 357.1     | 144.9           | 2.5×  |
| sincosf    | 212.8     | 98.0            | 2.2×  |
| exp10f     | 322.6     | 151.5           | 2.1×  |
| cosf       | 232.6     | 119.0           | 2.0×  |
| sinf       | 227.3     | 123.5           | 1.8×  |
| tanf       | 208.3     | 113.6           | 1.8×  |

- [x] **6.1** Replace `pcr_roundeven` / `pcr_roundevenf` with SSE4.1 `ROUNDSD` / `ROUNDSS`
  - Hot path of `sincos_rltl0`, `sincos_rltl`, `tanf_rltl`, `exp10f`, `exp2m1f` and others
  - Current Pascal implementation is ~30 lines with multiple branches; the instruction is 1 cycle
  - Use `{$IFDEF AVX2}` guard; keep existing bit-manipulation code as the `{$ELSE}` fallback
  - imm8 = 12 (0x0C): override MXCSR with round-to-nearest-even, suppress precision exception
  - **Primary fix for sinf / cosf / sincosf / tanf gap**

- [x] **6.2** Replace `pcr_fma` / `pcr_fmaf` with FMA3 `VFMADD213SD` / `VFMADD213SS`
  - Current implementation uses 80-bit `Extended`, forcing x87 mode on every call despite `{$FPUTYPE SSE64}`
  - `pcr_fma(x,y,z)` is called heavily in `muldd`/`polydd` (atan2f), `powf` (~20 calls), `compoundf`, and others
  - Hardware FMA is also *correctly rounded* — the Extended path is only an approximation for doubles
  - Use `{$IFDEF AVX2}` guard; keep Extended fallback for non-x86-64 targets
  - `VFMADD213SD xmm0, xmm1, xmm2`: xmm0 = xmm0×xmm1 + xmm2
  - **Primary fix for atan2f / powf gap; also expected to fix Bug B (pcr_powf rounding error)**

- [x] **6.3** Inline `lgamma_as_r7` / `lgamma_as_r8` into `pcr_lgammaf`
  - These helpers use `const c: array of Double` (open array) which FPC cannot inline
  - The caller passes compile-time constant arrays (`rn_sm`, `rd_sm`, `rn_md`, `rd_md`) that become
    constant loads when inlined — open arrays force pointer indirection instead
  - Replace by inlining the 7- and 8-term product expressions directly at each call site in `lgammaf`
  - **Primary fix for lgammaf 4.3× gap**

- [x] **6.4** Upgrade FPC build flag from `-O2` to `-O3` in `build.sh`
  - GCC `-O2` is far more aggressive than FPC `-O2`; FPC `-O3` closes some of the gap
  - Also consider adding `-CpCOREI` or `-CpCOREI7` for CPU-specific scheduling
  - Small across-the-board gain (~5–10%), no correctness risk

- [x] **6.5** Replace `pcr_fmax` / `pcr_fmin` with `MAXSD` / `MINSD` (or branchless bit-ops)
  - Current implementation branches through `IsNan()` on every call
  - C `fmax`/`fmin` compile to a single `MAXSD`/`MINSD` instruction
  - Note: `MAXSD` propagates the *second* operand on NaN input — verify call sites tolerate this
    before switching; use branchless bit-op fallback where NaN semantics must be exact

- [x] **6.6** Reduce `tgammaf` reduction loop overhead
  - For non-integer negative x, `tgammaf` uses `while j_tg > 0 do w_tg *= z_tg` — up to ~44 iterations
  - Uniformly distributed benchmark inputs frequently hit this path
  - Consider precomputed product table or logarithm-based reduction for large |ii|

---

## Phase 7 — x87-avoidance pass (post-port optimization, binary32)

Phase 7 mirrors `tasklist64.md` Phase 6 (Pillars A/B/C) but targets `src/pascoremath32.pas`.
The goal is the same — keep FPC's code generator on SSE and out of the x87 FPU — but the
shape of the work is different in three important ways, documented up front so a future
pass doesn't blanket-apply the 64-bit recipe:

1. **Monolithic file.** All 32-bit function bodies live inside `src/pascoremath32.pas`.
   There is no `_port_32.inc` / `_const_32.inc` split and no `tools/*_gen.py` script,
   so the "auto-generated reverts" footgun from 64-bit Phase 6 does not apply.
2. **Cast type depends on the surrounding expression.** 64-bit Pillar C uses `Double(...)`
   uniformly. In `pascoremath32.pas`, accurate paths declare coefficient arrays as
   `array[0..N] of Double` even when the function consumes/produces `Single` — the bare
   `1.0` in those expressions wants `Double(1.0)`, not `Single(1.0)`. Pick the cast that
   matches the variable being assigned, not the literal in isolation. Single-typed
   fast-path helpers (e.g. `ir: array[0..1] of Single`) want `Single(...)`.
3. **`tools/x87_audit.py` has 32-bit-blind spots.** Its [B] heuristic looks for the
   `Tb64u64`-style `cFoo[3].f` pattern (literal index + `.f` field). The 32-bit code
   uses plain `array[0..N] of Double`, so literal-indexed reads have no `.f` suffix
   and are silently skipped. As of 2026-04-26, baseline `--summary` reports
   `pascoremath32.pas: A=0 B=0 C=422` — the C count is real, the B count is a
   heuristic miss (actual: 580 literal-indexed reads across 73 distinct array names).
   Subtask 7.0 must extend the auditor before any [B] pass-gate becomes meaningful.

### Three pillars (same as 64-bit Phase 6, with the deltas above)

- **Pillar A** — unroll fixed-trip-count polynomial loops. Already baseline-clean for
  binary32 (only 6 `for` loops in the entire file; audit A=0). Spot-check during each
  subtask, no dedicated work expected.
- **Pillar B** — lift literal-indexed const-array reads to named scalars. Real work:
  ~580 reads across ~73 arrays. Apply only to the small fixed-trip-count coefficient
  tables; runtime-indexed lookup tables (`S_TABLE[0..127]`, `lix_asinh_acosh[0..128]`,
  `c_table[…]`, the `array[0..63/64/127/128/157]` families) are out of scope — same
  guard as `tasklist64.md:1183`.
- **Pillar C** — `Single(...)` / `Double(...)` typecast on every bare float literal,
  using the cast-selection rule above. ~422 hits per audit baseline.

### Subtasks

- [x] **7.0** Tooling & policy
  - Extended `tools/x87_audit.py` with `PLAIN_ARRAY_DECL_RE` /
    `PLAIN_ARRAY_READ_RE`: a pre-pass collects names declared as
    `array[0..N] of {Single,Double}`, and any `<name>[<int>]` read (no
    `.f` suffix) on those names is flagged as [B]. Declaration lines are
    suppressed so the `array[0..N]` token in the decl itself is not
    self-matched.
  - Cast-selection rule documented in the per-function checklist (above):
    cast matches the assignment target's type, not the literal.
  - Re-baselined `--summary src/pascoremath32.pas` after the rule lands.
  - **Baseline (2026-04-26, pre-rule auditor):** `A=0 B=0 C=422`.
  - **Baseline (2026-04-26, post-rule auditor):** `A=0 B=571 C=422`.
    The new [B] rule recovers the 32-bit literal-indexed reads that were
    invisible under the Tb64u64 `.f`-suffix heuristic. (User-provided
    manual estimate was 580; the 9-hit gap is in the noise — likely
    multi-dim or commented-out reads not covered by the simple regex.)
  - 64-bit files unaffected: `pascoremath64.pas: A=1 B=0 C=8`,
    `pascoremathtypes.pas: A=0 B=0 C=30` after the change — no false
    positives on the 64-bit code-base.

- [x] **7.1** exp family — `expf`, `exp2f`, `exp10f`, `expm1f`, `exp2m1f`,
  `exp10m1f` regions of `pascoremath32.pas`. Pillar B done: lifted
  `b`/`c` (exp2f), `c_fast`/`ch`/`b_small` (expm1f), `b_exp10`/`c_exp10`
  (exp10f), `c_e10` and `cp4_e10`..`cp9_e10` (exp10m1f) to named scalars.
  `expf` was already pre-lifted (`c_exp_*`/`b_exp_*`); `exp2m1f` already
  uses inline `c0v`..`c7v` named scalars per the `c_table` regions.
  Audit B count: 571 → 488 (-83 reads). Tests: all 42 functions pass at
  `--pct 1`. Bench (taskset -c 1, AVX2): expf 232 / exp2f 163 / exp10f
  144 / expm1f 138 / exp2m1f 239 / exp10m1f 209 Mops/s — all ≥ C
  reference except expm1f (138 vs C 161 — still the largest gap, may
  reward a future Pillar C pass).

- [x] **7.2** log family — `logf`, `log10f`, `log1pf`, `log2f`, `log10p1f`,
  `log2p1f` regions. Pillar B done: lifted `bcoef`/`ccoef` (log2f),
  `c[7]` and `tl[0]` (logf — `cL0..cL6`, `logf_tl0`), `b[7]`/`c[3]`
  (log1pf — `b1p1..b1p7`, `c1p0..c1p3`), `b10`/`c10`/`tl10[0]` (log10f),
  `c_l2p1[6]` (log2p1f), `c_l10[9]` (log10p1f) to named scalars. The
  large runtime-indexed tables `tr`, `tl`, `tl10`, `tr10`, `ix_l2p1`,
  `lix_l2p1`, `tr_l10`, `tl_l10`, `x0`, `lix` remain arrays (out of
  scope — runtime-indexed). Audit B count: 488 → 427 (-61 reads).
  Tests: all 42 functions pass at `--pct 1`.
  Bench (taskset -c 1, AVX2): log2f 390 / logf 364 / log10f 340
  (FASTER than C 271) / log1pf 298 (TIE) / log2p1f 201 / log10p1f 226
  Mops/s. Note `c_aln` is in `pcr_powf` (misc family, 7.7), not log.

- [ ] **7.3** trig family — `sinf`, `cosf`, `sincosf`, `tanf`, `sinpif`,
  `cospif`, `tanpif` regions. Pillar B candidates: `sn`, `cn`, the small
  Estrin/Horner blocks. Watch for the `cn_a2` / `cd_a2` / `cn_a2p` / `cd_a2p`
  cross-listing under 7.4.

- [ ] **7.4** inverse-trig family — `atanf`, `atan2f`, `asinf`, `acosf`,
  `atanpif`, `atan2pif`, `asinpif`, `acospif`, `atanhf` regions. Pillar B
  candidates: `cn_a2`, `cd_a2`, `cn_a2p`, `cd_a2p`, `c_atanh_s`,
  `c_atanh_acc`, `c_near_ac`. The `cn_a2p` / `cd_a2p` Estrin block is the
  exemplar that motivated this phase — full lift to `cn_a2p_<i>` /
  `cd_a2p_<i>` named scalars, drop the array if no other reader.

- [ ] **7.5** hyperbolic family — `sinhf`, `coshf`, `tanhf`, `asinhf`,
  `acoshf` regions. Pillar B candidates: `c_sinh`, `ch_sinh`, `cp_sinh`,
  `c_asinh`, `cm_asinh`, `cp_asinh`, `cm_acosh`, `cp_acosh`. Note:
  `lix_asinh_acosh[0..128]` is runtime-indexed → out of scope.

- [ ] **7.6** special-functions family — `erff`, `erfcf`, `tgammaf`,
  `lgammaf`, `lgammaf_as_r7`/`r8` regions. Pillar B candidates:
  `c_erf_small`, `c_tg`, `stir2` / `stir4` / `stir8`, `rn_sm` / `rd_sm`
  / `rn_md` / `rd_md` (already inlined per 6.3 — verify whether the
  underlying arrays still exist post-inline), `CF_P1C` / `CF_P2C` /
  `CF_Q1C` / `CF_Q2C` / `CF_ERR_E22`. **Cross-cutting catch:** if any
  of these are referenced from two functions (analogous to 64-bit
  `cAcoshC` / `cTanhExpCh` sharing — `tasklist64.md:1404–1411`), update
  every reader in the same commit.

- [ ] **7.7** miscellaneous — `hypotf`, `cbrtf`, `rsqrtf`, `powf`,
  `compoundf`, plus any helpers (`muldd_f`, `polydd_f`, etc.). Spot-fix
  Pillar C; B unlikely to apply outside `powf`. Verify
  `tools/x87_audit.py` reports zero hits before close.

### Subtask template (per function family)

1. Run `python3 tools/x87_audit.py src/pascoremath32.pas` to get a current count;
   note the baseline.
2. Apply Pillar A — unroll any audit-flagged loops in scope. Commit.
   Run `taskset -c 1 env Benchmark32 <fn>` + `TestHarness32 --pct 1`.
3. Apply Pillar B — lift literal-indexed reads to named scalars. Drop the array
   only after grepping `pascoremath32.pas` for other readers (cross-function
   sharing risk per 64-bit Phase 6.4/6.5). Commit. Re-run.
4. Apply Pillar C — typecast sweep using the Single/Double cast-selection rule.
   Commit. Re-run.
5. Verify the audit reports zero hits for the touched function(s). Record
   Mops/s delta in the subtask note.

### Bench expectation

The 64-bit Phase 6 gains came largely from Pillar B unblocking constant folding
(see `tasklist64.md:1413–1420`, +88 % sinh, +83 % acosh, +33 % expm1). For 32-bit,
expect:
- **Pillar B**: meaningful gains on the small-coefficient families (atan2f, the
  Estrin blocks, the `c_atanh_acc` / `c_*_sm` Horner chains). This is the main
  perf lever.
- **Pillar C**: flat-to-small-positive — primarily codegen-correctness insurance.
  A flat bench post-Pillar-C is not a failure signal.
- **Pillar A**: N/A (already clean).

### Sharp edges to remember (binary32-specific)

- **Cast type matches the variable, not the literal.** `Double(1.0)` inside an
  accurate-path expression assigning to a `Double` accumulator; `Single(1.0)`
  inside a fast-path expression assigning to a `Single`. Blanket-replacing with
  `Single(...)` would silently downgrade precision in the dd/Horner refinement
  paths.
- **No generators.** Edits stick — no regen step to keep in sync. Conversely,
  there is no scripted way to mass-rewrite; every lift is a hand edit.
- **Out-of-scope tables stay arrays.** `S_TABLE`, `lix_asinh_acosh`, `lix_l2p1`,
  `c_table`, the `array[0..63/64/127/128/157]` families are runtime-indexed and
  must remain arrays. Lifting them is wrong (and impossible without unrolling
  the indexing function).

---

## Closed bugs (detected 2026-04-11, `--pct 1` sampling)

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

**How to debug:**
- Add `--diag 10` to the TestHarness32 bivariate output (already supported as of commit `25d897a`)
  and focus on tiny subnormal x (bits 0x0000_0001 to roughly 0x007F_FFFF) with large positive y.
- Compare `pcr_compoundf` step-by-step with C for input `x=$000001AD, y=$55555703`.
- Grep for `subnormal` in `compoundf.c` and verify the corresponding Pascal guard conditions.

This bug has been fixed.

---

### Bug B — `pcr_powf`: 1-ULP rounding error for large x, tiny negative y
Fixed: https://github.com/joaopauloschuler/pas-core-math/commit/909abae63d2d88559bc125d0c96fec873e2d4a79
