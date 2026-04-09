# pas-core-math Task List

Port of the CORE-MATH binary32 (single-precision) library to Free Pascal.
Goal: bit-exact, correctly-rounded results matching the C reference for all 2^32 float inputs.

---

## Phase 0 — Infrastructure (prerequisite for everything)

- [ ] **0.1** Define `TUInt128` as a pure record in `CoreMathTypes.pas`:
  ```pascal
  type TUInt128 = record
    lo, hi: UInt64;
  end;
  ```

- [ ] **0.2** Implement `MulWide(a, b: UInt64): TUInt128` in `CoreMathTypes.pas`.
  - Primary path: x86-64 inline assembly using the `MUL` instruction (`rdx:rax = rax * src`).
  - Portable fallback: four 32-bit partial products for non-x86-64 targets (ARM, etc.).

- [ ] **0.3** Overload `+` for `TUInt128 + UInt64 → TUInt128` (implements the `p1 += p0>>64` pattern as `p1 := p1 + p0.hi`).
  This is a genuine addition with carry propagation — implement as:
  ```pascal
  operator+(const a: TUInt128; b: UInt64): TUInt128; inline;
  begin
    Result.lo := a.lo + b;
    Result.hi := a.hi + UInt64(Result.lo < b);  // carry
  end;
  ```
  Mark `inline` — it is a hot-path operation called inside `MulWide`. Also try marking `MulWide` itself as `inline`; trust the compiler to handle it.

- [ ] **0.4** Define type-punning records in `CoreMathTypes.pas`:
  ```pascal
  type Tb32u32 = record case Boolean of
    False: (f: Single);  True: (u: LongWord); end;
  type Tb64u64 = record case Boolean of
    False: (f: Double);  True: (u: UInt64);   end;
  ```

- [ ] **0.5** Implement builtin equivalents in `CoreMathBuiltins.pas`:
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

- [ ] **0.6** Write a hex-float conversion utility (script or small tool) to convert C99
  hex float literals (e.g. `0x1.62e42fefa39efp-1`) to Pascal `Double` / `UInt64`
  constants. All lookup tables must be converted using this tool, not by hand.

- [ ] **0.7** Set up a test harness that compiles and runs both the C reference and the
  Pascal implementation, then compares results bit-for-bit for all 2^32 `Single`
  inputs (exhaustive). For bivariate functions, agree on a sampling strategy.

---

## Phase 1 — Simple univariate (no u128, no FMA, ≤ 130 lines)

Port in this order. Each function lives in its own `.pas` unit named after the function.

- [ ] **1.01** `rsqrt`   — 89 lines
- [ ] **1.02** `tanh`    — 89 lines
- [ ] **1.03** `atanpi`  — 106 lines
- [ ] **1.04** `cospi`   — 114 lines
- [ ] **1.05** `acos`    — 115 lines
- [ ] **1.06** `cbrt`    — 117 lines
- [ ] **1.07** `sinpi`   — 117 lines
- [ ] **1.08** `atan`    — 118 lines
- [ ] **1.09** `asin`    — 120 lines
- [ ] **1.10** `acospi`  — 126 lines  *(uses FMA — verify `Math.FMA` correctness first)*
- [ ] **1.11** `log2`    — 126 lines
- [ ] **1.12** `asinpi`  — 128 lines  *(uses FMA)*
- [ ] **1.13** `tanpi`   — 130 lines
- [ ] **1.14** `cosh`    — 132 lines

---

## Phase 2 — Medium univariate

- [ ] **2.01** `log`       — 133 lines
- [ ] **2.02** `exp2`      — 138 lines
- [ ] **2.03** `log1p`     — 140 lines
- [ ] **2.04** `exp2m1`    — 146 lines
- [ ] **2.05** `expm1`     — 146 lines
- [ ] **2.06** `exp10`     — 150 lines
- [ ] **2.07** `log10`     — 150 lines
- [ ] **2.08** `erfc`      — 151 lines
- [ ] **2.09** `log2p1`    — 151 lines
- [ ] **2.10** `erf`       — 152 lines
- [ ] **2.11** `sinh`      — 156 lines
- [ ] **2.12** `exp`       — 115 lines *(listed here due to two-pass rounding logic)*
- [ ] **2.13** `atanh`     — 158 lines
- [ ] **2.14** `exp10m1`   — 158 lines
- [ ] **2.15** `log10p1`   — 162 lines
- [ ] **2.16** `asinh`     — 168 lines
- [ ] **2.17** `acosh`     — 173 lines

---

## Phase 3 — Longer special univariate

- [ ] **3.01** `tgamma`  — 205 lines
- [ ] **3.02** `lgamma`  — 259 lines

---

## Phase 4 — Bivariate and compound (FMA-heavy)

- [ ] **4.01** `hypot`    — 282 lines  *(uses FMA, fenv rounding modes)*
- [ ] **4.02** `atan2`    — 231 lines  *(uses FMA)*
- [ ] **4.03** `atan2pi`  — 190 lines  *(uses FMA)*
- [ ] **4.04** `compound` — 611 lines  *(uses FMA, fenv)*
- [ ] **4.05** `pow`      — 325 lines  *(uses FMA, fenv)*

---

## Phase 5 — u128 functions (hardest)

Depends on Phase 0.2–0.3 being fully correct. Validate `MulWide` independently before
starting this phase.

- [ ] **5.01** `sin`    — 222 lines
- [ ] **5.02** `cos`    — 206 lines
- [ ] **5.03** `sincos` — 245 lines
- [ ] **5.04** `tan`    — 199 lines

---

## Per-function porting checklist

Apply this checklist to every function before marking it done:

- [ ] Hex float literals converted via the conversion utility (0.6), not by hand
- [ ] Lookup tables moved to unit-level `const` (no `static` locals)
- [ ] All type-punning uses `Tb32u32` / `Tb64u64` records (no unsafe casts)
- [ ] `__builtin_expect` wrappers removed entirely
- [ ] `__attribute__((noinline))` replaced with `[noinline]`
- [ ] `CORE_MATH_SUPPORT_ERRNO` blocks omitted (out of scope for Pascal port)
- [ ] Exhaustive test passes (bit-exact match against C reference for all inputs)
- [ ] Function named `cr_<name>f` in C becomes `cr_<name>f` in Pascal for traceability

---

## Architectural notes and known pitfalls

1. **`rbig()` must be extracted into a shared unit.** The large-argument range-reduction
   helper `rbig()` is byte-for-byte identical in `sin`, `cos`, `tan`, and `sincos`. It must
   live in one shared unit (e.g. `CoreMathTrig.pas`) and be called by all four, not duplicated.
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
   wrapped in `CoreMathBuiltins.pas` before Phase 4 starts — it is not safe to defer until
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
