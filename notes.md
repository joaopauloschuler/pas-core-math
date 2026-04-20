# lgammaf

## Baseline
- Pascal: 39.4 Mops/s vs C: 153.8 Mops/s (26% of C speed)

## Profiling findings (assembly analysis)
FPC generates x87 code in several key places, causing expensive SSE↔x87 pipeline
transitions:

1. **`lgamma_as_ln` and `lgamma_as_sinpi` not inlined**: Function call overhead per
   iteration; FPC can only inline functions if `inline` is in the declaration.

2. **Float literal constants stored as 80-bit extended by FPC**: When a Pascal floating-
   point literal (e.g. `0.6931471805599453`, `0.4189385332046727`, `3.373...`) appears
   in an expression without an explicit type cast, FPC stores it as `extended` (80-bit).
   Any arithmetic or comparison involving an extended constant forces the use of x87
   instructions (`fldl`, `fldt`, `fcomip`, `faddp`, etc.) even with `-CfAVX2`.
   Fix: wrap every such literal in `Double(...)` or `Single(...)`.

3. **Float comparisons using x87**: Comparing `Single` variables against untyped
   floating-point literals (e.g. `ax_lg > 3.373...`) triggers x87 comparison
   (`flds`/`fldt`/`fcomip`). Fix: use `t_lg.u > hex_bits` integer comparison
   (for `ax_lg` whose bits are already in `t_lg.u`) or `Single(constant)` cast.

4. **`Floor(Double(x))` → x87 function call**: FPC's `Math.Floor` takes an `Extended`
   argument. The conversion to extended and the function call both use x87. This path
   is taken for ~59% of the benchmark inputs (all `|x| < 2^23 = 8388608`).
   Fix: replace with SSE-based floor using `Trunc` (→ `vcvttss2si`) + integer adjust.

## Round 1 fixes (50.3 Mops/s, +28%)
- Added `inline` to `lgamma_as_sinpi` and `lgamma_as_ln`
- Fixed `lgamma_as_ln`: `0.6931471805599453` → `Double(0.6931471805599453)`
- Replaced `Floor(Double(x))` with SSE Trunc-based floor (no x87)
- `ax_lg > 3.373...` → `t_lg.u > $4057E0D0` (integer compare, avoids x87)
- `ax_lg > 10.666...` → `t_lg.u > $412AAAAB` (integer compare, avoids x87)
- `x >= 4.085e36` → `x >= Single(4.085e36)` (SSE comiss)
- `+ 0.4189385332046727` → `+ Double(0.4189385332046727)` (SSE addsd)
- `1.1447298858494002 - f_lg - lz_lg` → `Double(1.1447298858494002) - ...`
- Near-gamma-zero h_lg corrections all wrapped in `Double(...)` casts
- Result: all x87 eliminated from `pcr_lgammaf` (verified in assembly)
- Tests: PASS (4294967296 values, 0 mismatches)

## Still to investigate
- Large path: `lgamma_as_ln` now inlined but instruction count remains high
- The `divsd` in the small/medium rational paths may be a bottleneck
- Can the hot path (large positive x) be further reduced in latency?

