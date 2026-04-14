// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
//
// FixedTest32: call every single-precision function with a fixed input and
// print results from three implementations side-by-side:
//   FPC  - FPC/Math built-in (or closest mathematical equivalent; N/A = no
//          direct built-in, column shows NaN)
//   C    - CORE-MATH reference C library (cr_*)
//   pcr  - pas-core-math Pascal port (pcr_*)
//
// This is a quick sanity check, not an exhaustive correctness test.
// Use TestHarness32 for bit-exact validation across all inputs.
//
// pcr_* functions are always called through function pointers (same pattern
// as TestHarness32/Benchmark32) to prevent FPC 3.2 internal error 200312205
// caused by register-allocator exhaustion from deeply-nested inline chains.
{$I ../pascoremath.inc}
program FixedTest32;

uses
  pascoremathtypes, pascoremath32, ccoremath32, SysUtils, StrUtils, Math;

const
  TEST_X:  Single = 0.5;  // fixed input for all univariate and bivariate tests
  ACOSH_X: Single = 1.5;  // acoshf domain requires x >= 1; tested separately

type
  TUniFuncC  = function(x: Single): Single; cdecl;
  TUniFuncP  = function(x: Single): Single;
  TBivarFuncC = function(x, y: Single): Single; cdecl;
  TBivarFuncP = function(x, y: Single): Single;

var
  s_fpc, c_fpc, s_c, c_c, s_pcr, c_pcr: Single;

// Format a Single as a signed float string: (+1.2345678e+00) or (-1.23...)
// FPC 3.2 Format does not support the '+' flag modifier, so we handle it
// explicitly; negative sign is included naturally by %.7e for negative values.
function FmtSingle(v: Single): string;
var
  d: Double;
  s: string;
begin
  d := Double(v);
  s := Format('%.7e', [d]);
  // NaN/Inf: Format already includes a sign for Inf; NaN has no sign.
  // Only prepend '+' for finite non-negative values.
  if IsNaN(d) or IsInfinite(d) then
    Result := '(' + s + ')'
  else if d >= 0.0 then
    Result := '(+' + s + ')'
  else
    Result := '(' + s + ')';
end;

// Print one row: function label + three (hex, float) result pairs + MATCH/ERROR.
// MATCH when C and pcr are bit-identical (NaN == any NaN); ERROR otherwise.
procedure PrintRow(const Name: string; FPCRes, CRes, PcrRes: Single);
var
  f, c, p: Tb32u32;
  verdict: string;
begin
  f.f := FPCRes; c.f := CRes; p.f := PcrRes;
  // Bit-exact comparison: both NaN => match; otherwise identical bits required.
  if ((c.u and $7F800000) = $7F800000) and ((c.u and $007FFFFF) <> 0) then
    verdict := IfThen(((p.u and $7F800000) = $7F800000) and ((p.u and $007FFFFF) <> 0), 'MATCH', 'ERROR')
  else
    verdict := IfThen(c.u = p.u, 'MATCH', 'ERROR');
  WriteLn(
    Format('%-22s', [Name]) +
    '  FPC=$' + IntToHex(f.u, 8) + FmtSingle(FPCRes) +
    '  C=$'   + IntToHex(c.u, 8) + FmtSingle(CRes)   +
    '  pcr=$' + IntToHex(p.u, 8) + FmtSingle(PcrRes) +
    '  ' + verdict
  );
end;

// Evaluate and print one univariate row; label is built from FuncName + x.
// pcr_* called via pointer to avoid FPC inline chain issues.
procedure PrintUni(const FuncName: string; FPCResult: Single;
                   pfC: TUniFuncC; pfP: TUniFuncP; x: Single);
begin
  PrintRow(Format('%s(%g)', [FuncName, x]), FPCResult, pfC(x), pfP(x));
end;

// Two-column row: C and pcr only (no FPC column), used for functions with no
// FPC equivalent. MATCH/ERROR compares C vs pcr bit-exactly.
// The FPC column (31 chars = '  FPC=$' + 8 hex digits + 16-char FmtSingle) is
// replaced by spaces so C and pcr columns stay aligned with 3-column rows.
procedure PrintRowCP(const Name: string; CRes, PcrRes: Single);
const
  FPC_COL_WIDTH = 31;  // '  FPC=$'(7) + hex(8) + FmtSingle normal float(16)
var
  c, p: Tb32u32;
  verdict: string;
begin
  c.f := CRes; p.f := PcrRes;
  if ((c.u and $7F800000) = $7F800000) and ((c.u and $007FFFFF) <> 0) then
    verdict := IfThen(((p.u and $7F800000) = $7F800000) and ((p.u and $007FFFFF) <> 0), 'MATCH', 'ERROR')
  else
    verdict := IfThen(c.u = p.u, 'MATCH', 'ERROR');
  WriteLn(
    Format('%-22s', [Name]) +
    StringOfChar(' ', FPC_COL_WIDTH) +
    '  C=$'   + IntToHex(c.u, 8) + FmtSingle(CRes)   +
    '  pcr=$' + IntToHex(p.u, 8) + FmtSingle(PcrRes) +
    '  ' + verdict
  );
end;

// Evaluate and print a C+pcr-only univariate row (no FPC equivalent).
procedure PrintUniCP(const FuncName: string; pfC: TUniFuncC; pfP: TUniFuncP; x: Single);
begin
  PrintRowCP(Format('%s(%g)', [FuncName, x]), pfC(x), pfP(x));
end;

// Evaluate and print one bivariate row; label is built from FuncName + x + y.
procedure PrintBivar(const FuncName: string; FPCResult: Single;
                     pfC: TBivarFuncC; pfP: TBivarFuncP; x, y: Single);
begin
  PrintRow(Format('%s(%g,%g)', [FuncName, x, y]), FPCResult, pfC(x, y), pfP(x, y));
end;

begin
  {$IFDEF AVX2}
  WriteLn('Compiled with AVX2.');
  {$ENDIF}

  SetExceptionMask([exInvalidOp, exDenormalized, exZeroDivide,
                    exOverflow, exUnderflow, exPrecision]);

  WriteLn('=== FixedTest32: FPC vs C (cr_*) vs pas-core-math (pcr_*) ===');
  WriteLn('Functions with no FPC equivalent show only C and pcr columns.');
  WriteLn;

  // -----------------------------------------------------------------------
  // Univariate functions
  // -----------------------------------------------------------------------
  WriteLn(Format('--- Univariate (x = %g, except acoshf which uses x = %g) ---',
    [TEST_X, ACOSH_X]));

  PrintUni('acosf',    ArcCos(TEST_X),              @cr_acosf,    @pcr_acosf,    TEST_X);
  // acoshf domain requires x >= 1; use ACOSH_X
  PrintUni('acoshf',   ArcCosh(ACOSH_X),            @cr_acoshf,   @pcr_acoshf,   ACOSH_X);
  // acospif(x) = acos(x)/pi; no dedicated FPC built-in
  PrintUni('acospif',  ArcCos(TEST_X)/Pi,           @cr_acospif,  @pcr_acospif,  TEST_X);
  PrintUni('asinf',    ArcSin(TEST_X),              @cr_asinf,    @pcr_asinf,    TEST_X);
  PrintUni('asinhf',   ArcSinh(TEST_X),             @cr_asinhf,   @pcr_asinhf,   TEST_X);
  // asinpif(x) = asin(x)/pi; no dedicated FPC built-in
  PrintUni('asinpif',  ArcSin(TEST_X)/Pi,           @cr_asinpif,  @pcr_asinpif,  TEST_X);
  PrintUni('atanf',    ArcTan(TEST_X),              @cr_atanf,    @pcr_atanf,    TEST_X);
  PrintUni('atanhf',   ArcTanh(TEST_X),             @cr_atanhf,   @pcr_atanhf,   TEST_X);
  // atanpif(x) = atan(x)/pi; no dedicated FPC built-in
  PrintUni('atanpif',  ArcTan(TEST_X)/Pi,           @cr_atanpif,  @pcr_atanpif,  TEST_X);
  // cbrtf(x) = x^(1/3); no dedicated FPC built-in (Power used here, valid for x>0)
  PrintUni('cbrtf',    Power(TEST_X, 1/3),          @cr_cbrtf,    @pcr_cbrtf,    TEST_X);
  PrintUni('cosf',     Cos(TEST_X),                 @cr_cosf,     @pcr_cosf,     TEST_X);
  PrintUni('coshf',    Cosh(TEST_X),                @cr_coshf,    @pcr_coshf,    TEST_X);
  // cospif(x) = cos(x*pi); no dedicated FPC built-in
  // Note: FPC Cos(TEST_X*Pi) differs from CORE-MATH cospif(TEST_X) at the bit level
  PrintUni('cospif',   Cos(TEST_X*Pi),              @cr_cospif,   @pcr_cospif,   TEST_X);
  // erff / erfcf: no FPC built-in
  PrintUniCP('erff',                                @cr_erff,     @pcr_erff,     TEST_X);
  PrintUniCP('erfcf',                               @cr_erfcf,    @pcr_erfcf,    TEST_X);
  PrintUni('expf',     Exp(TEST_X),                 @cr_expf,     @pcr_expf,     TEST_X);
  // exp10f(x) = 10^x; no dedicated FPC built-in
  PrintUni('exp10f',   Power(10, TEST_X),           @cr_exp10f,   @pcr_exp10f,   TEST_X);
  // exp10m1f(x) = 10^x - 1; no dedicated FPC built-in
  PrintUni('exp10m1f', Power(10, TEST_X)-1,         @cr_exp10m1f, @pcr_exp10m1f, TEST_X);
  // exp2f(x) = 2^x; no dedicated FPC built-in
  PrintUni('exp2f',    Power(2, TEST_X),            @cr_exp2f,    @pcr_exp2f,    TEST_X);
  // exp2m1f(x) = 2^x - 1; no dedicated FPC built-in
  PrintUni('exp2m1f',  Power(2, TEST_X)-1,          @cr_exp2m1f,  @pcr_exp2m1f,  TEST_X);
  // expm1f(x) = e^x - 1; no dedicated FPC built-in
  PrintUni('expm1f',   Exp(TEST_X)-1,               @cr_expm1f,   @pcr_expm1f,   TEST_X);
  // lgammaf: no FPC built-in
  PrintUniCP('lgammaf',                             @cr_lgammaf,  @pcr_lgammaf,  TEST_X);
  PrintUni('logf',     Ln(TEST_X),                  @cr_logf,     @pcr_logf,     TEST_X);
  PrintUni('log10f',   Log10(TEST_X),               @cr_log10f,   @pcr_log10f,   TEST_X);
  // log10p1f(x) = log10(1+x); no dedicated FPC built-in
  PrintUni('log10p1f', Log10(1+TEST_X),             @cr_log10p1f, @pcr_log10p1f, TEST_X);
  // log1pf(x) = ln(1+x); no dedicated FPC built-in
  PrintUni('log1pf',   Ln(1+TEST_X),                @cr_log1pf,   @pcr_log1pf,   TEST_X);
  PrintUni('log2f',    Log2(TEST_X),                @cr_log2f,    @pcr_log2f,    TEST_X);
  // log2p1f(x) = log2(1+x); no dedicated FPC built-in
  PrintUni('log2p1f',  Log2(1+TEST_X),              @cr_log2p1f,  @pcr_log2p1f,  TEST_X);
  // rsqrtf(x) = 1/sqrt(x); no dedicated FPC built-in
  PrintUni('rsqrtf',   1/Sqrt(TEST_X),              @cr_rsqrtf,   @pcr_rsqrtf,   TEST_X);
  PrintUni('sinf',     Sin(TEST_X),                 @cr_sinf,     @pcr_sinf,     TEST_X);
  PrintUni('sinhf',    Sinh(TEST_X),                @cr_sinhf,    @pcr_sinhf,    TEST_X);
  // sinpif(x) = sin(x*pi); no dedicated FPC built-in
  PrintUni('sinpif',   Sin(TEST_X*Pi),              @cr_sinpif,   @pcr_sinpif,   TEST_X);
  PrintUni('tanf',     Tan(TEST_X),                 @cr_tanf,     @pcr_tanf,     TEST_X);
  PrintUni('tanhf',    Tanh(TEST_X),                @cr_tanhf,    @pcr_tanhf,    TEST_X);
  // tanpif(TEST_X=0.5) = tan(pi/2) = +Inf; FPC Tan(TEST_X*Pi) returns a large
  // finite value because floating-point pi/2 is not exact, while CORE-MATH
  // correctly returns +Inf.
  PrintUni('tanpif',   Tan(TEST_X*Pi),              @cr_tanpif,   @pcr_tanpif,   TEST_X);
  // tgammaf: no FPC built-in
  PrintUniCP('tgammaf',                             @cr_tgammaf,  @pcr_tgammaf,  TEST_X);

  WriteLn;

  // -----------------------------------------------------------------------
  // Bivariate functions
  // -----------------------------------------------------------------------
  WriteLn(Format('--- Bivariate (x = %g, y = %g) ---', [TEST_X, TEST_X]));

  PrintBivar('atan2f',    ArcTan2(TEST_X, TEST_X),        @cr_atan2f,    @pcr_atan2f,    TEST_X, TEST_X);
  // atan2pif(y,x) = atan2(y,x)/pi; no dedicated FPC built-in
  PrintBivar('atan2pif',  ArcTan2(TEST_X, TEST_X)/Pi,     @cr_atan2pif,  @pcr_atan2pif,  TEST_X, TEST_X);
  PrintBivar('hypotf',    Hypot(TEST_X, TEST_X),          @cr_hypotf,    @pcr_hypotf,    TEST_X, TEST_X);
  PrintBivar('powf',      Power(TEST_X, TEST_X),          @cr_powf,      @pcr_powf,      TEST_X, TEST_X);
  // compoundf(x,n) = (1+x)^n; no dedicated FPC built-in
  PrintBivar('compoundf', Power(1+TEST_X, TEST_X),        @cr_compoundf, @pcr_compoundf, TEST_X, TEST_X);

  WriteLn;

  // -----------------------------------------------------------------------
  // sincosf — two output values
  // -----------------------------------------------------------------------
  WriteLn(Format('--- sincosf (x = %g) ---', [TEST_X]));
  SinCos(TEST_X, s_fpc, c_fpc);
  cr_sincosf(TEST_X, @s_c, @c_c);
  pcr_sincosf(TEST_X, s_pcr, c_pcr);
  PrintRow(Format('sincosf sin(%g)', [TEST_X]), s_fpc, s_c, s_pcr);
  PrintRow(Format('sincosf cos(%g)', [TEST_X]), c_fpc, c_c, c_pcr);
end.
