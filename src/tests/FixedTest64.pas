// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
//
// FixedTest64: call every double-precision function with a fixed input and
// print results from three implementations side-by-side:
//   FPC  - FPC/Math built-in (or closest mathematical equivalent; N/A = no
//          direct built-in)
//   C    - CORE-MATH reference C library (cr_*)
//   pcr  - pas-core-math Pascal port (pcr_*)
//
// Quick sanity check; for bit-exact validation use TestHarness64.
{$I ../pascoremath.inc}
program FixedTest64;

uses
  pascoremathtypes, pascoremath64, ccoremath64, SysUtils, StrUtils, Math;

const
  TEST_X:  Double = 0.5;
  ACOSH_X: Double = 1.5;  // acosh domain requires x >= 1

type
  TUniFuncC   = function(x: Double): Double; cdecl;
  TUniFuncP   = function(x: Double): Double;
  TBivarFuncC = function(x, y: Double): Double; cdecl;
  TBivarFuncP = function(x, y: Double): Double;

var
  s_fpc, c_fpc, s_c, c_c, s_pcr, c_pcr: Double;

function FmtDouble(v: Double): string;
var
  s: string;
begin
  s := Format('%.7e', [v]);
  if IsNaN(v) or IsInfinite(v) then
    Result := '(' + s + ')'
  else if v >= 0.0 then
    Result := '(+' + s + ')'
  else
    Result := '(' + s + ')';
end;

procedure PrintRow(const Name: string; FPCRes, CRes, PcrRes: Double);
var
  f, c, p: Tb64u64;
  verdict: string;
begin
  f.f := FPCRes; c.f := CRes; p.f := PcrRes;
  if ((c.u and $7FF0000000000000) = $7FF0000000000000) and ((c.u and $000FFFFFFFFFFFFF) <> 0) then
    verdict := IfThen(((p.u and $7FF0000000000000) = $7FF0000000000000) and ((p.u and $000FFFFFFFFFFFFF) <> 0), 'MATCH', 'ERROR')
  else
    verdict := IfThen(c.u = p.u, 'MATCH', 'ERROR');
  WriteLn(
    Format('%-22s', [Name]) +
    '  FPC=$' + IntToHex(f.u, 16) + FmtDouble(FPCRes) +
    '  C=$'   + IntToHex(c.u, 16) + FmtDouble(CRes)   +
    '  pcr=$' + IntToHex(p.u, 16) + FmtDouble(PcrRes) +
    '  ' + verdict
  );
end;

procedure PrintUni(const FuncName: string; FPCResult: Double;
                   pfC: TUniFuncC; pfP: TUniFuncP; x: Double);
begin
  PrintRow(Format('%s(%g)', [FuncName, x]), FPCResult, pfC(x), pfP(x));
end;

procedure PrintRowCP(const Name: string; CRes, PcrRes: Double);
const
  FPC_COL_WIDTH = 39;  // '  FPC=$'(7) + hex(16) + FmtDouble normal float(16)
var
  c, p: Tb64u64;
  verdict: string;
begin
  c.f := CRes; p.f := PcrRes;
  if ((c.u and $7FF0000000000000) = $7FF0000000000000) and ((c.u and $000FFFFFFFFFFFFF) <> 0) then
    verdict := IfThen(((p.u and $7FF0000000000000) = $7FF0000000000000) and ((p.u and $000FFFFFFFFFFFFF) <> 0), 'MATCH', 'ERROR')
  else
    verdict := IfThen(c.u = p.u, 'MATCH', 'ERROR');
  WriteLn(
    Format('%-22s', [Name]) +
    StringOfChar(' ', FPC_COL_WIDTH) +
    '  C=$'   + IntToHex(c.u, 16) + FmtDouble(CRes)   +
    '  pcr=$' + IntToHex(p.u, 16) + FmtDouble(PcrRes) +
    '  ' + verdict
  );
end;

procedure PrintUniCP(const FuncName: string; pfC: TUniFuncC; pfP: TUniFuncP; x: Double);
begin
  PrintRowCP(Format('%s(%g)', [FuncName, x]), pfC(x), pfP(x));
end;

procedure PrintBivar(const FuncName: string; FPCResult: Double;
                     pfC: TBivarFuncC; pfP: TBivarFuncP; x, y: Double);
begin
  PrintRow(Format('%s(%g,%g)', [FuncName, x, y]), FPCResult, pfC(x, y), pfP(x, y));
end;

// cdecl wrappers for bivariate C functions
function wrap_atan2_c(y, x: Double): Double; cdecl;   begin Result := cr_atan2(y, x);   end;
function wrap_atan2pi_c(y, x: Double): Double; cdecl; begin Result := cr_atan2pi(y, x); end;
function wrap_hypot_c(x, y: Double): Double; cdecl;   begin Result := cr_hypot(x, y);   end;
function wrap_pow_c(x, y: Double): Double; cdecl;     begin Result := cr_pow(x, y);     end;
// Pascal bivariate wrappers
function wrap_atan2_p(y, x: Double): Double;   begin Result := pcr_atan2(y, x);   end;
function wrap_atan2pi_p(y, x: Double): Double; begin Result := pcr_atan2pi(y, x); end;
function wrap_hypot_p(x, y: Double): Double;   begin Result := pcr_hypot(x, y);   end;
function wrap_pow_p(x, y: Double): Double;     begin Result := pcr_pow(x, y);     end;

begin
  {$IFDEF AVX2}
  WriteLn('Compiled with AVX2.');
  {$ENDIF}
  SetExceptionMask([exInvalidOp, exDenormalized, exZeroDivide,
                    exOverflow, exUnderflow, exPrecision]);

  WriteLn('=== FixedTest64: FPC vs C (cr_*) vs pas-core-math (pcr_*) ===');
  WriteLn('Functions with no FPC equivalent show only C and pcr columns.');
  WriteLn;
  WriteLn('--- Univariate (x = 0.5, except acosh which uses x = 1.5) ---');

  PrintUni('acos(0.5)',   ArcCos(TEST_X),  @cr_acos,   @pcr_acos,   TEST_X);
  PrintUni('acosh(1.5)',  ArcCosh(ACOSH_X),@cr_acosh,  @pcr_acosh,  ACOSH_X);
  PrintUniCP('acospi(0.5)',              @cr_acospi,  @pcr_acospi,  TEST_X);
  PrintUni('asin(0.5)',   ArcSin(TEST_X),  @cr_asin,   @pcr_asin,   TEST_X);
  PrintUni('asinh(0.5)',  ArcSinh(TEST_X), @cr_asinh,  @pcr_asinh,  TEST_X);
  PrintUniCP('asinpi(0.5)',              @cr_asinpi,  @pcr_asinpi,  TEST_X);
  PrintUni('atan(0.5)',   ArcTan(TEST_X),  @cr_atan,   @pcr_atan,   TEST_X);
  PrintUni('atanh(0.5)',  ArcTanh(TEST_X), @cr_atanh,  @pcr_atanh,  TEST_X);
  PrintUniCP('atanpi(0.5)',              @cr_atanpi,  @pcr_atanpi,  TEST_X);
  PrintUniCP('cbrt(0.5)',               @cr_cbrt,    @pcr_cbrt,    TEST_X);
  PrintUni('cos(0.5)',    Cos(TEST_X),     @cr_cos,    @pcr_cos,    TEST_X);
  PrintUni('cosh(0.5)',   Cosh(TEST_X),    @cr_cosh,   @pcr_cosh,   TEST_X);
  PrintUniCP('cospi(0.5)',              @cr_cospi,   @pcr_cospi,   TEST_X);
  PrintUniCP('erf(0.5)',                @cr_erf,     @pcr_erf,     TEST_X);
  PrintUniCP('erfc(0.5)',               @cr_erfc,    @pcr_erfc,    TEST_X);
  PrintUni('exp(0.5)',    Exp(TEST_X),     @cr_exp,    @pcr_exp,    TEST_X);
  PrintUniCP('exp10(0.5)',              @cr_exp10,   @pcr_exp10,   TEST_X);
  PrintUniCP('exp10m1(0.5)',            @cr_exp10m1, @pcr_exp10m1, TEST_X);
  PrintUniCP('exp2(0.5)',               @cr_exp2,    @pcr_exp2,    TEST_X);
  PrintUniCP('exp2m1(0.5)',             @cr_exp2m1,  @pcr_exp2m1,  TEST_X);
  PrintUniCP('expm1(0.5)',              @cr_expm1,   @pcr_expm1,   TEST_X);
  PrintUniCP('lgamma(0.5)',             @cr_lgamma,  @pcr_lgamma,  TEST_X);
  PrintUni('log(0.5)',    Ln(TEST_X),      @cr_log,    @pcr_log,    TEST_X);
  PrintUni('log10(0.5)',  Log10(TEST_X),   @cr_log10,  @pcr_log10,  TEST_X);
  PrintUniCP('log10p1(0.5)',            @cr_log10p1, @pcr_log10p1, TEST_X);
  PrintUniCP('log1p(0.5)',              @cr_log1p,   @pcr_log1p,   TEST_X);
  PrintUni('log2(0.5)',   Log2(TEST_X),    @cr_log2,   @pcr_log2,   TEST_X);
  PrintUniCP('log2p1(0.5)',             @cr_log2p1,  @pcr_log2p1,  TEST_X);
  PrintUniCP('rsqrt(0.5)',              @cr_rsqrt,   @pcr_rsqrt,   TEST_X);
  PrintUni('sin(0.5)',    Sin(TEST_X),     @cr_sin,    @pcr_sin,    TEST_X);
  PrintUni('sinh(0.5)',   Sinh(TEST_X),    @cr_sinh,   @pcr_sinh,   TEST_X);
  PrintUniCP('sinpi(0.5)',              @cr_sinpi,   @pcr_sinpi,   TEST_X);
  PrintUni('tan(0.5)',    Tan(TEST_X),     @cr_tan,    @pcr_tan,    TEST_X);
  PrintUni('tanh(0.5)',   Tanh(TEST_X),    @cr_tanh,   @pcr_tanh,   TEST_X);
  PrintUniCP('tanpi(0.5)',              @cr_tanpi,   @pcr_tanpi,   TEST_X);
  PrintUniCP('tgamma(0.5)',             @cr_tgamma,  @pcr_tgamma,  TEST_X);

  WriteLn;
  WriteLn('--- Bivariate (x = 0.5, y = 0.5) ---');

  PrintBivar('atan2(0.5,0.5)',   ArcTan2(TEST_X,TEST_X), @wrap_atan2_c,   @wrap_atan2_p,   TEST_X, TEST_X);
  PrintRowCP('atan2pi(0.5,0.5)', cr_atan2pi(TEST_X, TEST_X), pcr_atan2pi(TEST_X, TEST_X));
  PrintBivar('hypot(0.5,0.5)',   Hypot(TEST_X,TEST_X),   @wrap_hypot_c,   @wrap_hypot_p,   TEST_X, TEST_X);
  PrintBivar('pow(0.5,0.5)',     Power(TEST_X,TEST_X),   @wrap_pow_c,     @wrap_pow_p,     TEST_X, TEST_X);

  WriteLn;
  WriteLn('--- sincos(0.5) ---');
  cr_sincos(TEST_X, @s_c, @c_c);
  pcr_sincos(TEST_X, s_pcr, c_pcr);
  SinCos(TEST_X, s_fpc, c_fpc);
  PrintRow('sincos sin(0.5)', s_fpc, s_c, s_pcr);
  PrintRow('sincos cos(0.5)', c_fpc, c_c, c_pcr);
end.
