// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
{$I ../pascoremath.inc}
program TestHarness64;

uses
  pascoremathtypes, pascoremath64, ccoremath64, SysUtils, Math;

type
  TUniFuncC   = function(x: Double): Double; cdecl;
  TUniFuncP   = function(x: Double): Double;
  TBivarFuncC = function(x, y: Double): Double; cdecl;
  TBivarFuncP = function(x, y: Double): Double;

const
  // Number of random samples per univariate function
  SAMPLES_UNI = 100000000;  // 10^8; increase to 10^9 for final validation
  SAMPLES_BIV = 10000000;   // 10^7

var
  TotalPass, TotalFail: Int32;
  DiagMax: Int32;
  StartTick: UInt64;
  Filter: string = '';
  SamplePct: Double = 100.0;  // percentage of SAMPLES_* to actually run; --pct <n>

// Fast xorshift64 PRNG (period 2^64 - 1)
var RngState: UInt64 = $6E0BEEF1CAFE7A5E;  // arbitrary seed; override with --seed

function Xorshift64: UInt64; inline;
begin
  RngState := RngState xor (RngState shl 13);
  RngState := RngState xor (RngState shr 7);
  RngState := RngState xor (RngState shl 17);
  Result := RngState;
end;

procedure ReportResult(const FuncName: string; Tested, Mismatches: Int64; MaxUlp: Double);
begin
  if Mismatches = 0 then
  begin
    WriteLn(Format('%-16s  tested=%10d  mismatches=0  PASS', [FuncName, Tested]));
    Inc(TotalPass);
  end
  else
  begin
    WriteLn(Format('%-16s  tested=%10d  mismatches=%d  FAIL max_ulp=%.6g', [FuncName, Tested, Mismatches, MaxUlp]));
    Inc(TotalFail);
  end;
end;

// Bit-exact comparison: both NaN => match; otherwise bits must be identical.
function BitsMatch(a, b: Double): Boolean; inline;
var
  ua, ub: Tb64u64;
begin
  ua.f := a;
  ub.f := b;
  // If a is NaN, b must also be NaN (any NaN)
  if ((ua.u and $7FF0000000000000) = $7FF0000000000000) and ((ua.u and $000FFFFFFFFFFFFF) <> 0) then
  begin
    Result := ((ub.u and $7FF0000000000000) = $7FF0000000000000) and ((ub.u and $000FFFFFFFFFFFFF) <> 0);
    Exit;
  end;
  Result := ua.u = ub.u;
end;

// ---- special values to always include ----
const
  SPECIAL_COUNT = 16;
// Constants for special bit patterns (avoids large-literal warnings in var init)
  SP_POS_ZERO  = UInt64($0000000000000000);
  SP_ONE       = UInt64($3FF0000000000000);
  SP_HALF      = UInt64($3FE0000000000000);
  SP_TWO       = UInt64($4000000000000000);
  SP_INF       = UInt64($7FF0000000000000);
  SP_QNAN_POS  = UInt64($7FF8000000000000);
  SP_SNAN      = UInt64($7FF0000000000001);
  SP_MIN_SUB   = UInt64($0000000000000001);
  SP_MAX_SUB   = UInt64($000FFFFFFFFFFFFF);
  SP_MIN_NORM  = UInt64($0010000000000000);
  SP_MAX_NORM  = UInt64($7FEFFFFFFFFFFFFF);
  SP_PI        = UInt64($3FF921FB54442D18);

var
  SPECIAL_BITS: array[0..SPECIAL_COUNT-1] of UInt64;

// ---- univariate sampled test ----

procedure TestUni(const name: string; pfC: TUniFuncC; pfP: TUniFuncP);
var
  i: Int32;
  bits: UInt64;
  v, rc, rp: Tb64u64;
  cr, pr: Double;
  mismatches, tested: Int64;
  diagShown: Int32;
  maxUlp: Double;
begin
  if (Filter <> '') and (LowerCase(name) <> Filter) then Exit;
  mismatches := 0;
  tested := 0;
  diagShown := 0;
  maxUlp := 0;

  // Special values first
  for i := 0 to SPECIAL_COUNT - 1 do
  begin
    v.u := SPECIAL_BITS[i];
    cr := pfC(v.f);
    pr := pfP(v.f);
    if not BitsMatch(cr, pr) then
    begin
      Inc(mismatches);
      if diagShown < DiagMax then
      begin
        rc.f := cr; rp.f := pr;
        WriteLn(Format('  [%s] input=$%16.16x(%.17g)  C=$%16.16x  P=$%16.16x',
                       [name, SPECIAL_BITS[i], v.f, rc.u, rp.u]));
        Inc(diagShown);
      end;
    end;
    Inc(tested);
  end;

  // Random sampling
  for i := 1 to Int64(Round(SAMPLES_UNI * SamplePct / 100.0)) do
  begin
    bits := Xorshift64;
    v.u := bits;
    cr := pfC(v.f);
    pr := pfP(v.f);
    if not BitsMatch(cr, pr) then
    begin
      Inc(mismatches);
      if diagShown < DiagMax then
      begin
        rc.f := cr; rp.f := pr;
        WriteLn(Format('  [%s] input=$%16.16x(%.17g)  C=$%16.16x  P=$%16.16x',
                       [name, bits, v.f, rc.u, rp.u]));
        Inc(diagShown);
      end;
    end;
    Inc(tested);
  end;

  ReportResult(name, tested, mismatches, maxUlp);
end;

// ---- bivariate sampled test ----

procedure TestBivar(const name: string; pfC: TBivarFuncC; pfP: TBivarFuncP);
var
  i: Int32;
  vx, vy, rc, rp: Tb64u64;
  cr, pr: Double;
  mismatches: Int64;
  diagShown: Int32;
  maxUlp: Double;
begin
  if (Filter <> '') and (LowerCase(name) <> Filter) then Exit;
  mismatches := 0;
  diagShown := 0;
  maxUlp := 0;
  for i := 1 to Int64(Round(SAMPLES_BIV * SamplePct / 100.0)) do
  begin
    vx.u := Xorshift64;
    vy.u := Xorshift64;
    cr := pfC(vx.f, vy.f);
    pr := pfP(vx.f, vy.f);
    if not BitsMatch(cr, pr) then
    begin
      Inc(mismatches);
      if diagShown < DiagMax then
      begin
        rc.f := cr; rp.f := pr;
        WriteLn(Format('  [%s] x=$%16.16x y=$%16.16x  C=$%16.16x  P=$%16.16x',
                       [name, vx.u, vy.u, rc.u, rp.u]));
        Inc(diagShown);
      end;
    end;
  end;
  ReportResult(name, Int64(Round(SAMPLES_BIV * SamplePct / 100.0)), mismatches, maxUlp);
end;

// ---- sincos test ----

procedure TestSinCos;
var
  i: Int32;
  v: Tb64u64;
  cs, cc, ps, pc: Double;
  mismatches, tested: Int64;
begin
  if (Filter <> '') and (Filter <> 'sincos') then Exit;
  mismatches := 0;
  tested := 0;
  for i := 0 to SPECIAL_COUNT - 1 do
  begin
    v.u := SPECIAL_BITS[i];
    cr_sincos(v.f, @cs, @cc);
    pcr_sincos(v.f, ps, pc);
    if not BitsMatch(cs, ps) then Inc(mismatches);
    if not BitsMatch(cc, pc) then Inc(mismatches);
    Inc(tested);
  end;
  for i := 1 to Int64(Round(SAMPLES_UNI * SamplePct / 100.0)) do
  begin
    v.u := Xorshift64;
    cr_sincos(v.f, @cs, @cc);
    pcr_sincos(v.f, ps, pc);
    if not BitsMatch(cs, ps) then Inc(mismatches);
    if not BitsMatch(cc, pc) then Inc(mismatches);
    Inc(tested);
  end;
  ReportResult('sincos', tested, mismatches, 0);
end;

// cdecl wrappers for bivariate C functions
function wrap_atan2_c(y, x: Double): Double; cdecl;   begin Result := cr_atan2(y, x);    end;
function wrap_atan2pi_c(y, x: Double): Double; cdecl; begin Result := cr_atan2pi(y, x);  end;
function wrap_hypot_c(x, y: Double): Double; cdecl;   begin Result := cr_hypot(x, y);    end;
function wrap_pow_c(x, y: Double): Double; cdecl;     begin Result := cr_pow(x, y);      end;

// Pascal bivariate wrappers
function wrap_atan2_p(y, x: Double): Double;   begin Result := pcr_atan2(y, x);   end;
function wrap_atan2pi_p(y, x: Double): Double; begin Result := pcr_atan2pi(y, x); end;
function wrap_hypot_p(x, y: Double): Double;   begin Result := pcr_hypot(x, y);   end;
function wrap_pow_p(x, y: Double): Double;     begin Result := pcr_pow(x, y);     end;

procedure ParseArgs;
var
  i: Int32;
begin
  DiagMax := 3;
  i := 1;
  while i <= ParamCount do
  begin
    if (ParamStr(i) = '--diag') and (i < ParamCount) then
      DiagMax := StrToIntDef(ParamStr(i + 1), 3)
    else if (ParamStr(i) = '--func') and (i < ParamCount) then
      Filter := LowerCase(ParamStr(i + 1))
    else if (ParamStr(i) = '--seed') and (i < ParamCount) then
      RngState := StrToQWordDef(ParamStr(i + 1), RngState)
    else if (ParamStr(i) = '--pct') and (i < ParamCount) then
      SamplePct := StrToFloatDef(ParamStr(i + 1), 100.0);
    Inc(i);
  end;
end;

begin
  // Initialize special values array (high-bit constants are handled here to avoid
  // range-check warnings in typed constant declarations)
  SPECIAL_BITS[0]  := SP_POS_ZERO;
  SPECIAL_BITS[1]  := UInt64(High(Int64)) + 1;    // $8000000000000000 = -0
  SPECIAL_BITS[2]  := SP_ONE;
  SPECIAL_BITS[3]  := SP_ONE or (UInt64(1) shl 63); // $BFF0000000000000 = -1.0
  SPECIAL_BITS[4]  := SP_TWO;
  SPECIAL_BITS[5]  := SP_HALF;
  SPECIAL_BITS[6]  := SP_INF;
  SPECIAL_BITS[7]  := SP_INF  or (UInt64(1) shl 63); // $FFF0000000000000 = -Inf
  SPECIAL_BITS[8]  := SP_QNAN_POS;
  SPECIAL_BITS[9]  := SP_QNAN_POS or (UInt64(1) shl 63); // neg quiet NaN
  SPECIAL_BITS[10] := SP_SNAN;
  SPECIAL_BITS[11] := SP_MIN_SUB;
  SPECIAL_BITS[12] := SP_MAX_SUB;
  SPECIAL_BITS[13] := SP_MIN_NORM;
  SPECIAL_BITS[14] := SP_MAX_NORM;
  SPECIAL_BITS[15] := SP_PI;

  {$IFDEF AVX2}
  WriteLn('Compiled with AVX2.');
  {$ENDIF}

  SetExceptionMask([exInvalidOp, exDenormalized, exZeroDivide,
                    exOverflow, exUnderflow, exPrecision]);

  ParseArgs;
  StartTick := GetTickCount64;

  TotalPass := 0;
  TotalFail := 0;

  WriteLn('=== TestHarness64: comparing Pascal (pcr_*) vs C (cr_*) ===');
  WriteLn(Format('(random sampling: %d inputs per univariate function)', [SAMPLES_UNI]));
  if Filter <> '' then
    WriteLn(Format('(filter: %s)', [Filter]));
  WriteLn;

  // Univariate functions — add each after pcr_* is implemented in pascoremath64.pas
  TestUni('acos',    @cr_acos,    @pcr_acos);
  TestUni('acosh',   @cr_acosh,   @pcr_acosh);
  TestUni('acospi',  @cr_acospi,  @pcr_acospi);
  TestUni('asin',    @cr_asin,    @pcr_asin);
  TestUni('asinh',   @cr_asinh,   @pcr_asinh);
  TestUni('asinpi',  @cr_asinpi,  @pcr_asinpi);
  TestUni('atan',    @cr_atan,    @pcr_atan);
  TestUni('atanh',   @cr_atanh,   @pcr_atanh);
  TestUni('atanpi',  @cr_atanpi,  @pcr_atanpi);
  TestUni('cbrt',    @cr_cbrt,    @pcr_cbrt);
  TestUni('cos',     @cr_cos,     @pcr_cos);
  TestUni('cosh',    @cr_cosh,    @pcr_cosh);
  TestUni('cospi',   @cr_cospi,   @pcr_cospi);
  TestUni('erf',     @cr_erf,     @pcr_erf);
  TestUni('erfc',    @cr_erfc,    @pcr_erfc);
  TestUni('exp',     @cr_exp,     @pcr_exp);
  TestUni('exp10',   @cr_exp10,   @pcr_exp10);
  TestUni('exp10m1', @cr_exp10m1, @pcr_exp10m1);
  TestUni('exp2',    @cr_exp2,    @pcr_exp2);
  TestUni('exp2m1',  @cr_exp2m1,  @pcr_exp2m1);
  TestUni('expm1',   @cr_expm1,   @pcr_expm1);
  TestUni('lgamma',  @cr_lgamma,  @pcr_lgamma);
  TestUni('log',     @cr_log,     @pcr_log);
  TestUni('log10',   @cr_log10,   @pcr_log10);
  TestUni('log10p1', @cr_log10p1, @pcr_log10p1);
  TestUni('log1p',   @cr_log1p,   @pcr_log1p);
  TestUni('log2',    @cr_log2,    @pcr_log2);
  TestUni('log2p1',  @cr_log2p1,  @pcr_log2p1);
  TestUni('rsqrt',   @cr_rsqrt,   @pcr_rsqrt);
  TestUni('sin',     @cr_sin,     @pcr_sin);
  TestUni('sinh',    @cr_sinh,    @pcr_sinh);
  TestUni('sinpi',   @cr_sinpi,   @pcr_sinpi);
  TestUni('tan',     @cr_tan,     @pcr_tan);
  TestUni('tanh',    @cr_tanh,    @pcr_tanh);
  TestUni('tanpi',   @cr_tanpi,   @pcr_tanpi);
  TestUni('tgamma',  @cr_tgamma,  @pcr_tgamma);

  // Bivariate functions
  TestBivar('atan2',   @wrap_atan2_c,   @wrap_atan2_p);
  TestBivar('atan2pi', @wrap_atan2pi_c, @wrap_atan2pi_p);
  TestBivar('hypot',   @wrap_hypot_c,   @wrap_hypot_p);
  TestBivar('pow',     @wrap_pow_c,     @wrap_pow_p);

  // sincos
  TestSinCos;

  WriteLn;
  WriteLn(Format('Elapsed time: %.3f s', [(GetTickCount64 - StartTick) / 1000.0]));
  if (Filter <> '') and (TotalPass + TotalFail = 0) then
  begin
    WriteLn(Format('No function matched filter %s', [Filter]));
    Halt(1);
  end;
  WriteLn(Format('=== TOTAL: %d PASS, %d FAIL ===', [TotalPass, TotalFail]));
  if TotalFail = 0 then
    WriteLn('OVERALL: PASS')
  else
    WriteLn('OVERALL: FAIL');

  if TotalFail > 0 then
    Halt(1);
end.
