// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
{$I ../pascoremath.inc}
program Benchmark64;

uses
  pascoremathtypes, pascoremath64, ccoremath64, SysUtils, Math, DateUtils, StrUtils;

type
  TUniFuncC   = function(x: Double): Double; cdecl;
  TUniFuncP   = function(x: Double): Double;
  TBivarFuncC = function(x, y: Double): Double; cdecl;
  TBivarFuncP = function(x, y: Double): Double;

const
  BENCH_N        = 200000000;
  // Stride across the UInt64 space: keep increment odd to visit many exponent classes
  STRIDE: UInt64 = $7E3779B97F4A7C15;  // large odd increment to cover many exponent classes
  TIE_THRESHOLD  = 0.05;

var
  GlobalSink: UInt64 = 0;
  PWins, CWins, PTies: Int32;
  TotalSpeedup: Double = 0.0;
  BenchCount: Int32 = 0;
  Filter: string = '';

procedure BenchUni(const name: string; pfC: TUniFuncC; pfP: TUniFuncP);
var
  i: Int32;
  u: UInt64;
  v, r: Tb64u64;
  sink: UInt64;
  t0, t1: TDateTime;
  msC, msP: Int64;
  mopsC, mopsP: Double;
  cSink, pSink: UInt64;
begin
  if (Filter <> '') and (LowerCase(name) <> Filter) then Exit;
  // C version
  sink := 0;
  u := 0;
  t0 := Now;
  for i := 1 to BENCH_N do
  begin
    v.u := u;
    r.f := pfC(v.f);
    sink := sink xor r.u;
    Inc(u, STRIDE);
  end;
  t1 := Now;
  msC := MillisecondsBetween(t1, t0);
  cSink := sink;

  // Pascal version
  sink := 0;
  u := 0;
  t0 := Now;
  for i := 1 to BENCH_N do
  begin
    v.u := u;
    r.f := pfP(v.f);
    sink := sink xor r.u;
    Inc(u, STRIDE);
  end;
  t1 := Now;
  msP := MillisecondsBetween(t1, t0);
  pSink := sink;

  if msC > 0 then mopsC := BENCH_N / msC / 1000.0 else mopsC := 999.9;
  if msP > 0 then mopsP := BENCH_N / msP / 1000.0 else mopsP := 999.9;

  GlobalSink := GlobalSink xor cSink xor pSink;

  if mopsP > mopsC * (1.0 + TIE_THRESHOLD) then Inc(PWins)
  else if mopsC > mopsP * (1.0 + TIE_THRESHOLD) then Inc(CWins)
  else Inc(PTies);
  if mopsC > 0 then
  begin
    TotalSpeedup := TotalSpeedup + mopsP / mopsC;
    Inc(BenchCount);
  end;
  WriteLn(Format('%-16s  C: %6.1f Mops/s  Pascal: %6.1f Mops/s  sink=%s%s',
    [name, mopsC, mopsP,
     IfThen(cSink = pSink, 'MATCH', 'MISMATCH'),
     IfThen(mopsP > mopsC * (1.0 + TIE_THRESHOLD), '  FASTER! YAY!',
       IfThen(mopsC <= mopsP * (1.0 + TIE_THRESHOLD), '  TIE', ''))]));
end;

procedure BenchBivar(const name: string; pfC: TBivarFuncC; pfP: TBivarFuncP);
var
  i: Int32;
  ux, uy: UInt64;
  vx, vy, r: Tb64u64;
  sink: UInt64;
  t0, t1: TDateTime;
  msC, msP: Int64;
  mopsC, mopsP: Double;
  cSink, pSink: UInt64;
begin
  if (Filter <> '') and (LowerCase(name) <> Filter) then Exit;
  // C version
  sink := 0;
  ux := 0;
  uy := $5555555555555555;
  t0 := Now;
  for i := 1 to BENCH_N do
  begin
    vx.u := ux;
    vy.u := uy;
    r.f := pfC(vx.f, vy.f);
    sink := sink xor r.u;
    Inc(ux, STRIDE);
    Inc(uy, STRIDE + 1);
  end;
  t1 := Now;
  msC := MillisecondsBetween(t1, t0);
  cSink := sink;

  // Pascal version
  sink := 0;
  ux := 0;
  uy := $5555555555555555;
  t0 := Now;
  for i := 1 to BENCH_N do
  begin
    vx.u := ux;
    vy.u := uy;
    r.f := pfP(vx.f, vy.f);
    sink := sink xor r.u;
    Inc(ux, STRIDE);
    Inc(uy, STRIDE + 1);
  end;
  t1 := Now;
  msP := MillisecondsBetween(t1, t0);
  pSink := sink;

  if msC > 0 then mopsC := BENCH_N / msC / 1000.0 else mopsC := 999.9;
  if msP > 0 then mopsP := BENCH_N / msP / 1000.0 else mopsP := 999.9;

  GlobalSink := GlobalSink xor cSink xor pSink;

  if mopsP > mopsC * (1.0 + TIE_THRESHOLD) then Inc(PWins)
  else if mopsC > mopsP * (1.0 + TIE_THRESHOLD) then Inc(CWins)
  else Inc(PTies);
  if mopsC > 0 then
  begin
    TotalSpeedup := TotalSpeedup + mopsP / mopsC;
    Inc(BenchCount);
  end;
  WriteLn(Format('%-16s  C: %6.1f Mops/s  Pascal: %6.1f Mops/s  sink=%s%s',
    [name, mopsC, mopsP,
     IfThen(cSink = pSink, 'MATCH', 'MISMATCH'),
     IfThen(mopsP > mopsC * (1.0 + TIE_THRESHOLD), '  FASTER! YAY!',
       IfThen(mopsC <= mopsP * (1.0 + TIE_THRESHOLD), '  TIE', ''))]));
end;

procedure BenchSinCos;
var
  i: Int32;
  u: UInt64;
  v, rs, rc: Tb64u64;
  sink: UInt64;
  t0, t1: TDateTime;
  msC, msP: Int64;
  mopsC, mopsP: Double;
  cSink, pSink: UInt64;
  ps, pc: Double;
begin
  if (Filter <> '') and (Filter <> 'sincos') then Exit;
  // C version
  sink := 0;
  u := 0;
  t0 := Now;
  for i := 1 to BENCH_N do
  begin
    v.u := u;
    cr_sincos(v.f, @rs.f, @rc.f);
    sink := sink xor rs.u xor rc.u;
    Inc(u, STRIDE);
  end;
  t1 := Now;
  msC := MillisecondsBetween(t1, t0);
  cSink := sink;

  // Pascal version
  sink := 0;
  u := 0;
  t0 := Now;
  for i := 1 to BENCH_N do
  begin
    v.u := u;
    pcr_sincos(v.f, ps, pc);
    rs.f := ps;
    rc.f := pc;
    sink := sink xor rs.u xor rc.u;
    Inc(u, STRIDE);
  end;
  t1 := Now;
  msP := MillisecondsBetween(t1, t0);
  pSink := sink;

  if msC > 0 then mopsC := BENCH_N / msC / 1000.0 else mopsC := 999.9;
  if msP > 0 then mopsP := BENCH_N / msP / 1000.0 else mopsP := 999.9;

  GlobalSink := GlobalSink xor cSink xor pSink;

  if mopsP > mopsC * (1.0 + TIE_THRESHOLD) then Inc(PWins)
  else if mopsC > mopsP * (1.0 + TIE_THRESHOLD) then Inc(CWins)
  else Inc(PTies);
  if mopsC > 0 then
  begin
    TotalSpeedup := TotalSpeedup + mopsP / mopsC;
    Inc(BenchCount);
  end;
  WriteLn(Format('%-16s  C: %6.1f Mops/s  Pascal: %6.1f Mops/s  sink=%s%s',
    ['sincos', mopsC, mopsP,
     IfThen(cSink = pSink, 'MATCH', 'MISMATCH'),
     IfThen(mopsP > mopsC * (1.0 + TIE_THRESHOLD), '  FASTER! YAY!',
       IfThen(mopsC <= mopsP * (1.0 + TIE_THRESHOLD), '  TIE', ''))]));
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

begin
  {$IFDEF AVX2}
  WriteLn('Compiled with AVX2.');
  {$ENDIF}
  SetExceptionMask([exInvalidOp, exDenormalized, exZeroDivide,
                    exOverflow, exUnderflow, exPrecision]);

  if ParamCount >= 1 then Filter := LowerCase(ParamStr(1));

  if Filter = '' then
    WriteLn(Format('=== Benchmark64: %d calls per function ===', [BENCH_N]))
  else
    WriteLn(Format('=== Benchmark64: %d calls per function (filter=%s) ===', [BENCH_N, Filter]));
  WriteLn;

  BenchUni('acos',    @cr_acos,    @pcr_acos);
  BenchUni('acosh',   @cr_acosh,   @pcr_acosh);
  BenchUni('acospi',  @cr_acospi,  @pcr_acospi);
  BenchUni('asin',    @cr_asin,    @pcr_asin);
  BenchUni('asinh',   @cr_asinh,   @pcr_asinh);
  BenchUni('asinpi',  @cr_asinpi,  @pcr_asinpi);
  BenchUni('atan',    @cr_atan,    @pcr_atan);
  BenchUni('atanh',   @cr_atanh,   @pcr_atanh);
  BenchUni('atanpi',  @cr_atanpi,  @pcr_atanpi);
  BenchUni('cbrt',    @cr_cbrt,    @pcr_cbrt);
  BenchUni('cos',     @cr_cos,     @pcr_cos);
  BenchUni('cosh',    @cr_cosh,    @pcr_cosh);
  BenchUni('cospi',   @cr_cospi,   @pcr_cospi);
  BenchUni('erf',     @cr_erf,     @pcr_erf);
  BenchUni('erfc',    @cr_erfc,    @pcr_erfc);
  BenchUni('exp',     @cr_exp,     @pcr_exp);
  BenchUni('exp10',   @cr_exp10,   @pcr_exp10);
  BenchUni('exp10m1', @cr_exp10m1, @pcr_exp10m1);
  BenchUni('exp2',    @cr_exp2,    @pcr_exp2);
  BenchUni('exp2m1',  @cr_exp2m1,  @pcr_exp2m1);
  BenchUni('expm1',   @cr_expm1,   @pcr_expm1);
  BenchUni('lgamma',  @cr_lgamma,  @pcr_lgamma);
  BenchUni('log',     @cr_log,     @pcr_log);
  BenchUni('log10',   @cr_log10,   @pcr_log10);
  BenchUni('log10p1', @cr_log10p1, @pcr_log10p1);
  BenchUni('log1p',   @cr_log1p,   @pcr_log1p);
  BenchUni('log2',    @cr_log2,    @pcr_log2);
  BenchUni('log2p1',  @cr_log2p1,  @pcr_log2p1);
  BenchUni('rsqrt',   @cr_rsqrt,   @pcr_rsqrt);
  BenchUni('sin',     @cr_sin,     @pcr_sin);
  BenchUni('sinh',    @cr_sinh,    @pcr_sinh);
  BenchUni('sinpi',   @cr_sinpi,   @pcr_sinpi);
  BenchUni('tan',     @cr_tan,     @pcr_tan);
  BenchUni('tanh',    @cr_tanh,    @pcr_tanh);
  BenchUni('tanpi',   @cr_tanpi,   @pcr_tanpi);
  BenchUni('tgamma',  @cr_tgamma,  @pcr_tgamma);

  BenchBivar('atan2',   @wrap_atan2_c,   @wrap_atan2_p);
  BenchBivar('atan2pi', @wrap_atan2pi_c, @wrap_atan2pi_p);
  BenchBivar('hypot',   @wrap_hypot_c,   @wrap_hypot_p);
  BenchBivar('pow',     @wrap_pow_c,     @wrap_pow_p);

  BenchSinCos;

  WriteLn;
  if BenchCount = 0 then
    WriteLn(Format('No function matched filter %s', [Filter]))
  else if BenchCount = 1 then
    WriteLn(Format('GlobalSink = %u', [GlobalSink]))
  else
  begin
    WriteLn(Format('Pascal won: %d  |  C won: %d  |  Ties (<%d%%): %d', [PWins, CWins, Round(TIE_THRESHOLD * 100), PTies]));
    if TotalSpeedup / BenchCount >= 1.0 then
      WriteLn(Format('On average, Pascal is %.2fx faster than C (arithmetic mean over %d functions)',
        [TotalSpeedup / BenchCount, BenchCount]))
    else
      WriteLn(Format('On average, Pascal is %.2fx slower than C (arithmetic mean over %d functions)',
        [BenchCount / TotalSpeedup, BenchCount]));
    WriteLn(Format('GlobalSink = %u (prevents dead-code elimination)', [GlobalSink]));
  end;
end.
