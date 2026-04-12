// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//                                                                                                                                                                                                      
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
{$I ../pascoremath.inc}
program Benchmark;

uses
  pascoremathtypes, pascoremath, ccoremath, SysUtils, Math, DateUtils, StrUtils;

type
  TUniFuncC   = function(x: Single): Single; cdecl;
  TUniFuncP   = function(x: Single): Single;
  TBivarFuncC = function(x, y: Single): Single; cdecl;
  TBivarFuncP = function(x, y: Single): Single;

const
  BENCH_N = 50000000;
  STRIDE  = High(Cardinal) div BENCH_N;

var
  GlobalSink: LongWord = 0;

procedure BenchUni(const name: string; pfC: TUniFuncC; pfP: TUniFuncP);
var
  i: Integer;
  u: Cardinal;
  v, r: Tb32u32;
  sink: LongWord;
  t0, t1: TDateTime;
  msC, msP: Int64;
  mopsC, mopsP: Double;
  cSink, pSink: LongWord;
begin
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

  WriteLn(Format('%-16s  C: %6.1f Mops/s  Pascal: %6.1f Mops/s  sink=%s%s',
    [name, mopsC, mopsP,
     IfThen(cSink = pSink, 'MATCH', 'MISMATCH'),
     IfThen(mopsP > mopsC, '  FASTER! YAY!', '')]));
end;

procedure BenchBivar(const name: string; pfC: TBivarFuncC; pfP: TBivarFuncP);
var
  i: Integer;
  ux, uy: Cardinal;
  vx, vy, r: Tb32u32;
  sink: LongWord;
  t0, t1: TDateTime;
  msC, msP: Int64;
  mopsC, mopsP: Double;
  cSink, pSink: LongWord;
begin
  // C version
  sink := 0;
  ux := 0;
  uy := High(Cardinal) div 3;
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
  uy := High(Cardinal) div 3;
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

  WriteLn(Format('%-16s  C: %6.1f Mops/s  Pascal: %6.1f Mops/s  sink=%s%s',
    [name, mopsC, mopsP,
     IfThen(cSink = pSink, 'MATCH', 'MISMATCH'),
     IfThen(mopsP > mopsC, '  FASTER! YAY!', '')]));
end;

procedure BenchSinCos;
var
  i: Integer;
  u: Cardinal;
  v, rs, rc: Tb32u32;
  sink: LongWord;
  t0, t1: TDateTime;
  msC, msP: Int64;
  mopsC, mopsP: Double;
  cSink, pSink: LongWord;
  ps, pc: Single;
begin
  // C version
  sink := 0;
  u := 0;
  t0 := Now;
  for i := 1 to BENCH_N do
  begin
    v.u := u;
    cr_sincosf(v.f, @rs.f, @rc.f);
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
    pcr_sincosf(v.f, ps, pc);
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

  WriteLn(Format('%-16s  C: %6.1f Mops/s  Pascal: %6.1f Mops/s  sink=%s%s',
    ['sincosf', mopsC, mopsP,
     IfThen(cSink = pSink, 'MATCH', 'MISMATCH'),
     IfThen(mopsP > mopsC, '  FASTER! YAY!', '')]));
end;

// cdecl wrappers for bivariate C functions
function wrap_atan2_c(y, x: Single): Single; cdecl;    begin Result := cr_atan2f(y, x);     end;
function wrap_atan2pi_c(y, x: Single): Single; cdecl;  begin Result := cr_atan2pif(y, x);   end;
function wrap_hypot_c(x, y: Single): Single; cdecl;    begin Result := cr_hypotf(x, y);     end;
function wrap_pow_c(x, y: Single): Single; cdecl;      begin Result := cr_powf(x, y);       end;
function wrap_compound_c(x, y: Single): Single; cdecl; begin Result := cr_compoundf(x, y);  end;

// Pascal bivariate wrappers (register convention)
function wrap_atan2_p(y, x: Single): Single;    begin Result := pcr_atan2f(y, x);    end;
function wrap_atan2pi_p(y, x: Single): Single;  begin Result := pcr_atan2pif(y, x);  end;
function wrap_hypot_p(x, y: Single): Single;    begin Result := pcr_hypotf(x, y);    end;
function wrap_pow_p(x, y: Single): Single;      begin Result := pcr_powf(x, y);      end;
function wrap_compound_p(x, y: Single): Single; begin Result := pcr_compoundf(x, y); end;

begin
  // Mask all FP exceptions: we iterate over all bit patterns including NaN/Inf
  SetExceptionMask([exInvalidOp, exDenormalized, exZeroDivide,
                    exOverflow, exUnderflow, exPrecision]);

  WriteLn(Format('=== Benchmark: %d calls per function ===', [BENCH_N]));
  WriteLn;

  BenchUni('acosf',    @cr_acosf,    @pcr_acosf);
  BenchUni('acoshf',   @cr_acoshf,   @pcr_acoshf);
  BenchUni('acospif',  @cr_acospif,  @pcr_acospif);
  BenchUni('asinf',    @cr_asinf,    @pcr_asinf);
  BenchUni('asinhf',   @cr_asinhf,   @pcr_asinhf);
  BenchUni('asinpif',  @cr_asinpif,  @pcr_asinpif);
  BenchUni('atanf',    @cr_atanf,    @pcr_atanf);
  BenchUni('atanhf',   @cr_atanhf,   @pcr_atanhf);
  BenchUni('atanpif',  @cr_atanpif,  @pcr_atanpif);
  BenchUni('cbrtf',    @cr_cbrtf,    @pcr_cbrtf);
  BenchUni('cosf',     @cr_cosf,     @pcr_cosf);
  BenchUni('coshf',    @cr_coshf,    @pcr_coshf);
  BenchUni('cospif',   @cr_cospif,   @pcr_cospif);
  BenchUni('erff',     @cr_erff,     @pcr_erff);
  BenchUni('erfcf',    @cr_erfcf,    @pcr_erfcf);
  BenchUni('expf',     @cr_expf,     @pcr_expf);
  BenchUni('exp10f',   @cr_exp10f,   @pcr_exp10f);
  BenchUni('exp10m1f', @cr_exp10m1f, @pcr_exp10m1f);
  BenchUni('exp2f',    @cr_exp2f,    @pcr_exp2f);
  BenchUni('exp2m1f',  @cr_exp2m1f,  @pcr_exp2m1f);
  BenchUni('expm1f',   @cr_expm1f,   @pcr_expm1f);
  BenchUni('lgammaf',  @cr_lgammaf,  @pcr_lgammaf);
  BenchUni('logf',     @cr_logf,     @pcr_logf);
  BenchUni('log10f',   @cr_log10f,   @pcr_log10f);
  BenchUni('log10p1f', @cr_log10p1f, @pcr_log10p1f);
  BenchUni('log1pf',   @cr_log1pf,   @pcr_log1pf);
  BenchUni('log2f',    @cr_log2f,    @pcr_log2f);
  BenchUni('log2p1f',  @cr_log2p1f,  @pcr_log2p1f);
  BenchUni('rsqrtf',   @cr_rsqrtf,   @pcr_rsqrtf);
  BenchUni('sinf',     @cr_sinf,     @pcr_sinf);
  BenchUni('sinhf',    @cr_sinhf,    @pcr_sinhf);
  BenchUni('sinpif',   @cr_sinpif,   @pcr_sinpif);
  BenchUni('tanf',     @cr_tanf,     @pcr_tanf);
  BenchUni('tanhf',    @cr_tanhf,    @pcr_tanhf);
  BenchUni('tanpif',   @cr_tanpif,   @pcr_tanpif);
  BenchUni('tgammaf',  @cr_tgammaf,  @pcr_tgammaf);

  BenchBivar('atan2f',    @wrap_atan2_c,    @wrap_atan2_p);
  BenchBivar('atan2pif',  @wrap_atan2pi_c,  @wrap_atan2pi_p);
  BenchBivar('hypotf',    @wrap_hypot_c,    @wrap_hypot_p);
  BenchBivar('powf',      @wrap_pow_c,      @wrap_pow_p);
  BenchBivar('compoundf', @wrap_compound_c, @wrap_compound_p);

  BenchSinCos;

  WriteLn;
  WriteLn(Format('GlobalSink = %u (prevents dead-code elimination)', [GlobalSink]));
end.
