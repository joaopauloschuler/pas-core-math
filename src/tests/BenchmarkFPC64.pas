// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
{$I ../pascoremath.inc}
program BenchmarkFPC64;

uses
  pascoremathtypes, pascoremath64, SysUtils, Math, DateUtils, StrUtils;

type
  TUniFuncP   = function(x: Double): Double;
  TBivarFuncP = function(x, y: Double): Double;

const
  BENCH_N        = 50000000;
  STRIDE: UInt64 = $7E3779B97F4A7C15;
  TIE_THRESHOLD  = 0.05;

var
  GlobalSink: UInt64 = 0;
  PCMWins, FPCWins, PTies: Int32;
  TotalSpeedup: Double = 0.0;
  BenchCount: Int32 = 0;
  Filter: string = '';

procedure BenchUni(const name: string; pfFPC: TUniFuncP; pfPCM: TUniFuncP);
var
  i: Int32;
  u: UInt64;
  v, r: Tb64u64;
  sink: UInt64;
  t0, t1: TDateTime;
  msFPC, msPCM: Int64;
  mopsFPC, mopsPCM: Double;
begin
  if (Filter <> '') and (LowerCase(name) <> Filter) then Exit;
  // FPC version
  sink := 0;
  u := 0;
  t0 := Now;
  for i := 1 to BENCH_N do
  begin
    v.u := u;
    r.f := pfFPC(v.f);
    sink := sink xor r.u;
    Inc(u, STRIDE);
  end;
  t1 := Now;
  msFPC := MillisecondsBetween(t1, t0);
  GlobalSink := GlobalSink xor sink;

  // PCM version
  sink := 0;
  u := 0;
  t0 := Now;
  for i := 1 to BENCH_N do
  begin
    v.u := u;
    r.f := pfPCM(v.f);
    sink := sink xor r.u;
    Inc(u, STRIDE);
  end;
  t1 := Now;
  msPCM := MillisecondsBetween(t1, t0);
  GlobalSink := GlobalSink xor sink;

  if msFPC > 0 then mopsFPC := BENCH_N / msFPC / 1000.0 else mopsFPC := 999.9;
  if msPCM > 0 then mopsPCM := BENCH_N / msPCM / 1000.0 else mopsPCM := 999.9;

  if mopsPCM > mopsFPC * (1.0 + TIE_THRESHOLD) then Inc(PCMWins)
  else if mopsFPC > mopsPCM * (1.0 + TIE_THRESHOLD) then Inc(FPCWins)
  else Inc(PTies);
  if mopsFPC > 0 then
  begin
    TotalSpeedup := TotalSpeedup + mopsPCM / mopsFPC;
    Inc(BenchCount);
  end;
  WriteLn(Format('%-16s  FPC: %6.1f Mops/s  PCM: %6.1f Mops/s%s',
    [name, mopsFPC, mopsPCM,
     IfThen(mopsPCM > mopsFPC * (1.0 + TIE_THRESHOLD), '  FASTER! YAY!',
       IfThen(mopsFPC <= mopsPCM * (1.0 + TIE_THRESHOLD), '  TIE', ''))]));
end;

procedure BenchBivar(const name: string; pfFPC: TBivarFuncP; pfPCM: TBivarFuncP);
var
  i: Int32;
  ux, uy: UInt64;
  vx, vy, r: Tb64u64;
  sink: UInt64;
  t0, t1: TDateTime;
  msFPC, msPCM: Int64;
  mopsFPC, mopsPCM: Double;
begin
  if (Filter <> '') and (LowerCase(name) <> Filter) then Exit;
  // FPC version
  sink := 0;
  ux := 0;
  uy := $5555555555555555;
  t0 := Now;
  for i := 1 to BENCH_N do
  begin
    vx.u := ux;
    vy.u := uy;
    r.f := pfFPC(vx.f, vy.f);
    sink := sink xor r.u;
    Inc(ux, STRIDE);
    Inc(uy, STRIDE + 1);
  end;
  t1 := Now;
  msFPC := MillisecondsBetween(t1, t0);
  GlobalSink := GlobalSink xor sink;

  // PCM version
  sink := 0;
  ux := 0;
  uy := $5555555555555555;
  t0 := Now;
  for i := 1 to BENCH_N do
  begin
    vx.u := ux;
    vy.u := uy;
    r.f := pfPCM(vx.f, vy.f);
    sink := sink xor r.u;
    Inc(ux, STRIDE);
    Inc(uy, STRIDE + 1);
  end;
  t1 := Now;
  msPCM := MillisecondsBetween(t1, t0);
  GlobalSink := GlobalSink xor sink;

  if msFPC > 0 then mopsFPC := BENCH_N / msFPC / 1000.0 else mopsFPC := 999.9;
  if msPCM > 0 then mopsPCM := BENCH_N / msPCM / 1000.0 else mopsPCM := 999.9;

  if mopsPCM > mopsFPC * (1.0 + TIE_THRESHOLD) then Inc(PCMWins)
  else if mopsFPC > mopsPCM * (1.0 + TIE_THRESHOLD) then Inc(FPCWins)
  else Inc(PTies);
  if mopsFPC > 0 then
  begin
    TotalSpeedup := TotalSpeedup + mopsPCM / mopsFPC;
    Inc(BenchCount);
  end;
  WriteLn(Format('%-16s  FPC: %6.1f Mops/s  PCM: %6.1f Mops/s%s',
    [name, mopsFPC, mopsPCM,
     IfThen(mopsPCM > mopsFPC * (1.0 + TIE_THRESHOLD), '  FASTER! YAY!',
       IfThen(mopsFPC <= mopsPCM * (1.0 + TIE_THRESHOLD), '  TIE', ''))]));
end;

procedure BenchSinCos;
var
  i: Int32;
  u: UInt64;
  v, rs, rc: Tb64u64;
  sink: UInt64;
  t0, t1: TDateTime;
  msFPC, msPCM: Int64;
  mopsFPC, mopsPCM: Double;
  ps, pc: Double;
  fs, fc: Double;
begin
  if (Filter <> '') and (Filter <> 'sincos') then Exit;
  // FPC version (SinCos from Math unit)
  sink := 0;
  u := 0;
  t0 := Now;
  for i := 1 to BENCH_N do
  begin
    v.u := u;
    SinCos(v.f, fs, fc);
    rs.f := fs;
    rc.f := fc;
    sink := sink xor rs.u xor rc.u;
    Inc(u, STRIDE);
  end;
  t1 := Now;
  msFPC := MillisecondsBetween(t1, t0);
  GlobalSink := GlobalSink xor sink;

  // PCM version
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
  msPCM := MillisecondsBetween(t1, t0);
  GlobalSink := GlobalSink xor sink;

  if msFPC > 0 then mopsFPC := BENCH_N / msFPC / 1000.0 else mopsFPC := 999.9;
  if msPCM > 0 then mopsPCM := BENCH_N / msPCM / 1000.0 else mopsPCM := 999.9;

  if mopsPCM > mopsFPC * (1.0 + TIE_THRESHOLD) then Inc(PCMWins)
  else if mopsFPC > mopsPCM * (1.0 + TIE_THRESHOLD) then Inc(FPCWins)
  else Inc(PTies);
  if mopsFPC > 0 then
  begin
    TotalSpeedup := TotalSpeedup + mopsPCM / mopsFPC;
    Inc(BenchCount);
  end;
  WriteLn(Format('%-16s  FPC: %6.1f Mops/s  PCM: %6.1f Mops/s%s',
    ['sincos', mopsFPC, mopsPCM,
     IfThen(mopsPCM > mopsFPC * (1.0 + TIE_THRESHOLD), '  FASTER! YAY!',
       IfThen(mopsFPC <= mopsPCM * (1.0 + TIE_THRESHOLD), '  TIE', ''))]));
end;

// FPC Double wrappers
function fpc_sin(x: Double): Double;     begin Result := Sin(x);       end;
function fpc_cos(x: Double): Double;     begin Result := Cos(x);       end;
function fpc_tan(x: Double): Double;     begin Result := Tan(x);       end;
function fpc_asin(x: Double): Double;    begin Result := ArcSin(x);    end;
function fpc_acos(x: Double): Double;    begin Result := ArcCos(x);    end;
function fpc_atan(x: Double): Double;    begin Result := ArcTan(x);    end;
function fpc_sinh(x: Double): Double;    begin Result := Sinh(x);      end;
function fpc_cosh(x: Double): Double;    begin Result := Cosh(x);      end;
function fpc_tanh(x: Double): Double;    begin Result := Tanh(x);      end;
function fpc_asinh(x: Double): Double;   begin Result := ArcSinh(x);   end;
function fpc_acosh(x: Double): Double;   begin Result := ArcCosh(x);   end;
function fpc_atanh(x: Double): Double;   begin Result := ArcTanh(x);   end;
function fpc_exp(x: Double): Double;     begin Result := Exp(x);       end;
function fpc_log(x: Double): Double;     begin Result := Ln(x);        end;
function fpc_log2(x: Double): Double;    begin Result := Log2(x);      end;
function fpc_log10(x: Double): Double;   begin Result := Log10(x);     end;
// FPC bivariate wrappers
function fpc_atan2(y, x: Double): Double;  begin Result := ArcTan2(y, x); end;
function fpc_hypot(x, y: Double): Double;  begin Result := Hypot(x, y);   end;
function fpc_pow(x, y: Double): Double;    begin Result := Power(x, y);   end;
// PCM bivariate wrappers
function pcm_atan2(y, x: Double): Double;  begin Result := pcr_atan2(y, x);  end;
function pcm_hypot(x, y: Double): Double;  begin Result := pcr_hypot(x, y);  end;
function pcm_pow(x, y: Double): Double;    begin Result := pcr_pow(x, y);    end;

begin
  {$IFDEF AVX2}
  WriteLn('Compiled with AVX2.');
  {$ENDIF}
  SetExceptionMask([exInvalidOp, exDenormalized, exZeroDivide,
                    exOverflow, exUnderflow, exPrecision]);

  if ParamCount >= 1 then Filter := LowerCase(ParamStr(1));

  if Filter = '' then
    WriteLn(Format('=== FPC vs Pascal CORE-MATH (PCM) Benchmark64: %d calls per function ===', [BENCH_N]))
  else
    WriteLn(Format('=== FPC vs PCM Benchmark64: %d calls (filter=%s) ===', [BENCH_N, Filter]));
  WriteLn;

  BenchUni('sin',    @fpc_sin,    @pcr_sin);
  BenchUni('cos',    @fpc_cos,    @pcr_cos);
  BenchUni('tan',    @fpc_tan,    @pcr_tan);
  BenchUni('asin',   @fpc_asin,   @pcr_asin);
  BenchUni('acos',   @fpc_acos,   @pcr_acos);
  BenchUni('atan',   @fpc_atan,   @pcr_atan);
  BenchUni('sinh',   @fpc_sinh,   @pcr_sinh);
  BenchUni('cosh',   @fpc_cosh,   @pcr_cosh);
  BenchUni('tanh',   @fpc_tanh,   @pcr_tanh);
  BenchUni('asinh',  @fpc_asinh,  @pcr_asinh);
  BenchUni('acosh',  @fpc_acosh,  @pcr_acosh);
  BenchUni('atanh',  @fpc_atanh,  @pcr_atanh);
  BenchUni('exp',    @fpc_exp,    @pcr_exp);
  BenchUni('log',    @fpc_log,    @pcr_log);
  BenchUni('log2',   @fpc_log2,   @pcr_log2);
  BenchUni('log10',  @fpc_log10,  @pcr_log10);

  BenchBivar('atan2',  @fpc_atan2,  @pcm_atan2);
  BenchBivar('hypot',  @fpc_hypot,  @pcm_hypot);
  BenchBivar('pow',    @fpc_pow,    @pcm_pow);

  BenchSinCos;

  WriteLn;
  if BenchCount = 0 then
    WriteLn(Format('No function matched filter %s', [Filter]))
  else if BenchCount = 1 then
    WriteLn(Format('GlobalSink = %u', [GlobalSink]))
  else
  begin
    WriteLn(Format('PCM won: %d  |  FPC won: %d  |  Ties (<%d%%): %d', [PCMWins, FPCWins, Round(TIE_THRESHOLD * 100), PTies]));
    if TotalSpeedup / BenchCount >= 1.0 then
      WriteLn(Format('On average, PCM is %.2fx faster than FPC (arithmetic mean over %d functions)',
        [TotalSpeedup / BenchCount, BenchCount]))
    else
      WriteLn(Format('On average, PCM is %.2fx slower than FPC (arithmetic mean over %d functions)',
        [BenchCount / TotalSpeedup, BenchCount]));
    WriteLn(Format('GlobalSink = %u (prevents dead-code elimination)', [GlobalSink]));
  end;
end.
