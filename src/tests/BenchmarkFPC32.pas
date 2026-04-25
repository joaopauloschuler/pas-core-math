// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
{$I ../pascoremath.inc}
program BenchmarkFPC32;

uses
  pascoremathtypes, pascoremath32, SysUtils, Math, DateUtils, StrUtils;

type
  TUniFuncP   = function(x: Single): Single;
  TBivarFuncP = function(x, y: Single): Single;

const
  BENCH_N        = 50000000;
  STRIDE         = High(Cardinal) div BENCH_N;
  TIE_THRESHOLD  = 0.05;

var
  GlobalSink: UInt32 = 0;
  PCMWins, FPCWins, PTies: Int32;
  TotalSpeedup: Double = 0.0;
  BenchCount: Int32 = 0;
  Filter: string = '';

procedure BenchUni(const name: string; pfFPC: TUniFuncP; pfPCM: TUniFuncP);
var
  i: Int32;
  u: Cardinal;
  v, r: Tb32u32;
  sink: UInt32;
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
  ux, uy: Cardinal;
  vx, vy, r: Tb32u32;
  sink: UInt32;
  t0, t1: TDateTime;
  msFPC, msPCM: Int64;
  mopsFPC, mopsPCM: Double;
begin
  if (Filter <> '') and (LowerCase(name) <> Filter) then Exit;
  // FPC version
  sink := 0;
  ux := 0;
  uy := High(Cardinal) div 3;
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
  uy := High(Cardinal) div 3;
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
  u: Cardinal;
  v, rs, rc: Tb32u32;
  sink: UInt32;
  t0, t1: TDateTime;
  msFPC, msPCM: Int64;
  mopsFPC, mopsPCM: Double;
  ps, pc: Single;
  fs, fc: Double;
begin
  if (Filter <> '') and (Filter <> 'sincosf') then Exit;
  // FPC version
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
    pcr_sincosf(v.f, ps, pc);
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
    ['sincosf', mopsFPC, mopsPCM,
     IfThen(mopsPCM > mopsFPC * (1.0 + TIE_THRESHOLD), '  FASTER! YAY!',
       IfThen(mopsFPC <= mopsPCM * (1.0 + TIE_THRESHOLD), '  TIE', ''))]));
end;

// FPC wrappers (Single -> Double -> Single)
function fpc_sinf(x: Single): Single;    begin Result := Sin(x);        end;
function fpc_cosf(x: Single): Single;    begin Result := Cos(x);        end;
function fpc_tanf(x: Single): Single;    begin Result := Tan(x);        end;
function fpc_asinf(x: Single): Single;   begin Result := ArcSin(x);     end;
function fpc_acosf(x: Single): Single;   begin Result := ArcCos(x);     end;
function fpc_atanf(x: Single): Single;   begin Result := ArcTan(x);     end;
function fpc_sinhf(x: Single): Single;   begin Result := Sinh(x);       end;
function fpc_coshf(x: Single): Single;   begin Result := Cosh(x);       end;
function fpc_tanhf(x: Single): Single;   begin Result := Tanh(x);       end;
function fpc_asinhf(x: Single): Single;  begin Result := ArcSinh(x);    end;
function fpc_acoshf(x: Single): Single;  begin Result := ArcCosh(x);    end;
function fpc_atanhf(x: Single): Single;  begin Result := ArcTanh(x);    end;
function fpc_expf(x: Single): Single;    begin Result := Exp(x);        end;
function fpc_logf(x: Single): Single;    begin Result := Ln(x);         end;
function fpc_log2f(x: Single): Single;   begin Result := Log2(x);       end;
function fpc_log10f(x: Single): Single;  begin Result := Log10(x);      end;
// FPC bivariate wrappers
function fpc_atan2f(y, x: Single): Single;  begin Result := ArcTan2(y, x);  end;
function fpc_hypotf(x, y: Single): Single;  begin Result := Hypot(x, y);    end;
function fpc_powf(x, y: Single): Single;    begin Result := Power(x, y);    end;

// PCM bivariate wrappers
function pcm_atan2f(y, x: Single): Single;  begin Result := pcr_atan2f(y, x);  end;
function pcm_hypotf(x, y: Single): Single;  begin Result := pcr_hypotf(x, y);  end;
function pcm_powf(x, y: Single): Single;    begin Result := pcr_powf(x, y);    end;

begin
  PCMWins := 0; FPCWins := 0; PTies := 0;
  {$IFDEF AVX2}
  WriteLn('Compiled with AVX2.');
  {$ENDIF}
  SetExceptionMask([exInvalidOp, exDenormalized, exZeroDivide,
                    exOverflow, exUnderflow, exPrecision]);

  if ParamCount >= 1 then Filter := LowerCase(ParamStr(1));

  if Filter = '' then
    WriteLn(Format('=== FPC vs Pascal CORE-MATH (PCM) Benchmark: %d calls per function ===', [BENCH_N]))
  else
    WriteLn(Format('=== FPC vs Pascal CORE-MATH (PCM) Benchmark: %d calls per function (filter=%s) ===', [BENCH_N, Filter]));
  WriteLn;

  BenchUni('sinf',    @fpc_sinf,    @pcr_sinf);
  BenchUni('cosf',    @fpc_cosf,    @pcr_cosf);
  BenchUni('tanf',    @fpc_tanf,    @pcr_tanf);
  BenchUni('asinf',   @fpc_asinf,   @pcr_asinf);
  BenchUni('acosf',   @fpc_acosf,   @pcr_acosf);
  BenchUni('atanf',   @fpc_atanf,   @pcr_atanf);
  BenchUni('sinhf',   @fpc_sinhf,   @pcr_sinhf);
  BenchUni('coshf',   @fpc_coshf,   @pcr_coshf);
  BenchUni('tanhf',   @fpc_tanhf,   @pcr_tanhf);
  BenchUni('asinhf',  @fpc_asinhf,  @pcr_asinhf);
  BenchUni('acoshf',  @fpc_acoshf,  @pcr_acoshf);
  BenchUni('atanhf',  @fpc_atanhf,  @pcr_atanhf);
  BenchUni('expf',    @fpc_expf,    @pcr_expf);
  BenchUni('logf',    @fpc_logf,    @pcr_logf);
  BenchUni('log2f',   @fpc_log2f,   @pcr_log2f);
  BenchUni('log10f',  @fpc_log10f,  @pcr_log10f);
  BenchBivar('atan2f',  @fpc_atan2f,  @pcm_atan2f);
  BenchBivar('hypotf',  @fpc_hypotf,  @pcm_hypotf);
  BenchBivar('powf',    @fpc_powf,    @pcm_powf);

  BenchSinCos;

  WriteLn;
  if BenchCount = 0 then
    WriteLn(Format('No function matched filter %s', [Filter]))
  else if BenchCount = 1 then
    WriteLn(Format('GlobalSink = %u', [GlobalSink]))
  else
  begin
    WriteLn(Format('PCM won: %d  |  FPC won: %d  |  Ties (<%d%%): %d',
      [PCMWins, FPCWins, Round(TIE_THRESHOLD * 100), PTies]));
    if TotalSpeedup / BenchCount >= 1.0 then
      WriteLn(Format('On average, PCM is %.2fx faster than FPC (arithmetic mean over %d functions)',
        [TotalSpeedup / BenchCount, BenchCount]))
    else
      WriteLn(Format('On average, PCM is %.2fx slower than FPC (arithmetic mean over %d functions)',
        [BenchCount / TotalSpeedup, BenchCount]));
    WriteLn(Format('GlobalSink = %u (prevents dead-code elimination)', [GlobalSink]));
  end;
end.
