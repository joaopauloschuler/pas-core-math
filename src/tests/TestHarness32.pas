// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//                                                                                                                                                                                                      
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
{$I ../pascoremath.inc}
program TestHarness32;

uses
  pascoremathtypes, pascoremath32, ccoremath32, SysUtils, Math;

type
  TUniFuncC   = function(x: Single): Single; cdecl;
  TUniFuncP   = function(x: Single): Single;
  TBivarFuncC = function(x, y: Single): Single; cdecl;
  TBivarFuncP = function(x, y: Single): Single;

var
  TotalPass, TotalFail: Int32;
  Stride: Cardinal;
  DiagMax: Int32;  // max mismatches to print (0 = none)

procedure ReportResult(const FuncName: string; Tested, Mismatches: Int64; MaxError: Double);
begin
  if Mismatches = 0 then
  begin
    WriteLn(Format('%-16s  tested=%10d  mismatches=0  PASS', [FuncName, Tested]));
    Inc(TotalPass);
  end
  else
  begin
    WriteLn(Format('%-16s  tested=%10d  mismatches=%d  FAIL max_error=%.20f', [FuncName, Tested, Mismatches, MaxError]));
    Inc(TotalFail);
  end;
end;

// Bit-exact comparison: both NaN => match; otherwise must be identical bits
function BitsMatch(a, b: Single): Boolean; inline;
var
  ua, ub: Tb32u32;
begin
  ua.f := a;
  ub.f := b;
  // If a is NaN, b must also be NaN (any NaN)
  if ((ua.u and $7F800000) = $7F800000) and ((ua.u and $007FFFFF) <> 0) then
  begin
    Result := ((ub.u and $7F800000) = $7F800000) and ((ub.u and $007FFFFF) <> 0);
    Exit;
  end;
  Result := ua.u = ub.u;
end;

// ---- univariate exhaustive tests ----

procedure TestUni(const name: string; pfC: TUniFuncC; pfP: TUniFuncP);
var
  u: Cardinal;
  v, rc, rp: Tb32u32;
  cr, pr: Single;
  mismatches, tested: Int64;
  diagShown: Int32;
  error, max_error: Single;
begin
  mismatches := 0;
  tested := 0;
  diagShown := 0;
  max_error := 0;
  u := 0;
  repeat
    v.u := u;
    cr := pfC(v.f);
    pr := pfP(v.f);
    if not BitsMatch(cr, pr) then
    begin
      Inc(mismatches);
      if diagShown < DiagMax then
      begin
        rc.f := cr; rp.f := pr;
        WriteLn(Format('  [%s] input=$%8.8x  C=$%8.8x  P=$%8.8x', [name, u, rc.u, rp.u]));
        Inc(diagShown);
        error := abs(cr - pr);
        if error > max_error then max_error := error;
      end;
    end;
    Inc(tested);
    if u > High(Cardinal) - Stride then Break;
    Inc(u, Stride);
  until False;
  ReportResult(name, tested, mismatches, max_error);
end;

// ---- bivariate sampled tests ----

procedure TestBivar(const name: string; pfC: TBivarFuncC; pfP: TBivarFuncP);
const
  SAMPLES = 10000000;
  STRIDE  = High(Cardinal) div SAMPLES;
var
  i: Int32;
  ux, uy: Cardinal;
  vx, vy, rc, rp: Tb32u32;
  cr, pr: Single;
  mismatches: Int64;
  diagShown: Int32;
  error, max_error: Single;
begin
  mismatches := 0;
  diagShown := 0;
  max_error := 0;
  ux := 0;
  uy := High(Cardinal) div 3;
  for i := 1 to SAMPLES do
  begin
    vx.u := ux;
    vy.u := uy;
    cr := pfC(vx.f, vy.f);
    pr := pfP(vx.f, vy.f);
    if not BitsMatch(cr, pr) then
    begin
      Inc(mismatches);
      if diagShown < DiagMax then
      begin
        rc.f := cr; rp.f := pr;
        WriteLn(Format('  [%s] x=$%8.8x y=$%8.8x  C=$%8.8x  P=$%8.8x',
                       [name, ux, uy, rc.u, rp.u]));
        Inc(diagShown);
        error := abs(cr - pr);
        if error > max_error then max_error := error;
      end;
    end;
    Inc(ux, STRIDE);
    Inc(uy, STRIDE + 1);
  end;
  ReportResult(name, SAMPLES, mismatches, max_error);
end;

// ---- sincosf exhaustive test ----

procedure TestSinCos;
var
  u: Cardinal;
  v: Tb32u32;
  cs, cc, ps, pc: Single;
  mismatches, tested: Int64;
  error, max_error: Single;
begin
  mismatches := 0;
  tested := 0;
  max_error := 0;
  u := 0;
  repeat
    v.u := u;
    cr_sincosf(v.f, @cs, @cc);
    pcr_sincosf(v.f, ps, pc);
    if not BitsMatch(cs, ps) then
    begin
      Inc(mismatches);
      error := abs(cs - ps);
      if error > max_error then max_error := error;
    end;
    if not BitsMatch(cc, pc) then
    begin
      Inc(mismatches);
      error := abs(cc - pc);
      if error > max_error then max_error := error;
    end;
    Inc(tested);
    if u > High(Cardinal) - Stride then Break;
    Inc(u, Stride);
  until False;
  ReportResult('sincosf', tested, mismatches, max_error);
end;

// Pascal wrappers for bivariate C functions (adapts cdecl -> register)
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

function ParsePct: Cardinal;
var
  i: Int32;
  pct: Int32;
begin
  Result := 1; // default: 100% => stride 1
  DiagMax := 0;
  i := 1;
  while i <= ParamCount do
  begin
    if (ParamStr(i) = '--pct') and (i < ParamCount) then
    begin
      pct := StrToIntDef(ParamStr(i + 1), 100);
      if (pct < 1) or (pct > 100) then
      begin
        WriteLn('Error: --pct must be between 1 and 100');
        Halt(1);
      end;
      Result := 100 div pct;
    end
    else if (ParamStr(i) = '--diag') and (i < ParamCount) then
      DiagMax := StrToIntDef(ParamStr(i + 1), 3);
    Inc(i);
  end;
end;

begin
  // Mask all FP exceptions: we iterate over all bit patterns including NaN/Inf
  SetExceptionMask([exInvalidOp, exDenormalized, exZeroDivide,
                    exOverflow, exUnderflow, exPrecision]);

  Stride := ParsePct;

  TotalPass := 0;
  TotalFail := 0;

  WriteLn('=== TestHarness: comparing Pascal (pcr_*) vs C (cr_*) ===');
  if Stride > 1 then
    WriteLn(Format('(sampling mode: stride=%d, ~%d%% of inputs)', [Stride, 100 div Stride]));
  WriteLn;

  // Univariate exhaustive (36 functions)
  TestUni('acosf',    @cr_acosf,    @pcr_acosf);
  TestUni('acoshf',   @cr_acoshf,   @pcr_acoshf);
  TestUni('acospif',  @cr_acospif,  @pcr_acospif);
  TestUni('asinf',    @cr_asinf,    @pcr_asinf);
  TestUni('asinhf',   @cr_asinhf,   @pcr_asinhf);
  TestUni('asinpif',  @cr_asinpif,  @pcr_asinpif);
  TestUni('atanf',    @cr_atanf,    @pcr_atanf);
  TestUni('atanhf',   @cr_atanhf,   @pcr_atanhf);
  TestUni('atanpif',  @cr_atanpif,  @pcr_atanpif);
  TestUni('cbrtf',    @cr_cbrtf,    @pcr_cbrtf);
  TestUni('cosf',     @cr_cosf,     @pcr_cosf);
  TestUni('coshf',    @cr_coshf,    @pcr_coshf);
  TestUni('cospif',   @cr_cospif,   @pcr_cospif);
  TestUni('erff',     @cr_erff,     @pcr_erff);
  TestUni('erfcf',    @cr_erfcf,    @pcr_erfcf);
  TestUni('expf',     @cr_expf,     @pcr_expf);
  TestUni('exp10f',   @cr_exp10f,   @pcr_exp10f);
  TestUni('exp10m1f', @cr_exp10m1f, @pcr_exp10m1f);
  TestUni('exp2f',    @cr_exp2f,    @pcr_exp2f);
  TestUni('exp2m1f',  @cr_exp2m1f,  @pcr_exp2m1f);
  TestUni('expm1f',   @cr_expm1f,   @pcr_expm1f);
  TestUni('lgammaf',  @cr_lgammaf,  @pcr_lgammaf);
  TestUni('logf',     @cr_logf,     @pcr_logf);
  TestUni('log10f',   @cr_log10f,   @pcr_log10f);
  TestUni('log10p1f', @cr_log10p1f, @pcr_log10p1f);
  TestUni('log1pf',   @cr_log1pf,   @pcr_log1pf);
  TestUni('log2f',    @cr_log2f,    @pcr_log2f);
  TestUni('log2p1f',  @cr_log2p1f,  @pcr_log2p1f);
  TestUni('rsqrtf',   @cr_rsqrtf,   @pcr_rsqrtf);
  TestUni('sinf',     @cr_sinf,     @pcr_sinf);
  TestUni('sinhf',    @cr_sinhf,    @pcr_sinhf);
  TestUni('sinpif',   @cr_sinpif,   @pcr_sinpif);
  TestUni('tanf',     @cr_tanf,     @pcr_tanf);
  TestUni('tanhf',    @cr_tanhf,    @pcr_tanhf);
  TestUni('tanpif',   @cr_tanpif,   @pcr_tanpif);
  TestUni('tgammaf',  @cr_tgammaf,  @pcr_tgammaf);

  // Bivariate sampled (5 functions)
  TestBivar('atan2f',    @wrap_atan2_c,    @wrap_atan2_p);
  TestBivar('atan2pif',  @wrap_atan2pi_c,  @wrap_atan2pi_p);
  TestBivar('hypotf',    @wrap_hypot_c,    @wrap_hypot_p);
  TestBivar('powf',      @wrap_pow_c,      @wrap_pow_p);
  TestBivar('compoundf', @wrap_compound_c, @wrap_compound_p);

  // sincosf exhaustive
  TestSinCos;

  WriteLn;
  WriteLn(Format('=== TOTAL: %d PASS, %d FAIL ===', [TotalPass, TotalFail]));
  if TotalFail = 0 then
    WriteLn('OVERALL: PASS')
  else
    WriteLn('OVERALL: FAIL');

  if TotalFail > 0 then
    Halt(1);
end.
