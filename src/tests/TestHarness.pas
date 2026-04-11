{$I ../pascoremath.inc}
program TestHarness;

uses
  pascoremathtypes, pascoremath, ccoremath, SysUtils, Math;

type
  TUniFuncC   = function(x: Single): Single; cdecl;
  TUniFuncP   = function(x: Single): Single;
  TBivarFuncC = function(x, y: Single): Single; cdecl;
  TBivarFuncP = function(x, y: Single): Single;

var
  TotalPass, TotalFail: Integer;

procedure ReportResult(const FuncName: string; Tested, Mismatches: Int64);
begin
  if Mismatches = 0 then
  begin
    WriteLn(Format('%-16s  tested=%10d  mismatches=0  PASS', [FuncName, Tested]));
    Inc(TotalPass);
  end
  else
  begin
    WriteLn(Format('%-16s  tested=%10d  mismatches=%d  FAIL', [FuncName, Tested, Mismatches]));
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
  v: Tb32u32;
  cr, pr: Single;
  mismatches: Int64;
begin
  mismatches := 0;
  for u := 0 to High(Cardinal) do
  begin
    v.u := u;
    cr := pfC(v.f);
    pr := pfP(v.f);
    if not BitsMatch(cr, pr) then
      Inc(mismatches);
  end;
  ReportResult(name, Int64(High(Cardinal)) + 1, mismatches);
end;

// ---- bivariate sampled tests ----

procedure TestBivar(const name: string; pfC: TBivarFuncC; pfP: TBivarFuncP);
const
  SAMPLES = 10000000;
  STRIDE  = High(Cardinal) div SAMPLES;
var
  i: Integer;
  ux, uy: Cardinal;
  vx, vy: Tb32u32;
  cr, pr: Single;
  mismatches: Int64;
begin
  mismatches := 0;
  ux := 0;
  uy := High(Cardinal) div 3;
  for i := 1 to SAMPLES do
  begin
    vx.u := ux;
    vy.u := uy;
    cr := pfC(vx.f, vy.f);
    pr := pfP(vx.f, vy.f);
    if not BitsMatch(cr, pr) then
      Inc(mismatches);
    Inc(ux, STRIDE);
    Inc(uy, STRIDE + 1);
  end;
  ReportResult(name, SAMPLES, mismatches);
end;

// ---- sincosf exhaustive test ----

procedure TestSinCos;
var
  u: Cardinal;
  v: Tb32u32;
  cs, cc, ps, pc: Single;
  mismatches: Int64;
begin
  mismatches := 0;
  for u := 0 to High(Cardinal) do
  begin
    v.u := u;
    cr_sincosf(v.f, @cs, @cc);
    pcr_sincosf(v.f, ps, pc);
    if not BitsMatch(cs, ps) then
      Inc(mismatches)
    else if not BitsMatch(cc, pc) then
      Inc(mismatches);
  end;
  ReportResult('sincosf', Int64(High(Cardinal)) + 1, mismatches);
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

begin
  // Mask all FP exceptions: we iterate over all bit patterns including NaN/Inf
  SetExceptionMask([exInvalidOp, exDenormalized, exZeroDivide,
                    exOverflow, exUnderflow, exPrecision]);

  TotalPass := 0;
  TotalFail := 0;

  WriteLn('=== TestHarness: comparing Pascal (pcr_*) vs C (cr_*) ===');
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
