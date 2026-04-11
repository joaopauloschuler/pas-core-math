{$I pascoremath.inc}
unit pascoremath;

interface

uses pascoremathtypes, pascoremathhelperfuncs;

function pcr_rsqrtf(x: Single): Single;
function pcr_tanhf(x: Single): Single;

implementation

// 1.01: rsqrtf — reciprocal square root (1/sqrt(x))
// Port of cr_rsqrtf from CORE-MATH binary32/rsqrtf.c
function pcr_rsqrtf(x: Single): Single;
const
  tb: array[0..1] of Tb32u32 = ((u: $000C1740), (u: $005222E0));
var
  xd: Double;
  ix: Tb32u32;
  m: LongWord;
  e, k: Integer;
  r, dr: Tb32u32;
begin
  xd := x;
  ix.f := x;
  if (ix.u >= LongWord($FF shl 23)) or (ix.u = 0) then
  begin
    // +/-0: return +/-inf
    if (ix.u shl 1) = 0 then
    begin
      Result := 1.0 / x;
      Exit;
    end;
    // Negative finite or negative NaN
    if (ix.u shr 31) <> 0 then
    begin
      ix.u := ix.u and $7FFFFFFF;   // clear sign bit
      if ix.u > LongWord($FF shl 23) then
      begin
        Result := x + x;            // negative NaN -> propagate NaN
        Exit;
      end;
      // Negative finite: invalid operation
      pcr_feraiseexcept_invalid;
      Result := pcr_nanf('');
      Exit;
    end;
    // +inf: rsqrt(+inf) = 0
    if (ix.u shl 9) = 0 then
    begin
      Result := 0.0;
      Exit;
    end;
    // Positive NaN: propagate
    Result := x + x;
    Exit;
  end;
  // Check for three special inputs that need corrected rounding
  m := ix.u shl 8;
  if (ix.u = $002F7E2A) or (m = $BDF8A800) or (m = $55B7BD00) then
  begin
    if ix.u <> $0055B7BD then
    begin
      e := Integer(ix.u shr 23);
      k := 1;
      if ix.u = $002F7E2A then e := -1;
      if m = $55B7BD00 then k := 0;
      r := tb[k];
      e := (512 - e) div 2 - 578;
      r.u := r.u or LongWord(e shl 23);
      dr.u := LongWord((e - 25) shl 23);
      Result := r.f - dr.f;
      Exit;
    end;
  end;
  // General case: (1/x) * sqrt(x) computed in double
  Result := (1.0 / xd) * pcr_sqrt(xd);
end;

// 1.02: tanhf — hyperbolic tangent
// Port of cr_tanhf from CORE-MATH binary32/tanhf.c
function pcr_tanhf(x: Single): Single;
const
  // Numerator polynomial coefficients (hex floats converted to decimal doubles)
  cn: array[0..7] of Double = (
    1.0,
    0.14869591254532963,      // 0x1.30877b8b72d33p-3
    0.00551287098907202,      // 0x1.694aa09ae9e5ep-8
    7.653349704714027e-05,    // 0x1.4101377abb729p-14
    4.4724281332217524e-07,   // 0x1.e0392b1db0018p-22
    1.0666590627970085e-09,   // 0x1.2533756e546f7p-30
    8.352093632538344e-13,    // 0x1.d62e5abe6ae8ap-41
    9.376645859884988e-17     // 0x1.b06be534182dep-54
  );
  // Denominator polynomial coefficients
  cd: array[0..7] of Double = (
    1.0,
    0.4820292458786627,       // 0x1.ed99131b0ebeap-2
    0.03285595294862704,      // 0x1.0d27ed6c95a69p-5
    0.0007262056643542124,    // 0x1.7cbdaca0e9fccp-11
    6.510296665448557e-06,    // 0x1.b4e60b892578ep-18
    2.4619801106746077e-08,   // 0x1.a6f707c5c71abp-26
    3.5204157099784045e-11,   // 0x1.35a8b6e2cd94cp-35
    1.2726168760182741e-14    // 0x1.ca8230677aa01p-47
  );
  ir: array[0..1] of Single = (1.0, -1.0);
  // -0x1.555556p-2f = -1/3 as Single
  c_neg_third: Single = -0.3333333432674408;
  // 0x1p-25f = 2^-25 as Single
  c_two_neg25: Single = 2.9802322387695312e-08;
var
  z: Double;
  t: Tb32u32;
  ux: UInt32;
  e: Integer;
  x2: Single;
  z2, z4, z8: Double;
  n0, n2, n4, n6: Double;
  d0, d2, d4, d6: Double;
  r: Double;
begin
  z := x;
  t.f := x;
  ux := t.u;
  e := Integer((ux shr 23) and $FF);

  // Inf or NaN
  if e = $FF then
  begin
    if (ux shl 9) <> 0 then
    begin
      Result := x + x;   // NaN -> propagate
      Exit;
    end;
    Result := ir[ux shr 31];   // +-inf -> +-1
    Exit;
  end;

  // |x| < 2^-13: small argument approximations
  if e < 115 then
  begin
    if e < 102 then   // |x| < 2^-26
    begin
      if (ux shl 1) = 0 then
      begin
        Result := x;   // +/-0
        Exit;
      end;
      // tanh(x) ~ x - x*|x|  (keeps correct sign for underflow)
      Result := pcr_fmaf(-x, pcr_fabsf(x), x);
      Exit;
    end;
    // tanh(x) ~ x - x^3/3
    x2 := x * x;
    Result := pcr_fmaf(x, c_neg_third * x2, x);
    Exit;
  end;

  // |x| large enough that tanh(x) rounds to +-1 (with a tiny correction)
  if (ux shl 1) > LongWord($82205966) then   // 0x41102CB3u << 1
  begin
    Result := pcr_copysignf(1.0, x) - pcr_copysignf(c_two_neg25, x);
    Exit;
  end;

  // General case: rational minimax approximation in double
  z2 := z * z;
  z4 := z2 * z2;
  z8 := z4 * z4;

  n0 := cn[0] + z2 * cn[1];
  n2 := cn[2] + z2 * cn[3];
  n4 := cn[4] + z2 * cn[5];
  n6 := cn[6] + z2 * cn[7];
  n0 := n0 + z4 * n2;
  n4 := n4 + z4 * n6;
  n0 := n0 + z8 * n4;

  d0 := cd[0] + z2 * cd[1];
  d2 := cd[2] + z2 * cd[3];
  d4 := cd[4] + z2 * cd[5];
  d6 := cd[6] + z2 * cd[7];
  d0 := d0 + z4 * d2;
  d4 := d4 + z4 * d6;
  d0 := d0 + z8 * d4;

  r := z * n0 / d0;
  Result := r;
end;

end.
