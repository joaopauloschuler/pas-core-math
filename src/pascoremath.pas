{$I pascoremath.inc}
unit pascoremath;

interface

uses Math, pascoremathtypes, pascoremathhelperfuncs;

function pcr_rsqrtf(x: Single): Single;
function pcr_tanhf(x: Single): Single;
function pcr_atanpif(x: Single): Single;
function pcr_cospif(x: Single): Single;
function pcr_acosf(x: Single): Single;
function pcr_cbrtf(x: Single): Single;
function pcr_sinpif(x: Single): Single;
function pcr_atanf(x: Single): Single;
function pcr_asinf(x: Single): Single;
function pcr_acospif(x: Single): Single;
function pcr_log2f(x: Single): Single;
function pcr_asinpif(x: Single): Single;
function pcr_tanpif(x: Single): Single;
function pcr_coshf(x: Single): Single;
function pcr_logf(x: Single): Single;
function pcr_exp2f(x: Single): Single;
function pcr_log1pf(x: Single): Single;
function pcr_exp2m1f(x: Single): Single;
function pcr_expm1f(x: Single): Single;
function pcr_exp10f(x: Single): Single;
function pcr_log10f(x: Single): Single;
function pcr_erfcf(x: Single): Single;
function pcr_log2p1f(x: Single): Single;
function pcr_erff(x: Single): Single;
function pcr_sinhf(x: Single): Single;
function pcr_expf(x: Single): Single;
function pcr_atanhf(x: Single): Single;
function pcr_exp10m1f(x: Single): Single;
function pcr_log10p1f(x: Single): Single;
function pcr_asinhf(x: Single): Single;
function pcr_acoshf(x: Single): Single;
function pcr_tgammaf(x: Single): Single;
function pcr_lgammaf(x: Single): Single;
function pcr_hypotf(x, y: Single): Single;
function pcr_atan2f(y, x: Single): Single;
function pcr_atan2pif(y, x: Single): Single;
function pcr_powf(x0, y0: Single): Single;
function pcr_compoundf(x, n: Single): Single;
function pcr_sinf(x: Single): Single;
function pcr_cosf(x: Single): Single;
procedure pcr_sincosf(x: Single; out s, c: Single);
function pcr_tanf(x: Single): Single;

implementation

// Arithmetic right shift helpers (FPC shr is always logical)
function sar_i32(x: Integer; n: Integer): Integer; inline;
begin
  if x >= 0 then Result := x shr n
  else Result := not(not(x) shr n);
end;

function sar_i64(x: Int64; n: Integer): Int64; inline;
begin
  if x >= 0 then Result := x shr n
  else Result := not(not(x) shr n);
end;

// Shared polynomial evaluator degree-12 (used by acosf, asinf)
function pcr_poly12(z: Double; const c: array of Double): Double;
var z2, z4, c0, c2, c4, c6, c8, c10: Double;
begin
  z2 := z * z; z4 := z2 * z2;
  c0  := c[0]  + z * c[1];
  c2  := c[2]  + z * c[3];
  c4  := c[4]  + z * c[5];
  c6  := c[6]  + z * c[7];
  c8  := c[8]  + z * c[9];
  c10 := c[10] + z * c[11];
  c0 := c0 + c2 * z2;
  c4 := c4 + c6 * z2;
  c8 := c8 + z2 * c10;
  c0 := c0 + z4 * (c4 + z4 * c8);
  Result := c0;
end;

// Shared S[] table (sin(i*pi/64) for i=0..127) used by cospif and sinpif
const
  S_TABLE: array[0..127] of Double = (
    0.0, 0.049067674327418015, 0.0980171403295606, 0.14673047445536175,
    0.19509032201612828, 0.2429801799032639, 0.2902846772544624, 0.33688985339222005,
    0.3826834323650898, 0.4275550934302821, 0.47139673682599764, 0.5141027441932218,
    0.5555702330196022, 0.5956993044924334, 0.6343932841636455, 0.6715589548470184,
    0.7071067811865476, 0.7409511253549591, 0.773010453362737, 0.8032075314806449,
    0.8314696123025452, 0.8577286100002721, 0.881921264348355, 0.9039892931234433,
    0.9238795325112867, 0.9415440651830208, 0.9569403357322088, 0.970031253194544,
    0.9807852804032304, 0.989176509964781, 0.9951847266721969, 0.9987954562051724,
    1.0, 0.9987954562051724, 0.9951847266721969, 0.989176509964781,
    0.9807852804032304, 0.970031253194544, 0.9569403357322088, 0.9415440651830208,
    0.9238795325112867, 0.9039892931234433, 0.881921264348355, 0.8577286100002721,
    0.8314696123025452, 0.8032075314806449, 0.773010453362737, 0.7409511253549591,
    0.7071067811865476, 0.6715589548470184, 0.6343932841636455, 0.5956993044924334,
    0.5555702330196022, 0.5141027441932218, 0.47139673682599764, 0.4275550934302821,
    0.3826834323650898, 0.33688985339222005, 0.2902846772544624, 0.2429801799032639,
    0.19509032201612828, 0.14673047445536175, 0.0980171403295606, 0.049067674327418015,
    0.0, -0.049067674327418015, -0.0980171403295606, -0.14673047445536175,
    -0.19509032201612828, -0.2429801799032639, -0.2902846772544624, -0.33688985339222005,
    -0.3826834323650898, -0.4275550934302821, -0.47139673682599764, -0.5141027441932218,
    -0.5555702330196022, -0.5956993044924334, -0.6343932841636455, -0.6715589548470184,
    -0.7071067811865476, -0.7409511253549591, -0.773010453362737, -0.8032075314806449,
    -0.8314696123025452, -0.8577286100002721, -0.881921264348355, -0.9039892931234433,
    -0.9238795325112867, -0.9415440651830208, -0.9569403357322088, -0.970031253194544,
    -0.9807852804032304, -0.989176509964781, -0.9951847266721969, -0.9987954562051724,
    -1.0, -0.9987954562051724, -0.9951847266721969, -0.989176509964781,
    -0.9807852804032304, -0.970031253194544, -0.9569403357322088, -0.9415440651830208,
    -0.9238795325112867, -0.9039892931234433, -0.881921264348355, -0.8577286100002721,
    -0.8314696123025452, -0.8032075314806449, -0.773010453362737, -0.7409511253549591,
    -0.7071067811865476, -0.6715589548470184, -0.6343932841636455, -0.5956993044924334,
    -0.5555702330196022, -0.5141027441932218, -0.47139673682599764, -0.4275550934302821,
    -0.3826834323650898, -0.33688985339222005, -0.2902846772544624, -0.2429801799032639,
    -0.19509032201612828, -0.14673047445536175, -0.0980171403295606, -0.049067674327418015
  );

// ── 1.01 rsqrtf ──────────────────────────────────────────────────────────────
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
    if (ix.u shl 1) = 0 then begin Result := 1.0 / x; Exit; end;
    if (ix.u shr 31) <> 0 then
    begin
      ix.u := ix.u and $7FFFFFFF;
      if ix.u > LongWord($FF shl 23) then begin Result := x + x; Exit; end;
      pcr_feraiseexcept_invalid;
      Result := pcr_nanf('');
      Exit;
    end;
    if (ix.u shl 9) = 0 then begin Result := 0.0; Exit; end;
    Result := x + x;
    Exit;
  end;
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
  Result := (1.0 / xd) * pcr_sqrt(xd);
end;

// ── 1.02 tanhf ───────────────────────────────────────────────────────────────
function pcr_tanhf(x: Single): Single;
const
  cn: array[0..7] of Double = (
    1.0, 0.14869591254532963, 0.00551287098907202, 7.653349704714027e-05,
    4.4724281332217524e-07, 1.0666590627970085e-09, 8.352093632538344e-13, 9.376645859884988e-17);
  cd: array[0..7] of Double = (
    1.0, 0.4820292458786627, 0.03285595294862704, 0.0007262056643542124,
    6.510296665448557e-06, 2.4619801106746077e-08, 3.5204157099784045e-11, 1.2726168760182741e-14);
  ir: array[0..1] of Single = (1.0, -1.0);
  c_neg_third: Single = -0.3333333432674408;
  c_two_neg25: Single = 2.9802322387695312e-08;
var
  z: Double;
  t: Tb32u32;
  ux: LongWord;
  e: Integer;
  x2: Single;
  z2, z4, z8: Double;
  n0, n2, n4, n6: Double;
  d0, d2, d4, d6: Double;
  r: Double;
begin
  z := x; t.f := x; ux := t.u;
  e := Integer((ux shr 23) and $FF);
  if e = $FF then
  begin
    if (ux shl 9) <> 0 then begin Result := x + x; Exit; end;
    Result := ir[ux shr 31]; Exit;
  end;
  if e < 115 then
  begin
    if e < 102 then
    begin
      if (ux shl 1) = 0 then begin Result := x; Exit; end;
      Result := pcr_fmaf(-x, pcr_fabsf(x), x); Exit;
    end;
    x2 := x * x;
    Result := pcr_fmaf(x, c_neg_third * x2, x); Exit;
  end;
  if (ux shl 1) > LongWord($82205966) then
  begin
    Result := pcr_copysignf(1.0, x) - pcr_copysignf(c_two_neg25, x); Exit;
  end;
  z2 := z * z; z4 := z2 * z2; z8 := z4 * z4;
  n0 := cn[0] + z2*cn[1]; n2 := cn[2] + z2*cn[3];
  n4 := cn[4] + z2*cn[5]; n6 := cn[6] + z2*cn[7];
  n0 := n0 + z4*n2; n4 := n4 + z4*n6; n0 := n0 + z8*n4;
  d0 := cd[0] + z2*cd[1]; d2 := cd[2] + z2*cd[3];
  d4 := cd[4] + z2*cd[5]; d6 := cd[6] + z2*cd[7];
  d0 := d0 + z4*d2; d4 := d4 + z4*d6; d0 := d0 + z8*d4;
  r := z * n0 / d0;
  Result := r;
end;

// ── 2.01 atanpif ─────────────────────────────────────────────────────────────
function pcr_atanpif(x: Single): Single;
const
  cn: array[0..5] of Double = (
    0.31830988618379064, 0.7250620755086127, 0.5797844040060893,
    0.193473170705847, 0.02469825010811925, 0.0008063015432615248);
  cd: array[0..6] of Double = (
    1.0, 2.6111830231477096, 2.4918407653440666, 1.0590480183430666,
    0.19415473041607811, 0.012196596718179518, 0.00011321825378267113);
var
  t: Tb32u32;
  e: Integer;
  gt: Boolean;
  f: Single;
  z, z2, z4, z8: Double;
  cn0, cn2, cn4: Double;
  cd0, cd2, cd4, cd6: Double;
  r, sx: Double;
  ax: LongWord;
begin
  t.f := x;
  e := Integer((t.u shr 23) and $FF);
  gt := e >= 127;
  if e > 127 + 24 then  // |x| >= 2^25
  begin
    f := pcr_copysignf(0.5, x);
    if e = $FF then
    begin
      if (t.u shl 9) <> 0 then begin Result := x + x; Exit; end;
      Result := f; Exit;
    end;
    // |x| >= 0x1.45f306p+124
    if pcr_fabsf(x) >= Single(2.7078809278823703e+37) then
      Result := f - pcr_copysignf(Single(1.4901161193847656e-08), x)
    else
      Result := f - Single(0.31830987334251404) / x;
    Exit;
  end;
  z := x;
  if e < 127 - 13 then  // |x| < 2^-13
  begin
    sx := z * 0.3183098861837907;  // 0x1.45f306dc9c883p-2
    if e < 127 - 25 then  // |x| < 2^-25
    begin
      Result := sx; Exit;
    end;
    Result := sx - (0.3333333333333333 * sx) * (z * z); Exit;
  end;
  ax := t.u and $7FFFFFFF;
  if ax = $3FA267DD then
  begin
    Result := pcr_copysignf(Single(0.2875366806983948), x) - pcr_copysignf(Single(2.7755575615628914e-17), x);
    Exit;
  end;
  if ax = $3F693531 then
  begin
    Result := pcr_copysignf(Single(0.23518063127994537), x) + pcr_copysignf(Single(3.725290298461914e-09), x);
    Exit;
  end;
  if ax = $3F800000 then
  begin
    Result := pcr_copysignf(0.25, x); Exit;
  end;
  if gt then z := 1.0 / z;
  z2 := z*z; z4 := z2*z2; z8 := z4*z4;
  cn0 := cn[0] + z2*cn[1];
  cn2 := cn[2] + z2*cn[3];
  cn4 := cn[4] + z2*cn[5];
  cn0 := cn0 + z4*cn2;
  cn0 := cn0 + z8*cn4;
  cn0 := cn0 * z;
  cd0 := cd[0] + z2*cd[1];
  cd2 := cd[2] + z2*cd[3];
  cd4 := cd[4] + z2*cd[5];
  cd6 := cd[6];
  cd0 := cd0 + z4*cd2;
  cd4 := cd4 + z4*cd6;
  cd0 := cd0 + z8*cd4;
  r := cn0 / cd0;
  if gt then r := pcr_copysign(0.5, z) - r;
  Result := r;
end;

// ── 2.02 cospif ──────────────────────────────────────────────────────────────
function pcr_cospif(x: Single): Single;
const
  sn: array[0..2] of Double = (
    1.142904749427467e-11, -2.488163196168101e-34, 1.625023320396236e-57);
  cn: array[0..2] of Double = (
    -6.531156331319305e-23, 7.109333835435933e-46, -3.0954114513195225e-69);
var
  ix: Tb32u32;
  e, m_int, s, p: Integer;
  m_uw: LongWord;
  k: Integer;
  iq: LongWord;
  is_idx, ic_idx: LongWord;
  ts, tc, z, z2, fs, fc, r: Double;
  ax: LongWord;
begin
  ix.f := x;
  e := Integer((ix.u shr 23) and $FF);
  if e = $FF then
  begin
    if (ix.u shl 9) = 0 then
    begin
      pcr_feraiseexcept_invalid;
      Result := pcr_nanf('');
    end else
      Result := x + x;
    Exit;
  end;
  m_uw := (ix.u and $7FFFFF) or $800000;
  m_int := Integer(m_uw);
  s := 143 - e;
  p := e - 112;
  if p < 0 then  // |x| < 2^-15
  begin
    ax := ix.u and $7FFFFFFF;
    if ax >= $19F030 then
      Result := pcr_fmaf(Single(-4.934802055358887) * x, x, 1.0)
    else
      Result := pcr_fmaf(-x, x, 1.0);
    Exit;
  end;
  if p > 31 then
  begin
    if p > 63 then begin Result := 1.0; Exit; end;
    iq := LongWord(Integer(LongWord(m_int) shl (p - 32)));
    Result := S_TABLE[(iq + 32) and 127];
    Exit;
  end;
  k := Integer(LongWord(m_int) shl p);
  if k = 0 then
  begin
    iq := LongWord(m_int) shr (32 - p);
    Result := S_TABLE[(iq + 32) and 127];
    Exit;
  end;
  z := k; z2 := z * z;
  fs := sn[0] + z2 * (sn[1] + z2 * sn[2]);
  fc := cn[0] + z2 * (cn[1] + z2 * cn[2]);
  iq := LongWord(m_int) shr s; iq := (iq + 1) shr 1;
  is_idx := iq and 127;
  ic_idx := (iq + 32) and 127;
  ts := S_TABLE[ic_idx];
  tc := S_TABLE[is_idx];
  r := ts + (ts * z2) * fc - (tc * z) * fs;
  Result := r;
end;

// ── 2.03 acosf ───────────────────────────────────────────────────────────────
function pcr_acosf(x: Single): Single;
const
  pi2: Double = 1.5707963267948966;
  b: array[0..15] of Double = (
    0.9999999997220561, 0.1666675305523315, 0.07491953938381704, 0.047534405138862854,
    -0.024905344107261872, 0.6698889818036169, -5.003757071019054, 27.02642690834356,
    -103.66551324982036, 288.04495822181497, -580.9121849063603, 842.6925540871983,
    -857.2868238883075, 581.0567760763246, -235.92908248702702, 43.51567221246845);
  c1: array[0..11] of Double = (
    0.1666666666666473, 0.07500000000425495, 0.044642856775806136, 0.030381960865898193,
    0.022371723076598973, 0.01736016508415668, 0.01388117521087077, 0.012193412697105537,
    0.0064317722535114155, 0.019772599269663224, -0.016582844751635805, 0.03214361520381252);
  c2: array[0..11] of Double = (
    1.4142135623730947, 0.11785113019794026, 0.026516504277464867, 0.007891817376506467,
    0.0026853981502991025, 0.000988848836905083, 0.00038253952347123667, 0.00015842231966484147,
    5.141249514992934e-05, 5.100236375743145e-05, -1.66352623873716e-05, 2.1931983490736225e-05);
var
  t: Tb32u32;
  xs: Double;
  ax: LongWord;
  z, z2, z4, z8, z16, r: Double;
  ub_s, lb_s: Single;
  bx, z_sq: Double;
  s_val: Double;
begin
  t.f := x;
  xs := x;
  ax := t.u shl 1;
  if ax >= LongWord($7F shl 24) then
  begin
    // as_special inline
    if t.u = (LongWord($7F) shl 23) then begin Result := 0.0; Exit; end;  // x=1
    if t.u = (LongWord($17F) shl 23) then begin Result := Single(1.5707963705062866) + Single(-4.371138828673793e-08); Exit; end;  // x=-1
    if ax > (LongWord($FF) shl 24) then begin Result := x + x; Exit; end;  // nan
    pcr_feraiseexcept_invalid;
    Result := pcr_nanf('');
    Exit;
  end;
  if ax < $7EC2A1DC then  // |x| < 0x1.c2a1dcp-1
  begin
    if ax < $40000000 then  // |x| < 2^-63 (spurious underflow guard)
    begin
      Result := Single(1.5707963705062866) + Single(-4.371138828673793e-08);
      Exit;
    end;
    z := xs; z2 := z*z; z4 := z2*z2; z8 := z4*z4; z16 := z8*z8;
    r := z * ((((b[0] + z2*b[1]) + z4*(b[2] + z2*b[3])) + z8*((b[4] + z2*b[5]) + z4*(b[6] + z2*b[7]))) +
              z16*(((b[8] + z2*b[9]) + z4*(b[10] + z2*b[11])) + z8*((b[12] + z2*b[13]) + z4*(b[14] + z2*b[15]))));
    ub_s := Single(1.5707963270725467 - r);
    lb_s := Single(1.5707963265172467 - r);
    if ub_s = lb_s then begin Result := ub_s; Exit; end;
  end;
  // accurate path
  if ax < (LongWord($7E) shl 24) then
  begin
    if t.u = $328885A3 then begin Result := Single(1.5707963705062866) + Single(2.9802322387695312e-08); Exit; end;
    if t.u = $39826222 then begin Result := Single(1.570796012878418) + Single(2.9802322387695312e-08); Exit; end;
    z_sq := xs * xs;
    r := (pi2 - xs) - (xs * z_sq) * pcr_poly12(z_sq, c1);
  end else
  begin
    bx := pcr_fabs(xs);
    z := 1.0 - bx;
    s_val := pcr_copysign(pcr_sqrt(z), xs);
    if (t.u shr 31) = 0 then
      r := s_val * pcr_poly12(z, c2)
    else
      r := 3.141592653589793 + s_val * pcr_poly12(z, c2);
  end;
  Result := r;
end;

// ── 2.04 cbrtf ───────────────────────────────────────────────────────────────
function pcr_cbrtf(x: Single): Single;
const
  escale_u: array[0..2] of UInt64 = (
    UInt64($3FF0000000000000),  // 1.0
    UInt64($3FF428A2F98D728B),  // 2^(1/3)
    UInt64($3FF965FEA53D6E3D)); // 2^(2/3)
  c: array[0..7] of Double = (
    0.5685564078059381, 0.7024960185339382, -0.39381000363475277, 0.21397507019181075,
    -0.08593966563932363, 0.023134567971640832, -0.003702862366439682, 0.00026571366637555694);
var
  t: Tb32u32;
  u: LongWord;
  au: LongWord;
  sgn: LongWord;
  e: LongWord;
  nz: Integer;
  mant: LongWord;
  cvt1, cvt2: Tb64u64;
  et, it_: LongWord;
  isc: UInt64;
  z, r0, z2, z4, f, r: Double;
  ub, lb: Single;
  u0: Double;
  h: Double;
  m0, m1: Int64;
begin
  t.f := x;
  u := t.u;
  au := u shl 1;
  sgn := u shr 31;
  e := au shr 24;
  if (au < LongWord(1 shl 24)) or (au >= LongWord($FF shl 24)) then
  begin
    if au >= LongWord($FF shl 24) then begin Result := x + x; Exit; end;  // inf/nan
    if au = 0 then begin Result := x; Exit; end;  // +-0
    // subnormal
    nz := 24 - Integer(pcr_bsr32(au));  // = clz(au) - 7
    au := au shl nz;
    if nz > 1 then e := e - LongWord(nz - 1) else e := e;
  end;
  mant := au and $FFFFFF;
  cvt1.u := (UInt64(mant) shl 28) or (UInt64($3FF) shl 52);
  e := e + 899;
  et := e div 3;
  it_ := e mod 3;
  isc := escale_u[it_];
  isc := isc + UInt64(Int64(Int64(et) - 342) shl 52);
  isc := isc or (UInt64(sgn) shl 63);
  cvt2.u := isc;
  z := cvt1.f;
  r0 := -0.024975246527242426 / z;  // -0x1.9931c6c2d19d1p-6 / z
  z2 := z * z; z4 := z2 * z2;
  f := ((c[0] + z*c[1]) + z2*(c[2] + z*c[3])) + z4*((c[4] + z*c[5]) + z2*(c[6] + z*c[7])) + r0;
  r := f * cvt2.f;
  ub := Single(r);
  lb := Single(r - cvt2.f * 1.4182e-9);
  if ub = lb then
  begin
    cvt2.f := r;
    Result := ub;
    Exit;
  end;
  u0 := -13.34654827009379;  // -0x1.ab16ec65d138fp+3
  h := f*f*f - z;
  f := f - (f * r0 * u0) * h;
  r := f * cvt2.f;
  cvt1.f := r;
  ub := Single(r);
  m0 := Int64(cvt1.u shl 19);
  if m0 < 0 then m1 := Int64(-1) else m1 := Int64(0);
  if (m0 xor m1) < (Int64(1) shl 31) then
  begin
    cvt1.u := (cvt1.u + (UInt64(1) shl 31)) and UInt64($FFFFFFFF00000000);
    ub := Single(cvt1.f);
  end;
  Result := ub;
end;

// ── 2.05 sinpif ──────────────────────────────────────────────────────────────
function pcr_sinpif(x: Single): Single;
const
  sn: array[0..2] of Double = (
    1.142904749427467e-11, -2.488163196168101e-34, 1.625023320396236e-57);
  cn: array[0..2] of Double = (
    -6.531156331319305e-23, 7.109333835435933e-46, -3.0954114513195225e-69);
var
  ix: Tb32u32;
  e, m_int, sgn, s, si: Integer;
  iq: LongWord;
  is_idx, ic_idx: LongWord;
  ts, tc, z, z2, fs, fc, r: Double;
  k: Integer;
begin
  ix.f := x;
  e := Integer((ix.u shr 23) and $FF);
  if e = $FF then
  begin
    if (ix.u shl 9) = 0 then
    begin
      pcr_feraiseexcept_invalid;
      Result := pcr_nanf('');
    end else
      Result := x + x;
    Exit;
  end;
  m_int := Integer((ix.u and $7FFFFF) or $800000);
  sgn := -Integer(ix.u shr 31);  // 0 or -1
  m_int := (m_int xor sgn) - sgn;
  s := 143 - e;
  if s < 0 then  // |x| >= 2^17
  begin
    if s < -6 then begin Result := pcr_copysignf(0.0, x); Exit; end;  // |x| >= 2^23
    iq := LongWord(m_int) shl (-s - 1);
    iq := iq and 127;
    if (iq = 0) or (iq = 64) then begin Result := pcr_copysignf(0.0, x); Exit; end;
    Result := S_TABLE[iq];
    Exit;
  end else if s > 30 then  // |x| < 2^-14
  begin
    z := x; z2 := z * z;
    Result := z * (3.141592653589793 + z2 * (-5.16771278004997));
    Exit;
  end;
  si := 25 - s;
  if (si >= 0) and ((LongWord(m_int) shl si) = 0) then
  begin
    Result := pcr_copysignf(0.0, x); Exit;
  end;
  k := Integer(LongWord(m_int) shl (31 - s));
  z := k; z2 := z * z;
  fs := sn[0] + z2 * (sn[1] + z2 * sn[2]);
  fc := cn[0] + z2 * (cn[1] + z2 * cn[2]);
  iq := LongWord(sar_i32(m_int, s)); iq := (iq + 1) shr 1;
  is_idx := iq and 127;
  ic_idx := (iq + 32) and 127;
  ts := S_TABLE[is_idx];
  tc := S_TABLE[ic_idx];
  r := ts + (ts * z2) * fc + (tc * z) * fs;
  Result := r;
end;

// ── 2.06 atanf ───────────────────────────────────────────────────────────────
function pcr_atanf(x: Single): Single;
const
  pi2: Double = 1.5707963267948966;
  cn: array[0..6] of Double = (
    0.33000489885804146, 0.8269936260181494, 0.7536692267812706, 0.3041250206581639,
    0.052585465033265374, 0.0030928116297212196, 2.6680447001914062e-05);
  cd: array[0..6] of Double = (
    0.3300048988580414, 0.9369952589708292, 1.0, 0.4972028591750377,
    0.1155090060414157, 0.0109022453539874, 0.00027322693677761577);
  PI_OVER2_H: Double = 1.5625;
  PI_OVER2_L: Double = 0.008296326794896619;
var
  t: Tb32u32;
  ta: LongWord;
  e: Integer;
  gt: Boolean;
  z, z2, z4, z8: Double;
  cn0, cn2, cn4, cn6: Double;
  cd0, cd2, cd4, cd6: Double;
  r: Double;
begin
  t.f := x;
  e := Integer((t.u shr 23) and $FF);
  gt := e >= 127;
  ta := t.u and $7FFFFFFF;
  if ta >= $4C700518 then  // |x| >= 0x1.e00a3p+25
  begin
    if ta > $7F800000 then begin Result := x + x; Exit; end;  // nan
    Result := pcr_copysign(pi2, x); Exit;  // inf or large
  end;
  if e < 127 - 13 then  // |x| < 2^-13
  begin
    if e < 127 - 25 then  // |x| < 2^-25
    begin
      if (t.u shl 1) = 0 then begin Result := x; Exit; end;
      Result := pcr_fmaf(-x, pcr_fabsf(x), x); Exit;
    end;
    Result := pcr_fmaf(Single(-0.3333333333333333) * x, x * x, x); Exit;
  end;
  z := x;
  if gt then z := 1.0 / z;
  z2 := z*z; z4 := z2*z2; z8 := z4*z4;
  cn0 := cn[0] + z2*cn[1]; cn2 := cn[2] + z2*cn[3];
  cn4 := cn[4] + z2*cn[5]; cn6 := cn[6];
  cn0 := cn0 + z4*cn2; cn4 := cn4 + z4*cn6; cn0 := cn0 + z8*cn4;
  cn0 := cn0 * z;
  cd0 := cd[0] + z2*cd[1]; cd2 := cd[2] + z2*cd[3];
  cd4 := cd[4] + z2*cd[5]; cd6 := cd[6];
  cd0 := cd0 + z4*cd2; cd4 := cd4 + z4*cd6; cd0 := cd0 + z8*cd4;
  r := cn0 / cd0;
  if not gt then begin Result := r; Exit; end;
  r := (pcr_copysign(PI_OVER2_L, z) - r) + pcr_copysign(PI_OVER2_H, z);
  Result := r;
end;

// ── 2.07 asinf ───────────────────────────────────────────────────────────────
function pcr_asinf(x: Single): Single;
const
  pi2: Double = 1.5707963267948966;
  b: array[0..15] of Double = (
    1.000000000000001, 0.16666694674143204, 0.07497112542795417, 0.0458179575336707,
    0.005331008900413985, 0.34410258152367046, -2.680930042099564, 15.541270760972983,
    -63.17329833405016, 184.79515144873312, -390.0198166803775, 589.2790780950768,
    -621.89777643639, 435.8403729646551, -182.48552714860514, 34.63705332873756);
  c1: array[0..11] of Double = (
    0.1666666666666473, 0.07500000000425495, 0.044642856775806136, 0.030381960865898193,
    0.022371723076598973, 0.01736016508415668, 0.01388117521087077, 0.012193412697105537,
    0.0064317722535114155, 0.019772599269663224, -0.016582844751635805, 0.03214361520381252);
  c2: array[0..11] of Double = (
    1.4142135623730947, 0.11785113019794026, 0.026516504277464867, 0.007891817376506467,
    0.0026853981502991025, 0.000988848836905083, 0.00038253952347123667, 0.00015842231966484147,
    5.141249514992934e-05, 5.100236375743145e-05, -1.66352623873716e-05, 2.1931983490736225e-05);
var
  t: Tb32u32;
  xs: Double;
  ax: LongWord;
  z, z2, z4, z8, z16, r, c0: Double;
  ub_s, lb_s: Single;
  bx, s_val: Double;
begin
  t.f := x; xs := x;
  ax := t.u shl 1;
  if ax > LongWord($7F shl 24) then
  begin
    // as_special
    if ax > (LongWord($FF) shl 24) then begin Result := x + x; Exit; end;  // nan
    pcr_feraiseexcept_invalid;
    Result := pcr_nanf('');
    Exit;
  end;
  if ax < $7EC29000 then  // |x| < 0x1.c25p-1 approximately
  begin
    if ax < (115 shl 24) then  // |x| < 2^-12
    begin
      Result := pcr_fmaf(x, Single(2.9802322387695312e-08), x); Exit;
    end;
    z := xs; z2 := z*z; z4 := z2*z2; z8 := z4*z4; z16 := z8*z8;
    r := z * ((((b[0] + z2*b[1]) + z4*(b[2] + z2*b[3])) + z8*((b[4] + z2*b[5]) + z4*(b[6] + z2*b[7]))) +
              z16*(((b[8] + z2*b[9]) + z4*(b[10] + z2*b[11])) + z8*((b[12] + z2*b[13]) + z4*(b[14] + z2*b[15]))));
    ub_s := Single(r);
    lb_s := Single(r - z * 9.015999891115456e-10);
    if ub_s = lb_s then begin Result := ub_s; Exit; end;
  end;
  if ax < (LongWord($7E) shl 24) then
  begin
    z := xs; z2 := z * z;
    c0 := pcr_poly12(z2, c1);
    r := z + (z * z2) * c0;
  end else
  begin
    if ax = $7E55688A then begin Result := pcr_copysignf(Single(0.7299242615699768), x) + pcr_copysignf(Single(1.4901161193847656e-08), x); Exit; end;
    if ax = $7E107434 then begin Result := pcr_copysignf(Single(0.5611220598220825), x) + pcr_copysignf(Single(1.4901161193847656e-08), x); Exit; end;
    bx := pcr_fabs(xs);
    z := 1.0 - bx;
    s_val := pcr_sqrt(z);
    r := pi2 - s_val * pcr_poly12(z, c2);
    r := pcr_copysign(r, xs);
  end;
  Result := r;
end;

// ── 2.08 acospif ─────────────────────────────────────────────────────────────
function pcr_acospif(x: Single): Single;
const
  ch: array[0..15, 0..7] of Double = (
    (0.31830988618379064, 0.05305164769736625, 0.023873241441162645, 0.014210265648377423,
     0.00967069454939692, 0.007127294943345352, 0.0054123088854984595, 0.0054941925043860435),
    (0.49999999998995076, -0.06830988526542055, 0.028345020403629382, -0.016167029449478044,
     0.01070943653063092, -0.007652978565306432, 0.005335533266464206, -0.0026333050557125395),
    (0.49999999952078344, -0.06830986073948833, 0.028344464767354173, -0.016159946510885016,
     0.010654488243861775, -0.007393236218625197, 0.004642228742710111, -0.001826894552556536),
    (0.49999999483623414, -0.06830969111395047, 0.02834181593211837, -0.01613681384471718,
     0.010532431562201716, -0.007004045811776552, 0.00394774332469213, -0.00129186689983247),
    (0.499999972035478, -0.0683090629285235, 0.028334370301561714, -0.016087593321521987,
     0.010336417923139498, -0.0065337733449434485, 0.003318343434143512, -0.0009293620753800756),
    (0.4999998982483317, -0.06830742580615927, 0.02831876408813582, -0.01600473147816068,
     0.010071756089994829, -0.006025242185043103, 0.002774074164502572, -0.0006790543452815604),
    (0.49999971401979315, -0.06830400771647019, 0.028291536785706447, -0.01588402406007313,
     0.00975009246422934, -0.005509998564443676, 0.002314724387438669, -0.0005032260093581129),
    (0.4999993288811753, -0.06829787058456117, 0.02824956969145477, -0.015724378671429857,
     0.009385225179622859, -0.005008988611530038, 0.001932016650266531, -0.0003777698015384437),
    (0.49999862054189037, -0.06828798201035122, 0.028190346834476125, -0.015527129485538686,
     0.008990643644453681, -0.004534905112928049, 0.0016152454388246105, -0.00028696596361868227),
    (0.49999743715562245, -0.0682732854890257, 0.028112062536902023, -0.01529527652877646,
     0.008578305783998557, -0.00409455631822447, 0.001353776823971926, -0.00022037480993851409),
    (0.49999560172579033, -0.06825275950876084, 0.02801362144639405, -0.01503281810632914,
     0.008158178849581539, -0.0036907838870853128, 0.0011380495632999778, -0.00017094592147442263),
    (0.49999291770384074, -0.06822546192844604, 0.027894574076084967, -0.014744230356412454,
     0.007738204747834519, -0.003323878775111286, 0.0009598741584255129, -0.00013384370789719762),
    (0.49998917500047757, -0.06819055970415733, 0.027755020451573003, -0.014434092971678629,
     0.0073244756411289, -0.0029925745075133816, 0.00081241793688122, -0.00010570399192569352),
    (0.4999841558466798, -0.06814734587118984, 0.027595503682921574, -0.014106838410046394,
     0.006921494448252957, -0.0026947204203676355, 0.0006900639292580665, -8.415513850335474e-05),
    (0.4999776401448277, -0.06809524636661032, 0.0274169066792878, -0.013766596827913393,
     0.006532451900726738, -0.0024277253159630096, 0.0005882324316513122, -6.750456758710982e-05),
    (0.4999694101123286, -0.06803381928269747, 0.027220359059896233, -0.013417110954495652,
     0.006159485915362247, -0.002188841583759408, 0.0005032052395692806, -5.453042765077486e-05)
  );
var
  t: Tb32u32;
  ax_f: Single;
  az, z: Double;
  e, s, i: Integer;
  z2, z4, c0_v, c2_v, c4_v, c6_v, f, r: Double;
  o_val: Double;
begin
  t.f := x;
  ax_f := pcr_fabsf(x);
  az := ax_f;
  z := x;
  e := Integer((t.u shr 23) and $FF);
  if e >= 127 then
  begin
    if x = 1.0 then begin Result := 0.0; Exit; end;
    if x = -1.0 then begin Result := 1.0; Exit; end;
    if (e = $FF) and ((t.u shl 9) <> 0) then begin Result := x + x; Exit; end;
    pcr_feraiseexcept_invalid;
    Result := pcr_nanf('');
    Exit;
  end;
  s := 146 - e;
  i := 0;
  if s < 32 then
    i := Integer(((t.u and $7FFFFF) or $800000) shr s);
  z2 := z * z; z4 := z2 * z2;
  if i = 0 then
  begin
    c0_v := ch[0, 0] + z2 * ch[0, 1];
    c2_v := ch[0, 2] + z2 * ch[0, 3];
    c4_v := ch[0, 4] + z2 * ch[0, 5];
    c6_v := ch[0, 6] + z2 * ch[0, 7];
    c0_v := c0_v + c2_v * z4;
    c4_v := c4_v + c6_v * z4;
    c0_v := pcr_fma(c4_v * z4, z4, c0_v);
    Result := 0.5 - z * c0_v;
  end else
  begin
    f := pcr_sqrt(1.0 - az);
    c0_v := ch[i, 0] + az * ch[i, 1];
    c2_v := ch[i, 2] + az * ch[i, 3];
    c4_v := ch[i, 4] + az * ch[i, 5];
    c6_v := ch[i, 6] + az * ch[i, 7];
    c0_v := c0_v + c2_v * z2;
    c4_v := c4_v + c6_v * z2;
    c0_v := c0_v + c4_v * z4;
    if (t.u shr 31) = 0 then o_val := 0.0 else o_val := 1.0;
    r := o_val + c0_v * pcr_copysign(f, x);
    Result := r;
  end;
end;

// ── 2.09 log2f ───────────────────────────────────────────────────────────────
function pcr_log2f(x: Single): Single;
const
  ix_arr: array[0..64] of Double = (
    1.0, 0.9846153855323792, 0.969696968793869, 0.9552238807082176,
    0.9411764703691006, 0.9275362323969603, 0.9142857138067484, 0.9014084506779909,
    0.8888888880610466, 0.876712329685688, 0.8648648653179407, 0.8533333335071802,
    0.8421052638441324, 0.8311688303947449, 0.8205128200352192, 0.8101265821605921,
    0.8000000007450581, 0.7901234570890665, 0.7804878056049347, 0.7710843365639448,
    0.7619047611951828, 0.7529411762952805, 0.744186045601964, 0.7356321830302477,
    0.7272727265954018, 0.7191011235117912, 0.7111111115664244, 0.7032967042177916,
    0.6956521738320589, 0.6881720423698425, 0.6808510646224022, 0.673684211447835,
    0.666666666045785, 0.6597938146442175, 0.6530612241476774, 0.6464646458625793,
    0.6400000005960464, 0.6336633656173944, 0.627450980246067, 0.6213592234998941,
    0.615384615957737, 0.6095238104462624, 0.6037735845893621, 0.5981308408081532,
    0.5925925932824612, 0.5871559642255306, 0.5818181820213795, 0.5765765774995089,
    0.5714285708963871, 0.5663716811686754, 0.5614035092294216, 0.5565217398107052,
    0.551724137738347, 0.5470085479319096, 0.5423728805035353, 0.5378151256591082,
    0.533333333209157, 0.5289256200194359, 0.5245901644229889, 0.5203252024948597,
    0.5161290317773819, 0.5120000001043081, 0.5079365074634552, 0.5039370078593493,
    0.5);
  lix_arr: array[0..64] of Double = (
    0.0, -0.02236781168484005, -0.0443941207020679, -0.06608919028982063,
    -0.08746284158624303, -0.10852445598039796, -0.1292830177007496, -0.14974711954667003,
    -0.16992500278592682, -0.18982455736845097, -0.20945336487316665, -0.2288186902019652,
    -0.24792751226792284, -0.2667865420385158, -0.28540221970200735, -0.30378074838704267,
    -0.3219280935437479, -0.3398500023387814, -0.35755200327446923, -0.3750394328165031,
    -0.39231742412237475, -0.40939093647360536, -0.4262647564655919, -0.4429434975702343,
    -0.4594316199809117, -0.47573343113434957, -0.49185309540593974, -0.5077946383092384,
    -0.5235619562249647, -0.5391588124516459, -0.5545888499981193, -0.5698556063575141,
    -0.5849625020647706, -0.5999128417252602, -0.6147098448709913, -0.6293566214232241,
    -0.6438561884311103, -0.6582114843893249, -0.6724253423073993, -0.6865005267213509,
    -0.7004397167974777, -0.7142455154827492, -0.7279204553189823, -0.74146698715693,
    -0.7548875004839505, -0.7681843225095769, -0.7813597130208042, -0.7944158640407686,
    -0.8073549234012186, -0.820178963045007, -0.832890012989079, -0.8454900491808812,
    -0.8579809956314275, -0.8703647171481033, -0.8826430516291907, -0.8948177643576423,
    -0.9068905959444221, -0.9188632367707391, -0.9307373362192718, -0.9425145074386375,
    -0.9541963117304897, -0.9657842843681714, -0.9772799248435309, -0.9886846868141538,
    -1.0);
  bcoef: array[0..2] of Double = (
    1.4426950429725995, -0.7213691893530103, 0.4808376657770807);
  ccoef: array[0..5] of Double = (
    1.4426950408889683, -0.7213475204443797, 0.48089834631236367,
    -0.36067376567480197, 0.2885606699465501, -0.24038686869298112);
var
  t: Tb32u32;
  xd: Tb64u64;
  ux: LongWord;
  n: Integer;
  e: Integer;
  m: LongWord;
  j: Integer;
  z, z2, el, f: Double;
  lb_s, ub_s: Single;
  c0: Double;
  neg_inf_t: Tb32u32;
begin
  t.f := x;
  ux := t.u;
  if ux >= $7F800000 then  // special: <=−0, nan, +inf
  begin
    // as_special inline
    if (ux shl 1) = 0 then  // -0
    begin
      pcr_feraiseexcept_divbyzero;
      neg_inf_t.u := $FF800000;
      Result := neg_inf_t.f;
      Exit;
    end;
    if ux = $7F800000 then begin Result := x; Exit; end;  // +inf
    if (ux shl 1) > $FF000000 then begin Result := x + x; Exit; end;  // nan
    pcr_feraiseexcept_invalid;
    Result := pcr_nanf('');
    Exit;
  end;
  if ux < $800000 then  // subnormal
  begin
    if ux = 0 then  // +0
    begin
      pcr_feraiseexcept_divbyzero;
      neg_inf_t.u := $FF800000;
      Result := neg_inf_t.f;
      Exit;
    end;
    n := 23 - Integer(pcr_bsr32(ux));  // = clz(ux) - 8
    ux := ux shl n;
    ux := ux - LongWord(n shl 23);
  end;
  e := sar_i32(Integer(ux), 23) - $7F;
  m := ux and $7FFFFF;
  if m = 0 then begin Result := Single(e); Exit; end;
  j := Integer((m + (1 shl 16)) shr 17);
  xd.u := (UInt64(m) shl 29) or (UInt64($3FF) shl 52);
  z := xd.f * ix_arr[j] - 1.0;
  z2 := z * z;
  el := Double(e) - lix_arr[j];
  f := (el + z * bcoef[0]) + z2 * (bcoef[1] + z * bcoef[2]);
  lb_s := Single(f);
  ub_s := Single(f + 3.256559466535691e-10);
  if lb_s = ub_s then begin Result := lb_s; Exit; end;
  c0 := ccoef[0] + z * ccoef[1];
  c0 := c0 + z2 * ((ccoef[2] + z * ccoef[3]) + z2 * (ccoef[4] + z * ccoef[5]));
  Result := Single(el + z * c0);
end;

// ── 2.10 asinpif ─────────────────────────────────────────────────────────────
function pcr_asinpif(x: Single): Single;
const
  ch: array[0..15, 0..7] of Double = (
    (0.31830988618379064, 0.05305164769736625, 0.023873241441162645, 0.014210265648377423,
     0.00967069454939692, 0.007127294943345352, 0.0054123088854984595, 0.0054941925043860435),
    // row 1 differs from acospif by 1 ULP in slots 0 and 1
    (0.4999999999899508, -0.06830988526542066, 0.028345020403629382, -0.016167029449478044,
     0.01070943653063092, -0.007652978565306432, 0.005335533266464206, -0.0026333050557125395),
    (0.49999999952078344, -0.06830986073948833, 0.028344464767354173, -0.016159946510885016,
     0.010654488243861775, -0.007393236218625197, 0.004642228742710111, -0.001826894552556536),
    (0.49999999483623414, -0.06830969111395047, 0.02834181593211837, -0.01613681384471718,
     0.010532431562201716, -0.007004045811776552, 0.00394774332469213, -0.00129186689983247),
    (0.499999972035478, -0.0683090629285235, 0.028334370301561714, -0.016087593321521987,
     0.010336417923139498, -0.0065337733449434485, 0.003318343434143512, -0.0009293620753800756),
    (0.4999998982483317, -0.06830742580615927, 0.02831876408813582, -0.01600473147816068,
     0.010071756089994829, -0.006025242185043103, 0.002774074164502572, -0.0006790543452815604),
    (0.49999971401979315, -0.06830400771647019, 0.028291536785706447, -0.01588402406007313,
     0.00975009246422934, -0.005509998564443676, 0.002314724387438669, -0.0005032260093581129),
    (0.4999993288811753, -0.06829787058456117, 0.02824956969145477, -0.015724378671429857,
     0.009385225179622859, -0.005008988611530038, 0.001932016650266531, -0.0003777698015384437),
    (0.49999862054189037, -0.06828798201035122, 0.028190346834476125, -0.015527129485538686,
     0.008990643644453681, -0.004534905112928049, 0.0016152454388246105, -0.00028696596361868227),
    (0.49999743715562245, -0.0682732854890257, 0.028112062536902023, -0.01529527652877646,
     0.008578305783998557, -0.00409455631822447, 0.001353776823971926, -0.00022037480993851409),
    (0.49999560172579033, -0.06825275950876084, 0.02801362144639405, -0.01503281810632914,
     0.008158178849581539, -0.0036907838870853128, 0.0011380495632999778, -0.00017094592147442263),
    (0.49999291770384074, -0.06822546192844604, 0.027894574076084967, -0.014744230356412454,
     0.007738204747834519, -0.003323878775111286, 0.0009598741584255129, -0.00013384370789719762),
    (0.49998917500047757, -0.06819055970415733, 0.027755020451573003, -0.014434092971678629,
     0.0073244756411289, -0.0029925745075133816, 0.00081241793688122, -0.00010570399192569352),
    (0.4999841558466798, -0.06814734587118984, 0.027595503682921574, -0.014106838410046394,
     0.006921494448252957, -0.0026947204203676355, 0.0006900639292580665, -8.415513850335474e-05),
    (0.4999776401448277, -0.06809524636661032, 0.0274169066792878, -0.013766596827913393,
     0.006532451900726738, -0.0024277253159630096, 0.0005882324316513122, -6.750456758710982e-05),
    (0.4999694101123286, -0.06803381928269747, 0.027220359059896233, -0.013417110954495652,
     0.006159485915362247, -0.002188841583759408, 0.0005032052395692806, -5.453042765077486e-05)
  );
var
  t: Tb32u32;
  ax_f: Single;
  az, z: Double;
  e, s, i: Integer;
  z2, z4, c0_v, c2_v, c4_v, c6_v, f, r: Double;
begin
  t.f := x;
  ax_f := pcr_fabsf(x);
  az := ax_f;
  z := x;
  e := Integer((t.u shr 23) and $FF);
  if e >= 127 then
  begin
    if ax_f = 1.0 then begin Result := pcr_copysignf(0.5, x); Exit; end;
    if (e = $FF) and ((t.u shl 9) <> 0) then begin Result := x + x; Exit; end;
    pcr_feraiseexcept_invalid;
    Result := pcr_nanf('');
    Exit;
  end;
  s := 146 - e;
  i := 0;
  if s < 32 then
    i := Integer(((t.u and $7FFFFF) or $800000) shr s);
  z2 := z * z; z4 := z2 * z2;
  if i = 0 then
  begin
    c0_v := ch[0, 0] + z2 * ch[0, 1];
    c2_v := ch[0, 2] + z2 * ch[0, 3];
    c4_v := ch[0, 4] + z2 * ch[0, 5];
    c6_v := ch[0, 6] + z2 * ch[0, 7];
    c0_v := c0_v + c2_v * z4;
    c4_v := c4_v + c6_v * z4;
    c0_v := c0_v + c4_v * (z4 * z4);
    Result := z * c0_v;
  end else
  begin
    f := pcr_sqrt(1.0 - az);
    c0_v := ch[i, 0] + az * ch[i, 1];
    c2_v := ch[i, 2] + az * ch[i, 3];
    c4_v := ch[i, 4] + az * ch[i, 5];
    c6_v := ch[i, 6] + az * ch[i, 7];
    c0_v := c0_v + c2_v * z2;
    c4_v := c4_v + c6_v * z2;
    c0_v := c0_v + c4_v * z4;
    r := pcr_fma(-c0_v, pcr_copysign(f, x), pcr_copysign(0.5, x));
    Result := r;
  end;
end;

// ── 2.11 tanpif ──────────────────────────────────────────────────────────────
function pcr_tanpif(x: Single): Single;
const
  cn: array[0..3] of Double = (
    0.7853981633974484, -0.2805387264887832, 0.02201158908691473, -0.00023103959012326923);
  cd: array[0..3] of Double = (
    1.0, -0.6470611340915767, 0.0973140255480054, -0.0032269805489163333);
var
  ix: Tb32u32;
  e_bits: LongWord;
  x4, nx4, dx4, ni, zf: Single;
  k: Integer;
  a: LongWord;
  z, z2, z4, r: Double;
  pos_inf_t, neg_inf_t: Tb32u32;
begin
  ix.f := x;
  e_bits := ix.u and ($FF shl 23);
  if e_bits > (150 shl 23) then  // |x| > 2^23
  begin
    if e_bits = ($FF shl 23) then  // nan or inf
    begin
      if (ix.u shl 9) = 0 then  // inf
      begin
        pcr_feraiseexcept_invalid;
        Result := pcr_nanf('');
      end else
        Result := x + x;
    end else
      Result := pcr_copysignf(0.0, x);
    Exit;
  end;
  x4 := 4.0 * x;
  nx4 := pcr_roundevenf(x4);
  dx4 := x4 - nx4;
  ni := pcr_roundevenf(x);
  zf := x - ni;
  if dx4 = 0.0 then  // 4*x is integer
  begin
    k := Trunc(x4);
    if (k and 1) <> 0 then begin Result := pcr_copysignf(1.0, zf); Exit; end;  // x = 1/4 mod 1/2
    k := k and 6;
    if k = 0 then begin Result := pcr_copysignf(0.0, x); Exit; end;   // x = 0 mod 2
    if k = 4 then begin Result := -pcr_copysignf(0.0, x); Exit; end;  // x = 1 mod 2
    pos_inf_t.u := $7F800000; neg_inf_t.u := $FF800000;
    if k = 2 then begin Result := pos_inf_t.f; Exit; end;  // x = 1/2 mod 2 → +inf
    Result := neg_inf_t.f; Exit;                            // x = -1/2 mod 2 → -inf
  end;
  ix.f := zf;
  a := ix.u and $7FFFFFFF;
  if a = $3E933802 then
  begin
    Result := pcr_copysignf(Single(1.2687946557998657), zf) + pcr_copysignf(Single(2.9802322387695312e-08), zf);
    Exit;
  end;
  if a = $38F26685 then
  begin
    Result := pcr_copysignf(Single(0.000363122730050236), zf) + pcr_copysignf(Single(7.275957614183426e-12), zf);
    Exit;
  end;
  z := zf; z2 := z * z; z4 := z2 * z2;
  r := (z - z*z2) * ((cn[0] + z2*cn[1]) + z4*(cn[2] + z2*cn[3])) /
       (((cd[0] + z2*cd[1]) + z4*(cd[2] + z2*cd[3])) * (0.25 - z2));
  Result := r;
end;

// ── 2.12 coshf ───────────────────────────────────────────────────────────────
function pcr_coshf(x: Single): Single;
const
  c_arr: array[0..3] of Double = (
    1.0, 0.021660849391257477, 0.0002345984913513542, 1.6938658699950235e-06);
  ch_arr: array[0..6] of Double = (
    1.0, 0.02166084939249829, 0.0002345961982022468, 1.6938509724129055e-06,
    9.172562701702629e-09, 3.973729405780548e-11, 1.4345723178374038e-13);
  // tb[k] = 2^(k/32) * 0.5, stored as uint64 bit patterns
  tb_arr: array[0..31] of UInt64 = (
    UInt64($3FE0000000000000), UInt64($3FE059B0D3158574), UInt64($3FE0B5586CF9890F), UInt64($3FE11301D0125B51),
    UInt64($3FE172B83C7D517B), UInt64($3FE1D4873168B9AA), UInt64($3FE2387A6E756238), UInt64($3FE29E9DF51FDEE1),
    UInt64($3FE306FE0A31B715), UInt64($3FE371A7373AA9CB), UInt64($3FE3DEA64C123422), UInt64($3FE44E086061892D),
    UInt64($3FE4BFDAD5362A27), UInt64($3FE5342B569D4F82), UInt64($3FE5AB07DD485429), UInt64($3FE6247EB03A5585),
    UInt64($3FE6A09E667F3BCD), UInt64($3FE71F75E8EC5F74), UInt64($3FE7A11473EB0187), UInt64($3FE82589994CCE13),
    UInt64($3FE8ACE5422AA0DB), UInt64($3FE93737B0CDC5E5), UInt64($3FE9C49182A3F090), UInt64($3FEA5503B23E255D),
    UInt64($3FEAE89F995AD3AD), UInt64($3FEB7F76F2FB5E47), UInt64($3FEC199BDD85529C), UInt64($3FECB720DCEF9069),
    UInt64($3FED5818DCFBA487), UInt64($3FEDFC97337B9B5F), UInt64($3FEEA4AFA2A490DA), UInt64($3FEF50765B6E4540));
  iln2: Double = 46.16624130844683;
  cp_arr: array[0..3] of Double = (
    0.4999999999999984, 0.04166666666748819, 0.0013888887416776143, 2.4812354013894482e-05);
var
  t: Tb32u32;
  z: Double;
  ax: LongWord;
  a_d, ia, h, h2: Double;
  ja: Tb64u64;
  jp: Int64;
  jm: Int64;
  jp_idx, jm_idx: Int64;
  jp_shr5: Int64;
  jm_shr5: Int64;
  sp, sm: Tb64u64;
  te, to_, rp, rm, r: Double;
  ub_s, lb_s: Single;
  iln2h: Double;
  iln2l: Double;
  z2: Double;
begin
  t.f := x;
  z := x;
  ax := t.u shl 1;
  if ax > $8565A9F8 then  // |x| >~ 89.4
  begin
    if ax >= $FF000000 then begin Result := x * x; Exit; end;  // inf or nan
    Result := 2.0 * Single(3.4028234663852886e+38);  // 2.0f*0x1.fffffep127f => overflow to +inf
    Exit;
  end;
  if ax < $7C000000 then  // |x| < 0.125
  begin
    if ax < $74000000 then  // |x| < 0x1p-11
    begin
      if ax < $66000000 then  // |x| < ~2^-25
        Result := pcr_fmaf(pcr_fabsf(x), Single(2.9802322387695312e-08), 1.0)
      else
        Result := Single((0.5 * z) * z + 1.0);
      Exit;
    end;
    z2 := z * z;
    Result := Single(1.0 + z2 * ((cp_arr[0] + z2*cp_arr[1]) + (z2*z2)*(cp_arr[2] + z2*cp_arr[3])));
    Exit;
  end;
  a_d := iln2 * z;
  ia := pcr_roundeven(a_d);
  h := a_d - ia;
  h2 := h * h;
  ja.f := ia + 6755399441055744.0;  // ia + 0x1.8p52
  jp := Int64(ja.u);
  jm := -jp;
  jp_idx := jp and 31;
  jp_shr5 := jp shr 5;  // jp > 0, so logical = arithmetic
  jm_idx := jm and 31;
  jm_shr5 := sar_i64(jm, 5);
  sp.u := tb_arr[jp_idx] + (UInt64(jp_shr5) shl 52);
  sm.u := tb_arr[jm_idx] + (UInt64(jm_shr5) shl 52);
  te := c_arr[0] + h2 * c_arr[2];
  to_ := c_arr[1] + h2 * c_arr[3];
  rp := sp.f * (te + h * to_);
  rm := sm.f * (te - h * to_);
  r := rp + rm;
  ub_s := Single(r);
  lb_s := Single(r - 1.45e-10 * r);
  if ub_s <> lb_s then
  begin
    iln2h := 46.16624128818512;  // 0x1.7154765p+5
    iln2l := 2.026170940661134e-08;  // 0x1.5c17f0bbbe88p-26
    h := (iln2h * z - ia) + iln2l * z;
    h2 := h * h;
    te := ch_arr[0] + h2*ch_arr[2] + (h2*h2)*(ch_arr[4] + h2*ch_arr[6]);
    to_ := ch_arr[1] + h2*(ch_arr[3] + h2*ch_arr[5]);
    r := sp.f*(te + h*to_) + sm.f*(te - h*to_);
    ub_s := Single(r);
  end;
  Result := ub_s;
end;

// ── 2.01 logf ────────────────────────────────────────────────────────────────
function pcr_logf(x: Single): Single;
const
  tr: array[0..64] of Double = (
    1.0, 0.9846153855323792, 0.969696968793869, 0.9552238807082176,
    0.9411764703691006, 0.9275362323969603, 0.9142857138067484, 0.9014084506779909,
    0.8888888880610466, 0.876712329685688, 0.8648648653179407, 0.8533333335071802,
    0.8421052638441324, 0.8311688303947449, 0.8205128200352192, 0.8101265821605921,
    0.8000000007450581, 0.7901234570890665, 0.7804878056049347, 0.7710843365639448,
    0.7619047611951828, 0.7529411762952805, 0.744186045601964, 0.7356321830302477,
    0.7272727265954018, 0.7191011235117912, 0.7111111115664244, 0.7032967042177916,
    0.6956521738320589, 0.6881720423698425, 0.6808510646224022, 0.673684211447835,
    0.666666666045785, 0.6597938146442175, 0.6530612241476774, 0.6464646458625793,
    0.6400000005960464, 0.6336633656173944, 0.627450980246067, 0.6213592234998941,
    0.615384615957737, 0.6095238104462624, 0.6037735845893621, 0.5981308408081532,
    0.5925925932824612, 0.5871559642255306, 0.5818181820213795, 0.5765765774995089,
    0.5714285708963871, 0.5663716811686754, 0.5614035092294216, 0.5565217398107052,
    0.551724137738347, 0.5470085479319096, 0.5423728805035353, 0.5378151256591082,
    0.533333333209157, 0.5289256200194359, 0.5245901644229889, 0.5203252024948597,
    0.5161290317773819, 0.5120000001043081, 0.5079365074634552, 0.5039370078593493,
    0.5
  );
  tl: array[0..64] of Double = (
    -3.5e-14, 0.01550418560460768, 0.030771659598041262, 0.04580953591484388,
    0.060624622049230484, 0.07522342068457975, 0.08961215921352109, 0.10379679371071239,
    0.11778303658767103, 0.13157635674094637, 0.14518200932059394, 0.15860502997287676,
    0.17185025611171698, 0.18492233942529956, 0.19782574391196148, 0.2105647692528338,
    0.2231435503828522, 0.23556607093438212, 0.2478361629732237, 0.25995752545552514,
    0.27193371641492936, 0.28376817336344023, 0.29546421411616175, 0.3070250364881339,
    0.3184537320498222, 0.32975328648884833, 0.34092658633027395, 0.3519764218474708,
    0.36290549380574877, 0.37371641072487166, 0.3844116977461438, 0.39499380687295393,
    0.40546510903945193, 0.41582789482353383, 0.426084395834734, 0.4362367677062056,
    0.4462871016970619, 0.456237434616602, 0.46608973015739485, 0.47584590454978676,
    0.4855078148503432, 0.4950772652844173, 0.5045560112762293, 0.5139457516260683,
    0.5232481426003596, 0.53246479729783, 0.5415972820834634, 0.5506471163519167,
    0.5596157888667103, 0.5685047357891911, 0.5773153642198814, 0.5860490437811823,
    0.5947071080959038, 0.6032908497500271, 0.6118015426775647, 0.6202404104794182,
    0.6286086596551698, 0.6369074618877882, 0.6451379604422272, 0.6533012734679021,
    0.6613984831766526, 0.6694306537388675, 0.6773988245230937, 0.6853040031279882,
    0.6931471805599103
  );
  b: array[0..2] of Double = (
    1.0000000014444432, -0.5000150201101822, 0.3332912677830366
  );
  c: array[0..6] of Double = (
    -0.5, 0.3333333333337377, -0.25000000000064687, 0.1999999917184808,
    -0.16666665689564816, 0.14291141156288886, -0.12505426680916787
  );
var
  t: Tb32u32;
  tz: Tb64u64;
  ux, m: LongWord;
  j: LongWord;
  e, n: Integer;
  z, z2, r, el, dr, f: Double;
  ub, lb: Single;
begin
  t.f := x;
  ux := t.u;
  if (ux < LongWord(1 shl 23)) or (ux >= $7F800000) then
  begin
    if (ux = 0) or (ux >= $7F800000) then
    begin
      if (ux shl 1) = 0 then // +/-0 -> -inf
        begin t.u := $FF800000; Result := t.f; Exit; end;
      if ux = $7F800000 then begin Result := x; Exit; end; // +inf
      if (ux shl 1) > $FF000000 then begin Result := x + x; Exit; end; // nan
      // x < 0 -> NaN
      t.u := $FFC00000; Result := t.f; Exit;
    end;
    // subnormal
    n := 23 - Integer(pcr_bsr32(ux));
    ux := ux shl n;
    ux := ux - LongWord(n shl 23);
  end;
  if ux = LongWord(127 shl 23) then begin Result := 0.0; Exit; end;
  m := ux and LongWord((1 shl 23) - 1);
  j := (m + LongWord(1 shl (23-7))) shr (23-6);
  e := sar_i32(Integer(ux), 23) - 127;
  tz.u := (UInt64(m) shl 29) or (UInt64($3FF) shl 52);
  z := tz.f * tr[j] - 1.0;
  z2 := z * z;
  r := ((e * 0.6931471805599453 + tl[j]) + z * b[0]) + z2 * (b[1] + z * b[2]);
  ub := Single(r);
  lb := Single(r + 2.2572521629626863e-10); // 0x1.f06p-33
  if ub <> lb then
  begin
    f := z2 * ((c[0] + z*c[1]) + z2*((c[2] + z*c[3]) + z2*(c[4] + z*c[5] + z2*c[6])));
    if pcr_fabsf(x - 1.0) < 9.765625e-04 then // 0x1p-10
    begin
      Result := Single(z + f); Exit;
    end;
    f := f - 1.8641886737243033e-15 * e; // 0x1.0ca86c3898dp-49
    f := f + z;
    f := f + tl[j] - tl[0];
    el := e * 0.6931471805599472; // 0x1.62e42fefa3ap-1
    r := el + f;
    ub := Single(r);
    tz.f := r;
    if (tz.u and ((UInt64(1) shl 28) - 1)) = 0 then
    begin
      dr := (el - r) + f;
      r := r + dr * 64.0;
      ub := Single(r);
    end;
  end;
  Result := ub;
end;

// ── 2.02 exp2f ───────────────────────────────────────────────────────────────
function pcr_exp2f(x: Single): Single;
const
  tb: array[0..63] of Tb64u64 = (
    (u: $3FF0000000000000), (u: $3FF02C9A3E778061), (u: $3FF059B0D3158574), (u: $3FF0874518759BC8),
    (u: $3FF0B5586CF9890F), (u: $3FF0E3EC32D3D1A2), (u: $3FF11301D0125B51), (u: $3FF1429AAEA92DE0),
    (u: $3FF172B83C7D517B), (u: $3FF1A35BEB6FCB75), (u: $3FF1D4873168B9AA), (u: $3FF2063B88628CD6),
    (u: $3FF2387A6E756238), (u: $3FF26B4565E27CDD), (u: $3FF29E9DF51FDEE1), (u: $3FF2D285A6E4030B),
    (u: $3FF306FE0A31B715), (u: $3FF33C08B26416FF), (u: $3FF371A7373AA9CB), (u: $3FF3A7DB34E59FF7),
    (u: $3FF3DEA64C123422), (u: $3FF4160A21F72E2A), (u: $3FF44E086061892D), (u: $3FF486A2B5C13CD0),
    (u: $3FF4BFDAD5362A27), (u: $3FF4F9B2769D2CA7), (u: $3FF5342B569D4F82), (u: $3FF56F4736B527DA),
    (u: $3FF5AB07DD485429), (u: $3FF5E76F15AD2148), (u: $3FF6247EB03A5585), (u: $3FF6623882552225),
    (u: $3FF6A09E667F3BCD), (u: $3FF6DFB23C651A2F), (u: $3FF71F75E8EC5F74), (u: $3FF75FEB564267C9),
    (u: $3FF7A11473EB0187), (u: $3FF7E2F336CF4E62), (u: $3FF82589994CCE13), (u: $3FF868D99B4492ED),
    (u: $3FF8ACE5422AA0DB), (u: $3FF8F1AE99157736), (u: $3FF93737B0CDC5E5), (u: $3FF97D829FDE4E50),
    (u: $3FF9C49182A3F090), (u: $3FFA0C667B5DE565), (u: $3FFA5503B23E255D), (u: $3FFA9E6B5579FDBF),
    (u: $3FFAE89F995AD3AD), (u: $3FFB33A2B84F15FB), (u: $3FFB7F76F2FB5E47), (u: $3FFBCC1E904BC1D2),
    (u: $3FFC199BDD85529C), (u: $3FFC67F12E57D14B), (u: $3FFCB720DCEF9069), (u: $3FFD072D4A07897C),
    (u: $3FFD5818DCFBA487), (u: $3FFDA9E603DB3285), (u: $3FFDFC97337B9B5F), (u: $3FFE502EE78B3FF6),
    (u: $3FFEA4AFA2A490DA), (u: $3FFEFA1BEE615A27), (u: $3FFF50765B6E4540), (u: $3FFFA7C1819E90D8)
  );
  b: array[0..3] of Double = (
    1.0, 0.6931471805202393, 0.2402288551437867, 0.05550459682799693
  );
  c: array[0..5] of Double = (
    0.6931471805599453, 0.24022650695910072, 0.05550410866402609,
    0.009618129107500536, 0.001333362331326638, 0.00015403602972146417
  );
var
  t: Tb32u32;
  u: Tb32u32;
  sv: Tb64u64;
  ux: LongWord;
  k, m_int, msk: Integer;
  offd, xd, h, h2, r: Double;
  ub, lb: Single;
begin
  t.f := x;
  // check if x is an exact integer (low 16 bits of mantissa zero)
  if (t.u and $FFFF) = 0 then
  begin
    k := Integer((t.u shr 23) and $FF) - 127;
    if (k >= 0) and (k < 9) and ((t.u shl (9 + k)) = 0) then
    begin
      msk := sar_i32(Integer(t.u), 31); // -1 if negative, 0 if positive
      m_int := Integer((t.u and $7FFFFF) or (1 shl 23)) shr (23 - k);
      m_int := (m_int xor msk) - msk + 127;
      if (m_int > 0) and (m_int < 255) then
      begin
        t.u := LongWord(m_int) shl 23;
        Result := t.f; Exit;
      end
      else if (m_int <= 0) and (m_int > -23) then
      begin
        t.u := LongWord(1) shl LongWord(22 + m_int);
        Result := t.f; Exit;
      end;
    end;
  end;
  ux := t.u shl 1;
  if (ux >= $86000000) or (ux < $65000000) then
  begin
    if ux < $65000000 then begin Result := 1.0 + x; Exit; end; // |x| < 0x1p-26
    // x in [-128, -149) falls through to main computation (produces subnormal result)
    if not ((t.u >= $C3000000) and (t.u < $C3150000)) then
    begin
      // as_special
      if ux >= $FF000000 then // inf or nan
      begin
        if ux > $FF000000 then begin Result := x + x; Exit; end; // nan
        if t.u shr 31 <> 0 then begin Result := 0.0; Exit; end; // -inf -> 0
        Result := x; Exit; // +inf
      end;
      if t.u >= $C3150000 then // x <= -149
      begin
        // underflow path
        xd := x;
        h := 1.401298464324817e-45 + (xd + 149.0) * 7.006492321624085e-46; // 0x1p-149 + (z+149)*0x1p-150
        h := pcr_fmax(h, 3.503246160812043e-46); // 0x1p-151
        Result := Single(h); Exit;
      end;
      // x >= 128 -> overflow
      Result := 1.7014118346046923e+38 * 1.7014118346046923e+38; Exit; // 0x1p127 * 0x1p127 = overflow
    end;
    // x in [-128, -149): fall through to main computation
  end;
  offd := 105553116266496.0; // 0x1.8p46
  xd := x;
  h := xd - ((xd + offd) - offd);
  h2 := h * h;
  u.f := x + 196608.0; // 0x1.8p17f
  sv := tb[u.u and $3F];
  sv.u := sv.u + (UInt64(u.u shr 6) shl 52);
  r := sv.f * ((b[0] + h * b[1]) + h2 * (b[2] + h * b[3]));
  ub := Single(r);
  lb := Single(r - r * 1.4438228390645236e-10); // r * eps where eps = 0x1.3d8p-33
  if ub <> lb then
  begin
    if ux <= $79E7526E then
    begin
      if t.u = $3B429D37 then begin Result := Single(1.0020605325698853) - Single(2.9802322387695312e-08); Exit; end;
      if t.u = $BCF3A937 then begin Result := Single(0.9795942902565002) - Single(1.4901161193847656e-08); Exit; end;
      if t.u = $B8D3D026 then begin Result := Single(0.9999299645423889) + Single(1.4901161193847656e-08); Exit; end;
    end;
    r := sv.f + (sv.f * h) * ((c[0] + h*c[1]) + h2*((c[2] + h*c[3]) + h2*(c[4] + h*c[5])));
    ub := Single(r);
  end;
  Result := ub;
end;

// ── 2.03 log1pf ──────────────────────────────────────────────────────────────
function pcr_log1pf(x: Single): Single;
const
  x0: array[0..31] of Double = (
    0.9846153259277344, 0.9552239179611206, 0.9275362491607666, 0.9014084339141846,
    0.8767123222351074, 0.8533333539962769, 0.8311688899993896, 0.810126543045044,
    0.790123462677002, 0.7710843086242676, 0.7529411315917969, 0.7356321811676025,
    0.719101071357727, 0.7032966613769531, 0.6881721019744873, 0.6736842393875122,
    0.6597938537597656, 0.6464647054672241, 0.6336634159088135, 0.6213592290878296,
    0.6095237731933594, 0.5981308221817017, 0.5871559381484985, 0.5765765905380249,
    0.5663716793060303, 0.5565217733383179, 0.5470085144042969, 0.5378150939941406,
    0.5289256572723389, 0.5203251838684082, 0.5119999647140503, 0.5039370059967041
  );
  lix: array[0..31] of Double = (
    0.015504246151250851, 0.04580949692638588, 0.07522340262177525, 0.10379681231873428,
    0.13157636524993893, 0.15860500597289098, 0.18492226772413786, 0.21056481754676373,
    0.2355660638728254, 0.2599575617004688, 0.2837682327459302, 0.3070250390308412,
    0.3297533590262705, 0.3519764827724638, 0.37371632412254996, 0.3949937654106705,
    0.41582783554970704, 0.43623667551594897, 0.45623735526113346, 0.4758458955673773,
    0.4950773264131371, 0.5139457827778414, 0.5324648417209502, 0.5506470937489146,
    0.5685047390885981, 0.5860489835469292, 0.6032909110533698, 0.6202404693671432,
    0.6369073914671951, 0.6533013092762884, 0.6694307228711411, 0.6853040068348488
  );
  b: array[0..7] of Double = (
    1.0, -0.5, 0.333333333333704, -0.2500000000005929,
    0.1999999921853749, -0.16666665744658113, 0.1429098594542405, -0.125052714602758
  );
  c: array[0..3] of Double = (
    0.9999999964978914, -0.49999999241150506, 0.33339251544971726, -0.2500690494115668
  );
  ln2: Double = 0.6931471805599453;   // 0x1.62e42fefa39efp-1
  ln2h: Double = 0.693145751953125;   // 0x1.62e4p-1
  ln2l: Double = 1.4286068203094173e-06; // 0x1.7f7d1cf79abcap-20
  lix0_offset: Double = 1.063904520037795e-11; // 0x1.7654p-37
var
  t: Tb32u32;
  tp: Tb64u64;
  r: Tb64u64;
  xd: Tb64u64;
  ux, ax: LongWord;
  e: Integer;
  m52: UInt64;
  j: LongWord;
  zd, z2, z4, f: Double;
  Lh, Ll, rh, rl: Double;
  fh_s: Single;
  Fdbl: Double;
  fl_s, tfl: Single;
  ub, lb: Single;
begin
  zd := x;
  t.f := x;
  ux := t.u;
  if ux >= $BF800000 then // x <= -1 (including -inf, -nan)
  begin
    if t.u = $BF800000 then // x = -1 -> -inf
      begin t.u := $FF800000; Result := t.f; Exit; end;
    if t.u = $7F800000 then begin Result := x; Exit; end; // +inf
    ax := t.u shl 1;
    if ax > $FF000000 then begin Result := x + x; Exit; end; // nan
    // x < -1 -> NaN
    t.u := $FFC00000; Result := t.f; Exit;
  end;
  ax := ux and $7FFFFFFF;
  if ax >= $7F800000 then // +inf or +nan
  begin
    if ax > $7F800000 then begin Result := x + x; Exit; end;
    Result := x; Exit;
  end;
  if ax < $3C880000 then // |x| < 0x1.1p-6
  begin
    if ax < $33000000 then // |x| < 0x1p-25
    begin
      if ax = 0 then begin Result := x; Exit; end;
      Result := pcr_fmaf(x, -x, x);
      Exit;
    end;
    z2 := zd * zd;
    z4 := z2 * z2;
    f := z2 * ((b[1] + zd*b[2]) + z2*(b[3] + zd*b[4]) + z4*(b[5] + zd*(b[6] + zd*b[7])));
    r.f := zd + f;
    if (r.u and $FFFFFFF) = 0 then
      r.f := r.f + 16384.0 * (f + (zd - r.f)); // 0x1p14
    Result := Single(r.f);
    Exit;
  end;
  // larger path
  tp.f := zd + 1.0;
  e := Integer(tp.u shr 52) - $3FF;
  m52 := tp.u and ($FFFFFFFFFFFFFFFF shr 12); // ~0ull>>12
  j := LongWord((tp.u shr (52-5)) and 31);
  xd.u := m52 or (UInt64($3FF) shl 52);
  zd := xd.f * x0[j] - 1.0; // z is exact for x<0x1.0cp+30
  z2 := zd * zd;
  rh := (ln2 * e + lix[j]) + zd * ((c[0] + zd*c[1]) + z2*(c[2] + zd*c[3]));
  ub := Single(rh);
  lb := Single(rh - 2.1555e-11); // eps
  if ub <> lb then
  begin
    z4 := z2 * z2;
    f := z2 * ((b[1] + zd*b[2]) + z2*(b[3] + zd*b[4]) + z4*(b[5] + zd*(b[6] + zd*b[7])));
    Lh := ln2h * e;
    Ll := ln2l * e;
    Ll := Ll + zd;
    rl := lix[j] - lix0_offset; // subtract offset 0x1.7654p-37
    rh := Lh + rl;
    rl := ((Lh - rh) + rl) + (Ll + f);
    fh_s := Single(rh + rl);
    Fdbl := (rh - Double(fh_s)) + rl;
    fl_s := Single(Fdbl);
    tfl := fl_s * 2.0;
    if (fh_s + tfl) - fh_s = tfl then
      fl_s := fl_s + pcr_copysignf(0.5, Single(Fdbl - Double(fl_s))) * pcr_fabsf(fl_s);
    ub := fh_s + fl_s;
  end;
  Result := ub;
end;

// ── 2.04 exp2m1f ─────────────────────────────────────────────────────────────
function pcr_exp2m1f(x: Single): Single;
const
  tb_e2m1: array[0..15] of Double = (
    1.0, 1.0442737824274138, 1.0905077326652577, 1.1387886347566916,
    1.189207115002721, 1.241857812073484, 1.2968395546510096, 1.3542555469368927,
    1.4142135623730951, 1.4768261459394993, 1.5422108254079405, 1.6104903319492543,
    1.681792830507429, 1.7562521603732995, 1.8340080864093424, 1.9152065613971474
  );
  // c_table[0]: 1-coeff (just ln2)
  // c_table[1..]: various polynomial sets stored inline below
var
  t: Tb32u32;
  su: Tb64u64;
  zd: Double;
  ux, ax: LongWord;
  z2, r: Double;
  ia: Double;
  i64a: Int64;
  j, e64: Int64;
  s, c0, c2, c4, w: Double;
  c0v, c1v, c2v, c3v, c4v, c5v, c6v, c7v: Double;
begin
  t.f := x;
  zd := x;
  ux := t.u;
  ax := ux and $7FFFFFFF;
  if ux >= $C1C80000 then // x <= -25
  begin
    if ax > ($FF shl 23) then begin Result := x + x; Exit; end; // nan
    if ux = $FF800000 then begin Result := -1.0; Exit; end; // -inf
    Result := -1.0 + Single(1.4901161193847656e-08); Exit; // -1 + 2^-26
  end;
  if ax >= $43000000 then // x >= 128
  begin
    if ax > ($FF shl 23) then begin Result := x + x; Exit; end; // nan
    if ux = $7F800000 then begin Result := x; Exit; end; // +inf
    Result := 3.4028234663852886e+38 + 1.0141204801825835e+31; Exit; // overflow
  end;
  if ax < $3DF95F1F then // |x| < 8.44e-2/log(2) (small path)
  begin
    z2 := zd * zd;
    if ax < $3D67A4CC then // |x| < 3.92e-2/log(2)
    begin
      if ax < $3CAA2FEE then // |x| < 1.44e-2/log(2)
      begin
        if ax < $3BAC1405 then // |x| < 3.64e-3/log(2)
        begin
          if ax < $3A358876 then // |x| < 4.8e-4/log(2)
          begin
            if ax < $37D32EF6 then // |x| < 1.745e-5/log(2)
            begin
              if ax < $331FDD82 then // |x| < 2.58e-8/log(2)
              begin
                if ax < $2538AA3B then // |x| < 0x1.715476p-53
                begin
                  r := 0.6931471805599453; // c[0] only
                end
                else
                begin
                  r := 0.6931471805599454 + zd * 0.24022650695910072; // c[0]+z*c[1]
                end;
              end
              else
              begin
                if ux = $B3D85005 then begin Result := Single(-6.981959899121648e-08) - Single(9.926167350636332e-24); Exit; end;
                if ux = $3338428D then begin Result := Single(2.973696133778958e-08) + Single(8.271806125530277e-25); Exit; end;
                c0v := 0.6931471805599453; c1v := 0.24022650696367256; c2v := 0.05550410866482101;
                r := c0v + zd * (c1v + zd * c2v);
              end;
            end
            else
            begin
              if ux = $388BCA4F then begin Result := Single(4.620431354851461e-05) - Single(5.082197683525802e-21); Exit; end;
              c0v := 0.6931471805599453; c1v := 0.24022650695910072;
              c2v := 0.05550410930422927; c3v := 0.009618129107686644;
              r := (c0v + zd * c1v) + z2 * (c2v + zd * c3v);
            end;
          end
          else
          begin
            c0v := 0.6931471805599453; c1v := 0.2402265069590641;
            c2v := 0.055504108664832436; c3v := 0.009618134417479019; c4v := 0.001333355815169557;
            r := (c0v + zd * c1v) + z2 * (c2v + zd * (c3v + zd * c4v));
          end;
        end
        else
        begin
          c0v := 0.6931471805599454; c1v := 0.24022650695910067;
          c2v := 0.05550410866322344; c3v := 0.009618129107951784;
          c4v := 0.0013333656890870747; c5v := 0.0001540353035411431;
          r := (c0v + zd * c1v) + z2 * ((c2v + zd * c3v) + z2 * (c4v + zd * c5v));
        end;
      end
      else
      begin
        c0v := 0.6931471805599453; c1v := 0.24022650695910544;
        c2v := 0.05550410866481867; c3v := 0.009618129095800282;
        c4v := 0.0013333558164648996; c5v := 0.0001540427006814203; c6v := 1.5252733783448092e-05;
        r := (c0v + zd * c1v) + z2 * ((c2v + zd * c3v) + z2 * (c4v + zd * (c5v + zd * c6v)));
      end;
    end
    else
    begin
      c0v := 0.6931471805599453; c1v := 0.24022650695910078;
      c2v := 0.05550410866490448; c3v := 0.009618129107593989;
      c4v := 0.0013333557866797964; c5v := 0.0001540353074233086;
      c6v := 1.5255751829253785e-05; c7v := 1.3215486693701843e-06;
      r := ((c0v + zd * c1v) + z2 * (c2v + zd * c3v)) +
           (z2*z2) * ((c4v + zd * c5v) + z2 * (c6v + zd * c7v));
    end;
    r := r * zd;
    Result := Single(r);
    Exit;
  end;
  // main table path
  c0v := 0.043321698784995886; c1v := 0.0009383847928200837;
  c2v := 1.3550807712983854e-05; c3v := 1.4676119301623784e-07;
  c4v := 1.271309415715539e-09; c5v := 9.382438953978075e-12;
  ia := 16.0 * zd;
  i64a := Int64(Trunc(ia)); // floor via trunc (ia > 0 here since |x|>=0.0844 and we're in large path)
  // Actually need proper floor:
  if ia < 0.0 then
  begin
    if ia <> i64a then i64a := i64a - 1;
  end;
  ia := Double(i64a);
  j := i64a and $F;
  e64 := i64a - j;
  e64 := sar_i64(e64, 4);
  s := tb_e2m1[j];
  su.u := UInt64(e64 + $3FF) shl 52;
  s := s * su.f;
  r := ia - ia; // h = a - ia but a=16*z, ia=floor(a)
  // h = 16*zd - ia
  r := 16.0 * zd - Double(i64a);  // h
  z2 := r * r;
  c0 := c0v + r * c1v;
  c2 := c2v + r * c3v;
  c4 := c4v + r * c5v;
  c0 := c0 + z2 * (c2 + z2 * c4);
  w := s * r;
  Result := Single((s - 1.0) + w * c0);
end;

// ── 2.05 expm1f ──────────────────────────────────────────────────────────────
function pcr_expm1f(x: Single): Single;
const
  c_fast: array[0..3] of Double = (
    1.0, 0.021660849391257477, 0.0002345984913513542, 1.6938658699950235e-06
  );
  ch: array[0..5] of Double = (
    0.02166084939249829, 0.0002345961982022468, 1.6938509724129055e-06,
    9.172562701702629e-09, 3.973729405780548e-11, 1.4345723178374038e-13
  );
  td: array[0..31] of Tb64u64 = (
    (u: $3FF0000000000000), (u: $3FF059B0D3158574), (u: $3FF0B5586CF9890F), (u: $3FF11301D0125B51),
    (u: $3FF172B83C7D517B), (u: $3FF1D4873168B9AA), (u: $3FF2387A6E756238), (u: $3FF29E9DF51FDEE1),
    (u: $3FF306FE0A31B715), (u: $3FF371A7373AA9CB), (u: $3FF3DEA64C123422), (u: $3FF44E086061892D),
    (u: $3FF4BFDAD5362A27), (u: $3FF5342B569D4F82), (u: $3FF5AB07DD485429), (u: $3FF6247EB03A5585),
    (u: $3FF6A09E667F3BCD), (u: $3FF71F75E8EC5F74), (u: $3FF7A11473EB0187), (u: $3FF82589994CCE13),
    (u: $3FF8ACE5422AA0DB), (u: $3FF93737B0CDC5E5), (u: $3FF9C49182A3F090), (u: $3FFA5503B23E255D),
    (u: $3FFAE89F995AD3AD), (u: $3FFB7F76F2FB5E47), (u: $3FFC199BDD85529C), (u: $3FFCB720DCEF9069),
    (u: $3FFD5818DCFBA487), (u: $3FFDFC97337B9B5F), (u: $3FFEA4AFA2A490DA), (u: $3FFF50765B6E4540)
  );
  b_small: array[0..7] of Double = (
    0.49999999999999656, 0.16666666666667135, 0.041666666668544565, 0.008333333332479211,
    0.0013888886118215516, 0.00019841274040338812, 2.4816724201894197e-05, 2.755731951095977e-06
  );
  iln2: Double = 46.16624130844683;  // 0x1.71547652b82fep+5
  big: Double = 6755399441055744.0;  // 0x1.8p52
  iln2h: Double = 46.16624128818512; // 0x1.7154765p+5
  iln2l: Double = 2.026170940661134e-08; // 0x1.5c17f0bbbe88p-26
var
  t: Tb32u32;
  u: Tb64u64;
  sv: Tb64u64;
  ux, ax: LongWord;
  zd, a, ia, h, h2, r: Double;
  c2d, c0d: Double;
  ub, lb: Single;
begin
  t.f := x;
  ux := t.u;
  ax := ux shl 1;
  zd := x;
  if ax < $7C400000 then // |x| < 0.15625
  begin
    if ax < $676A09E8 then // |x| < 0x1.6a09e8p-24
    begin
      if ax = 0 then begin Result := x; Exit; end;
      Result := pcr_fmaf(pcr_fabsf(x), Single(2.9802322387695312e-08), x); // fmaf(|x|,0x1p-25,x)
      Exit;
    end;
    // Horner polynomial for small x
    h2 := zd * zd;
    r := zd + h2 * ((b_small[0] + zd*b_small[1]) + h2*(b_small[2] + zd*b_small[3]) +
         (h2*h2)*((b_small[4] + zd*b_small[5]) + h2*(b_small[6] + zd*b_small[7])));
    Result := Single(r); Exit;
  end;
  if ax >= $8562E430 then // |x| > 88.72
  begin
    if ax > ($FF shl 24) then begin Result := x + x; Exit; end; // nan
    if (ux shr 31) <> 0 then // x < 0
    begin
      if ax = ($FF shl 24) then begin Result := -1.0; Exit; end; // -inf
      Result := -1.0 + Single(1.4901161193847656e-08); Exit;
    end;
    if ax = ($FF shl 24) then begin Result := x * x; Exit; end; // +inf
    t.u := $7F7FFFFF;  // 0x1.fffffep127
    Result := t.f * Single(zd); Exit;
  end;
  a := iln2 * zd;
  ia := pcr_roundeven(a);
  h := a - ia;
  h2 := h * h;
  u.f := ia + big;
  c2d := c_fast[2] + h * c_fast[3];
  c0d := c_fast[0] + h * c_fast[1];
  sv.u := td[u.u and $1F].u + (UInt64(u.u shr 5) shl 52);
  r := (c0d + h2 * c2d) * sv.f - 1.0;
  ub := Single(r);
  lb := Single(r - sv.f * 1.433306806575274e-10); // sv.f * 0x1.3b3p-33
  if ub <> lb then
  begin
    if ux > $C18AA123 then // x < -17.32
      begin Result := -1.0 + Single(1.4901161193847656e-08); Exit; end;
    h := (iln2h * zd - ia) + iln2l * zd;
    h2 := h * h;
    r := (sv.f - 1.0) + (sv.f * h) * ((ch[0] + h*ch[1]) + h2*((ch[2] + h*ch[3]) + h2*(ch[4] + h*ch[5])));
    ub := Single(r);
  end;
  Result := ub;
end;

// ── 2.06 exp10f ──────────────────────────────────────────────────────────────
function pcr_exp10f(x: Single): Single;
const
  c_exp10: array[0..5] of Double = (
    0.6931471805599453, 0.24022650695910072, 0.05550410866402609,
    0.009618129107500536, 0.001333362331326638, 0.00015403602972146417
  );
  b_exp10: array[0..3] of Double = (
    1.0, 0.021660849391257477, 0.0002345984913513542, 1.6938658699950235e-06
  );
  tb_e10: array[0..31] of UInt64 = (
    $3FF0000000000000, $3FF059B0D3158574, $3FF0B5586CF9890F, $3FF11301D0125B51,
    $3FF172B83C7D517B, $3FF1D4873168B9AA, $3FF2387A6E756238, $3FF29E9DF51FDEE1,
    $3FF306FE0A31B715, $3FF371A7373AA9CB, $3FF3DEA64C123422, $3FF44E086061892D,
    $3FF4BFDAD5362A27, $3FF5342B569D4F82, $3FF5AB07DD485429, $3FF6247EB03A5585,
    $3FF6A09E667F3BCD, $3FF71F75E8EC5F74, $3FF7A11473EB0187, $3FF82589994CCE13,
    $3FF8ACE5422AA0DB, $3FF93737B0CDC5E5, $3FF9C49182A3F090, $3FFA5503B23E255D,
    $3FFAE89F995AD3AD, $3FFB7F76F2FB5E47, $3FFC199BDD85529C, $3FFCB720DCEF9069,
    $3FFD5818DCFBA487, $3FFDFC97337B9B5F, $3FFEA4AFA2A490DA, $3FFF50765B6E4540
  );
  ex10: array[0..9] of Single = (
    10.0, 100.0, 1000.0, 10000.0, 100000.0,
    1000000.0, 10000000.0, 100000000.0, 1000000000.0, 10000000000.0
  );
  iln102:  Double = 106.30169903639559;    // 0x1.a934f0979a371p+6
  iln102h: Double = 3.3219280913472176;    // 0x1.a934f09p+1
  iln102l: Double = 3.5401447880558664e-09; // 0x1.e68dc57f2496p-29
var
  t: Tb32u32;
  sv: Tb64u64;
  ux: LongWord;
  zd, a, ia, h, h2, r: Double;
  ja: Int64;
  ub, lb: Single;
  k, bt, msk, cnt: LongWord;
  pv: LongWord;
begin
  t.f := x;
  ux := t.u shl 1;
  zd := x;
  if (ux > $84344134) or (ux < $72ADF1C6) then
  begin
    if ux < $72ADF1C6 then // |x| < 0x1.adf1c6p-13
    begin
      Result := Single(1.0 + zd*(2.302585092994046 + zd*(2.650949055239199 + zd*2.034678592293476)));
      Exit;
    end;
    if ux >= LongWord($FF shl 24) then // inf or nan
    begin
      if ux > LongWord($FF shl 24) then begin Result := x + x; Exit; end;
      if t.u shr 31 <> 0 then begin Result := 0.0; Exit; end;
      Result := x; Exit;
    end;
    if t.u > $C23369F4 then // x < -0x1.66d3e8p+5
    begin
      h := 1.401298464324817e-45 + (zd + 44.8534693539332) * 2.3275063689815626e-45;
      h := pcr_fmax(h, 3.503246160812043e-46);
      Result := Single(h); Exit;
    end;
    if t.u < $80000000 then // x > 0x1.344134p+5
    begin
      Result := 1.7014118346046923e+38 * 1.7014118346046923e+38; Exit;
    end;
  end;
  // check for integer power of 10
  if (t.u shl 12) = 0 then
  begin
    k := (t.u shr 20) - 1016;
    if k <= 26 then
    begin
      bt := LongWord(1) shl k;
      msk := $7551101;
      if (bt and msk) <> 0 then
      begin
        // popcount(msk & (bt-1))
        pv := msk and (bt - 1);
        cnt := 0;
        while pv <> 0 do begin Inc(cnt); pv := pv and (pv - 1); end;
        Result := ex10[cnt]; Exit;
      end;
    end;
  end;
  a := iln102 * zd;
  ia := pcr_roundeven(a);
  h := a - ia;
  ja := Int64(Trunc(ia));
  sv.u := tb_e10[ja and $1F] + (UInt64(ja shr 5) shl 52);
  h2 := h * h;
  r := ((b_exp10[0] + h*b_exp10[1]) + h2*(b_exp10[2] + h*b_exp10[3])) * sv.f;
  ub := Single(r);
  lb := Single(r - r * 1.45e-10);
  if ub <> lb then
  begin
    h := (iln102h * zd - ia * 0.03125) + iln102l * zd;
    h2 := h * h;
    r := sv.f + (sv.f * h) * ((c_exp10[0] + h*c_exp10[1]) + h2*((c_exp10[2] + h*c_exp10[3]) + h2*(c_exp10[4] + h*c_exp10[5])));
    ub := Single(r);
  end;
  Result := ub;
end;

// ── 2.07 log10f ──────────────────────────────────────────────────────────────
function pcr_log10f(x: Single): Single;
const
  tr10: array[0..64] of Double = (
    1.0, 0.9846153855323792, 0.969696968793869, 0.9552238807082176,
    0.9411764703691006, 0.9275362323969603, 0.9142857138067484, 0.9014084506779909,
    0.8888888880610466, 0.876712329685688, 0.8648648653179407, 0.8533333335071802,
    0.8421052638441324, 0.8311688303947449, 0.8205128200352192, 0.8101265821605921,
    0.8000000007450581, 0.7901234570890665, 0.7804878056049347, 0.7710843365639448,
    0.7619047611951828, 0.7529411762952805, 0.744186045601964, 0.7356321830302477,
    0.7272727265954018, 0.7191011235117912, 0.7111111115664244, 0.7032967042177916,
    0.6956521738320589, 0.6881720423698425, 0.6808510646224022, 0.673684211447835,
    0.666666666045785, 0.6597938146442175, 0.6530612241476774, 0.6464646458625793,
    0.6400000005960464, 0.6336633656173944, 0.627450980246067, 0.6213592234998941,
    0.615384615957737, 0.6095238104462624, 0.6037735845893621, 0.5981308408081532,
    0.5925925932824612, 0.5871559642255306, 0.5818181820213795, 0.5765765774995089,
    0.5714285708963871, 0.5663716811686754, 0.5614035092294216, 0.5565217398107052,
    0.551724137738347, 0.5470085479319096, 0.5423728805035353, 0.5378151256591082,
    0.533333333209157, 0.5289256200194359, 0.5245901644229889, 0.5203252024948597,
    0.5161290317773819, 0.5120000001043081, 0.5079365074634552, 0.5039370078593493,
    0.5
  );
  tl10: array[0..64] of Double = (
    -1.5987211554602254e-14, 0.006733382254484161, 0.01336396196243377, 0.019894828666364744,
    0.026328938823450224, 0.032669116513199134, 0.03891806625786707, 0.04507837474781176,
    0.051152522851833554, 0.05714288568152596, 0.06305174551955964, 0.06888128931931946,
    0.07463361794297847, 0.08031075159304697, 0.0859146289593699, 0.09144711736973643,
    0.09691001260357217, 0.10230504473043137, 0.10763387799534528, 0.1128981188345579,
    0.11809931248244676, 0.12323895183150664, 0.12831847779052916, 0.13333927915294033,
    0.1383026985707337, 0.14321003271156815, 0.1480625351773498, 0.15286141776840695,
    0.15760785341221065, 0.16230297497450022, 0.1669478791102102, 0.17154363071088186,
    0.17609125946013351, 0.18059176014330572, 0.1850461019361051, 0.18945522101811502,
    0.1938200256116286, 0.1981414002916851, 0.20242019787913146, 0.20665725058223308,
    0.21085336491040893, 0.215009324428774, 0.21912589150838047, 0.2232038039288199,
    0.22724378099746123, 0.2312465232741803, 0.2352127110226463, 0.23914300410757447,
    0.24303804909074672, 0.24689846968911106, 0.2507248769986597, 0.25451786583884395,
    0.2582780153946909, 0.2620058870291598, 0.2657020340047624, 0.2693669877246184,
    0.2730012721648387, 0.2766053961808713, 0.2801798562863768, 0.2837251380874764,
    0.28724171158280015, 0.29073003893567584, 0.294190571538128, 0.29762374698469335,
    0.3010299956639652
  );
  // st lookup table for exact powers of 10 (uint32 bits of float)
  st_u: array[0..15] of LongWord = (
    $501502F9, $41200000, $42C80000, $00000000,
    $447A0000, $00000000, $461C4000, $47C35000,
    $00000000, $49742400, $00000000, $4B189680,
    $4CBEBC20, $00000000, $4E6E6B28, $3F800000
  );
  b10: array[0..2] of Double = (
    0.4342944825305097, -0.2171537639402152, 0.14474655973900713
  );
  c10: array[0..6] of Double = (
    0.4342944819032518, -0.2171472409516272, 0.14476482730105739, -0.10857362030408772,
    0.08685889777743865, -0.07238812530018697, 0.062026410488936715
  );
  ln10:  Double = 0.3010299956639812;    // 0x1.34413509f79ffp-2
  ln10h: Double = 0.30102999566398125;   // 0x1.34413509f7ap-2
  ln10l: Double = -5.8314879359043e-17;  // -0x1.0cee0ed4ca7e9p-54
var
  t: Tb32u32;
  tz: Tb64u64;
  ux: LongWord;
  m64: Int64;
  j: Int64;
  e, n: Integer;
  z, z2, r, f, el, dr: Double;
  ub, lb: Single;
  je, st_u_val: LongWord;
  je_idx: Integer;
begin
  t.f := x;
  ux := t.u;
  if ux >= $7F800000 then // x <= 0, nan, inf
  begin
    // as_special
    if (ux shl 1) = 0 then // x = +/-0
    begin
      pcr_feraiseexcept_divbyzero;
      t.u := $FF800000; Result := t.f; Exit; // -inf
    end;
    if ux = $7F800000 then begin Result := x; Exit; end; // +inf
    if (ux shl 1) > $FF000000 then begin Result := x + x; Exit; end; // nan
    // x < 0 -> NaN
    pcr_feraiseexcept_invalid;
    t.u := $FFC00000; Result := t.f; Exit;
  end;
  // check for exact power of 10
  st_u_val := st_u[(ux shr 24) and $F];
  if ux = st_u_val then
  begin
    je := (LongWord(Integer(ux) shr 23) - 126);
    je_idx := Integer((je * $4D104D4) shr 28);
    Result := Single(je_idx); Exit;
  end;
  if ux < $00800000 then // subnormal
  begin
    if ux = 0 then // +0
    begin
      pcr_feraiseexcept_divbyzero;
      t.u := $FF800000; Result := t.f; Exit;
    end;
    n := 23 - Integer(pcr_bsr32(ux));
    ux := ux shl n;
    ux := ux - LongWord(n shl 23);
  end;
  e := sar_i32(Integer(ux), 23) - 127;
  m64 := Int64(ux and LongWord((1 shl 23) - 1));
  j := (m64 + (1 shl (23-7))) shr (23-6);
  tz.u := (UInt64(m64) shl 29) or (UInt64($3FF) shl 52);
  z := tz.f * tr10[j] - 1.0;
  z2 := z * z;
  r := ((e * ln10 + tl10[j]) + z * b10[0]) + z2 * (b10[1] + z * b10[2]);
  ub := Single(r);
  lb := Single(r + 9.802997302799099e-11); // 0x1.af23fp-34
  if ub <> lb then
  begin
    f := z * ((c10[0] + z*c10[1]) + z2*((c10[2] + z*c10[3]) + z2*(c10[4] + z*c10[5] + z2*c10[6])));
    f := f + ln10l * e;
    f := f + tl10[j] - tl10[0];
    el := e * ln10h;
    r := el + f;
    ub := Single(r);
    tz.f := r;
    if (tz.u and ($FFFFFFF)) = 0 then
    begin
      dr := (el - r) + f;
      r := r + dr * 32.0;
      ub := Single(r);
    end;
  end;
  Result := ub;
end;

// ── 2.08 erfcf ───────────────────────────────────────────────────────────────
function pcr_erfcf(x: Single): Single;
const
  E_tbl: array[0..127] of Double = (
    1.0, 1.0054299011128027, 1.0108892860517005, 1.016378314910953,
    1.0218971486541166, 1.0274459491187637, 1.0330248790212284, 1.0386341019613787,
    1.0442737824274138, 1.0499440858006872, 1.0556451783605572, 1.061377227289262,
    1.0671404006768237, 1.0729348675259756, 1.0787607977571199, 1.0846183622133092,
    1.0905077326652577, 1.0964290818163769, 1.102382583307841, 1.1083684117236787,
    1.1143867425958924, 1.1204377524096067, 1.1265216186082418, 1.1326385195987192,
    1.1387886347566916, 1.1449721444318042, 1.1511892299529827, 1.1574400736337511,
    1.1637248587775775, 1.1700437696832502, 1.1763969916502812, 1.182784710984341,
    1.189207115002721, 1.1956643920398273, 1.202156731452703, 1.2086843236265816,
    1.215247359980469, 1.2218460329727576, 1.22848053610687, 1.2351510639369334,
    1.241857812073484, 1.2486009771892048, 1.255380757024691, 1.2621973503942507,
    1.2690509571917332, 1.275941778396392, 1.2828700160787783, 1.2898358734066657,
    1.2968395546510096, 1.3038812651919358, 1.3109612115247644, 1.318079601266064,
    1.3252366431597413, 1.3324325470831615, 1.339667524053303, 1.3469417862329458,
    1.3542555469368927, 1.3616090206382248, 1.3690024229745905, 1.3764359707545302,
    1.383909881963832, 1.3914243757719262, 1.3989796725383112, 1.4065759938190154,
    1.4142135623730951, 1.4218926021691656, 1.42961333839197, 1.4373759974489824,
    1.4451808069770467, 1.4530279958490526, 1.460917794180647, 1.4688504333369818,
    1.4768261459394993, 1.4848451658727524, 1.4929077282912648, 1.5010140696264256,
    1.5091644275934228, 1.5173590411982147, 1.5255981507445384, 1.533881997840956,
    1.5422108254079407, 1.550584877685, 1.559004400237837, 1.567469639965553,
    1.5759808451078865, 1.5845382652524937, 1.593142151342267, 1.6017927556826934,
    1.6104903319492543, 1.6192351351948637, 1.6280274218573478, 1.6368674497669644,
    1.645755478153965, 1.6546917676561943, 1.6636765803267364, 1.6727101796415966,
    1.681792830507429, 1.6909247992693053, 1.7001063537185235, 1.709337763100463,
    1.718619298122478, 1.7279512309618377, 1.7373338352737062, 1.746767386199169,
    1.7562521603732995, 1.7657884359332727, 1.7753764925265212, 1.785016611318935,
    1.7947090750031072, 1.804454167806624, 1.8142521755003989, 1.8241033854070534,
    1.8340080864093424, 1.843966568958626, 1.8539791250833855, 1.864046048397789,
    1.8741676341103, 1.8843441790323345, 1.8945759815869656, 1.9048633418176741,
    1.9152065613971474, 1.925605943636125, 1.9360617934922943, 1.9465744175792332,
    1.9571441241754002, 1.9677712232331759, 1.978456026387951, 1.9891988469672663
  );
  ch_e: array[0..3] of Double = (
    -0.4999999999998181, 0.16666666666681407, -0.04166669845578799, 0.008333328785338493
  );
  ct0: array[0..15] of Double = (
    0.8777023949849978, 3.7, 0.4634594459136497, -1.4411533473283251,
    2.4804529237471646, -3.2268872806885516, 3.1822789433641323, -2.3097311555393483,
    1.120893497742982, -0.24054063081368343, -0.0923431569772441, 0.07445755843308839,
    0.001303398643906896, -0.016107485379421634, 0.0011106794945719596, 0.003590883161041724
  );
  ct1: array[0..15] of Double = (
    4.304476145626969, 2.95, 0.1277870732084466, -0.2050575647371633,
    0.11677187850253941, -0.05604779654508989, 0.021876078119428374, -0.006440300708301159,
    0.0011412012141349293, 3.8348710612507015e-05, -8.722199940010998e-05, 1.5843323715748103e-05,
    4.583168277792045e-06, -2.13078125377227e-06, -2.1255636530549933e-07, 2.2754647679221864e-07
  );
  c_sm: array[0..4] of Double = (  // small |x| polynomial
    1.1283791670955126, -0.37612638903148427, 0.11283791635934358,
    -0.02686604912025618, 0.005206760160490499
  );
  iln2_e: Double = 1.4426950408889634;  // 0x1.71547652b82fep+0
  ln2h_e: Double = 0.005415212348111709; // 0x1.62e42fefap-8
  ln2l_e: Double = 1.2864023133262396e-14; // 0x1.cf79abd6f5dc8p-47
var
  t: Tb32u32;
  jt: Tb64u64;
  S: Tb64u64;
  at_u: LongWord;
  sgn: LongWord;
  i: Int64;
  ax_f: Single;
  axd, x2, d, d2, e0, f, z, z2, z4, z8, s_val, r, y: Double;
  j_idx: Int64;
  c8, c9, c10, c11, c12: Double;
  c0, c1, c2, c3, c4, c5, c6, c7: Double;
begin
  ax_f := pcr_fabsf(x);
  t.f := x;
  at_u := t.u and $7FFFFFFF;
  sgn := t.u shr 31;
  if at_u > $40051000 then i := 1 else i := 0; // i selects polynomial set
  // x < -0x1.ea8f94p+1 = -3.8325... => erfc rounds to 2
  if t.u > $C07547CA then
  begin
    if t.u >= $FF800000 then // -Inf or NaN
    begin
      if t.u = $FF800000 then begin Result := 2.0; Exit; end; // -Inf
      Result := x + x; Exit; // NaN
    end;
    Result := 2.0 - Single(2.9802322387695312e-08); Exit; // 2 - 2^-25
  end;
  // |x| >= 0x1.41bbf8p+3 = 10.054... => underflow
  if at_u >= $4120DDFC then
  begin
    if at_u >= $7F800000 then
    begin
      if at_u = $7F800000 then begin Result := 0.0; Exit; end; // +Inf
      Result := x + x; Exit; // NaN
    end;
    // 0x1p-149 * 0.25 rounds to 0 or 2^-149
    Result := Single(1.401298464324817e-45) * 0.25; Exit;
  end;
  // small |x| <= 0x1.7p-4 = 0.08984375
  if at_u <= $3DB80000 then
  begin
    if t.u = $B76C9F62 then // x = -0x1.d93ec4p-17
      begin Result := Single(1.0000158548355103) + Single(2.9802322387695312e-08); Exit; end;
    if at_u <= $32E2DFC4 then // |x| <= 0x1.c5bf88p-26
    begin
      if at_u = 0 then begin Result := 1.0; Exit; end;
      if sgn <> 0 then Result := 1.0 + Single(2.9802322387695312e-08)   // 1 + 2^-25 (sgn=1 => x<0)
      else Result := 1.0 - Single(2.9802322387695312e-08);               // 1 - 2^-25 (sgn=0 => x>0)
      Exit;
    end;
    axd := ax_f;
    x2 := axd * axd;
    f := Double(x) * (c_sm[0] + x2*(c_sm[1] + x2*(c_sm[2] + x2*(c_sm[3] + x2*c_sm[4]))));
    Result := Single(1.0 - f); Exit;
  end;
  // main path: -3.8325... <= x <= 10.054..., |x| > 0.0898...
  axd := ax_f;
  x2 := axd * axd;
  jt.f := x2 * iln2_e - 1024.00390625; // 0x1.00004p+10
  j_idx := sar_i64(Int64(jt.u shl 12), 48); // sign-extend 16-bit field
  S.u := UInt64(sar_i64(j_idx, 7) + (Int64($3FF) or Int64(sgn shl 11))) shl 52;
  d := (x2 + ln2h_e * j_idx) + ln2l_e * j_idx;
  d2 := d * d;
  e0 := E_tbl[j_idx and 127];
  f := d + d2 * ((ch_e[0] + d*ch_e[1]) + d2*(ch_e[2] + d*ch_e[3]));
  // select polynomial set
  if i = 0 then
  begin
    z := (axd - ct0[0]) / (axd + ct0[1]);
    c0 := ct0[3]; c1 := ct0[4]; c2 := ct0[5]; c3 := ct0[6];
    c4 := ct0[7]; c5 := ct0[8]; c6 := ct0[9]; c7 := ct0[10];
    c8 := ct0[11]; c9 := ct0[12]; c10 := ct0[13]; c11 := ct0[14]; c12 := ct0[15];
    s_val := ct0[2];
  end
  else
  begin
    z := (axd - ct1[0]) / (axd + ct1[1]);
    c0 := ct1[3]; c1 := ct1[4]; c2 := ct1[5]; c3 := ct1[6];
    c4 := ct1[7]; c5 := ct1[8]; c6 := ct1[9]; c7 := ct1[10];
    c8 := ct1[11]; c9 := ct1[12]; c10 := ct1[13]; c11 := ct1[14]; c12 := ct1[15];
    s_val := ct1[2];
  end;
  z2 := z * z; z4 := z2 * z2; z8 := z4 * z4;
  r := (((c0 + z*c1) + z2*(c2 + z*c3)) + z4*((c4 + z*c5) + z2*(c6 + z*c7))) +
       z8*(((c8 + z*c9) + z2*(c10 + z*c11)) + z4*c12);
  r := s_val + z * r;
  y := S.f * (e0 - f * e0) * r;
  if sgn <> 0 then y := 2.0 + y   // off[1] + r = 2 + r
  else y := 0.0 + y;               // off[0] + r = r
  Result := Single(y);
end;
// ── 3.01 log2p1f ────────────────────────────────────────────────────────────
function pcr_log2p1f(x: Single): Single;
const
  // reciprocal 1/(1+j/64) rounded to 24 bits (65 entries)
  ix_l2p1: array[0..64] of Double = (
    1, 0.98461538553237915, 0.96969699859619141, 0.95522385835647583,
    0.94117647409439087, 0.9275362491607666, 0.91428571939468384, 0.90140843391418457,
    0.8888888955116272, 0.87671232223510742, 0.86486488580703735, 0.85333335399627686,
    0.84210526943206787, 0.83116883039474487, 0.82051283121109009, 0.81012660264968872,
    0.80000001192092896, 0.79012346267700195, 0.7804877758026123, 0.77108430862426758,
    0.76190477609634399, 0.75294119119644165, 0.74418604373931885, 0.73563218116760254,
    0.72727274894714355, 0.71910113096237183, 0.71111112833023071, 0.7032967209815979,
    0.69565218687057495, 0.68817204236984253, 0.6808510422706604, 0.67368423938751221,
    0.66666668653488159, 0.65979379415512085, 0.65306121110916138, 0.64646464586257935,
    0.63999998569488525, 0.6336633563041687, 0.62745100259780884, 0.62135922908782959,
    0.61538463830947876, 0.60952383279800415, 0.60377359390258789, 0.59813082218170166,
    0.59259259700775146, 0.58715593814849854, 0.58181816339492798, 0.5765765905380249,
    0.57142859697341919, 0.56637167930603027, 0.56140351295471191, 0.5565217137336731,
    0.5517241358757019, 0.54700857400894165, 0.54237288236618042, 0.5378151535987854,
    0.53333336114883423, 0.52892559766769409, 0.5245901346206665, 0.5203251838684082,
    0.5161290168762207, 0.51200002431869507, 0.50793653726577759, 0.5039370059967041,
    0.5);
  // log of reciprocal biased by 0x1.dp-45
  lix_l2p1: array[0..64] of Double = (
    5.1514348342607263e-14, -0.022367811684788536, -0.04439407636273985, -0.066089224048082804,
    -0.087462835875830064, -0.1085244299058286, -0.12928300888322822, -0.14974714637691999,
    -0.16992499069334521, -0.18982456962888145, -0.20945333069492272, -0.22881865556185382,
    -0.24792750269461833, -0.2667865420384643, -0.2854022000515945, -0.30378071189946143,
    -0.32192807338947965, -0.33984999213565759, -0.35755205836261156, -0.37503948509145263,
    -0.39231739590641984, -0.40939090792174687, -0.42626476007650427, -0.44294350122313464,
    -0.45943157564158366, -0.47573341618658727, -0.49185306139564766, -0.5077946039210548,
    -0.52356192918467237, -0.53915881245159436, -0.55458889736047823, -0.56985554652463255,
    -0.58496245772544264, -0.59991288652635399, -0.61470987367467511, -0.62935662142317261,
    -0.64385622202142057, -0.65821150559318919, -0.67242529091409553, -0.68650051374702237,
    -0.70043966439646321, -0.71424546257787935, -0.72792043306531651, -0.74146703208398779,
    -0.7548874914145014, -0.76818438658314137, -0.78135975920750045, -0.79441583141607897,
    -0.80735485756406, -0.82017896778959409, -0.83289000341577457, -0.84549011678143371,
    -0.85798100050197845, -0.8703646483717884, -0.88264304667456084, -0.89481768940909856,
    -0.90689052036605922, -0.918863297737195, -0.93073741817970457, -0.9425145590837678,
    -0.95419635338248709, -0.96578421613769971, -0.97727984019577085, -0.98868469214657217,
    -0.99999999999994849);
  b_l2p1: array[0..2] of Double = (
    1.4426950429726688, -0.72136918934353911, 0.48083766343529977);
  c_l2p1: array[0..5] of Double = (
    1.4426950408889683, -0.72134752044437966, 0.48089834631236422,
    -0.36067376567497511, 0.28856066993666346, -0.24038686635219694);
  g_l2p1: array[0..7] of Double = (
    1.9259629797116934e-08, -0.72134752044448125,
    0.48089834696963674, -0.36067376023288011,
    0.2885389476641852, -0.2404491020710918,
    0.2062755131304331, -0.18051311004316711);
var
  tv: Tb32u32;
  tp: Tb64u64;
  xd: Tb64u64;
  ux, ax: UInt32;
  zz: Double;
  e_exp: Integer;
  m64: UInt64;
  j32: Integer;
  dd, d2, el_v: Double;
  f_v, lb_v: Double;
  ub_v: Single;
  c0_v, c2_v, c4_v: Double;
  z2_v, z4_v: Double;
begin
  tv.f := x;
  ux := tv.u;
  zz := x;
  if ux >= $BF800000 then begin  // x <= -1
    if ux = $BF800000 then begin Result := Single(-1.0/0.0); Exit; end;
    ax := ux shl 1;
    if ax > $FF000000 then begin Result := x + x; Exit; end;
    Result := Single(0.0/0.0); Exit;
  end;
  ax := ux and ($7FFFFFFF);
  if ax >= $7F800000 then begin  // +inf or +nan
    if ax > $7F800000 then begin Result := x + x; Exit; end;
    Result := x; Exit;
  end;
  if ax < $3CC00000 then begin  // |x| < 0.0234375
    if ax <= $0058B90B then begin  // |x| <= 0x1p-126*ln(2) approx
      if ax = 0 then begin Result := x; Exit; end;
      Result := Single(zz * 1.4426950408889634); Exit;
    end else begin
      z2_v := zz*zz; z4_v := z2_v*z2_v;
      f_v := zz*((g_l2p1[0] + zz*g_l2p1[1]) + z2_v*(g_l2p1[2] + zz*g_l2p1[3]) +
                 z4_v*((g_l2p1[4] + zz*g_l2p1[5]) + z2_v*(g_l2p1[6] + zz*g_l2p1[7])));
      f_v := f_v + zz * 1.4426950216293335;
      Result := Single(f_v); Exit;
    end;
  end;
  tp.f := zz + 1.0;
  e_exp := Integer((tp.u shr 52) - $3FF);
  m64 := tp.u and (not UInt64(0) shr 12);
  if m64 = 0 then begin Result := Single(e_exp); Exit; end;
  j32 := Integer((m64 + (UInt64(1) shl (52-7))) shr (52-6));
  xd.u := m64 or (UInt64($3FF) shl 52);
  dd := xd.f * ix_l2p1[j32] - 1.0;
  d2 := dd * dd;
  el_v := Double(e_exp) - lix_l2p1[j32];
  f_v := (el_v + dd*b_l2p1[0]) + d2*(b_l2p1[1] + dd*b_l2p1[2]);
  lb_v := f_v;
  ub_v := Single(f_v + 3.2565594665356912e-10);
  if Single(lb_v) = ub_v then begin Result := Single(lb_v); Exit; end;
  // check two hard cases
  if ux = $4EBD09E3 then begin Result := Single(30.562536239624023) + Single(4.76837158203125e-07); Exit; end;
  if ux = $BD6D142E then begin Result := Single(-0.086018145084381104) + Single(1.862645149230957e-09); Exit; end;
  c0_v := c_l2p1[0] + dd*c_l2p1[1];
  c2_v := c_l2p1[2] + dd*c_l2p1[3];
  c4_v := c_l2p1[4] + dd*c_l2p1[5];
  c0_v := c0_v + d2*(c2_v + d2*c4_v);
  f_v := Double(e_exp) + (5.1514348342607263e-14 - lix_l2p1[j32]) + dd*c0_v;
  Result := Single(f_v);
end;

// ── 3.02 erff ──────────────────────────────────────────────────────────────
function pcr_erff(x: Single): Single;
const
  C_erf: array[0..55, 0..7] of Double = (
  (0.49261347321793791, 0.90579313677511597, -0.42459053286087173, -0.16924650407296887, 0.18119731439766232, 0.016799455750244524, -0.050923531014944678, 0.0028189848078158763),
  (0.54752844539954459, 0.85091390489333729, -0.45204801197240407, -0.12353763072512959, 0.18349734265047696, -0.001931893142801651, -0.04857222055475334, 0.0078289045929204205),
  (0.59891738659435068, 0.79313897153445778, -0.47092626434679835, -0.077971344208993376, 0.1801231551170506, -0.01938784142771208, -0.044180683107513059, 0.012105413837424785),
  (0.64663270800670813, 0.73353363653916459, -0.48138144897751484, -0.033906828253933478, 0.1715861542963612, -0.03486931103952922, -0.038117640929976632, 0.015442001235256559),
  (0.69059246870032676, 0.67312830249850619, -0.4838109674200044, 0.0074499877204104309, 0.15859297905578715, -0.04783047088142809, -0.030825391703067067, 0.017710052746510808),
  (0.73077729241084721, 0.61289029396304384, -0.47882054215834935, 0.045088934384124937, 0.14199398096472582, -0.057899791962494104, -0.022784650917334715, 0.018862586408709678),
  (0.76722566123234159, 0.55370023825077219, -0.46718457602429592, 0.078224577927713312, 0.12272719925514923, -0.064887795654489863, -0.014479299072080593, 0.018931025659311967),
  (0.80002789416423326, 0.49633368578049469, -0.44980240273920602, 0.10631105639249756, 0.10176194005833371, -0.068782013024985961, -0.0063639259802374245, 0.018015881670963955),
  (0.82931915059331518, 0.44144831841852705, -0.42765305846892365, 0.12904316078625924, 0.080045743478965578, -0.069730667314353956, 0.0011634788440777004, 0.016272717737345847),
  (0.85527181044291711, 0.38957677365552829, -0.40175104783348586, 0.14634492083187892, 0.058457922398875678, -0.068017363620575269, 0.0077819357929112525, 0.013895038628831281),
  (0.87808757515181002, 0.34112482093974111, -0.37310527290420986, 0.15834765451127528, 0.037772057744762118, -0.064029567074999405, 0.013260088620566789, 0.011095990722926325),
  (0.89798960920701787, 0.29637438049741344, -0.34268287745155002, 0.16535992453533385, 0.018628926693587164, -0.058223852582126694, 0.017460853750747827, 0.0080905546235608729),
  (0.91521500387173227, 0.25549068578278905, -0.31137927329914933, 0.16783209762723653, 0.0015204136487839254, -0.051090828748879519, 0.020338670246721237, 0.00507974412321688),
  (0.93000779678611956, 0.21853276425487753, -0.27999510420282292, 0.1663182300879168, -0.013215908292348482, -0.043122314960366406, 0.021930485566038334, 0.0022379234421607075),
  (0.94261272724427847, 0.18546634817174804, -0.24922040535687803, 0.16143783040824267, -0.025392568264442042, -0.034782843630157233, 0.022342010346704643, -0.00029606380825509344),
  (0.95326985102802309, 0.15617832355458153, -0.21962576749951612, 0.15383971584616857, -0.034959956502336915, -0.0264869400348791, 0.021730977365108452, -0.002423714217684135),
  (0.96221008419949117, 0.13049187367280624, -0.19165993945759699, 0.14416973282841652, -0.041987997666416713, -0.018582972719691608, 0.020289153123537188, -0.0040878050403306979),
  (0.96965169513984228, 0.10818156296067474, -0.16565301828397275, 0.13304360184528274, -0.046643832650907782, -0.011343735055122334, 0.018224699511709484, -0.0052699597172998196),
  (0.9757977205900914, 0.088987726356669475, -0.14182418888117229, 0.12102562523474721, -0.049167564302606825, -0.0049633676443465342, 0.015746207347497786, -0.0059855767905814643),
  (0.980834245966551, 0.072629665503254739, -0.12029288348981136, 0.10861350368635604, -0.04984792967583504, 0.00044019984475663525, 0.01304936946272983, -0.0062769854965132564),
  (0.9849294635082112, 0.058817295576895752, -0.10109222677268306, 0.096229077985516318, -0.048999455515149566, 0.0048183998057526531, 0.010306874047771723, -0.0062056917467217644),
  (0.98823340391229519, 0.047261024719377487, -0.084183700281167465, 0.084214469178176918, -0.046942279330100319, 0.0081820309209528885, 0.0076617369156954471, -0.0058445784272914587),
  (0.99087822750686272, 0.037679774086816623, -0.069472083472263149, 0.072832844573072042, -0.043985418995012925, 0.010589390996004632, 0.0052239404602482855, -0.0052706412849319812),
  (0.99297895874393471, 0.029807154664801713, -0.05681988857942559, 0.062272890182024226, -0.040413887400609159, 0.012133720215707683, 0.0030700270046512936, -0.0045588405007917481),
  (0.99463455163389425, 0.023395903761608296, -0.046060685530295575, 0.05265601517182094, -0.036479704989582112, 0.012930961556340355, 0.0012450878786196438, -0.0037773205965664566),
  (0.99592918229475302, 0.018220748170203777, -0.037010894720361932, 0.044045337210834937, -0.032396582564123359, 0.013108620909804976, -0.00023348258575818946, -0.0029841309437353475),
  (0.99693367664795451, 0.014079902896030661, -0.029479796688224556, 0.036455581912286546, -0.028337840157266567, 0.012796265600024895, -0.0013711206608561846, -0.0022254020682695714),
  (0.99770699512421801, 0.010795436020481278, -0.023277658918860803, 0.029863156022595917, -0.024436996993197537, 0.01211796242092659, -0.0021907201087155168, -0.0015348170953254647),
  (0.99829771087140551, 0.0082127346369952215, -0.018222004975576418, 0.024215804147522949, -0.02079040738204356, 0.011186744874436792, -0.0027272679424444036, -0.0009341362483086961),
  (0.99874543240701885, 0.0061992973446308129, -0.014142147067230512, 0.019441416216616599, -0.017461317417300278, 0.010101027053041023, -0.0030228888430427073, -0.00043449277633907422),
  (0.99908213518143396, 0.004643059175781691, -0.010882169943077325, 0.015455704144156798, -0.014484764137264841, 0.0089427550774859709, -0.0031225730415550134, -3.8175435227474089e-05),
  (0.99933337859119276, 0.0034504285962259303, -0.008302593809551756, 0.012168601370730614, -0.011872817852758623, 0.0077770068423028979, -0.0030707423369344036, 0.00025936263705818297),
  (0.99951939529124134, 0.0025441865040150146, -0.0062809604317091167, 0.009489351875802042, -0.0096197656435485761, 0.0066527131726320902, -0.0029087038558813071, 0.00046748657616679741),
  (0.99965604806931274, 0.0018613666092925042, -0.0047115842297264727, 0.007330342851164416, -0.0077069373257350306, 0.0056041714098085246, -0.0026729573762086668, 0.00059851791323776017),
  (0.99975565607993799, 0.0013512072527286062, -0.0035046938117459419, 0.005609797298487618, -0.0061069746975710324, 0.0046530473160883998, -0.0023942621671011019, 0.00066613791814829721),
  (0.99982769701468499, 0.00097323807446023569, -0.0025851636352860526, 0.004253481245918434, -0.0047874335624023094, 0.0038106040477139265, -0.0020973327307732058, 0.00068411850479664117),
  (0.99987939500763068, 0.00069554188547860368, -0.0018910045011603846, 0.0031955983631085285, -0.0037136816121056512, 0.003079949511136694, -0.0018010186419316083, 0.00066536561528653744),
  (0.99991620598062658, 0.00049321305113704235, -0.0013717487984997085, 0.0023790465467580144, -0.0028511118775535536, 0.0024581482524300242, -0.0015188199114472608, 0.00062134644950509856),
  (0.99994221297601882, 0.00034701872224903227, -0.00098683449142581868, 0.0017552008158370879, -0.002166731175259081, 0.0019380966981796372, -0.0012596119625042153, 0.00056171077626291158),
  (0.99996044405336681, 0.00024225811122990498, -0.00070406263579420188, 0.0012833686530412333, -0.0016302073633598833, 0.0015101056693959851, -0.0010284637718540465, 0.00049420242397801438),
  (0.99997312476875266, 0.00016780728923366131, -0.00049817788994454171, 0.00093004131071250547, -0.0012144706095426547, 0.0011631716338346787, -0.00082746577951179538, 0.00042471469692076022),
  (0.99998187630799307, 0.00011533215079327046, -0.00034960058212242002, 0.00066804045939168133, -0.00089596513866240348, 0.00088594574414336289, -0.00065650463138913176, 0.00035746461922862524),
  (0.99998786917170335, 7.8649693548360708e-05, -0.00024332248944276653, 0.00047563606991809409, -0.00065464190824693531, 0.00066742866651560353, -0.00051394552630094003, 0.0002952337204821917),
  (0.99999194103908062, 5.3217044645068047e-05, -0.00016796629718520607, 0.000335690068744768, -0.00047377200005840116, 0.00049743023257662819, -0.00039720193370637141, 0.00023963641078976976),
  (0.99999468616400899, 3.5728233590487455e-05, -0.00011500025189011397, 0.00023486196262317391, -0.00033964744765398142, 0.00036683758307549138, -0.0003031874090858513, 0.00019138655591156029),
  (0.99999652244808557, 2.3800134528024529e-05, -7.8094191437426663e-05, 0.00016289766557146129, -0.00024122249662388096, 0.00026773529294131634, -0.00022865497713078851, 0.0001505422894208187),
  (0.9999977412324873, 1.5730928349533599e-05, -5.260029168296787e-05, 0.00011201117407626997, -0.00016973517832097345, 0.00019341750017729041, -0.00017043644158216189, 0.00011671731569916309),
  (0.99999854387686993, 1.031659468378863e-05, -3.5140900653086202e-05, 7.6360263658514977e-05, -0.00011833738196773359, 0.00013832664342989577, -0.00012559762473224564, 8.9253526770131036e-05),
  (0.99999906835661023, 6.7131361956839387e-06, -2.3286191187812652e-05, 5.1611605046464615e-05, -8.1751767521628702e-05, 9.7947126025213962e-05, -9.1526669002494113e-05, 6.7354663319306958e-05),
  (0.99999940840719481, 4.3343264924529069e-06, -1.5305590433498795e-05, 3.4587135307086044e-05, -5.5966011335462916e-05, 6.867587301443768e-05, -6.5971875888997571e-05, 5.0184038564375413e-05),
  (0.99999962716703439, 2.776673790460282e-06, -9.9786714398456092e-06, 2.2981675722138639e-05, -3.7968947211198758e-05, 4.768587565057889e-05, -4.7043770337430147e-05, 3.6931273166747702e-05),
  (0.9999997668041688, 1.7649612458156294e-06, -6.4531395590728805e-06, 1.5141207255800815e-05, -2.5528952211500767e-05, 3.2793742702145623e-05, -3.3193703352104246e-05, 2.6853795581592947e-05),
  (0.99999985524310697, 1.1131471046973298e-06, -4.1395157986145951e-06, 9.8915005461396886e-06, -1.7012154759834782e-05, 2.2338138923730696e-05, -2.3178752814561998e-05, 1.9298911294543391e-05),
  (0.99999991082014172, 6.9658965845718095e-07, -2.6339796482601615e-06, 6.4076271416638353e-06, -1.1236415488351412e-05, 1.5072796640709194e-05, -1.6020238500711028e-05, 1.3711750607540815e-05),
  (0.99999994547438986, 4.3252235464692901e-07, -1.6625078022828812e-06, 4.1160021235915782e-06, -7.3562640776110176e-06, 1.0075459881387709e-05, -1.0961006204555806e-05, 9.6336299933705873e-06),
  (0.99999996691449189, 2.6646929256788535e-07, -1.0408956752451789e-06, 2.6218427223451289e-06, -4.7738154447972423e-06, 6.6725369652585968e-06, -7.4248432503650242e-06, 6.6944788653075042e-06)
  );
  c_erf_small: array[0..7] of Double = (
    1.1283791670955126, -0.3761263890317818,
    0.1128379167034242, -0.026866170388309935,
    0.0052239723351509325, -0.00085477344060515487,
    0.00012018447509482211, -1.3721145267025539e-05);
var
  tv_e: Tb32u32;
  ux_e: UInt32;
  ss, zz_e, vf: Double;
  ax_f: Single;
  i_erf: Integer;
  z2_e, z4_e, z8_e, c0_e, c2_e, c4_e, c6_e: Double;
  row_e: Integer;
begin
  ax_f := pcr_fabsf(x);
  tv_e.f := ax_f;
  ux_e := tv_e.u;
  ss := x;
  zz_e := ax_f;
  if ux_e > $407AD444 then begin  // |x| > 3.919...
    if ux_e > ($FF shl 23) then begin Result := x + x; Exit; end;  // nan
    if ux_e = ($FF shl 23) then begin Result := pcr_copysignf(1.0, x); Exit; end;  // +-inf
    Result := pcr_copysignf(1.0, x) - pcr_copysignf(2.9802322387695312e-08, x); Exit;
  end;
  i_erf := Trunc(16.0 * ax_f);
  if ux_e < $3EE00000 then begin  // |x| < 0.4375
    z2_e := ss*ss; z4_e := z2_e*z2_e; z8_e := z4_e*z4_e;
    c0_e := c_erf_small[0] + z2_e*c_erf_small[1];
    c2_e := c_erf_small[2] + z2_e*c_erf_small[3];
    c4_e := c_erf_small[4] + z2_e*c_erf_small[5];
    c6_e := c_erf_small[6] + z2_e*c_erf_small[7];
    c0_e := c0_e + z4_e*c2_e;
    c4_e := c4_e + z4_e*c6_e;
    c0_e := c0_e + z8_e*c4_e;
    Result := Single(ss * c0_e); Exit;
  end;
  zz_e := (zz_e - 0.03125) - 0.0625 * Floor(16.0 * Double(ax_f));
  row_e := i_erf - 7;
  z2_e := zz_e*zz_e; z4_e := z2_e*z2_e;
  c0_e := C_erf[row_e,0] + zz_e*C_erf[row_e,1];
  c2_e := C_erf[row_e,2] + zz_e*C_erf[row_e,3];
  c4_e := C_erf[row_e,4] + zz_e*C_erf[row_e,5];
  c6_e := C_erf[row_e,6] + zz_e*C_erf[row_e,7];
  c0_e := c0_e + z2_e*c2_e;
  c4_e := c4_e + z2_e*c6_e;
  c0_e := c0_e + z4_e*c4_e;
  Result := Single(pcr_copysign(c0_e, ss));
end;

// ── 3.03 sinhf ────────────────────────────────────────────────────────────
function pcr_sinhf(x: Single): Single;
const
  c_sinh: array[0..3] of Double = (1, 0.021660849391257477, 0.0002345984913513542, 1.6938658699950235e-06);
  ch_sinh: array[0..6] of Double = (
    1, 0.02166084939249829, 0.0002345961982022468, 1.6938509724129055e-06,
    9.1725627017026289e-09, 3.973729405780548e-11, 1.4345723178374038e-13);
  cp_sinh: array[0..3] of Double = (0.16666666666666666, 0.0083333333333572308, 0.00019841269076590929, 2.7565149135114762e-06);
  tb_sinh: array[0..31] of UInt64 = (
    UInt64($3FE0000000000000), UInt64($3FE059B0D3158574), UInt64($3FE0B5586CF9890F), UInt64($3FE11301D0125B51),
    UInt64($3FE172B83C7D517B), UInt64($3FE1D4873168B9AA), UInt64($3FE2387A6E756238), UInt64($3FE29E9DF51FDEE1),
    UInt64($3FE306FE0A31B715), UInt64($3FE371A7373AA9CB), UInt64($3FE3DEA64C123422), UInt64($3FE44E086061892D),
    UInt64($3FE4BFDAD5362A27), UInt64($3FE5342B569D4F82), UInt64($3FE5AB07DD485429), UInt64($3FE6247EB03A5585),
    UInt64($3FE6A09E667F3BCD), UInt64($3FE71F75E8EC5F74), UInt64($3FE7A11473EB0187), UInt64($3FE82589994CCE13),
    UInt64($3FE8ACE5422AA0DB), UInt64($3FE93737B0CDC5E5), UInt64($3FE9C49182A3F090), UInt64($3FEA5503B23E255D),
    UInt64($3FEAE89F995AD3AD), UInt64($3FEB7F76F2FB5E47), UInt64($3FEC199BDD85529C), UInt64($3FECB720DCEF9069),
    UInt64($3FED5818DCFBA487), UInt64($3FEDFC97337B9B5F), UInt64($3FEEA4AFA2A490DA), UInt64($3FEF50765B6E4540));
var
  ts: Tb32u32;
  ux_s: UInt32;
  zs: Double;
  a_s, ia_s, h_s, h2_s: Double;
  ja_s: Tb64u64;
  sp_s, sm_s: Tb64u64;
  jp_s, jm_s: Int64;
  te_s, to_s, rp_s, rm_s, r_s: Double;
  ub_s: Single;
  lb_s: Single;
  z2_s, z4_s: Double;
begin
  ts.f := x;
  ux_s := ts.u shl 1;
  zs := x;
  if ux_s > $8565A9F8 then begin  // |x| > 0x1.65a9f8p+6
    if ux_s >= $FF000000 then begin
      if (ux_s shl 8) <> 0 then begin Result := x + x; Exit; end;  // nan
      Result := x; Exit;  // +-inf
    end;
    Result := pcr_copysignf(Single(2.0) * 3.4028234663852886e+38, x); Exit;
  end;
  if ux_s < $7C000000 then begin  // |x| < 0.125
    if ux_s <= $74250BFE then begin  // |x| <= 0x1.250bfep-11
      if ux_s < $66000000 then begin  // |x| < 0x1p-24
        Result := pcr_fmaf(x, pcr_fabsf(x), x); Exit;
      end;
      if ux_s = $74250BFE then begin
        Result := pcr_copysignf(1.0, x)*Single(0.00055894249817356467) + pcr_copysignf(1.0, x)*Single(1.4551915228366852e-11); Exit;
      end;
      Result := Single((x*0.1666666716337204)*(x*x) + x); Exit;
    end;
    z2_s := zs*zs; z4_s := z2_s*z2_s;
    Result := Single(zs + (z2_s*zs)*((cp_sinh[0] + z2_s*cp_sinh[1]) + z4_s*(cp_sinh[2] + z2_s*cp_sinh[3]))); Exit;
  end;
  a_s := 46.166241308446828 * zs;
  ia_s := pcr_roundeven(a_s);
  h_s := a_s - ia_s;
  h2_s := h_s * h_s;
  ja_s.f := ia_s + 6755399441055744;
  jp_s := Int64(ja_s.u);
  jm_s := -jp_s;
  sp_s.u := tb_sinh[jp_s and 31] + (UInt64(jp_s shr 5) shl 52);
  sm_s.u := tb_sinh[jm_s and 31] + (UInt64(jm_s shr 5) shl 52);
  te_s := c_sinh[0] + h2_s*c_sinh[2];
  to_s := (c_sinh[1] + h2_s*c_sinh[3]);
  rp_s := sp_s.f*(te_s + h_s*to_s);
  rm_s := sm_s.f*(te_s - h_s*to_s);
  r_s := rp_s - rm_s;
  ub_s := Single(r_s);
  lb_s := Single(r_s - 1.52e-10*r_s);
  if ub_s <> lb_s then begin
    h_s := (46.16624128818512*zs - ia_s) + 2.026170940661134e-08*zs;
    h2_s := h_s*h_s;
    te_s := ch_sinh[0] + h2_s*ch_sinh[2] + (h2_s*h2_s)*(ch_sinh[4] + h2_s*ch_sinh[6]);
    to_s := ch_sinh[1] + h2_s*(ch_sinh[3] + h2_s*ch_sinh[5]);
    r_s := sp_s.f*(te_s + h_s*to_s) - sm_s.f*(te_s - h_s*to_s);
    ub_s := Single(r_s);
  end;
  Result := ub_s;
end;

// ── 3.04 expf ────────────────────────────────────────────────────────────
function pcr_expf(x: Single): Single;
const
  c_exp: array[0..5] of Double = (
    0.69314718055994529, 0.24022650695910072, 0.055504108664026088,
    0.0096181291075005358, 0.001333362331326638, 0.00015403602972146417);
  b_exp: array[0..3] of Double = (1, 0.69314718052023927, 0.2402288551437867, 0.055504596827996931);
  tb_exp: array[0..63] of UInt64 = (
    UInt64($3FF0000000000000), UInt64($3FF02C9A3E778061), UInt64($3FF059B0D3158574), UInt64($3FF0874518759BC8),
    UInt64($3FF0B5586CF9890F), UInt64($3FF0E3EC32D3D1A2), UInt64($3FF11301D0125B51), UInt64($3FF1429AAEA92DE0),
    UInt64($3FF172B83C7D517B), UInt64($3FF1A35BEB6FCB75), UInt64($3FF1D4873168B9AA), UInt64($3FF2063B88628CD6),
    UInt64($3FF2387A6E756238), UInt64($3FF26B4565E27CDD), UInt64($3FF29E9DF51FDEE1), UInt64($3FF2D285A6E4030B),
    UInt64($3FF306FE0A31B715), UInt64($3FF33C08B26416FF), UInt64($3FF371A7373AA9CB), UInt64($3FF3A7DB34E59FF7),
    UInt64($3FF3DEA64C123422), UInt64($3FF4160A21F72E2A), UInt64($3FF44E086061892D), UInt64($3FF486A2B5C13CD0),
    UInt64($3FF4BFDAD5362A27), UInt64($3FF4F9B2769D2CA7), UInt64($3FF5342B569D4F82), UInt64($3FF56F4736B527DA),
    UInt64($3FF5AB07DD485429), UInt64($3FF5E76F15AD2148), UInt64($3FF6247EB03A5585), UInt64($3FF6623882552225),
    UInt64($3FF6A09E667F3BCD), UInt64($3FF6DFB23C651A2F), UInt64($3FF71F75E8EC5F74), UInt64($3FF75FEB564267C9),
    UInt64($3FF7A11473EB0187), UInt64($3FF7E2F336CF4E62), UInt64($3FF82589994CCE13), UInt64($3FF868D99B4492ED),
    UInt64($3FF8ACE5422AA0DB), UInt64($3FF8F1AE99157736), UInt64($3FF93737B0CDC5E5), UInt64($3FF97D829FDE4E50),
    UInt64($3FF9C49182A3F090), UInt64($3FFA0C667B5DE565), UInt64($3FFA5503B23E255D), UInt64($3FFA9E6B5579FDBF),
    UInt64($3FFAE89F995AD3AD), UInt64($3FFB33A2B84F15FB), UInt64($3FFB7F76F2FB5E47), UInt64($3FFBCC1E904BC1D2),
    UInt64($3FFC199BDD85529C), UInt64($3FFC67F12E57D14B), UInt64($3FFCB720DCEF9069), UInt64($3FFD072D4A07897C),
    UInt64($3FFD5818DCFBA487), UInt64($3FFDA9E603DB3285), UInt64($3FFDFC97337B9B5F), UInt64($3FFE502EE78B3FF6),
    UInt64($3FFEA4AFA2A490DA), UInt64($3FFEFA1BEE615A27), UInt64($3FFF50765B6E4540), UInt64($3FFFA7C1819E90D8));
var
  te: Tb32u32;
  ux_exp: UInt32;
  z_exp: Double;
  a_exp: Double;
  u_exp: Tb64u64;
  sv: Tb64u64;
  ia_exp, h_exp, h2_exp, r_exp: Double;
  ub_exp: Single;
  lb_exp: Single;
  w_exp, s_exp: Double;
begin
  te.f := x;
  ux_exp := te.u shl 1;
  z_exp := x;
  a_exp := 1.4426950408889634 * z_exp;
  u_exp.f := a_exp + 105553116266496;
  if (ux_exp > $8562E42E) or (ux_exp < $6F93813E) then begin
    if ux_exp < $6F93813E then begin  // |x| < 0x1.93813ep-16
      Result := Single(1.0 + z_exp*(1.0 + z_exp*0.5)); Exit;
    end;
    if ux_exp >= ($FF shl 24) then begin
      if ux_exp > ($FF shl 24) then begin Result := x + x; Exit; end;  // nan
      if (te.u shr 31) <> 0 then Result := 0.0 else Result := x; Exit;  // +-inf
    end;
    if te.u > $C2CE8EC0 then begin  // x < -0x1.9d1d8p+6
      Result := Single(pcr_fmax(1.4012984643248171e-45 + (z_exp + 103.27892990343184)*1.0108231726433641e-45, 3.5032461608120427e-46)); Exit;
    end;
    if ((te.u shr 31) = 0) and (te.u > $42B17217) then begin  // x > 0x1.62e42ep+6
      Result := Single(3.4028234663852886e+38 * 3.4028234663852886e+38); Exit;
    end;
  end;
  ia_exp := 105553116266496 - u_exp.f;
  h_exp := a_exp + ia_exp;
  sv.u := tb_exp[u_exp.u and $3F] + ((u_exp.u shr 6) shl 52);
  h2_exp := h_exp * h_exp;
  r_exp := ((b_exp[0] + h_exp*b_exp[1]) + h2_exp*(b_exp[2] + h_exp*b_exp[3])) * sv.f;
  ub_exp := Single(r_exp);
  lb_exp := Single(r_exp - r_exp*1.45e-10);
  if ub_exp <> lb_exp then begin
    h_exp := (1.442695040255785*z_exp + ia_exp) + 6.3317841895660438e-10*z_exp;
    s_exp := sv.f;
    h2_exp := h_exp * h_exp;
    w_exp := s_exp * h_exp;
    r_exp := s_exp + w_exp*((c_exp[0] + h_exp*c_exp[1]) + h2_exp*((c_exp[2] + h_exp*c_exp[3]) + h2_exp*(c_exp[4] + h_exp*c_exp[5])));
    ub_exp := Single(r_exp);
  end;
  Result := ub_exp;
end;

// ── 3.05 atanhf ────────────────────────────────────────────────────────────
function pcr_atanhf(x: Single): Single;
const
  tr_atanh: array[0..63] of Double = (
    0.9922480620443821, 0.97709923610091209, 0.96240601502358913, 0.94814814813435078,
    0.93430656939744949, 0.92086330987513065, 0.90780141763389111, 0.89510489441454411,
    0.88275862112641335, 0.87074830010533333, 0.85906040295958519, 0.84768211841583252,
    0.83660130761563778, 0.82580645196139812, 0.81528662331402302, 0.80503144674003124,
    0.79503105580806732, 0.78527607396245003, 0.77575757540762424, 0.76646706648170948,
    0.75739645026624203, 0.74853801168501377, 0.73988439328968525, 0.73142857104539871,
    0.72316384129226208, 0.71508379839360714, 0.70718231983482838, 0.69945355132222176,
    0.69189189188182354, 0.68449197895824909, 0.67724867723882198, 0.67015706747770309,
    0.66321243532001972, 0.65641025640070438, 0.64974619261920452, 0.6432160809636116,
    0.63681592047214508, 0.63054187223315239, 0.62439024448394775, 0.61835748702287674,
    0.61244019120931625, 0.60663507133722305, 0.60093896649777889, 0.59534883685410023,
    0.58986175060272217, 0.58447488583624363, 0.57918552123010159, 0.57399103231728077,
    0.56888888962566853, 0.56387665122747421, 0.55895196460187435, 0.55411255359649658,
    0.54935622401535511, 0.54468085058033466, 0.54008438810706139, 0.53556485287845135,
    0.53112033195793629, 0.52674897201359272, 0.52244897931814194, 0.5182186234742403,
    0.51405622437596321, 0.50996015965938568, 0.50592885352671146, 0.50196078419685364);
  tl_atanh: array[0..63] of Double = (
    0.0038910702064755593, 0.011583529917253579, 0.019159432158344258, 0.0266222572666821,
    0.033975330925150045, 0.041221834321774946, 0.048364813673108971, 0.05540718355577081,
    0.062351739003096061, 0.069201160978450194, 0.075958020852849917, 0.082624786913314874,
    0.089203828481750638, 0.095697426288811963, 0.10210777126004227, 0.1084369690266159,
    0.11468705059063057, 0.12085996822529385, 0.12695760521603641, 0.1329817738483913,
    0.13893422513065432, 0.14481664629879729, 0.15063066513628579, 0.15637785526388293,
    0.16205973468362792, 0.16767777130253891, 0.17323338410238578, 0.17872794488290744,
    0.18416278058662977, 0.1895391762128262, 0.19485837557728855, 0.20012158250006382,
    0.20532996242715057, 0.21048464732934077, 0.21558673262018843, 0.22063727996588015,
    0.22563732201152162, 0.23058785731370252, 0.23548985714373422, 0.24034426540047171,
    0.24515199416811606, 0.24991393457449784, 0.25463095142604886, 0.25930388240233709,
    0.26393354527608248, 0.26852073295571777, 0.2730662180496442, 0.27757075296261952,
    0.28203506849484128, 0.28645987745028084, 0.29084587023204084, 0.29519372376674946,
    0.29950409405906614, 0.30377762555610432, 0.3080149386805166, 0.31221664463895504,
    0.31638333477824293, 0.32051558893193816, 0.3246139735744894, 0.32867903636145596,
    0.33271131678186228, 0.33671133731504505, 0.34067961262950625, 0.34461664073581982);
  ln2n_atanh: array[0..23] of Double = (
    0.34657359016697264, 0.69314718044694534, 1.039720770726918, 1.3862943610068905,
    1.7328679512868632, 2.0794415415668359, 2.4260151318468086, 2.7725887221267813,
    3.119162312406754, 3.4657359026867267, 3.8123094929666994, 4.1588830832466721,
    4.5054566735266448, 4.8520302638066175, 5.1986038540865902, 5.5451774443665629,
    5.8917510346465347, 6.2383246249265074, 6.5848982152064801, 6.9314718054864528,
    7.2780453957664255, 7.6246189860463982, 7.9711925763263709, 8.3177661666063436);
  b_atanh: array[0..2] of Double = (0.49999999981938592, -0.25000751197527959, 0.16667568182579121);
  c_atanh_s: array[0..3] of Double = (0.33333333333333076, 0.20000000005949592, 0.14285692689117752, 0.11136190102118979);
  c_atanh_acc: array[0..6] of Double = (
    0.5, -0.2500000000000015, 0.16666666666666946, -0.12499999980249463,
    0.099999999758626501, -0.083339906466727953, 0.071435144331263092);
var
  ta: Tb32u32;
  ux_a, ax_a: UInt32;
  sgn_a: Double;
  e_a: UInt32;
  md_a, mn_a: UInt32;
  nz_a: Integer;
  jn_a, jd_a: Integer;
  tn_a, td_a: Tb64u64;
  zn_a, zd_a, zn2_a, zd2_a: Double;
  rn_a, rd_a, r_a: Double;
  ub_a: Single;
  lb_a: Single;
  zz_a: Double;
  z2_a, z4_a: Double;
  zn4_a, zd4_a: Double;
  fn_a, fd_a, en_a: Double;
begin
  ta.f := x;
  ux_a := ta.u;
  ax_a := ux_a shl 1;
  if (ax_a < $7A300000) or (ax_a >= $7F000000) then begin
    if ax_a >= $7F000000 then begin  // NaN or |x|>=1
      if ax_a = $7F000000 then begin  // +-1
        Result := pcr_copysignf(Single(1.0/0.0), x); Exit;
      end;
      if ax_a > $FF000000 then begin Result := x + x; Exit; end;  // nan
      Result := Single(0.0/0.0); Exit;  // |x|>1
    end;
    if ax_a < $73713744 then begin  // |x| < 0.000352...
      if ax_a = 0 then begin Result := x; Exit; end;
      Result := pcr_fmaf(x, 2.9802322387695312e-08, x); Exit;
    end;
    zz_a := x; z2_a := zz_a*zz_a; z4_a := z2_a*z2_a;
    r_a := c_atanh_s[0] + z2_a*c_atanh_s[1] + z4_a*(c_atanh_s[2] + z2_a*c_atanh_s[3]);
    Result := Single(zz_a + (zz_a*z2_a)*r_a); Exit;
  end;
  if (ux_a shr 31) <> 0 then sgn_a := -1.0 else sgn_a := 1.0;
  e_a := ax_a shr 24;
  md_a := ((ux_a shl 8) or ($80000000)) shr (126 - e_a);
  mn_a := UInt32(-Int32(md_a));
  nz_a := 32 - pcr_bsr32(mn_a);  // = clz(mn)+1
  mn_a := mn_a shl nz_a;
  jn_a := Integer(mn_a shr 26);
  jd_a := Integer(md_a shr 26);
  tn_a.u := UInt64(Int64(Int64(mn_a) shl 20) or (Int64(1023) shl 52));
  td_a.u := UInt64(Int64(Int64(md_a) shl 20) or (Int64(1023) shl 52));
  zn_a := tn_a.f * tr_atanh[jn_a] - 1.0;
  zd_a := td_a.f * tr_atanh[jd_a] - 1.0;
  zn2_a := zn_a * zn_a;
  zd2_a := zd_a * zd_a;
  rn_a := ((tl_atanh[jn_a] - ln2n_atanh[nz_a-1]) + zn_a*b_atanh[0]) + zn2_a*(b_atanh[1] + zn_a*b_atanh[2]);
  rd_a := (tl_atanh[jd_a] + zd_a*b_atanh[0]) + zd2_a*(b_atanh[1] + zd_a*b_atanh[2]);
  r_a := sgn_a*(rd_a - rn_a);
  ub_a := Single(r_a);
  lb_a := Single(r_a + sgn_a*0.226e-9);
  if ub_a <> lb_a then begin
    zn4_a := zn2_a*zn2_a; zd4_a := zd2_a*zd2_a;
    fn_a := zn_a*(((c_atanh_acc[0] + zn_a*c_atanh_acc[1]) + zn2_a*(c_atanh_acc[2] + zn_a*c_atanh_acc[3])) +
                  zn4_a*((c_atanh_acc[4] + zn_a*c_atanh_acc[5]) + zn2_a*c_atanh_acc[6]));
    fn_a := fn_a + 9.3209433686215166e-16 * Double(nz_a);
    fn_a := fn_a + tl_atanh[jn_a];
    en_a := Double(nz_a) * 0.34657359027997359;
    fd_a := zd_a*(((c_atanh_acc[0] + zd_a*c_atanh_acc[1]) + zd2_a*(c_atanh_acc[2] + zd_a*c_atanh_acc[3])) +
                  zd4_a*((c_atanh_acc[4] + zd_a*c_atanh_acc[5]) + zd2_a*c_atanh_acc[6]));
    fd_a := fd_a + tl_atanh[jd_a];
    r_a := fd_a - fn_a + en_a;
    ub_a := Single(sgn_a * r_a);
  end;
  Result := ub_a;
end;

// ── 3.06 exp10m1f ────────────────────────────────────────────────────────────
function pcr_exp10m1f(x: Single): Single;
const
  c_e10: array[0..5] of Double = (
    0.043321698784995886, 0.00093838479282008368, 1.3550807712983854e-05,
    1.4676119301623784e-07, 1.2713094157155389e-09, 9.3824389539780747e-12);
  tb_e10: array[0..15] of Double = (
    1, 1.0442737824274138, 1.0905077326652577, 1.1387886347566916, 1.189207115002721, 1.241857812073484, 1.2968395546510096, 1.3542555469368927,
    1.4142135623730951, 1.4768261459394993, 1.5422108254079405, 1.6104903319492543, 1.681792830507429, 1.7562521603732995, 1.8340080864093424, 1.9152065613971474);
  cp4_e10: array[0..3] of Double = (2.3025850929940459, 2.6509490552391992, 2.0346785922934552, 1.1712551489193503);
  cp5_e10: array[0..4] of Double = (2.3025850929940459, 2.6509490552387951, 2.0346785922938739, 1.1712557955234444, 0.53938292940865262);
  cp6_e10: array[0..5] of Double = (
    2.3025850929940459, 2.6509490552391983, 2.0346785922348913, 1.1712551489516381, 0.53938692370821983, 0.20699584816918598);
  cp7_e10: array[0..6] of Double = (
    2.3025850929940459, 2.6509490552392512, 2.0346785922933694, 1.1712551474718793, 0.53938292993265813, 0.20700578860031116, 0.068089364982424655);
  cp8_e10: array[0..7] of Double = (
    2.3025850929940455, 2.6509490552391997, 2.0346785922965154, 1.1712551489080671, 0.53938291788367909, 0.20699585338612078, 0.068102837768701199, 0.019597694483460711);
  cp9_e10: array[0..8] of Double = (
    2.3025850929940455, 2.6509490552391819, 2.0346785922935298, 1.1712551489623777, 0.53938292914310315, 0.20699580881200672, 0.068089378992517491, 0.019609449708105794, 0.0050139289122738354);
var
  t10: Tb32u32;
  ux10, ax10: UInt32;
  z10: Double;
  a10, ia10, h10: Double;
  i10, j10: Int64;
  e10: Int64;
  s10: Double;
  su: Tb64u64;
  h2_10, c0_10, c2_10, c4_10, w10: Double;
  rr: Double;
  k10: UInt32;
begin
  t10.f := x;
  ux10 := t10.u;
  ax10 := ux10 and $7FFFFFFF;
  z10 := x;
  if ux10 > $C0F0D2F1 then begin  // x < -7.52575
    if ax10 > ($FF shl 23) then begin Result := x + x; Exit; end;  // nan
    if ux10 = $FF800000 then Result := -1.0
    else Result := -1.0 + Single(1.4901161193847656e-08); Exit;
  end else if ax10 > $421A209A then begin  // x > 38.5318
    if ax10 >= ($FF shl 23) then begin Result := x + x; Exit; end;
    Result := Single(3.4028234663852886e+38 + 3.4028234663852886e+38); Exit;
  end else if ax10 < $3D89C604 then begin  // |x| < 0.1549/log(10)
    rr := 0.0;
    if ax10 < $3D1622FB then begin
      if ax10 < $3C8B76A3 then begin
        if ax10 < $3BCCED04 then begin
          if ax10 < $3ACF33EB then begin
            if ax10 < $395A966B then begin
              if ax10 < $36FE4A4B then begin
                if ax10 < $32407F39 then begin
                  if ax10 < $245E5BD9 then begin
                    rr := 2.3025850929940459;
                  end else begin
                    if ux10 = $2C994B7B then begin Result := Single(1.003213657979618e-11) - Single(8.0779356694631609e-28); Exit; end;
                    rr := 2.3025850929940459 + z10 * 2.6509490552391992;
                  end;
                end else begin
                  if ux10 = $B6FA215B then begin Result := Single(-1.7164389646495692e-05) + Single(3.3881317890172014e-21); Exit; end;
                  rr := 2.3025850929940459 + z10*(2.6509490552896504 + z10*2.0346785922934552);
                end;
              end else begin
                rr := (cp4_e10[0] + z10*cp4_e10[1]) + z10*z10*(cp4_e10[2] + z10*cp4_e10[3]);
              end;
            end else begin
              rr := (cp5_e10[0] + z10*cp5_e10[1]) + z10*z10*(cp5_e10[2] + z10*(cp5_e10[3] + z10*cp5_e10[4]));
            end;
          end else begin
            rr := (cp6_e10[0] + z10*cp6_e10[1]) + z10*z10*((cp6_e10[2] + z10*cp6_e10[3]) + z10*z10*(cp6_e10[4] + z10*cp6_e10[5]));
          end;
        end else begin
          rr := (cp7_e10[0] + z10*cp7_e10[1]) + z10*z10*((cp7_e10[2] + z10*cp7_e10[3]) + z10*z10*(cp7_e10[4] + z10*(cp7_e10[5] + z10*cp7_e10[6])));
        end;
      end else begin
        rr := (cp8_e10[0] + z10*cp8_e10[1]) + z10*z10*((cp8_e10[2] + z10*cp8_e10[3]) + (z10*z10)*((cp8_e10[4] + z10*cp8_e10[5]) + z10*z10*(cp8_e10[6] + z10*cp8_e10[7])));
      end;
    end else begin
      rr := ((cp9_e10[0] + z10*cp9_e10[1]) + z10*z10*(cp9_e10[2] + z10*cp9_e10[3])) + (z10*z10*z10*z10)*((cp9_e10[4] + z10*cp9_e10[5]) + z10*z10*(cp9_e10[6] + z10*(cp9_e10[7] + z10*cp9_e10[8])));
    end;
    Result := Single(rr * z10); Exit;
  end else begin
    // check exact powers of 10
    if (ux10 shl 11) = 0 then begin
      k10 := (ux10 shr 21) - $1FC;
      if k10 <= $B then begin
        if k10 = 0 then begin Result := 10.0-1.0; Exit; end;
        if k10 = 4 then begin Result := 100.0-1.0; Exit; end;
        if k10 = 6 then begin Result := 1000.0-1.0; Exit; end;
        if k10 = 8 then begin Result := 10000.0-1.0; Exit; end;
        if k10 = 9 then begin Result := 100000.0-1.0; Exit; end;
        if k10 = 10 then begin Result := 1000000.0-1.0; Exit; end;
        if k10 = 11 then begin Result := 10000000.0-1.0; Exit; end;
      end;
    end;
    a10 := 53.150849461555481 * z10;
    ia10 := Floor(a10);
    h10 := (a10 - ia10) + 5.6642316608893863e-08 * z10;
    i10 := Trunc(ia10);
    j10 := i10 and $F;
    e10 := i10 - j10;
    e10 := e10 shr 4;
    s10 := tb_e10[j10];
    su.u := (UInt64(e10 + $3FF)) shl 52;
    s10 := s10 * su.f;
    h2_10 := h10 * h10;
    c0_10 := c_e10[0] + h10*c_e10[1];
    c2_10 := c_e10[2] + h10*c_e10[3];
    c4_10 := c_e10[4] + h10*c_e10[5];
    c0_10 := c0_10 + h2_10*(c2_10 + h2_10*c4_10);
    w10 := s10 * h10;
    Result := Single((s10 - 1.0) + w10*c0_10);
  end;
end;


// ── 3.07 log10p1f ────────────────────────────────────────────────────────────
function pcr_log10p1f(x: Single): Single;
const
  tr_l10: array[0..63] of Double = (
    0.9922480583190918, 0.97709923982620239, 0.96240603923797607, 0.94814813137054443,
    0.9343065619468689, 0.92086333036422729, 0.90780138969421387, 0.89510488510131836,
    0.88275861740112305, 0.87074828147888184, 0.85906040668487549, 0.84768211841583252,
    0.83660131692886353, 0.82580643892288208, 0.81528663635253906, 0.805031418800354,
    0.79503107070922852, 0.78527605533599854, 0.7757575511932373, 0.76646709442138672,
    0.75739645957946777, 0.74853801727294922, 0.73988437652587891, 0.73142856359481812,
    0.72316384315490723, 0.7150837779045105, 0.70718234777450562, 0.69945353269577026,
    0.69189190864562988, 0.68449199199676514, 0.67724865674972534, 0.67015707492828369,
    0.66321241855621338, 0.65641027688980103, 0.64974617958068848, 0.64321607351303101,
    0.63681590557098389, 0.63054186105728149, 0.62439024448394775, 0.61835747957229614,
    0.61244016885757446, 0.60663509368896484, 0.60093897581100464, 0.59534883499145508,
    0.58986175060272217, 0.58447486162185669, 0.57918554544448853, 0.57399106025695801,
    0.56888890266418457, 0.56387662887573242, 0.5589519739151001, 0.55411255359649658,
    0.54935622215270996, 0.54468083381652832, 0.5400843620300293, 0.5355648398399353,
    0.53112035989761353, 0.52674895524978638, 0.52244895696640015, 0.51821863651275635,
    0.51405620574951172, 0.50996017456054688, 0.50592887401580811, 0.50196081399917603);
  tl_l10: array[0..63] of Double = (
    0.0033797422691004653, 0.010061324592103849, 0.016641660398421521, 0.023123806531881508,
    0.029510600946365428, 0.035804820696601423, 0.04200915675927927, 0.048126072670659374,
    0.054158034204826393, 0.060107373998456282, 0.065976296741911244, 0.071766978049616167,
    0.077481456113724118, 0.083121735195996221, 0.088689676289720187, 0.094187169637755655,
    0.099615898294463157, 0.10497764486764279, 0.11027398831780565, 0.11550648532083185,
    0.12067672930309649, 0.1257861375083863, 0.13083614318801212, 0.13582808368965787,
    0.1407632959048486, 0.14564307407262186, 0.15046858843573072, 0.15524113201422163,
    0.15996174823881779, 0.16463162839464418, 0.1692518476704411, 0.17382339315055523,
    0.17834735028039525, 0.18282462816481015, 0.1872562654118731, 0.19164311141307008,
    0.19598609788417384, 0.20028607574785418, 0.20454389100326453, 0.20876038167368596,
    0.21293633243952884, 0.21707246947094125, 0.2211696275214583, 0.22522849188545682,
    0.22924976460497623, 0.23323416319093454, 0.23718228522931573, 0.24109487155885392,
    0.24497253794716645, 0.24881590534170478, 0.2526255058159062, 0.25640201064859103,
    0.26014595218693398, 0.26385790637563561, 0.26753839739443214, 0.27118794242299332,
    0.27480705007439105, 0.27839631709550916, 0.28195613352428522, 0.28548697269100148,
    0.28898939362644532, 0.29246373889003263, 0.29591053413566176, 0.29933018510220022);
  b_l10: array[0..3] of Double = (0.43429448180522795, -0.21714724073915409, 0.14477135206183736, -0.1085812344983514);
  c_l10: array[0..8] of Double = (
    0.43429448190325182, -0.21714724095162236, 0.14476482730107551, -0.10857362051434018, 0.086858896436697877,
    -0.072382301461477352, 0.062041939634640365, -0.054407774191098719, 0.048375922777279659);
  st_l10: array[0..15] of UInt32 = (
    $43928000, $41100000, $42C60000, $00000000, $4479C000, $00000000, $461C3C00, $47C34F80,
    $00000000, $497423F0, $00000000, $4B18967F, $00000000, $00000000, $00000000, $00000000);
var
  tl: Tb32u32;
  ux_l, ax_l: UInt32;
  zl: Double;
  tz_l: Tb64u64;
  m64_l: UInt64;
  e_l, j_l: Integer;
  vv_l, v2_l, v4_l: Double;
  f_l, lj_v: Double;
  ub_l: Single;
  lb_l: Single;
  ie_l, je_l, jeX: UInt32;
begin
  tl.f := x;
  ux_l := tl.u;
  zl := x;
  if ux_l >= $BF800000 then begin  // x <= -1
    if tl.u = $BF800000 then begin  // x = -1
      pcr_feraiseexcept_divbyzero;
      tl.u := $FF800000;  // -inf
      Result := tl.f; Exit;
    end;
    if (ux_l shl 1) > $FF000000 then begin Result := x + x; Exit; end;
    pcr_feraiseexcept_invalid;
    tl.u := $FFC00000; Result := tl.f; Exit;
  end;
  if Int32(ux_l) >= Int32($7F800000) then begin  // +inf or +nan (signed compare, excludes negatives)
    if ux_l > $7F800000 then begin Result := x + x; Exit; end;
    Result := x; Exit;
  end;
  if ux_l = st_l10[(ux_l shr 24) and $F] then begin
    ie_l := ux_l;
    ie_l := ie_l shr 23;
    je_l := ie_l - 126;
    jeX := (je_l * $9A209A8) shr 29;
    Result := Single(jeX); Exit;
  end;
  tz_l.f := zl + 1.0;
  m64_l := tz_l.u and (not UInt64(0) shr 12);
  e_l := Integer((tz_l.u shr 52) - 1023);
  j_l := Integer(m64_l shr 46);
  tz_l.u := m64_l or (UInt64($3FF) shl 52);
  if m64_l = 0 then begin
    if (ux_l = 0) or (ux_l = $80000000) then begin Result := x; Exit; end;
  end;
  vv_l := tz_l.f * tr_l10[j_l] - 1.0;
  v2_l := vv_l * vv_l;
  f_l := (Double(e_l)*0.3010299956639812 + tl_l10[j_l]) + vv_l*((b_l10[0] + vv_l*b_l10[1]) + v2_l*(b_l10[2] + vv_l*b_l10[3]));
  ub_l := Single(f_l);
  lb_l := Single(f_l + 3.0431213104975541e-13);
  if ub_l <> lb_l then begin
    if ux_l = $399A7C00 then begin Result := Single(0.00012794842768926173) + Single(3.637978807091713e-12); Exit; end;
    if ux_l = $B051E173 then begin Result := Single(-3.3160182932867599e-10) + Single(6.9388939039072284e-18); Exit; end;
    if ux_l = $BB0EA6D9 then begin Result := Single(-0.0009463561000302434) + Single(1.4551915228366852e-11); Exit; end;
    if ux_l = $BD86FFB9 then begin Result := Single(-0.029614737257361412) + Single(-4.6566128730773926e-10); Exit; end;
    lj_v := tl_l10[j_l] + 1.5315526624704034e-13;
    ax_l := ux_l and $7FFFFFFF;
    if ax_l < $3D100000 then begin  // |x| < 0x1.2p-5
      if ax_l < $33000000 then begin  // |x| < 0x1p-25
        ub_l := Single(zl*0.43429448176175356 + zl*(1.414982685385801e-10 + zl*c_l10[1]));
        Result := ub_l; Exit;
      end;
      e_l := 0;
      vv_l := x;
      v2_l := vv_l * vv_l;
      lj_v := 0.0;
    end;
    v4_l := v2_l * v2_l;
    f_l := vv_l*(((c_l10[0] + vv_l*c_l10[1]) + v2_l*(c_l10[2] + vv_l*c_l10[3])) + v4_l*((c_l10[4] + vv_l*c_l10[5]) + v2_l*((c_l10[6] + vv_l*c_l10[7]) + v2_l*c_l10[8])));
    f_l := f_l + Double(e_l)*-8.5323443170571066e-14;
    f_l := f_l + lj_v;
    f_l := f_l + Double(e_l)*0.30102999566406652;
    ub_l := Single(f_l);
  end;
  Result := ub_l;
end;

// Shared ix[] and lix[] tables for asinhf and acoshf (129 entries each)
const
  ix_asinh_acosh: array[0..128] of Double = (
    1, 0.99224806201527826, 0.98461538461560849, 0.97709923663933296,
    0.96969696969608776, 0.96240601503814105, 0.9552238805918023, 0.94814814814890269,
    0.9411764705873793, 0.93430656933924183, 0.92753623188764323, 0.92086330935126171,
    0.91428571428696159, 0.90780141843424644, 0.90140845070709474, 0.89510489509848412,
    0.88888888889050577, 0.88275862068985589, 0.87671232876891736, 0.87074829931952991,
    0.86486486486683134, 0.8590604026830988, 0.85333333333255723, 0.84768211920163594,
    0.84210526316019241, 0.83660130719363224, 0.83116883116599638, 0.82580645161215216,
    0.8205128205154324, 0.81528662420168985, 0.8101265822770074, 0.80503144653630443,
    0.80000000000291038, 0.79503105589537881, 0.79012345678347629, 0.78527607361320406,
    0.78048780487733893, 0.77575757575687021, 0.77108433734974824, 0.76646706587052904,
    0.76190476190822665, 0.75739644969871733, 0.75294117646990344, 0.74853801169956569,
    0.74418604651873466, 0.73988439305685461, 0.7356321839033626, 0.73142857142374851,
    0.72727272727934178, 0.72316384180157911, 0.71910112359910272, 0.71508379888837226,
    0.711111111115315, 0.70718232044600882, 0.7032967033010209, 0.69945355191885028,
    0.6956521739193704, 0.69189189189637545, 0.6881720430101268, 0.68449197860900313,
    0.68085106383659877, 0.6772486772533739, 0.67368421053106431, 0.6701570680597797,
    0.66666666667151731, 0.66321243523270823, 0.65979381442593876, 0.6564102564152563,
    0.65306122449692339, 0.64974619289569091, 0.64646464645920787, 0.64321608039608691,
    0.63999999999941792, 0.6368159203993855, 0.63366336633043829, 0.63054187192756217,
    0.6274509803915862, 0.62439024390187114, 0.62135922329616733, 0.61835748792509548,
    0.61538461539021228, 0.61244019138393924, 0.6095238095294917, 0.60663507108984049,
    0.60377358490950428, 0.60093896713806316, 0.59813084112829529, 0.59534883720334619,
    0.59259259259852115, 0.58986175115569495, 0.58715596330875997, 0.58447488585079554,
    0.58181818181765266, 0.57918552035698667, 0.5765765765827382, 0.57399103138595819,
    0.57142857143480796, 0.56888888888352085, 0.56637168141605798, 0.5638766519841738,
    0.56140350877831224, 0.55895196506753564, 0.55652173912676517, 0.55411255410581362,
    0.55172413792752195, 0.54935622317134403, 0.5470085470151389, 0.54468085106054787,
    0.54237288136209827, 0.54008438817982096, 0.53781512605200987, 0.53556485356239136,
    0.53333333334012423, 0.53112033194338437, 0.52892561983026098, 0.52674897119868547,
    0.52459016392822377, 0.52244897959462833, 0.52032520325155929, 0.51821862348879222,
    0.5161290322575951, 0.51405622489983216, 0.51200000000244472, 0.50996015936834738,
    0.5079365079436684, 0.50592885375954211, 0.50393700787390117, 0.50196078431326896,
    0.5);
  lix_asinh_acosh: array[0..128] of Double = (
    0, 0.0077821404422823226, 0.015504186535737881, 0.023167059283467056,
    0.030771658667663182, 0.038318864301568167, 0.045809536036751169, 0.053244514518016477,
    0.060624621817344335, 0.067950661912600477, 0.075223421233722179, 0.082443669212438828,
    0.089612158688322896, 0.096729626464576515, 0.10379679367846033, 0.11081436634745238,
    0.11778303565456447, 0.12470347850072987, 0.1315763577866729, 0.1384023228593465,
    0.14518200984222415, 0.15191604202754727, 0.15860503017754807, 0.16524957289962727,
    0.17185025692393074, 0.17840765746792978, 0.18492233849742259, 0.19139485300053896,
    0.19782574332673664, 0.20421554143130569, 0.21056476910916863, 0.2168739383062987,
    0.22314355131057179, 0.22937410107143966, 0.23556607132117974, 0.24171993689533061,
    0.24783616390549076, 0.25391520998187295, 0.2599575244364713, 0.26596354849418208,
    0.27193371547909428, 0.27786845101061858, 0.28376817313155411, 0.28963329257815412,
    0.2954642128842862, 0.30126133058725674, 0.30702503530127834, 0.31275571001049074,
    0.31845373110943964, 0.32411946866296587, 0.32975328636746576, 0.33535554191317973,
    0.34092658696468148, 0.34646676734052423, 0.35197642315103911, 0.35745588891282254,
    0.36290549368027353, 0.36832556115222748, 0.37371640979449355, 0.37907835293587894,
    0.38441169890032761, 0.38971675113309029, 0.39499380823382041, 0.40024316413156019,
    0.40546510810088843, 0.4106599249859505, 0.41582789515439755, 0.42096929463651261,
    0.42608439529998615, 0.43117346481484703, 0.4362367667833309, 0.44127456081408384,
    0.44628710262932902, 0.45127464413729851, 0.45623743349136464, 0.46117571511205202,
    0.46608972992550873, 0.4709797152197005, 0.47584590487769463, 0.48068852934188655,
    0.48550781577260588, 0.49030398805110553, 0.4950772667885292, 0.49982786955679037,
    0.50455601074602885, 0.50926190178662467, 0.5139457510908656, 0.5186077642180501,
    0.52324814375454343, 0.5278670896147033, 0.53246479885924003, 0.53704146588653812,
    0.5415972824336539, 0.54613243760677588, 0.55064711794197574, 0.55514150754777758,
    0.55961578792450872, 0.56407013829423902, 0.56850473535244139, 0.57291975355860225,
    0.57731536502345493, 0.58169173963098453, 0.58604904501017208, 0.59038744661434084,
    0.59470710775305924, 0.5990081896544962, 0.60329085142603345, 0.60755525023056722,
    0.61180154109462426, 0.6160298772263143, 0.62024040974890171, 0.62443328800086584,
    0.62860865940964117, 0.63276666958388439, 0.636907462245482, 0.6410311794109268,
    0.64513796138540813, 0.64922794661976657, 0.65330127201365518, 0.65735807269483126,
    0.66139848224627451, 0.66542263254463574, 0.66943065393785439, 0.67342267520079802,
    0.67739882357770898, 0.68135922479880817, 0.68530400309914674, 0.6892332812397185,
    0.69314718055994529);

// ── 3.08 asinhf ────────────────────────────────────────────────────────────
function pcr_asinhf(x: Single): Single;
const
  c_asinh: array[0..7] of Double = (
    0.1666666666666666, -0.074999999999870018, 0.044642857099780067, -0.03038193899998537, 0.022371820451468214, -0.017341279402218638, 0.013747204759994313, -0.0093574477267578029);
  cm_asinh: array[0..2] of Double = (1.0000000000932958, -0.50000378550500935, 0.33332252602066714);
  cp_asinh: array[0..5] of Double = (
    1, -0.5, 0.33333333331462334, -0.24999999997581948, 0.20000326978745125, -0.16666993701509006);
var
  tsi: Tb32u32;
  xs_a: Double;
  x2_a, x4_a, x8_a: Double;
  f_a: Double;
  xd_a: Double;
  tp_a: Tb64u64;
  m_a: UInt64;
  j_a, e_a: Integer;
  ww_a: Tb64u64;
  z_a, z2_a: Double;
  r_a: Tb64u64;
  c0_a, c2_a, c4_a: Double;
  Lh_a, Ll_a, hh_a: Double;
  res_f: Single;
begin
  tsi.f := x;
  tsi.u := tsi.u and $7FFFFFFF;  // clear sign
  xs_a := x;
  if tsi.u <= $3E815667 then begin  // |x| <= 0.252...
    if tsi.u <= $39DDB3D7 then begin  // |x| <= 0.000422...
      if tsi.u = 0 then begin Result := x; Exit; end;
      Result := pcr_fmaf(x, -2.9802322387695312e-08, x); Exit;
    end;
    x2_a := xs_a*xs_a; x4_a := x2_a*x2_a; x8_a := x4_a*x4_a;
    f_a := x2_a*(((c_asinh[0] + x2_a*c_asinh[1]) + x4_a*(c_asinh[2] + x2_a*c_asinh[3])) +
                  x8_a*((c_asinh[4] + x2_a*c_asinh[5]) + x4_a*(c_asinh[6] + x2_a*c_asinh[7])));
    Result := Single(xs_a - xs_a*f_a); Exit;
  end else begin
    if tsi.u >= $7F800000 then begin Result := x + x; Exit; end;
    xd_a := pcr_fabs(xs_a);
    x2_a := xd_a * xd_a;
    tp_a.f := xd_a + Sqrt(x2_a + 1.0);
    m_a := tp_a.u and (not UInt64(0) shr 12);
    j_a := Integer((m_a + (UInt64(1) shl (52-8))) shr (52-7));
    e_a := Integer((tp_a.u shr 52) - $3FF);
    ww_a.u := m_a or (UInt64($3FF) shl 52);
    z_a := ww_a.f * ix_asinh_acosh[j_a] - 1.0;
    z2_a := z_a * z_a;
    r_a.f := ((lix_asinh_acosh[128] * Double(e_a) + lix_asinh_acosh[j_a]) + z_a*cm_asinh[0]) + z2_a*(cm_asinh[1] + z_a*cm_asinh[2]);
    if ((r_a.u + 259000) and $0FFFFFFF) < 260000 then begin  // accurate path
      z2_a := z_a * z_a;
      c0_a := cp_asinh[0] + z_a*cp_asinh[1];
      c2_a := cp_asinh[2] + z_a*cp_asinh[3];
      c4_a := cp_asinh[4] + z_a*cp_asinh[5];
      c0_a := c0_a + z2_a*(c2_a + z2_a*c4_a);
      Lh_a := 0.693145751953125 * Double(e_a);
      Ll_a := 1.4286068203094173e-06 * Double(e_a);
      r_a.f := (z_a * c0_a + (Ll_a + lix_asinh_acosh[j_a])) + Lh_a;
      if (r_a.u and $0FFFFFFF) = 0 then begin
        hh_a := (z_a * c0_a + (Ll_a + lix_asinh_acosh[j_a])) + (Lh_a - r_a.f);
        r_a.f := r_a.f + 64.0*hh_a;
      end;
    end;
    Result := Single(pcr_copysign(r_a.f, xs_a));
  end;
end;

// ── 3.09 acoshf ────────────────────────────────────────────────────────────
function pcr_acoshf(x: Single): Single;
const
  c_near_ac: array[0..7] of Double = (
    -0.08333333333328993, 0.018749999994215116, -0.0055803568456194285, 0.0018988638601201038,
    -0.0006990184157670764, 0.00027017646829586764, -0.00010420707647055027, 3.116411306139845e-05);
  cm_acosh: array[0..2] of Double = (
    1.0000000000932958, -0.5000037855050093, 0.33332252602066714);
  cp_acosh: array[0..5] of Double = (
    1.0, -0.5, 0.33333333331462334, -0.24999999997581948, 0.20000326978745125, -0.16666993701509006);
var
  tv_ac: Tb32u32;
  tp_ac, ww_ac, r_ac: Tb64u64;
  ux_ac: UInt32;
  zf_ac: Single;
  zz_ac, a_ac, x2_ac: Double;
  z2_ac, z4_ac, f_ac: Double;
  m64_ac: UInt64;
  j_ac, e_ac: Integer;
  c0_ac, c2_ac, c4_ac: Double;
  Lh_ac, Ll_ac, hh_ac: Double;
begin
  tv_ac.f := x;
  ux_ac := tv_ac.u;
  if ux_ac <= $3F800000 then begin
    // as_special: x=1 → 0, x<1 → FE_INVALID, nan → propagate
    if ux_ac = $3F800000 then begin Result := 0.0; Exit; end;
    if (ux_ac shl 1) > $FF000000 then begin Result := x + x; Exit; end; // nan
    pcr_feraiseexcept_invalid;
    Result := pcr_nanf(nil);
    Exit;
  end
  else if ux_ac <= $3F99DB23 then begin
    // near-1 path: acosh(x) ≈ sqrt(2*(x-1)) * (1 + poly(x-1))
    zf_ac := x - 1.0;
    zz_ac := Double(zf_ac);
    a_ac  := Sqrt(2.0 * zz_ac);
    z2_ac := zz_ac * zz_ac;
    z4_ac := z2_ac * z2_ac;
    f_ac  := ((c_near_ac[0] + zz_ac*c_near_ac[1]) + z2_ac*(c_near_ac[2] + zz_ac*c_near_ac[3]))
           + z4_ac*((c_near_ac[4] + zz_ac*c_near_ac[5]) + z2_ac*(c_near_ac[6] + zz_ac*c_near_ac[7]));
    Result := Single(a_ac + (a_ac * zz_ac) * f_ac);
    Exit;
  end
  else if ux_ac < $7F800000 then begin
    // main path: acosh(x) = log(x + sqrt(x^2 - 1))
    zz_ac := Double(x);
    x2_ac := zz_ac * zz_ac;
    tp_ac.f := zz_ac + Sqrt(x2_ac - 1.0);
    m64_ac := tp_ac.u and (not UInt64(0) shr 12);
    j_ac   := Integer((m64_ac + (UInt64(1) shl (52-8))) shr (52-7));
    e_ac   := Integer((tp_ac.u shr 52) - $3FF);
    ww_ac.u := m64_ac or (UInt64($3FF) shl 52);
    zz_ac  := ww_ac.f * ix_asinh_acosh[j_ac] - 1.0;
    z2_ac  := zz_ac * zz_ac;
    r_ac.f := ((lix_asinh_acosh[128]*Double(e_ac) + lix_asinh_acosh[j_ac]) + zz_ac*cm_acosh[0])
            + z2_ac*(cm_acosh[1] + zz_ac*cm_acosh[2]);
    if ((r_ac.u + 259000) and $0FFFFFFF) < 260000 then begin  // accurate path
      z2_ac := zz_ac * zz_ac;
      c0_ac := cp_acosh[0] + zz_ac*cp_acosh[1];
      c2_ac := cp_acosh[2] + zz_ac*cp_acosh[3];
      c4_ac := cp_acosh[4] + zz_ac*cp_acosh[5];
      c0_ac := c0_ac + z2_ac*(c2_ac + z2_ac*c4_ac);
      Lh_ac := 0.693145751953125 * Double(e_ac);
      Ll_ac := 1.4286068203094173e-06 * Double(e_ac);
      r_ac.f := (zz_ac * c0_ac + (Ll_ac + lix_asinh_acosh[j_ac])) + Lh_ac;
      if (r_ac.u and $0FFFFFFF) = 0 then begin
        hh_ac := (zz_ac * c0_ac + (Ll_ac + lix_asinh_acosh[j_ac])) + (Lh_ac - r_ac.f);
        r_ac.f := r_ac.f + 64.0*hh_ac;
      end;
    end;
    Result := Single(r_ac.f);
    Exit;
  end
  else begin
    // as_special: +inf → +inf, nan → propagate, negative finite → FE_INVALID
    if ux_ac = $7F800000 then begin Result := x; Exit; end;  // +inf
    if (ux_ac shl 1) > $FF000000 then begin Result := x + x; Exit; end;  // nan
    pcr_feraiseexcept_invalid;
    Result := pcr_nanf(nil);
  end;
end;


// ── lgammaf helpers ───────────────────────────────────────────────────────────
function lgamma_as_r7(x: Double; const c: array of Double): Double;
begin
  Result := (((x-c[0])*(x-c[1]))*((x-c[2])*(x-c[3])))*(((x-c[4])*(x-c[5]))*((x-c[6])));
end;

function lgamma_as_r8(x: Double; const c: array of Double): Double;
begin
  Result := (((x-c[0])*(x-c[1]))*((x-c[2])*(x-c[3])))*(((x-c[4])*(x-c[5]))*((x-c[6])*(x-c[7])));
end;

function lgamma_as_sinpi(x0: Double): Double;
const
  c_sp: array[0..7] of Double = (
    4.0, -3.7392088021786822, 1.2780132969493234, -0.22899788751887007,
    0.025330969591019943, -0.0019036696909283705, 0.00010351520901665411, -4.131827870010798e-06);
var
  x_sp, x2_sp, x4_sp, x8_sp: Double;
begin
  x_sp := x0 - 0.5;
  x2_sp := x_sp * x_sp; x4_sp := x2_sp * x2_sp; x8_sp := x4_sp * x4_sp;
  Result := (0.25 - x2_sp) * ((c_sp[0] + x2_sp*c_sp[1]) + x4_sp*(c_sp[2] + x2_sp*c_sp[3])
    + x8_sp*((c_sp[4] + x2_sp*c_sp[5]) + x4_sp*(c_sp[6] + x2_sp*c_sp[7])));
end;

function lgamma_as_ln(x: Double): Double;
const
  c_aln: array[0..7] of Double = (
    0.9999999999999756, -0.4999999999895039, 0.33333333159684553, -0.2499998558809608,
    0.19999326364239384, -0.16648082067582767, 0.13983949497072432, -0.09792103500684109);
  il_aln: array[0..15] of Double = (
    9.372730068330667e-18, 0.06062462181643487, 0.11778303565638353, 0.17185025692665928,
    0.2231435513142097, 0.2719337154836418, 0.3184537311185346, 0.3629054936893685,
    0.40546510810816444, 0.4462871026284195, 0.48550781578170077, 0.5232481437645479,
    0.5596157879354228, 0.5947071077466928, 0.6286086594223742, 0.661398482245365);
  ix_aln: array[0..15] of Double = (
    1.0, 0.9411764705882353, 0.8888888888888888, 0.8421052631578947,
    0.8, 0.7619047619047619, 0.7272727272727273, 0.6956521739130435,
    0.6666666666666666, 0.64, 0.6153846153846154, 0.5925925925925926,
    0.5714285714285714, 0.5517241379310345, 0.5333333333333333, 0.5161290322580645);
var
  t_aln: Tb64u64;
  e_aln, i_aln: Integer;
  z_aln, z2_aln, z4_aln: Double;
begin
  t_aln.f := x;
  e_aln := Integer(Int64(t_aln.u shr 52) - $3FF);
  i_aln := Integer((t_aln.u shr 48) and $F);
  t_aln.u := (t_aln.u and $000FFFFFFFFFFFFF) or $3FF0000000000000;
  z_aln := ix_aln[i_aln] * t_aln.f - 1.0;
  z2_aln := z_aln * z_aln; z4_aln := z2_aln * z2_aln;
  Result := Double(e_aln) * 0.6931471805599453 + il_aln[i_aln]
    + z_aln * ((c_aln[0] + z_aln*c_aln[1]) + z2_aln*(c_aln[2] + z_aln*c_aln[3])
               + z4_aln*((c_aln[4] + z_aln*c_aln[5]) + z2_aln*(c_aln[6] + z_aln*c_aln[7])));
end;

// ── 4.01 tgammaf ─────────────────────────────────────────────────────────────
function pcr_tgammaf(x: Single): Single;
const
  tb_xu: array[0..9] of UInt32 = (
    $27DE86A9, $27E05475, $B63BEFB3, $3C7BB570, $41E886D1, $C067D177,
    $BD99DA31, $BF54C45A, $41EE77FE, $3F843A64);
  tb_f: array[0..9] of Single = (
    1.61908237795328e+14, 1.60606292279296e+14, -357083.5625, 64.52886962890625,
    3.801414727062912e+29, 0.24537095427513123,
    -13.968554496765137, -6.604792594909668, 4.6287780950070885e+30, 0.9819810390472412);
  tb_df: array[0..9] of Single = (
    4194304.0, 4194304.0, 0.0078125, 1.9073486328125e-06,
    9.44473296573929e+21, 3.725290298461914e-09,
    -2.384185791015625e-07, 1.1920928955078125e-07,
    -1.3234889800848443e-23, 1.4901161193847656e-08);
  c_tg: array[0..15] of Double = (
    1.7877108988969403, 1.5591939012079508, 1.0510493266811867, 0.47065801829337245,
    0.18881863831977497, 0.058831548411746724, 0.017825943652294146, 0.004228758148929772,
    0.001097917853717287, 0.00019456568933897892, 5.197131759674315e-05, 4.914144140610218e-06,
    2.437173471710669e-06, -1.4461519623063317e-07, 1.8260131876052383e-07, -4.919948895618967e-08);
var
  t_tg: Tb32u32;
  rt_tg: Tb64u64;
  ax_tg: UInt32;
  z_tg, d_tg, m_tg, ii_tg, step_tg: Double;
  d2_tg, d4_tg, d8_tg, f_tg: Double;
  w_tg, x0_tg, t0_tg: Double;
  fx_tg: Single;
  k_tg: Int32;
  jm_tg, j_tg, lp_tg, idx_tg: Integer;
begin
  t_tg.f := x;
  ax_tg := t_tg.u shl 1;
  if ax_tg >= $FF000000 then begin         // x = NaN or +/-Inf
    if ax_tg = $FF000000 then begin        // x = +/-Inf
      if (t_tg.u shr 31) <> 0 then begin  // x = -Inf: Invalid
        Result := x / x; Exit;
      end;
      Result := x; Exit;                   // x = +Inf
    end;
    Result := x + x; Exit;                 // x = NaN: propagate
  end;
  z_tg := x;
  if ax_tg < $6D000000 then begin          // |x| < 2^-18
    d_tg := (0.9890559953279725 - 0.9074790760808863 * z_tg) * z_tg - 0.5772156649015329;
    f_tg := 1.0 / z_tg + d_tg;
    rt_tg.f := f_tg;
    if ((rt_tg.u + 2) and $0FFFFFFF) < 4 then
      for idx_tg := 0 to 9 do
        if t_tg.u = tb_xu[idx_tg] then begin
          Result := tb_f[idx_tg] + tb_df[idx_tg]; Exit;
        end;
    Result := Single(f_tg); Exit;
  end;
  // Safe single-precision floor: |x|>=2^23 means x is already exact integer
  if (t_tg.u and $7FFFFFFF) >= $4B800000 then
    fx_tg := x
  else
    fx_tg := Single(Floor(Double(x)));
  if x >= 35.04010009765625 then begin     // overflow: 0x1.18522p+5
    Result := Single(1.7014118346046923e+38) * Single(1.7014118346046923e+38); Exit;
  end;
  if x <= -2147483648.0 then
    k_tg := Low(Int32)
  else
    k_tg := Int32(Trunc(fx_tg));
  if fx_tg = x then begin                  // x is an integer
    if x = 0.0 then begin
      Result := 1.0 / x; Exit;            // +/-Inf, raises DivByZero
    end;
    if x < 0.0 then begin                  // negative integer: undefined
      Result := Single(0.0/0.0); Exit;
    end;
    t0_tg := 1.0; x0_tg := 1.0;
    lp_tg := 1;
    while lp_tg < k_tg do begin
      t0_tg := t0_tg * x0_tg;
      x0_tg := x0_tg + 1.0;
      Inc(lp_tg);
    end;
    Result := Single(t0_tg); Exit;
  end;
  if x < -42.0 then begin                  // negative non-integer, |gamma| < 2^-151
    if (k_tg and 1) = 0 then
      Result := Single(5.877471754111438e-39) * Single(5.877471754111438e-39)
    else
      Result := Single(5.877471754111438e-39) * Single(-5.877471754111438e-39);
    Exit;
  end;
  // Main polynomial path: gamma(x) via polynomial around 2.875
  m_tg := z_tg - 2.875;
  ii_tg := Double(Round(m_tg));
  if ii_tg < 0.0 then step_tg := -1.0 else step_tg := 1.0;
  d_tg := m_tg - ii_tg;
  d2_tg := d_tg * d_tg; d4_tg := d2_tg * d2_tg; d8_tg := d4_tg * d4_tg;
  f_tg :=   (c_tg[0] + d_tg*c_tg[1]) + d2_tg*(c_tg[2] + d_tg*c_tg[3])
    + d4_tg*((c_tg[4] + d_tg*c_tg[5]) + d2_tg*(c_tg[6] + d_tg*c_tg[7]))
    + d8_tg*((c_tg[8] + d_tg*c_tg[9]) + d2_tg*(c_tg[10] + d_tg*c_tg[11])
      + d4_tg*((c_tg[12] + d_tg*c_tg[13]) + d2_tg*(c_tg[14] + d_tg*c_tg[15])));
  jm_tg := Trunc(Abs(ii_tg));
  w_tg := 1.0;
  if jm_tg <> 0 then begin
    z_tg := z_tg - 0.5 - step_tg * 0.5;
    w_tg := z_tg;
    j_tg := jm_tg - 1;
    while j_tg > 0 do begin
      z_tg := z_tg - step_tg;
      w_tg := w_tg * z_tg;
      Dec(j_tg);
    end;
  end;
  if ii_tg <= -0.5 then w_tg := 1.0 / w_tg;
  f_tg := f_tg * w_tg;
  rt_tg.f := f_tg;
  if ((rt_tg.u + 2) and $0FFFFFFF) < 8 then
    for idx_tg := 0 to 9 do
      if t_tg.u = tb_xu[idx_tg] then begin
        Result := tb_f[idx_tg] + tb_df[idx_tg]; Exit;
      end;
  Result := Single(f_tg);
end;

// ── 4.02 lgammaf ─────────────────────────────────────────────────────────────
function pcr_lgammaf(x: Single): Single;
const
  tb_lg_xu: array[0..26] of UInt32 = (
    $1B7679FF, $1E88452D, $2AD345DE, $39EEFE83, $3B7C53AA, $42468B59,
    $449ACF07, $46541516, $46B16323, $4F3F94C0, $50522F52, $65FCA09F,
    $716E5DD5, $77AC5674, $7943DEF8, $8212E5B3, $9B7679FF, $9E88452D,
    $A77A8E47, $AA6C2DFF, $ABA1BF3A, $B0D6F2CA, $B0E17820,
    $C02C060F, $C134EB14, $C33139A3, $C6F7E151);
  tb_lg_f: array[0..26] of Single = (
    49.94450759887695, 45.68510437011719, 28.611061096191406,
    7.693094253540039, 5.557419300079346, 143.14707946777344,
    7578.81298828125, 115584.2109375, 205035.484375,
    67147280384.0, 315532181504.0, 7.808042139034897e+24,
    8.054993914128835e+31, 5.378053227661709e+35, 5.030267977202206e+36,
    85.11940002441406, 49.94450759887695, 45.68510437011719,
    33.29256057739258, 29.192766189575195, 27.491884231567383,
    20.27604866027832, 20.228261947631836,
    -0.08386971801519394, -16.91703224182129, -742.2763061523438, -297142.9375);
  tb_lg_df: array[0..26] of Single = (
    -9.5367431640625e-07, -9.5367431640625e-07, 4.76837158203125e-07,
    -1.1920928955078125e-07, 1.1920928955078125e-07, 3.814697265625e-06,
    -0.0001220703125, 0.001953125, -0.00390625,
    -1024.0, -8192.0, -1.4411518807585587e+17,
    1.2089258196146292e+24, -9.903520314283042e+27, 7.922816251426434e+28,
    -1.9073486328125e-06, -9.5367431640625e-07, -9.5367431640625e-07,
    -9.5367431640625e-07, -4.76837158203125e-07, -4.76837158203125e-07,
    -4.76837158203125e-07, 4.76837158203125e-07,
    1.862645149230957e-09, -4.76837158203125e-07, -1.52587890625e-05, -0.0078125);
  // Rational approx for |x| < 0x1.52p-1 (0.66015625)
  rn_sm: array[0..7] of Double = (
    -21.02242974712385, -5.2778355372066486, 1.0000000000070894, -2.5841418677193797,
    -1.6593264187887429, -1.291958597504499, -1.1075438471412518, -1.0223486918257203);
  c0_sm: Double = 4.246153550603663;
  rd_sm: array[0..7] of Double = (
    -41.71231062445876, -8.19318858118602, -3.5633472115484848, -2.1293437187570254,
    -1.491491016168074, -1.2138533890414456, -1.0760801136819742, -1.0134457501866234);
  // Rational approx for 0x1.52p-1 <= |x| <= 0x1.afc1ap+1 (medium range)
  rn_md: array[0..6] of Double = (
    -44.80915068895768, -9.412835290371454, -3.443662014713007, -1.4777007076373339,
    -0.6141091163391442, -0.23157674609865453, -0.048874810545154704);
  c0_md: Double = 4.949273735344536;
  rd_md: array[0..7] of Double = (
    -82.27591560063709, -14.56281571764002, -5.201990312146981, -2.293271883973741,
    -1.0742628098109357, -0.44083977020065385, -0.15523770855500355, -0.027342864372157546);
  // Stirling series: 2-term (ax > 73.9)
  stir2: array[0..1] of Double = (0.08333333332119137, -0.002777571142708969);
  // Stirling series: 4-term (ax > 10.67)
  stir4: array[0..3] of Double = (
    0.08333333333290953, -0.0027777773522524095, 0.0007935128869641844, -0.0005771890268644819);
  // Stirling series: 8-term (ax <= 10.67)
  stir8: array[0..7] of Double = (
    0.08333333333309598, -0.002777777717605562, 0.0007936447627791816, -0.0005949176485969618,
    0.0008316045241367893, -0.0017144163881980096, 0.003738059988828462, -0.00511695529630956);
  // Near gamma-zero expansion: around x ≈ -2.7477 (8 coefficients)
  c_nz1: array[0..7] of Double = (
    -1.9143501856115988, 9.57518947570967, -20.095134916873302, 62.627282710722085,
    -194.766146075866, 646.9063894581668, -2194.0017863504936, 7596.5102418865035);
  // Near gamma-zero expansion: around x ≈ -2.457 (7 coefficients)
  c_nz2: array[0..6] of Double = (
    1.5156034480216574, 4.858320951634164, 1.411291143065722, 8.721782533183026,
    5.8004160810161975, 24.828757532203195, 23.995763760665245);
  // Near gamma-zero expansion: around x ≈ -3.143 (7 coefficients)
  c_nz3: array[0..6] of Double = (
    7.781884658131351, 25.831338372388036, 112.26898629717154, 588.8907416229549,
    3277.1937988413574, 19024.282530470795, 113553.17206745526);
var
  t_lg: Tb32u32;
  rt_lg: Tb64u64;
  ax_lg: Single;
  fx_lg: Single;
  z_lg, s_lg, f_lg: Double;
  lz_lg, iz_lg, iz2_lg, iz4_lg, iz8_lg: Double;
  p_lg, lp_lg: Double;
  h_lg, h2_lg, h4_lg: Double;
  r_lg: Single;
  tl_lg: UInt64;
  a_lg, b_lg, mi_lg: Integer;
begin
  ax_lg := Abs(x);
  t_lg.f := ax_lg;
  // Safe floorf: float32 values with |x| >= 2^23 are already exact integers
  if ax_lg >= 8388608.0 then
    fx_lg := x
  else
    fx_lg := Single(Floor(Double(x)));
  // NaN or Inf (t_lg.u = bits of |x|)
  if t_lg.u >= $7F800000 then begin
    if t_lg.u = $7F800000 then begin Result := x * x; Exit; end;  // +/-Inf → +Inf
    Result := x + x; Exit;  // NaN
  end;
  // Integer input
  if fx_lg = x then begin
    if x <= 0.0 then begin
      pcr_feraiseexcept_divbyzero;
      Result := Single(1.0/0.0); Exit;   // lgamma(<=0 integer) = +Inf
    end;
    if (x = 1.0) or (x = 2.0) then begin
      Result := 0.0; Exit;               // lgamma(1) = lgamma(2) = 0
    end;
    // positive integer > 2: fall through to main computation
  end;
  z_lg := ax_lg;
  s_lg := x;
  if ax_lg < 0.66015625 then begin        // |x| < 0x1.52p-1
    f_lg := (c0_sm * s_lg) * lgamma_as_r8(s_lg, rn_sm)
          / lgamma_as_r8(s_lg, rd_sm) - lgamma_as_ln(z_lg);
  end else begin
    if ax_lg > 3.373096466064453 then begin   // ax > 0x1.afc1ap+1
      if x >= 4.085003425410169e+36 then begin  // overflow threshold
        Result := pcr_fmaf(x, 83.30038452148438, 1.0812689350146765e+31); Exit;
      end;
      lz_lg := lgamma_as_ln(z_lg);
      f_lg := (z_lg - 0.5) * (lz_lg - 1.0) + 0.4189385332046727;
      if ax_lg < 1048576.0 then begin     // ax < 2^20: add Stirling correction
        iz_lg := 1.0 / z_lg; iz2_lg := iz_lg * iz_lg;
        if ax_lg > 1198.0 then begin
          f_lg := f_lg + iz_lg * (1.0/12.0);
        end else if ax_lg > 73.90081787109375 then begin  // 0x1.279a7p+6
          f_lg := f_lg + iz_lg * (stir2[0] + iz2_lg * stir2[1]);
        end else if ax_lg > 10.666666984558105 then begin  // 0x1.555556p+3
          iz4_lg := iz2_lg * iz2_lg;
          f_lg := f_lg + iz_lg * ((stir4[0] + iz2_lg*stir4[1])
                                 + iz4_lg*(stir4[2] + iz2_lg*stir4[3]));
        end else begin
          iz4_lg := iz2_lg * iz2_lg; iz8_lg := iz4_lg * iz4_lg;
          p_lg := ((stir8[0] + iz2_lg*stir8[1]) + iz4_lg*(stir8[2] + iz2_lg*stir8[3]))
                + iz8_lg*((stir8[4] + iz2_lg*stir8[5]) + iz4_lg*(stir8[6] + iz2_lg*stir8[7]));
          f_lg := f_lg + iz_lg * p_lg;
        end;
      end;
      if x < 0.0 then begin              // reflection for negative x
        f_lg := 1.1447298858494002 - f_lg - lz_lg;
        lp_lg := lgamma_as_ln(lgamma_as_sinpi(Double(x) - Double(fx_lg)));
        f_lg := f_lg - lp_lg;
      end;
    end else begin                        // medium x: 0x1.52p-1 <= ax <= 0x1.afc1ap+1
      f_lg := (z_lg - 1.0) * (z_lg - 2.0) * c0_md
            * lgamma_as_r7(z_lg, rn_md) / lgamma_as_r8(z_lg, rd_md);
      if x < 0.0 then begin
        // Near gamma-zeros: use local Taylor expansion for accuracy
        if (t_lg.u < $40301B93) and (t_lg.u > $402F95C2) then begin
          // near x ≈ -2.7477
          h_lg := (s_lg + 2.7476826467274127) - 9.055340329338315e-17;
          h2_lg := h_lg * h_lg; h4_lg := h2_lg * h2_lg;
          f_lg := h_lg * ((c_nz1[0] + h_lg*c_nz1[1]) + h2_lg*(c_nz1[2] + h_lg*c_nz1[3])
            + h4_lg*((c_nz1[4] + h_lg*c_nz1[5]) + h2_lg*(c_nz1[6] + h_lg*c_nz1[7])));
        end else if (t_lg.u > $401CECCB) and (t_lg.u < $401D95CA) then begin
          // near x ≈ -2.457
          h_lg := (s_lg + 2.4570247382208006) + 3.7075610815513266e-17;
          h2_lg := h_lg * h_lg; h4_lg := h2_lg * h2_lg;
          f_lg := h_lg * ((c_nz2[0] + h_lg*c_nz2[1]) + h2_lg*(c_nz2[2] + h_lg*c_nz2[3])
            + h4_lg*((c_nz2[4] + h_lg*c_nz2[5]) + h2_lg*(c_nz2[6])));
        end else if (t_lg.u > $40492009) and (t_lg.u < $404940EF) then begin
          // near x ≈ -3.143
          h_lg := (s_lg + 3.14358088834998) + 2.1818179852331714e-16;
          h2_lg := h_lg * h_lg; h4_lg := h2_lg * h2_lg;
          f_lg := h_lg * ((c_nz3[0] + h_lg*c_nz3[1]) + h2_lg*(c_nz3[2] + h_lg*c_nz3[3])
            + h4_lg*((c_nz3[4] + h_lg*c_nz3[5]) + h2_lg*(c_nz3[6])));
        end else begin
          f_lg := 1.1447298858494002 - f_lg;
          lp_lg := lgamma_as_ln(lgamma_as_sinpi(Double(x) - Double(fx_lg)) * z_lg);
          f_lg := f_lg - lp_lg;
        end;
      end;
    end;
  end;
  // Table lookup for exceptional cases
  rt_lg.f := f_lg;
  tl_lg := (rt_lg.u + 5) and $0FFFFFFF;
  r_lg := Single(f_lg);
  if tl_lg <= 31 then begin
    t_lg.f := x;
    a_lg := 0; b_lg := 27;
    while a_lg + 1 < b_lg do begin
      mi_lg := (a_lg + b_lg) div 2;
      if t_lg.u < tb_lg_xu[mi_lg] then b_lg := mi_lg
      else a_lg := mi_lg;
    end;
    if t_lg.u = tb_lg_xu[a_lg] then begin
      Result := tb_lg_f[a_lg] + tb_lg_df[a_lg]; Exit;
    end;
  end;
  Result := r_lg;
end;

{ ── muldd / polydd: double-double helpers shared by pcr_atan2f and pcr_atan2pif ── }

function muldd(xh, xl, ch, cl: Double; out l: Double): Double;
var
  ahlh, alhh, ahhh, ahhl: Double;
begin
  ahlh := ch * xl;
  alhh := cl * xh;
  ahhh := ch * xh;
  ahhl := pcr_fma(ch, xh, -ahhh);
  ahhl := ahhl + alhh + ahlh;
  ch := ahhh + ahhl;
  l := (ahhh - ch) + ahhl;
  Result := ch;
end;

{ c is flat: c[k*2]=c[k][0], c[k*2+1]=c[k][1] }
function polydd(xh, xl: Double; n: Integer; const c: array of Double; out l: Double): Double;
var
  i: Integer;
  ch, cl, th, tl: Double;
begin
  i := n - 1;
  ch := c[i*2];
  cl := c[i*2+1];
  Dec(i);
  while i >= 0 do begin
    ch := muldd(xh, xl, ch, cl, cl);
    th := ch + c[i*2];
    tl := (c[i*2] - th) + ch;
    ch := th;
    cl := cl + tl + c[i*2+1];
    Dec(i);
  end;
  l := cl;
  Result := ch;
end;

{ ── pcr_hypotf ────────────────────────────────────────────────────────────── }

function pcr_hypotf(x, y: Single): Single;
var
  ax, ay, at_h, c_h: Single;
  tx_h, ty_h: Tb32u32;
  snan_x, snan_y: Integer;
  xd, yd, x2, y2, r2, r_h, cd_h, ir2, dr2, rs, dz_h, dr_h, rh, rl: Double;
  t_h: Tb64u64;
begin
  ax := pcr_fabsf(x);
  ay := pcr_fabsf(y);
  tx_h.f := ax; ty_h.f := ay;
  if (tx_h.u >= $7F800000) or (ty_h.u >= $7F800000) then begin
    snan_x := 0;
    if (tx_h.u > $7F800000) and ((tx_h.u shr 22) and 1 = 0) then snan_x := 1;
    snan_y := 0;
    if (ty_h.u > $7F800000) and ((ty_h.u shr 22) and 1 = 0) then snan_y := 1;
    if (snan_x <> 0) or (snan_y <> 0) then begin Result := x + y; Exit; end;
    if tx_h.u = $7F800000 then begin Result := ax; Exit; end;
    if ty_h.u = $7F800000 then begin Result := ay; Exit; end;
    Result := ax + ay; Exit;
  end;
  at_h := pcr_fmaxf(ax, ay);
  ay   := pcr_fminf(ax, ay);
  xd := at_h; yd := ay;
  x2 := xd * xd; y2 := yd * yd; r2 := x2 + y2;
  if yd < xd * 0.00024414061044808477 then begin  { 0x1.fffffep-13 }
    c_h := pcr_fmaf(Single(0.0001220703125), ay, at_h);  { fmaf(0x1p-13f, ay, at) }
    Result := c_h; Exit;
  end;
  r_h := pcr_sqrt(r2);
  t_h.f := r_h;
  c_h := Single(r_h);
  if t_h.u > UInt64($47EFFFFFE0000000) then begin Result := c_h; Exit; end;
  if ((t_h.u + 1) and $0FFFFFFF) > 2 then begin Result := c_h; Exit; end;
  cd_h := c_h;
  if (cd_h*cd_h - x2) - y2 = 0.0 then begin Result := c_h; Exit; end;
  ir2 := 0.5 / r2;
  dr2 := (x2 - r2) + y2;
  rs := r_h * ir2;
  dz_h := dr2 - pcr_fma(r_h, r_h, -r2);
  dr_h := rs * dz_h;
  rh := r_h + dr_h;
  rl := dr_h + (r_h - rh);
  t_h.f := rh;
  if (t_h.u and $0FFFFFFF) = 0 then begin
    if rl > 0.0 then t_h.u := t_h.u + 1;
    if rl < 0.0 then t_h.u := t_h.u - 1;
  end;
  Result := t_h.f;
end;

{ ── pcr_atan2f_tiny: Taylor approx for tiny y/x ───────────────────────────── }

function pcr_atan2f_tiny(y, x: Single): Single;
const
  c_third = -0.3333333333333333;  { -0x1.5555555555555p-2 }
var
  dy, dx, z_t, e_t, zz_t, cz_t: Double;
  t_t: Tb64u64;
begin
  dy := y; dx := x;
  z_t := dy / dx;
  e_t := pcr_fma(-z_t, dx, dy);
  zz_t := z_t * z_t;
  cz_t := c_third * z_t;
  e_t := e_t / dx + cz_t * zz_t;
  t_t.f := z_t;
  if (t_t.u and $000000000FFFFFFF) = 0 then begin
    if z_t * e_t > 0.0 then t_t.u := t_t.u + 1
    else                     t_t.u := t_t.u - 1;
  end;
  Result := t_t.f;
end;

{ ── pcr_atan2f ────────────────────────────────────────────────────────────── }

function pcr_atan2f(y, x: Single): Single;
const
  { numerator poly coeffs }
  cn_a2: array[0..6] of Double = (
    1.0, 2.506848521335565, 2.2855336234350774, 0.9227540611112051,
    0.15965700667756133, 0.0093982071883745, 8.116266383809054e-05);
  { denominator poly coeffs }
  cd_a2: array[0..6] of Double = (
    1.0, 2.840181854668896, 3.03226090832491, 1.5083284691366383,
    0.35061013533424623, 0.03311601651598859, 0.0008307046818566012);
  m_a2: array[0..1] of Double = (0.0, 1.0);
  { pi=0x1.921fb54442d18p+1, pi2=0x1.921fb54442d18p+0, pi2l=0x1.1a62633145c07p-54 }
  off_a2: array[0..7] of Double = (
    0.0, 1.5707963267948966, 3.141592653589793, 1.5707963267948966,
    0.0, -1.5707963267948966, -3.141592653589793, -1.5707963267948966);
  offl_a2: array[0..7] of Double = (
    0.0, 6.123233995736766e-17, 1.2246467991473532e-16, 6.123233995736766e-17,
    0.0, -6.123233995736766e-17, -1.2246467991473532e-16, -6.123233995736766e-17);
  sgn_a2: array[0..1] of Double = (1.0, -1.0);
  { c[32][2] high-precision Taylor table, stored flat: c[k][0]=c_a2[k*2], c[k][1]=c_a2[k*2+1] }
  c_a2: array[0..63] of Double = (
    1.0, -9.999371390980276e-27,
    -0.3333333333333333, -1.8503696081891694e-17,
    0.2, -1.1109570448589454e-17,
    -0.14285714285714285, -6.906377603268134e-18,
    0.11111111111111104, -6.04986207547041e-19,
    -0.0909090909090874, -4.813186107306423e-18,
    0.07692307692296789, 6.0635784048294865e-18,
    -0.06666666666423005, 5.357368392422845e-19,
    0.05882352937094788, -3.0234229312063074e-18,
    -0.052631578418319884, 7.121906314230013e-19,
    0.047619042181978474, 3.2606056478649954e-18,
    -0.04347821570541953, 2.1474249882470664e-18,
    0.039999692056834395, -1.4496021955960109e-18,
    -0.037035291887394496, 8.394614713005368e-19,
    0.03447445331733282, 3.57321674047701e-19,
    -0.03222458593057924, -6.05919930630252e-19,
    0.030187892413408284, 8.550359777848272e-19,
    -0.028231467664296316, 8.645756239029875e-19,
    0.026160406443643782, -1.1013132206776929e-18,
    -0.02372363694711401, 5.82053126354201e-19,
    0.0206884675894418, 8.78953163493978e-19,
    -0.016981732349686016, 9.51948710384928e-19,
    0.012817485672589527, -7.742239346647683e-19,
    -0.008687313963361743, -5.63071564570556e-19,
    0.005163618439611815, 3.2806307955849715e-19,
    -0.00262691892479603, -1.373777510443008e-19,
    0.0011135875725313154, -6.618746076431752e-20,
    -0.0003808088264492978, 1.5620808802045265e-20,
    0.00010056717515235685, -5.194627331100907e-21,
    -1.919552508092125e-05, 8.202890655529531e-22,
    2.35167210652331e-06, -1.3883775330103416e-22,
    -1.3863591848022874e-07, -1.0947593470915086e-23);
var
  tx_a2, ty_a2: Tb32u32;
  ux, uy, ax_a2, ay_a2: UInt32;
  yinf_a2, xinf_a2, gt_a2: UInt32;
  i_a2: UInt32;
  zx_a2, zy_a2, z_a2, r_a2: Double;
  d_a2: Integer;
  z2_a2, z4_a2, z8_a2: Double;
  cn0_a2, cn2_a2, cn4_a2, cn6_a2: Double;
  cd0_a2, cd2_a2, cd4_a2, cd6_a2: Double;
  res_a2: Tb64u64;
  zh_a2, zl_a2, z2l_a2, z2h_a2: Double;
  pl_a2, ph_a2, sh_a2, sl_a2: Double;
  rf_a2: Single;
  th_a2, dh_a2, tm_a2: Double;
  tth_a2: Tb64u64;
begin
  tx_a2.f := x; ty_a2.f := y;
  ux := tx_a2.u; uy := ty_a2.u;
  ax_a2 := ux and $7FFFFFFF; ay_a2 := uy and $7FFFFFFF;
  if (ay_a2 >= $7F800000) or (ax_a2 >= $7F800000) then begin
    if ay_a2 > $7F800000 then begin Result := x + y; Exit; end;  { y nan }
    if ax_a2 > $7F800000 then begin Result := x + y; Exit; end;  { x nan }
    if ay_a2 = $7F800000 then yinf_a2 := 1 else yinf_a2 := 0;
    if ax_a2 = $7F800000 then xinf_a2 := 1 else xinf_a2 := 0;
    if (yinf_a2 and xinf_a2) <> 0 then begin
      if (ux shr 31) <> 0 then
        Result := Single(2.356194490192345 * sgn_a2[uy shr 31])   { +/-3pi/4 }
      else
        Result := Single(0.7853981633974483 * sgn_a2[uy shr 31]); { +/-pi/4 }
      Exit;
    end;
    if xinf_a2 <> 0 then begin
      if (ux shr 31) <> 0 then Result := Single(3.141592653589793 * sgn_a2[uy shr 31])
      else                      Result := Single(0.0 * sgn_a2[uy shr 31]);
      Exit;
    end;
    if yinf_a2 <> 0 then begin
      Result := Single(1.5707963267948966 * sgn_a2[uy shr 31]); Exit;
    end;
  end;
  if ay_a2 = 0 then begin
    if ax_a2 = 0 then begin
      i_a2 := (uy shr 31)*4 + (ux shr 31)*2;
      if (ux shr 31) <> 0 then Result := Single(off_a2[i_a2] + offl_a2[i_a2])
      else                      Result := Single(off_a2[i_a2]);
      Exit;
    end;
    if (ux shr 31) = 0 then begin
      Result := Single(0.0 * sgn_a2[uy shr 31]); Exit;
    end;
  end;
  if ay_a2 > ax_a2 then gt_a2 := 1 else gt_a2 := 0;
  i_a2 := (uy shr 31)*4 + (ux shr 31)*2 + gt_a2;
  zx_a2 := x; zy_a2 := y;
  { z = x/y if |y|>|x|, z = y/x otherwise }
  z_a2 := (m_a2[gt_a2]*zx_a2 + m_a2[1-gt_a2]*zy_a2) /
          (m_a2[gt_a2]*zy_a2 + m_a2[1-gt_a2]*zx_a2);
  d_a2 := Integer(ax_a2) - Integer(ay_a2);
  if (d_a2 < 226492416) and (d_a2 > -226492416) then begin  { 27<<23 }
    z2_a2 := z_a2 * z_a2;
    z4_a2 := z2_a2 * z2_a2;
    z8_a2 := z4_a2 * z4_a2;
    cn0_a2 := cn_a2[0] + z2_a2*cn_a2[1];
    cn2_a2 := cn_a2[2] + z2_a2*cn_a2[3];
    cn4_a2 := cn_a2[4] + z2_a2*cn_a2[5];
    cn6_a2 := cn_a2[6];
    cn0_a2 := cn0_a2 + z4_a2*cn2_a2;
    cn4_a2 := cn4_a2 + z4_a2*cn6_a2;
    cn0_a2 := cn0_a2 + z8_a2*cn4_a2;
    cd0_a2 := cd_a2[0] + z2_a2*cd_a2[1];
    cd2_a2 := cd_a2[2] + z2_a2*cd_a2[3];
    cd4_a2 := cd_a2[4] + z2_a2*cd_a2[5];
    cd6_a2 := cd_a2[6];
    cd0_a2 := cd0_a2 + z4_a2*cd2_a2;
    cd4_a2 := cd4_a2 + z4_a2*cd6_a2;
    cd0_a2 := cd0_a2 + z8_a2*cd4_a2;
    r_a2 := cn0_a2 / cd0_a2;
  end else
    r_a2 := 1.0;
  z_a2 := z_a2 * sgn_a2[gt_a2];
  r_a2 := z_a2 * r_a2 + off_a2[i_a2];
  res_a2.f := r_a2;
  if ((res_a2.u + 8) and $0FFFFFFF) <= 16 then begin
    { check tiny y/x }
    if (ay_a2 < ax_a2) and ((ax_a2 - ay_a2) shr 23 >= 25) then begin
      Result := pcr_atan2f_tiny(y, x); Exit;
    end;
    if gt_a2 = 0 then begin
      zh_a2 := zy_a2 / zx_a2;
      zl_a2 := pcr_fma(zh_a2, -zx_a2, zy_a2) / zx_a2;
    end else begin
      zh_a2 := zx_a2 / zy_a2;
      zl_a2 := pcr_fma(zh_a2, -zy_a2, zx_a2) / zy_a2;
    end;
    z2h_a2 := muldd(zh_a2, zl_a2, zh_a2, zl_a2, z2l_a2);
    ph_a2 := polydd(z2h_a2, z2l_a2, 32, c_a2, pl_a2);
    zh_a2 := zh_a2 * sgn_a2[gt_a2];
    zl_a2 := zl_a2 * sgn_a2[gt_a2];
    ph_a2 := muldd(zh_a2, zl_a2, ph_a2, pl_a2, pl_a2);
    sh_a2 := ph_a2 + off_a2[i_a2];
    sl_a2 := ((off_a2[i_a2] - sh_a2) + ph_a2) + pl_a2 + offl_a2[i_a2];
    rf_a2 := Single(sh_a2);
    th_a2 := rf_a2;
    dh_a2 := sh_a2 - th_a2;
    tm_a2 := dh_a2 + sl_a2;
    tth_a2.f := th_a2;
    if th_a2 + th_a2*8.673617379884035e-19 = th_a2 - th_a2*8.673617379884035e-19 then begin
      { 0x1p-60 = 8.673617379884035e-19 }
      tth_a2.u := tth_a2.u and $7FF0000000000000;
      tth_a2.u := tth_a2.u - $0180000000000000;  { subtract 24<<52 }
      if pcr_fabs(tm_a2) > tth_a2.f then
        tm_a2 := tm_a2 * 1.25
      else
        tm_a2 := tm_a2 * 0.75;
    end;
    r_a2 := th_a2 + tm_a2;
  end;
  Result := Single(r_a2);
end;

{ ── pcr_atan2pif ──────────────────────────────────────────────────────────── }

function pcr_atan2pif(y, x: Single): Single;
const
  cn_a2p: array[0..6] of Double = (
    0.3183098861837907, 0.7979546675063276, 0.7275079475448462,
    0.29372174016793834, 0.05082040362397926, 0.0029915422604631704,
    2.583487828867586e-05);
  cd_a2p: array[0..6] of Double = (
    1.0, 2.840181854668896, 3.03226090832491, 1.5083284691366383,
    0.35061013533424623, 0.03311601651598859, 0.0008307046818566012);
  m_a2p: array[0..1] of Double = (0.0, 1.0);
  off_a2p: array[0..7] of Double = (
    0.0, 0.5, 1.0, 0.5, 0.0, -0.5, -1.0, -0.5);
  sgnf_a2p: array[0..1] of Single = (1.0, -1.0);
  sgn_a2p: array[0..1] of Double = (1.0, -1.0);
  s_a2p: array[0..3] of Single = (0.25, 0.75, -0.25, -0.75);
  off2_a2p: array[0..3] of Double = (0.25, 0.75, -0.25, -0.75);
  { c[32][2] flat: c_a2p[k*2]=c[k][0], c_a2p[k*2+1]=c[k][1] }
  c_a2p: array[0..63] of Double = (
    0.3183098861837907, -1.9678676678365386e-17,
    -0.1061032953945969, 6.559565574705387e-18,
    0.06366197723675814, -1.1625142324443452e-18,
    -0.04547284088339867, 1.6330933027767818e-19,
    0.035367765131532274, -9.869723873189076e-19,
    -0.02893726238034349, -9.662744851500827e-20,
    0.024485375860256887, -1.149537084009613e-18,
    -0.02122065907814378, 1.7059737274086658e-18,
    0.018724110938995286, 8.81510938898954e-19,
    -0.016753151736008654, -1.3202448604516062e-18,
    0.015157611897126696, 6.75954372683546e-19,
    -0.01383954589266639, 8.078526646412535e-19,
    0.012732297425997631, 1.7919213416430043e-19,
    -0.011788699545460008, 6.66037210651152e-19,
    0.010973559311688615, 4.52760332181136e-19,
    -0.01025740427988246, 5.715852605373118e-19,
    0.00960910459824051, -3.6666946471079842e-19,
    -0.008986355259023526, -3.20433976550847e-19,
    0.008327115997597955, 8.618694896373386e-19,
    -0.007551468176501432, -8.492307901900414e-20,
    0.006585343763712261, 2.224514552737167e-19,
    -0.005405453291432151, -3.141319316546961e-19,
    0.00407993240560434, 1.6362318030205914e-19,
    -0.0027652579189205322, 1.658177124739989e-19,
    0.0016436307978093597, -5.343077563297589e-21,
    -0.0008361742639658701, -1.949710672957039e-20,
    0.00035446593346812673, -5.638324403413975e-21,
    -0.00012121521420485887, 5.168177557928168e-21,
    3.201152607657205e-05, -1.1919997370531597e-21,
    -6.110125403746142e-06, 1.333824885359168e-22,
    7.485604805690301e-07, -9.591512195984102e-24,
    -4.4129183432426896e-08, 1.786021469414917e-24);
var
  tx_a2p, ty_a2p: Tb32u32;
  ux_p, uy_p, ax_a2p, ay_a2p: UInt32;
  yinf_a2p, xinf_a2p, gt_a2p: UInt32;
  i_a2p: UInt32;
  zx_a2p, zy_a2p, z_a2p, r_a2p, z2_a2p: Double;
  z4_a2p, z8_a2p: Double;
  cn0_p, cn2_p, cn4_p, cn6_p: Double;
  cd0_p, cd2_p, cd4_p, cd6_p: Double;
  res_a2p: Tb64u64;
  zh_a2p, zl_a2p, z2l_a2p, z2h_a2p: Double;
  pl_a2p, ph_a2p, sh_a2p, sl_a2p: Double;
  rf_a2p: Single;
  th_a2p, dh_a2p, tm_a2p: Double;
  d_a2p: Tb64u64;
begin
  tx_a2p.f := x; ty_a2p.f := y;
  ux_p := tx_a2p.u; uy_p := ty_a2p.u;
  ax_a2p := ux_p and $7FFFFFFF; ay_a2p := uy_p and $7FFFFFFF;
  if (ay_a2p >= $7F800000) or (ax_a2p >= $7F800000) then begin
    if ay_a2p > $7F800000 then begin Result := x + y; Exit; end;
    if ax_a2p > $7F800000 then begin Result := x + y; Exit; end;
    if ay_a2p = $7F800000 then yinf_a2p := 1 else yinf_a2p := 0;
    if ax_a2p = $7F800000 then xinf_a2p := 1 else xinf_a2p := 0;
    if (yinf_a2p and xinf_a2p) <> 0 then begin
      if (ux_p shr 31) <> 0 then Result := 0.75 * sgnf_a2p[uy_p shr 31]
      else                        Result := 0.25 * sgnf_a2p[uy_p shr 31];
      Exit;
    end;
    if xinf_a2p <> 0 then begin
      if (ux_p shr 31) <> 0 then Result := sgnf_a2p[uy_p shr 31]
      else                        Result := 0.0 * sgnf_a2p[uy_p shr 31];
      Exit;
    end;
    if yinf_a2p <> 0 then begin
      Result := 0.5 * sgnf_a2p[uy_p shr 31]; Exit;
    end;
  end;
  if ay_a2p = 0 then begin
    if (ay_a2p or ax_a2p) = 0 then begin
      i_a2p := (uy_p shr 31)*4 + (ux_p shr 31)*2;
      Result := Single(off_a2p[i_a2p]); Exit;
    end;
    if (ux_p shr 31) = 0 then begin
      Result := 0.0 * sgnf_a2p[uy_p shr 31]; Exit;
    end;
  end;
  if ax_a2p = ay_a2p then begin
    i_a2p := (uy_p shr 31)*2 + (ux_p shr 31);
    Result := s_a2p[i_a2p]; Exit;
  end;
  if ay_a2p > ax_a2p then gt_a2p := 1 else gt_a2p := 0;
  i_a2p := (uy_p shr 31)*4 + (ux_p shr 31)*2 + gt_a2p;
  zx_a2p := x; zy_a2p := y;
  z_a2p := (m_a2p[gt_a2p]*zx_a2p + m_a2p[1-gt_a2p]*zy_a2p) /
           (m_a2p[gt_a2p]*zy_a2p + m_a2p[1-gt_a2p]*zx_a2p);
  r_a2p := cn_a2p[0];
  z2_a2p := z_a2p * z_a2p;
  z_a2p := z_a2p * sgn_a2p[gt_a2p];
  if z2_a2p > 5.551115123125783e-17 then begin  { 0x1p-54 }
    z4_a2p := z2_a2p * z2_a2p;
    z8_a2p := z4_a2p * z4_a2p;
    cn0_p := r_a2p    + z2_a2p*cn_a2p[1];
    cn2_p := cn_a2p[2] + z2_a2p*cn_a2p[3];
    cn4_p := cn_a2p[4] + z2_a2p*cn_a2p[5];
    cn6_p := cn_a2p[6];
    cn0_p := cn0_p + z4_a2p*cn2_p;
    cn4_p := cn4_p + z4_a2p*cn6_p;
    cn0_p := cn0_p + z8_a2p*cn4_p;
    cd0_p := cd_a2p[0] + z2_a2p*cd_a2p[1];
    cd2_p := cd_a2p[2] + z2_a2p*cd_a2p[3];
    cd4_p := cd_a2p[4] + z2_a2p*cd_a2p[5];
    cd6_p := cd_a2p[6];
    cd0_p := cd0_p + z4_a2p*cd2_p;
    cd4_p := cd4_p + z4_a2p*cd6_p;
    cd0_p := cd0_p + z8_a2p*cd4_p;
    r_a2p := cn0_p / cd0_p;
  end;
  r_a2p := z_a2p * r_a2p + off_a2p[i_a2p];
  res_a2p.f := r_a2p;
  if ((res_a2p.u shl 1) > UInt64($6D40000000000000)) and
     (((res_a2p.u + 8) and $0FFFFFFF) <= 16) then begin
    if ax_a2p = ay_a2p then begin
      r_a2p := off2_a2p[(uy_p shr 31)*2 + (ux_p shr 31)];
    end else begin
      if gt_a2p = 0 then begin
        zh_a2p := zy_a2p / zx_a2p;
        zl_a2p := pcr_fma(zh_a2p, -zx_a2p, zy_a2p) / zx_a2p;
      end else begin
        zh_a2p := zx_a2p / zy_a2p;
        zl_a2p := pcr_fma(zh_a2p, -zy_a2p, zx_a2p) / zy_a2p;
      end;
      z2h_a2p := muldd(zh_a2p, zl_a2p, zh_a2p, zl_a2p, z2l_a2p);
      ph_a2p := polydd(z2h_a2p, z2l_a2p, 32, c_a2p, pl_a2p);
      zh_a2p := zh_a2p * sgn_a2p[gt_a2p];
      zl_a2p := zl_a2p * sgn_a2p[gt_a2p];
      ph_a2p := muldd(zh_a2p, zl_a2p, ph_a2p, pl_a2p, pl_a2p);
      sh_a2p := ph_a2p + off_a2p[i_a2p];
      sl_a2p := ((off_a2p[i_a2p] - sh_a2p) + ph_a2p) + pl_a2p;
      rf_a2p := Single(sh_a2p);
      th_a2p := rf_a2p;
      dh_a2p := sh_a2p - th_a2p;
      tm_a2p := dh_a2p + sl_a2p;
      r_a2p := th_a2p + tm_a2p;
      d_a2p.f := r_a2p - th_a2p;
      if (d_a2p.u shl 12) = 0 then begin
        if pcr_fabs(d_a2p.f) > pcr_fabs(tm_a2p) then
          r_a2p := r_a2p - d_a2p.f * 9.765625e-04  { 0x1p-10 }
        else if pcr_fabs(d_a2p.f) < pcr_fabs(tm_a2p) then
          r_a2p := r_a2p + d_a2p.f * 9.765625e-04;
      end;
    end;
  end;
  Result := Single(r_a2p);
end;


{ ── pcr_powf helpers ───────────────────────────────────────────────────── }

function mulddd_pf(xh, xl, ch: Double; out l: Double): Double;
var ahlh, ahhh, ahhl: Double;
begin
  ahlh := ch * xl;
  ahhh := ch * xh;
  ahhl := pcr_fma(ch, xh, -ahhh);
  ahhl := ahhl + ahlh;
  ch := ahhh + ahhl;
  l := (ahhh - ch) + ahhl;
  Result := ch;
end;

function isint_pf(y0: Single): Boolean;
var wy: Tb32u32; ey, s: Integer;
begin
  wy.f := y0;
  ey := Integer((wy.u shr 23) and $FF) - 127;
  s := ey + 9;
  if ey >= 0 then begin
    if s >= 32 then begin Result := True; Exit; end;
    Result := (wy.u shl s) = 0;
    Exit;
  end;
  Result := (wy.u shl 1) = 0;
end;

function isodd_pf(y0: Single): Boolean;
var wy: Tb32u32; ey, s: Integer; oddb: UInt32;
begin
  wy.f := y0;
  ey := Integer((wy.u shr 23) and $FF) - 127;
  s := ey + 9;
  oddb := 0;
  if ey >= 0 then begin
    if (s < 32) and ((wy.u shl s) = 0) then oddb := (wy.u shr (32-s)) and 1;
    if s = 32 then oddb := wy.u and 1;
  end;
  Result := oddb <> 0;
end;

function is_signalingf_pf(x: Single): Boolean;
var u: Tb32u32;
begin
  u.f := x;
  u.u := u.u xor $00400000;
  Result := (u.u and $7FFFFFFF) > $7FC00000;
end;

function is_exact_pf(x0, y0: Single): Integer;
const
  xmax_ie: array[0..15] of UInt32 = (0, $FFFFFF, 4095, 255, 63, 27, 15, 9,
                                      7, 5, 5, 3, 3, 3, 3, 3);
var
  vw, ww: Tb32u32;
  m_ie, n_ie: UInt32;
  t_ie, y_int_ie: Integer;
  e_ie, f_ie: Int32;
  my64: UInt64;
  my32, n0_ie: UInt32;
  ez_ie: Int32;
  dm_ie: Single;
  sf_ie: Single;
begin
  vw.f := x0; ww.f := y0;
  { Early exit if |x|<>1 and low 16 bits of y.u are non-zero }
  if ((vw.u shl 1) <> $7F000000) and ((ww.u shl 16) <> 0) then begin
    Result := 0; Exit;
  end;
  { |x| = 1 }
  if (vw.u shl 1) = $7F000000 then begin
    Result := 1; Exit;
  end;
  { y >= 0 and y is an integer }
  if (y0 >= 0) and isint_pf(y0) then begin
    m_ie := vw.u and $7FFFFF;
    e_ie := Int32((vw.u shl 1) shr 24) - $96;
    if e_ie >= -149 then m_ie := m_ie or $800000
    else Inc(e_ie);
    t_ie := pcr_bsf32(m_ie);
    m_ie := m_ie shr t_ie;
    e_ie := e_ie + t_ie;
    if (y0 = 0) or (y0 = 1) then begin Result := 1; Exit; end;
    if m_ie = 1 then begin
      Result := Ord((-149 <= Int32(Trunc(y0)) * e_ie) and
                    (Int32(Trunc(y0)) * e_ie < 128));
      Exit;
    end;
    if (y0 < 0) or (y0 > 15) then begin Result := 0; Exit; end;
    y_int_ie := Integer(Trunc(y0));
    if m_ie > xmax_ie[y_int_ie] then begin Result := 0; Exit; end;
    my64 := UInt64(m_ie) * UInt64(m_ie);
    t_ie := 2;
    while t_ie < y_int_ie do begin my64 := my64 * UInt64(m_ie); Inc(t_ie); end;
    t_ie := pcr_bsr32(m_ie) + 1;  { 32 - clz(m_ie) }
    ez_ie := e_ie * y_int_ie + t_ie;
    if (ez_ie <= -149) or (128 < ez_ie) then begin Result := 0; Exit; end;
    Result := Ord(e_ie * y_int_ie >= -149);
    Exit;
  end;
  { Decompose |y| = n * 2^f with n odd }
  n_ie := ww.u and $7FFFFF;
  f_ie := Int32((ww.u shl 1) shr 24) - $96;
  if f_ie >= -149 then n_ie := n_ie or $800000
  else Inc(f_ie);
  t_ie := pcr_bsf32(n_ie);
  n_ie := n_ie shr t_ie;
  f_ie := f_ie + t_ie;
  { Decompose |x| = m * 2^e with m odd }
  m_ie := vw.u and $7FFFFF;
  e_ie := Int32((vw.u shl 1) shr 24) - $96;
  if e_ie >= -149 then m_ie := m_ie or $800000
  else Inc(e_ie);
  t_ie := pcr_bsf32(m_ie);
  m_ie := m_ie shr t_ie;
  e_ie := e_ie + t_ie;
  if y0 < 0 then begin
    if m_ie <> 1 then begin Result := 0; Exit; end;
    if f_ie >= 0 then begin
      if e_ie >= 0 then ez_ie := -(Int32(e_ie) shl f_ie) * Int32(n_ie)
      else ez_ie := (Int32(-e_ie) shl f_ie) * Int32(n_ie);
    end else begin
      t_ie := pcr_bsf32(UInt32(e_ie));
      if (-f_ie) > t_ie then begin Result := 0; Exit; end;
      ez_ie := sar_i32(-e_ie, -f_ie) * Int32(n_ie);
    end;
    Result := Ord((-149 <= ez_ie) and (ez_ie < 128));
    Exit;
  end;
  { y > 0, not integer, y = n*2^f with n odd and f < 0 }
  while f_ie <> 0 do begin
    Inc(f_ie);
    if (e_ie and 1) <> 0 then begin Result := 0; Exit; end;
    e_ie := e_ie div 2;
    dm_ie := Single(m_ie);
    sf_ie := Single(Round(pcr_sqrtf(dm_ie)));
    if sf_ie * sf_ie <> dm_ie then begin Result := 0; Exit; end;
    m_ie := UInt32(Round(sf_ie));
  end;
  if m_ie > 1 then begin
    if n_ie > 15 then begin Result := 0; Exit; end;
    if m_ie > xmax_ie[n_ie] then begin Result := 0; Exit; end;
  end;
  my32 := m_ie; n0_ie := n_ie;
  while n0_ie > 1 do begin Dec(n0_ie); my32 := my32 * m_ie; end;
  t_ie := pcr_bsr32(my32) + 1;  { 32 - clz(my32) }
  Result := Ord((-149 <= e_ie * Int32(n_ie)) and
                (e_ie * Int32(n_ie) + t_ie <= 128));
end;

{ ── pcr_powf_accurate2 ──────────────────────────────────────────────────── }

function pcr_powf_accurate2(x0, y0: Single; is_exact_val: Integer): Single;
const
  o_a2: array[0..1] of Double = (1.0, 2.0);
  { ch[][2] flat: 13 pairs }
  ch_a2: array[0..25] of Double = (
    2.8853900817779268,   4.071054748191002e-17,    { 0x1.71547652b82fep+1, 0x1.777d0ffda2b89p-55 }
    0.9617966939259756,   5.057761609823965e-17,    { 0x1.ec709dc3a03fdp-1, 0x1.d27f04ff73b3ap-55 }
    0.5770780163555853,   5.2552074957089844e-17,   { 0x1.2776c50ef9bfep-1, 0x1.e4b514251d0ecp-55 }
    0.4121985831111324,   1.2966716928545934e-17,   { 0x1.a61762a7aded9p-2, 0x1.de632dc7f6998p-57 }
    0.3205988979753255,   2.2721007160373763e-17,   { 0x1.484b13d7c02aep-2, 0x1.a320ec342ddb3p-56 }
    0.26230818925246924, -1.3180082845309455e-17,   { 0x1.0c9a84993fd48p-2, -0x1.e6425ce9a74a4p-57 }
    0.22195308322397042, -7.03730411210251e-18,     { 0x1.c68f568d8beafp-3, -0x1.03a175487feabp-57 }
    0.19235933776779732,  1.3451966316321192e-17,   { 0x1.89f3b14657dfbp-3, 0x1.f04a3acf0bcf7p-57 }
    0.16972889712828465, -4.892356301506132e-18,    { 0x1.5b9ad2f2d12ap-3, -0x1.68fdff6815a6fp-58 }
    0.15185945002402543,  1.1034764462533595e-18,   { 0x1.3702165b88acbp-3, 0x1.45b052ace6c8ep-60 }
    0.13749878152691522, -1.0236531880234422e-17,   { 0x1.1998f60f2f005p-3, -0x1.79a94f62fb524p-57 }
    0.12347055433477307, -2.2899638320795317e-18,   { 0x1.f9bc428e30809p-4, -0x1.51f063387e47p-59 }
    0.13806280141791244,  8.12205453916993e-18      { 0x1.1ac0ab871296ap-3, 0x1.2ba6a2e1a625bp-57 }
  );
  { ce[][2] flat: 18 pairs }
  ce_a2: array[0..35] of Double = (
    1.0,                   6.210306603644812e-30,   { 0x1p+0, 0x1.f7d70599926c4p-98 }
    0.6931471805599453,    2.3190468138467075e-17,  { 0x1.62e42fefa39efp-1, 0x1.abc9e3b39856bp-56 }
    0.24022650695910072,  -9.493931257207092e-18,   { 0x1.ebfbdff82c58fp-3, -0x1.5e43a540c283dp-57 }
    0.05550410866482158,  -3.1658222912778202e-18,  { 0x1.c6b08d704a0cp-5, -0x1.d3316277451e6p-59 }
    0.009618129107628477,  2.8324649708472296e-19,  { 0x1.3b2ab6fba4e77p-7, 0x1.4e66003ba7f85p-62 }
    0.0013333558146428443, 1.392811665254468e-20,   { 0x1.5d87fe78a6731p-10, 0x1.07183d46a9697p-66 }
    0.0001540353039338161, 1.1765991431639673e-20,  { 0x1.430912f86c787p-13, 0x1.bc81afca4c93p-67 }
    1.5252733804059841e-05,-8.0442948550066915e-22, { 0x1.ffcbfc588b0c7p-17, -0x1.e63f6f0116f4cp-71 }
    1.3215486790144314e-06,-6.293372509344185e-23,  { 0x1.62c0223a5c826p-20, -0x1.30542d98ea4a5p-74 }
    1.0178086009239703e-07,-1.3006186993534895e-24, { 0x1.b5253d395e7c6p-24, -0x1.9285a132ce05ep-80 }
    7.054911620796934e-09, -1.659037271615654e-25,  { 0x1.e4cf5158b7b01p-28, -0x1.9ac1facae1b88p-83 }
    4.4455382718682324e-10,-5.296718261089315e-28,  { 0x1.e8cac7351a7a8p-32, -0x1.4fb82adebd76bp-91 }
    2.5678436021925767e-11, 6.132224146744111e-28,  { 0x1.c3bd65182746dp-36, 0x1.84ad0689d30ep-91 }
    1.3691488868804127e-12,-2.8921457069199386e-29, { 0x1.8161931d765c3p-40, -0x1.254c6535279cep-95 }
    6.778715106350512e-14, -6.174675424236808e-30,  { 0x1.314943a26c9e2p-44, -0x1.f4f2fdc14fb82p-98 }
    3.132431569193972e-15,  1.9136583918642863e-31, { 0x1.c36e53b459602p-49, 0x1.f0d06a5a63c41p-103 }
    1.3594238037092167e-16,-8.23812497323582e-33,   { 0x1.397637b3876a4p-53, -0x1.5632c551ae458p-107 }
    5.542771640076379e-18, -4.787321441390354e-35   { 0x1.98fbfefdddb51p-58, -0x1.fd134923d52b4p-115 }
  );
var
  x_a2, y_a2: Double;
  t_a2: Tb64u64;
  e_a2, k_a2: Integer;
  xm_a2, xp_a2, zh_a2, zl_a2, z2l_a2, z2h_a2: Double;
  ey_a2, eh_a2, el_a2, ee_a2: Double;
  r_a2: Tb64u64;
  ty_a2: Tb32u32;
  et_a2: Integer;
  kk_a2: UInt32;
  isintflag_a2: Boolean;
  ll_a2, lh_a2: Tb64u64;
  res_a2: Single;
begin
  x_a2 := x0; y_a2 := y0;
  t_a2.f := x_a2;
  e_a2 := Integer((t_a2.u shr 52) and $7FF) - $3FF;
  t_a2.u := t_a2.u and ($FFFFFFFFFFFFFFFF shr 12);  { clear exponent+sign }
  k_a2 := Ord(t_a2.u > UInt64($6a09e667f3bcd));
  e_a2 := e_a2 + k_a2;
  t_a2.u := t_a2.u or (UInt64($3FF) shl 52);
  x_a2 := t_a2.f;
  xm_a2 := x_a2 - o_a2[k_a2];
  xp_a2 := x_a2 + o_a2[k_a2];
  zh_a2 := xm_a2 / xp_a2;
  zl_a2 := pcr_fma(zh_a2, -xp_a2, xm_a2) / xp_a2;
  z2h_a2 := muldd(zh_a2, zl_a2, zh_a2, zl_a2, z2l_a2);
  z2h_a2 := polydd(z2h_a2, z2l_a2, 13, ch_a2, z2l_a2);
  zh_a2 := muldd(zh_a2, zl_a2, z2h_a2, z2l_a2, zl_a2);
  zh_a2 := mulddd_pf(zh_a2, zl_a2, y_a2, zl_a2);
  ey_a2 := Double(e_a2) * y_a2;
  eh_a2 := ey_a2 + zh_a2;
  el_a2 := ((ey_a2 - eh_a2) + zh_a2) + zl_a2;
  ee_a2 := pcr_roundeven(eh_a2);
  eh_a2 := eh_a2 - ee_a2;
  eh_a2 := polydd(eh_a2, el_a2, 18, ce_a2, el_a2);
  r_a2.u := UInt64(Int64($3FF) + Int64(Trunc(ee_a2))) shl 52;
  { Check if y is an odd integer }
  ty_a2.f := y0;
  et_a2 := Integer((ty_a2.u shr 23) and $FF) - $7F;
  if (8 + et_a2) >= 0 then kk_a2 := ty_a2.u shl (8 + et_a2)
  else kk_a2 := ty_a2.u shr (-8 - et_a2);
  isintflag_a2 := (((kk_a2 shl 1) or UInt32(sar_i32(et_a2, 31))) = 0) or (et_a2 >= 23);
  ll_a2.f := el_a2;
  lh_a2.f := eh_a2;
  { Adjust for borderline rounding }
  if ((ll_a2.u shr 23) and UInt64($1FFFFFFF)) = UInt64($1FFFFFFF) then begin
    if eh_a2 < 1 then begin
      if el_a2 >= 5.551115123125783e-17 then begin    { 0x1p-54 }
        el_a2 := el_a2 - 1.1102230246251565e-16;     { 0x1p-53 }
        eh_a2 := eh_a2 + 1.1102230246251565e-16;
      end else if el_a2 <= -5.551115123125783e-17 then begin
        el_a2 := el_a2 + 1.1102230246251565e-16;
        eh_a2 := eh_a2 - 1.1102230246251565e-16;
      end;
    end else begin
      if el_a2 >= 1.1102230246251565e-16 then begin   { 0x1p-53 }
        el_a2 := el_a2 - 2.220446049250313e-16;      { 0x1p-52 }
        eh_a2 := eh_a2 + 2.220446049250313e-16;
      end else if el_a2 <= -1.1102230246251565e-16 then begin
        el_a2 := el_a2 + 2.220446049250313e-16;
        eh_a2 := eh_a2 - 2.220446049250313e-16;
      end;
    end;
  end else if ((ll_a2.u shr 23) and UInt64($1FFFFFFF)) = 0 then begin
    if el_a2 > 0 then begin
      if eh_a2 < 1 then begin
        if el_a2 >= 1.1102230246251565e-16 then begin  { 0x1p-53 }
          el_a2 := el_a2 - 1.1102230246251565e-16;
          eh_a2 := eh_a2 + 1.1102230246251565e-16;
        end;
      end else begin
        if el_a2 >= 2.220446049250313e-16 then begin   { 0x1p-52 }
          el_a2 := el_a2 - 2.220446049250313e-16;
          eh_a2 := eh_a2 + 2.220446049250313e-16;
        end;
      end;
    end else begin
      if eh_a2 < 1 then begin
        if el_a2 <= -1.1102230246251565e-16 then begin
          el_a2 := el_a2 + 1.1102230246251565e-16;
          eh_a2 := eh_a2 - 1.1102230246251565e-16;
        end;
      end else begin
        if el_a2 <= -2.220446049250313e-16 then begin
          el_a2 := el_a2 + 2.220446049250313e-16;
          eh_a2 := eh_a2 - 2.220446049250313e-16;
        end;
      end;
    end;
  end;
  ll_a2.f := el_a2;
  lh_a2.f := eh_a2;
  if (lh_a2.u and $FFFFFFF) = 0 then begin
    if pcr_fabs(ll_a2.f) > 4.0389678347315804e-28 then begin  { 0x1p-91 }
      if el_a2 < 0 then begin
        lh_a2.u := lh_a2.u - 1;
        eh_a2 := lh_a2.f;
      end else begin
        lh_a2.u := lh_a2.u + 1;
        eh_a2 := lh_a2.f;
      end;
    end;
  end;
  eh_a2 := eh_a2 * r_a2.f;
  el_a2 := el_a2 * r_a2.f;
  if isintflag_a2 and (kk_a2 <> 0) then
    eh_a2 := pcr_copysign(eh_a2, Double(x0));
  res_a2 := Single(eh_a2);
  Result := res_a2;
end;

{ ── pcr_powf ─────────────────────────────────────────────────────────────── }

function pcr_powf(x0, y0: Single): Single;
const
  ix_pf: array[0..32] of Double = (
    1.0,                   { 0x1p+0 }
    0.9696969696960878,    { 0x1.f07c1f07cp-1 }
    0.9411764705873793,    { 0x1.e1e1e1e1ep-1 }
    0.9142857142869616,    { 0x1.d41d41d42p-1 }
    0.8888888888905058,    { 0x1.c71c71c72p-1 }
    0.8648648648668313,    { 0x1.bacf914c2p-1 }
    0.8421052631601924,    { 0x1.af286bca2p-1 }
    0.8205128205154324,    { 0x1.a41a41a42p-1 }
    0.8000000000029104,    { 0x1.99999999ap-1 }
    0.7804878048773389,    { 0x1.8f9c18f9cp-1 }
    0.7619047619082266,    { 0x1.861861862p-1 }
    0.7441860465114587,    { 0x1.7d05f417dp-1 }
    0.7272727272720658,    { 0x1.745d1745dp-1 }
    0.711111111108039,     { 0x1.6c16c16c1p-1 }
    0.6956521739120944,    { 0x1.642c8590bp-1 }
    0.6808510638293228,    { 0x1.5c9882b93p-1 }
    0.6666666666642413,    { 0x1.555555555p-1 }
    0.6530612244896474,    { 0x1.4e5e0a72fp-1 }
    0.6399999999994179,    { 0x1.47ae147aep-1 }
    0.6274509803915862,    { 0x1.414141414p-1 }
    0.6153846153829363,    { 0x1.3b13b13b1p-1 }
    0.6037735849022283,    { 0x1.3521cfb2bp-1 }
    0.5925925925912452,    { 0x1.2f684bda1p-1 }
    0.5818181818176527,    { 0x1.29e4129e4p-1 }
    0.571428571427532,     { 0x1.249249249p-1 }
    0.5614035087710363,    { 0x1.1f7047dc1p-1 }
    0.551724137927522,     { 0x1.1a7b9611ap-1 }
    0.5423728813548223,    { 0x1.15b1e5f75p-1 }
    0.5333333333328483,    { 0x1.111111111p-1 }
    0.5245901639354997,    { 0x1.0c9714fbdp-1 }
    0.5161290322575951,    { 0x1.084210842p-1 }
    0.5079365079363924,    { 0x1.041041041p-1 }
    0.5                    { 0x1p-1 }
  );
  { lix[][2] flat: 33 pairs; lix[j][0]=lix[j*2], lix[j][1]=lix[j*2+1] }
  lix_pf: array[0..65] of Double = (
     0.0,                   0.0,                     { j=0 }
    -0.04443359375,         3.9474390234438854e-05,  { j=1 }
    -0.08740234375,        -6.049750165153175e-05,   { j=2 }
    -0.12890625,           -3.767669429982701e-04,   { j=3 }
    -0.169921875,          -3.1264396881159154e-06,  { j=4 }
    -0.208984375,          -4.689906256694731e-04,   { j=5 }
    -0.248046875,           1.193615603508767e-04,   { j=6 }
    -0.28515625,           -2.459688576559096e-04,   { j=7 }
    -0.322265625,           3.375301178861461e-04,   { j=8 }
    -0.357421875,          -1.3012961939581666e-04,  { j=9 }
    -0.392578125,           2.607022278003286e-04,   { j=10 }
    -0.42578125,           -4.8350470242596976e-04,  { j=11 }
    -0.458984375,          -4.472436386093797e-04,   { j=12 }
    -0.4921875,             3.3440366409270264e-04,  { j=13 }
     0.4765625,            -1.2445605898105753e-04,  { j=14 }
     0.4453125,             9.86483213785352e-05,    { j=15 }
     0.4140625,             9.749992735953246e-04,   { j=16 }
     0.384765625,           5.245308844637542e-04,   { j=17 }
     0.35546875,            6.750602239631808e-04,   { j=18 }
     0.328125,             -5.503419728077132e-04,   { j=19 }
     0.298828125,           7.321568549714691e-04,   { j=20 }
     0.271484375,           5.951704286000489e-04,   { j=21 }
     0.2451171875,         -4.689666748853101e-06,   { j=22 }
     0.21875,              -1.0971352597172757e-04,  { j=23 }
     0.1923828125,          2.6226543977164556e-04,  { j=24 }
     0.1669921875,          1.1779833296210864e-04,  { j=25 }
     0.1416015625,          4.174423632430148e-04,   { j=26 }
     0.1171875,             1.6945063520646333e-04,  { j=27 }
     0.09326171875,        -1.5231435983065282e-04,  { j=28 }
     0.0693359375,         -7.32750599339986e-05,    { j=29 }
     0.0458984375,         -9.47478881873323e-05,    { j=30 }
     0.022705078125,        1.4998374755498776e-05,  { j=31 }
     0.0,                   0.0                      { j=32 }
  );
  { c[]: log2 polynomial, 8 coefficients }
  c_pf: array[0..7] of Double = (
     1.4426950408889634,   { 0x1.71547652b82fep+0 }
    -0.7213475204444817,   { -0x1.71547652b82fep-1 }
     0.4808983469635712,   { 0x1.ec709dc3a2d0bp-2 }
    -0.36067376022317404,  { -0x1.71547652bc4a9p-2 }
     0.28853899623008594,  { 0x1.2776c441b72ep-2 }
    -0.24044915938489397,  { -0x1.ec709bdf453ecp-3 }
     0.20617758474822143,  { 0x1.a6406efd4b877p-3 }
    -0.1804151705675918    { -0x1.717d824a520f7p-3 }
  );
  { ce[]: exp2 polynomial, 6 coefficients }
  ce_pf: array[0..5] of Double = (
    0.043321698784995886,  { 0x1.62e42fefa398bp-5 }
    0.0009383847928200837, { 0x1.ebfbdff84555ap-11 }
    1.3550807712983854e-05,{ 0x1.c6b08d4ad86d3p-17 }
    1.4676119301623784e-07,{ 0x1.3b2ad1b1716a2p-23 }
    1.271309415715539e-09, { 0x1.5d7472718ce9dp-30 }
    9.382438953978075e-12  { 0x1.4a1d7f457ac56p-37 }
  );
  { tb[]: 2^(j/16) for j=0..15 }
  tb_pf: array[0..15] of Double = (
    1.0,                   { 0x1p+0 }
    1.0442737824274138,    { 0x1.0b5586cf9890fp+0 }
    1.0905077326652577,    { 0x1.172b83c7d517bp+0 }
    1.1387886347566916,    { 0x1.2387a6e756238p+0 }
    1.189207115002721,     { 0x1.306fe0a31b715p+0 }
    1.241857812073484,     { 0x1.3dea64c123422p+0 }
    1.2968395546510096,    { 0x1.4bfdad5362a27p+0 }
    1.3542555469368927,    { 0x1.5ab07dd485429p+0 }
    1.4142135623730951,    { 0x1.6a09e667f3bcdp+0 }
    1.4768261459394993,    { 0x1.7a11473eb0187p+0 }
    1.5422108254079407,    { 0x1.8ace5422aa0dbp+0 }
    1.6104903319492543,    { 0x1.9c49182a3f09p+0 }
    1.681792830507429,     { 0x1.ae89f995ad3adp+0 }
    1.7562521603732995,    { 0x1.c199bdd85529cp+0 }
    1.8340080864093424,    { 0x1.d5818dcfba487p+0 }
    1.9152065613971474     { 0x1.ea4afa2a490dap+0 }
  );
var
  x_pf, y_pf: Double;
  tx_pf, ty_pf: Tb64u64;
  m_pf: UInt64;
  e_pf, j_pf, k_pf: Integer;
  xd_pf: Tb64u64;
  z_pf, z2_pf, z4_pf: Double;
  c6_pf, c4_pf, c2_pf, c0_pf: Double;
  l_pf, zt_pf: Double;
  ia_pf: Double;
  h_pf, h2_pf: Double;
  il_pf, jl_pf, el_pf: Int64;
  s_pf: Double;
  su_pf: Tb64u64;
  w_pf: Double;
  rr_pf: Tb64u64;
  off_pf: UInt64;
  et_pf: Integer;
  kk_pf: UInt64;
  res_pf: Single;
begin
  x_pf := x0; y_pf := y0;
  tx_pf.f := x_pf; ty_pf.f := y_pf;
  { |x| = 1 }
  if (tx_pf.u shl 1) = (UInt64($3FF) shl 53) then begin
    if (tx_pf.u shr 63) <> 0 then begin  { x = -1 }
      if (ty_pf.u shl 1) > (UInt64($7FF) shl 53) then begin
        Result := y0 + y0; Exit;  { y = NaN }
      end;
      if isint_pf(y0) then begin
        if isodd_pf(y0) then Result := x0
        else Result := -x0;
        Exit;
      end;
      Result := (x_pf - x_pf) / (x_pf - x_pf); Exit;  { NaN }
    end;
    { x = +1 }
    if is_signalingf_pf(y0) then Result := x0 + y0
    else Result := x0;
    Exit;
  end;
  { y = 0 }
  if (ty_pf.u shl 1) = 0 then begin
    if is_signalingf_pf(x0) then Result := x0 + y0
    else Result := 1.0;
    Exit;
  end;
  { y = Inf or NaN }
  if (ty_pf.u shl 1) >= (UInt64($7FF) shl 53) then begin
    if (tx_pf.u shl 1) > (UInt64($7FF) shl 53) then begin
      Result := x0 + y0; Exit;  { x = NaN }
    end;
    if (ty_pf.u shl 1) = (UInt64($7FF) shl 53) then begin  { y = +-Inf }
      if (((tx_pf.u shl 1) < (UInt64($3FF) shl 53)) xor ((ty_pf.u shr 63) <> 0)) then
        Result := 0.0
      else
        Result := Single(1.0/0.0);
      Exit;
    end;
    Result := x0 + y0; Exit;  { y = NaN }
  end;
  { x = Inf, NaN, or negative }
  if tx_pf.u >= (UInt64($7FF) shl 52) then begin
    if (tx_pf.u shl 1) = (UInt64($7FF) shl 53) then begin  { x = +-Inf }
      if not isodd_pf(y0) then x0 := pcr_fabsf(x0);
      if (ty_pf.u shr 63) <> 0 then Result := 1.0/x0
      else Result := x0;
      Exit;
    end;
    if (tx_pf.u shl 1) > (UInt64($7FF) shl 53) then begin
      Result := x0 + x0; Exit;  { x = NaN }
    end;
    if tx_pf.u > (UInt64($7FF) shl 52) then begin  { x <= 0 }
      if (not isint_pf(y0)) and (x_pf <> 0) then begin
        Result := (x_pf - x_pf) / (x_pf - x_pf); Exit;  { NaN }
      end;
    end;
  end;
  { x = +0 or -0 }
  if (tx_pf.u shl 1) = 0 then begin
    if (ty_pf.u shr 63) <> 0 then begin  { y < 0 }
      if isodd_pf(y0) then Result := 1.0 / pcr_copysignf(0.0, x0)
      else Result := 1.0 / 0.0;
    end else begin  { y > 0 }
      if isodd_pf(y0) then Result := pcr_copysignf(1.0, x0) * 0.0
      else Result := 0.0;
    end;
    Exit;
  end;
  { Main path: x > 0 finite, y finite nonzero }
  m_pf := tx_pf.u and (not UInt64(0)) shr 12;  { 52-bit mantissa }
  e_pf := Integer((tx_pf.u shr 52) and $7FF) - $3FF;
  j_pf := Integer((Int64(m_pf) + (Int64(1) shl 46)) shr 47);
  k_pf := Ord(j_pf > 13);
  e_pf := e_pf + k_pf;
  xd_pf.u := m_pf or (UInt64($3FF) shl 52);
  z_pf := pcr_fma(xd_pf.f, ix_pf[j_pf], -1.0);
  z2_pf := z_pf * z_pf; z4_pf := z2_pf * z2_pf;
  c6_pf := c_pf[6] + z_pf * c_pf[7];
  c4_pf := c_pf[4] + z_pf * c_pf[5];
  c2_pf := c_pf[2] + z_pf * c_pf[3];
  c0_pf := c_pf[0] + z_pf * c_pf[1];
  c0_pf := c0_pf + z2_pf * c2_pf;
  c4_pf := c4_pf + z2_pf * c6_pf;
  c0_pf := c0_pf + z4_pf * c4_pf;
  l_pf := z_pf * c0_pf - lix_pf[j_pf*2+1];
  y_pf := y_pf * 16.0;
  zt_pf := (Double(e_pf) - lix_pf[j_pf*2]) * y_pf;
  z_pf := l_pf * y_pf + zt_pf;
  { Overflow check }
  if z_pf > 2048.0 then begin
    if isodd_pf(y0) then
      Result := pcr_copysignf(Single(1.7014118346046923e+38), x0) * Single(1.7014118346046923e+38)
    else
      Result := Single(1.7014118346046923e+38) * Single(1.7014118346046923e+38);
    Exit;
  end;
  { Underflow check }
  if z_pf < -2400.0 then begin
    if isodd_pf(y0) then
      Result := pcr_copysignf(Single(1.1754943508222875e-38), x0) * Single(1.1754943508222875e-38)
    else
      Result := Single(1.1754943508222875e-38) * Single(1.1754943508222875e-38);
    Exit;
  end;
  { Near 1: return 1 + z }
  if pcr_fabs(z_pf) < 1.4901161193847656e-08 then begin  { 0x1p-26 }
    Result := 1.0 + z_pf; Exit;
  end;
  ia_pf := Floor(z_pf);
  h_pf := pcr_fma(l_pf, y_pf, zt_pf - ia_pf);
  il_pf := Int64(Trunc(ia_pf));
  jl_pf := il_pf and $F;
  el_pf := sar_i64(il_pf - jl_pf, 4);
  s_pf := tb_pf[jl_pf];
  su_pf.u := UInt64(Int64($3FF) + el_pf) shl 52;
  s_pf := s_pf * su_pf.f;
  h2_pf := h_pf * h_pf;
  c0_pf := ce_pf[0] + h_pf * ce_pf[1];
  c2_pf := ce_pf[2] + h_pf * ce_pf[3];
  c4_pf := ce_pf[4] + h_pf * ce_pf[5];
  c0_pf := c0_pf + h2_pf * (c2_pf + h2_pf * c4_pf);
  w_pf := s_pf * h_pf;
  rr_pf.f := s_pf + w_pf * c0_pf;
  off_pf := 468;
  if ((rr_pf.u + off_pf) and $FFFFFFF) <= 2 * off_pf then begin
    Result := pcr_powf_accurate2(x0, y0, is_exact_pf(x0, y0));
    Exit;
  end;
  { Sign correction for odd integer y }
  et_pf := Integer((ty_pf.u shr 52) and $7FF) - $3FF;
  if et_pf >= -11 then kk_pf := ty_pf.u shl (11 + et_pf)
  else kk_pf := ty_pf.u shr (-11 - et_pf);
  if ((kk_pf shl 1) = 0) and (kk_pf <> 0) then
    rr_pf.f := pcr_copysign(rr_pf.f, x_pf);
  res_pf := Single(rr_pf.f);
  Result := res_pf;
end;



// ── compoundf: (1+x)^y correctly rounded ─────────────────────────────────────
// Constants: renamed to avoid case-insensitive clash with helper function names

const
  CF_INVLOG2: Double = 1.4426950408889634;
  // P1 polynomial coefficients (degree 1..7, index 0 unused)
  CF_P1C: array[0..7] of Double = (
    0.0,
    1.4426950408889634,
    -0.721347520444768,
    0.4808983469640691,
    -0.36067375082452474,
    0.28853899226737745,
    -0.24052620964966426,
    0.2061866781489112
  );
  // P2 polynomial (18 entries, pairs for degrees 1-5, singles for 6-13)
  CF_P2C: array[0..17] of Double = (
    1.4426950408889634,     // [0] deg1 hi
    2.0355273740931317e-17, // [1] deg1 lo
   -0.7213475204444817,    // [2] deg2 hi
   -1.0177636800051583e-17,// [3] deg2 lo
    0.4808983469629878,    // [4] deg3 hi
    2.5288808125554186e-17, // [5] deg3 lo
   -0.36067376022224085,   // [6] deg4 hi
   -5.096856780843964e-18, // [7] deg4 lo
    0.28853900817779266,   // [8] deg5 hi
    2.6289836446187457e-17, // [9] deg5 lo
   -0.24044917348149364,   // [10] deg6
    0.20609929155556583,   // [11] deg7
   -0.180336880114796,     // [12] deg8
    0.16029944899207862,   // [13] deg9
   -0.14426947904986057,   // [14] deg10
    0.13115406763151405,   // [15] deg11
   -0.12030649210107214,   // [16] deg12
    0.11105797113583867    // [17] deg13
  );
  // Q1 polynomial (degree 0..4)
  CF_Q1C: array[0..4] of Double = (
    1.0,
    0.6931471805351095,
    0.24022650695393627,
    0.055504515574106836,
    0.009618187397453178
  );
  // Q2 polynomial (12 entries)
  CF_Q2C: array[0..11] of Double = (
    1.0,                      // [0] deg0
    0.6931471805599453,       // [1] deg1 hi
    2.3190455425771328e-17,   // [2] deg1 lo
    0.24022650695910072,      // [3] deg2 hi
   -9.493917425934395e-18,    // [4] deg2 lo
    0.05550410866482158,      // [5] deg3 hi
   -2.4715450854778426e-18,   // [6] deg3 lo
    0.009618129107628477,     // [7] deg4
    0.0013333558146326069,    // [8] deg5
    0.00015403530393530196,   // [9] deg6
    1.5252789714172188e-05,   // [10] deg7
    1.321548082021043e-06     // [11] deg8
  );
  // inv[i] approximates 1/t for the i-th interval (46 entries)
  CF_INV: array[0..45] of Double = (
    1.40625, 1.375, 1.34375, 1.3125, 1.296875, 1.265625, 1.25,
    1.21875, 1.203125, 1.171875, 1.15625, 1.125, 1.109375, 1.09375,
    1.078125, 1.0625, 1.046875, 1.03125, 1.0, 1.0, 0.9765625,
    0.9609375, 0.9453125, 0.9375, 0.921875, 0.90625, 0.89453125,
    0.8828125, 0.87109375, 0.859375, 0.84765625, 0.8359375, 0.828125,
    0.81640625, 0.8046875, 0.796875, 0.78515625, 0.7734375, 0.765625,
    0.7578125, 0.75, 0.7421875, 0.73046875, 0.72265625, 0.71484375,
    0.70703125
  );
  // log2inv[i][0..1]: double-double approx of -log2(inv[i])
  CF_LOG2INV: array[0..45, 0..1] of Double = (
    (-0.4918530963296747,      1.0820682119194486e-17),
    (-0.45943161863729726,     3.8053583859449705e-19),
    (-0.42626475470209796,     1.9932012137193316e-17),
    (-0.3923174227787603,      1.6328502208352762e-17),
    (-0.37503943134692475,    -1.099000777384843e-17),
    (-0.33985000288462475,     2.0897960245560436e-17),
    (-0.32192809488736235,     3.717019964142682e-19),
    (-0.28540221886224837,     2.726283638197372e-17),
    (-0.2667865406949014,      1.148454798555715e-17),
    (-0.22881869049588088,     5.967894054218645e-18),
    (-0.20945336562894978,     1.747801539116594e-18),
    (-0.16992500144231237,     1.0448980122780218e-17),
    (-0.14974711950468206,    -3.3957331682262494e-18),
    (-0.12928301694496647,     1.147571414337692e-17),
    (-0.10852445677816905,    -5.4046572138033075e-18),
    (-0.0874628412503394,     -6.765321226991275e-18),
    (-0.06608919045777244,     4.130247852756734e-18),
    (-0.044394119358453436,   -1.3338680039226223e-18),
    (0.0, 0.0),
    (0.0, 0.0),
    (0.034215715337912955,     1.1151059892428047e-18),
    (0.057485494660760125,     1.1745696149950948e-19),
    (0.08113676272540549,      7.610716771889941e-19),
    (0.09310940439148147,      5.596192057804377e-18),
    (0.11735695063815874,      5.45905252946375e-18),
    (0.14201900487242788,     -4.898294009682521e-18),
    (0.16079621190305607,     -7.518564749957147e-18),
    (0.1798210375848123,      -7.144809625324702e-18),
    (0.19910010007969528,     -1.3263604826229526e-17),
    (0.2186402864753404,       7.522378350087652e-19),
    (0.2384487675555207,      -7.16757233710703e-18),
    (0.25853301359885306,     -4.3007535189465375e-18),
    (0.2720795454368008,       2.476475356878588e-17),
    (0.29264086791911725,     -3.485682515565736e-18),
    (0.31349947281678164,     -2.4630201066282264e-17),
    (0.3275746580285044,       2.6214744450027748e-17),
    (0.34894830882107136,      2.32325257219613e-17),
    (0.37064337992039037,      1.0829515961374715e-17),
    (0.38529015588479176,      2.2208024293925304e-17),
    (0.4000871578128723,       2.4103897311490816e-17),
    (0.4150374992788438,       5.224490061390109e-18),
    (0.43014439166905216,     -3.494516357745965e-18),
    (0.4531055401123633,       2.1370790227232135e-17),
    (0.46861853948368787,      2.119503535530862e-18),
    (0.4843001617159575,       1.9794762178834054e-17),
    (0.5001541129167947,      -4.028869598543938e-17)
  );
  // exp2_T[i] = (i-16)/32
  CF_EXP2T: array[0..32] of Double = (
    -0.5, -0.46875, -0.4375, -0.40625, -0.375, -0.34375, -0.3125,
    -0.28125, -0.25, -0.21875, -0.1875, -0.15625, -0.125, -0.09375,
    -0.0625, -0.03125, 0.0, 0.03125, 0.0625, 0.09375, 0.125, 0.15625,
    0.1875, 0.21875, 0.25, 0.28125, 0.3125, 0.34375, 0.375, 0.40625,
    0.4375, 0.46875, 0.5
  );
  // exp2_U[i][0..1]: double-double approx of 2^exp2_T[i]
  CF_EXP2U: array[0..32, 0..1] of Double = (
    (0.7071067811865476,   -4.833646656726457e-17),
    (0.7225904034885233,   -1.5118790674969937e-17),
    (0.7384130729697497,   -1.741997278446398e-17),
    (0.7545822137967114,   -5.082276638771475e-17),
    (0.7711054127039704,    3.9749174048488104e-17),
    (0.7879904225539432,   -5.068458235639152e-18),
    (0.8052451659746271,    1.2353596284898944e-17),
    (0.8228777390769825,   -5.062839956837386e-17),
    (0.8408964152537145,    4.099505010290748e-17),
    (0.859309649061239,    -9.256902091315555e-18),
    (0.8781260801866497,    1.4800703477244367e-17),
    (0.8973545375015536,    9.113729213956043e-18),
    (0.9170040432046712,    1.6415536121228136e-17),
    (0.93708381705515,     -3.061381706502071e-17),
    (0.9576032806985737,   -5.3099730280979813e-17),
    (0.9785720620877001,    4.480383895518334e-17),
    (1.0, 0.0),
    (1.0218971486541166,    5.109225028973444e-17),
    (1.0442737824274138,    8.551889705537965e-17),
    (1.0671404006768237,   -7.899853966841582e-17),
    (1.0905077326652577,   -3.046782079812471e-17),
    (1.1143867425958924,    1.0410278456845571e-16),
    (1.1387886347566916,    8.912812676025408e-17),
    (1.1637248587775775,    3.8292048369240935e-17),
    (1.189207115002721,     3.982015231465646e-17),
    (1.215247359980469,    -7.712630692681488e-17),
    (1.241857812073484,     4.658027591836937e-17),
    (1.2690509571917332,    2.667932131342186e-18),
    (1.2968395546510096,    2.5382502794888315e-17),
    (1.3252366431597413,   -2.8587312100388614e-17),
    (1.3542555469368927,    7.70094837980299e-17),
    (1.383909881963832,    -6.770511658794786e-17),
    (1.4142135623730951,   -9.667293313452913e-17)
  );
  // scale[e+29] = 2^(-e) for e = -29..128 (158 entries)
  CF_SCALE: array[0..157] of Double = (
    536870912.0, 268435456.0, 134217728.0, 67108864.0, 33554432.0,
    16777216.0, 8388608.0, 4194304.0, 2097152.0, 1048576.0, 524288.0,
    262144.0, 131072.0, 65536.0, 32768.0, 16384.0, 8192.0, 4096.0,
    2048.0, 1024.0, 512.0, 256.0, 128.0, 64.0, 32.0, 16.0, 8.0, 4.0,
    2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125,
    0.00390625, 0.001953125, 0.0009765625, 0.00048828125, 0.000244140625,
    0.0001220703125, 6.103515625e-05, 3.0517578125e-05, 1.52587890625e-05,
    7.62939453125e-06, 3.814697265625e-06, 1.9073486328125e-06,
    9.5367431640625e-07, 4.76837158203125e-07, 2.384185791015625e-07,
    1.1920928955078125e-07, 5.960464477539063e-08, 2.9802322387695312e-08,
    1.4901161193847656e-08, 7.450580596923828e-09, 3.725290298461914e-09,
    1.862645149230957e-09, 9.313225746154785e-10, 4.656612873077393e-10,
    2.3283064365386963e-10, 1.1641532182693481e-10, 5.820766091346741e-11,
    2.9103830456733704e-11, 1.4551915228366852e-11, 7.275957614183426e-12,
    3.637978807091713e-12, 1.8189894035458565e-12, 9.094947017729282e-13,
    4.547473508864641e-13, 2.2737367544323206e-13, 1.1368683772161603e-13,
    5.684341886080802e-14, 2.842170943040401e-14, 1.4210854715202004e-14,
    7.105427357601002e-15, 3.552713678800501e-15, 1.7763568394002505e-15,
    8.881784197001252e-16, 4.440892098500626e-16, 2.220446049250313e-16,
    1.1102230246251565e-16, 5.551115123125783e-17, 2.7755575615628914e-17,
    1.3877787807814457e-17, 6.938893903907228e-18, 3.469446951953614e-18,
    1.734723475976807e-18, 8.673617379884035e-19, 4.336808689942018e-19,
    2.168404344971009e-19, 1.0842021724855044e-19, 5.421010862427522e-20,
    2.710505431213761e-20, 1.3552527156068805e-20, 6.776263578034403e-21,
    3.3881317890172014e-21, 1.6940658945086007e-21, 8.470329472543003e-22,
    4.235164736271502e-22, 2.117582368135751e-22, 1.0587911840678754e-22,
    5.293955920339377e-23, 2.6469779601696886e-23, 1.3234889800848443e-23,
    6.617444900424222e-24, 3.308722450212111e-24, 1.6543612251060553e-24,
    8.271806125530277e-25, 4.1359030627651384e-25, 2.0679515313825692e-25,
    1.0339757656912846e-25, 5.169878828456423e-26, 2.5849394142282115e-26,
    1.2924697071141057e-26, 6.462348535570529e-27, 3.2311742677852644e-27,
    1.6155871338926322e-27, 8.077935669463161e-28, 4.0389678347315804e-28,
    2.0194839173657902e-28, 1.0097419586828951e-28, 5.048709793414476e-29,
    2.524354896707238e-29, 1.262177448353619e-29, 6.310887241768095e-30,
    3.1554436208840472e-30, 1.5777218104420236e-30, 7.888609052210118e-31,
    3.944304526105059e-31, 1.9721522630525295e-31, 9.860761315262648e-32,
    4.930380657631324e-32, 2.465190328815662e-32, 1.232595164407831e-32,
    6.162975822039155e-33, 3.0814879110195774e-33, 1.5407439555097887e-33,
    7.703719777548943e-34, 3.851859888774472e-34, 1.925929944387236e-34,
    9.62964972193618e-35, 4.81482486096809e-35, 2.407412430484045e-35,
    1.2037062152420224e-35, 6.018531076210112e-36, 3.009265538105056e-36,
    1.504632769052528e-36, 7.52316384526264e-37, 3.76158192263132e-37,
    1.88079096131566e-37, 9.4039548065783e-38, 4.70197740328915e-38,
    2.350988701644575e-38, 1.1754943508222875e-38, 5.877471754111438e-39,
    2.938735877055719e-39
  );
  // xmax[y] for 1<=y<=15: largest odd m such that m^y fits in 25 bits
  CF_XMAX: array[0..15] of UInt64 = (
    0, 16777215, 5791, 321, 75, 31, 17, 11, 7, 5, 5, 3, 3, 3, 3, 3
  );

// ── BSF/BSR 64-bit helpers ────────────────────────────────────────────────────
function cf_bsf64(x: UInt64): Integer;
{$IFDEF CPUX86_64}
var r: QWord;
begin
  asm
    bsf rax, x
    mov r, rax
  end;
  Result := Integer(r);
end;
{$ELSE}
var i: Integer;
begin
  Result := 0;
  for i := 0 to 63 do
    if (x and (UInt64(1) shl i)) <> 0 then begin Result := i; Exit; end;
end;
{$ENDIF}

function cf_bsr64(x: UInt64): Integer;
{$IFDEF CPUX86_64}
var r: QWord;
begin
  asm
    bsr rax, x
    mov r, rax
  end;
  Result := Integer(r);
end;
{$ELSE}
var i: Integer;
begin
  Result := 0;
  for i := 63 downto 0 do
    if (x and (UInt64(1) shl i)) <> 0 then begin Result := i; Exit; end;
end;
{$ENDIF}

// ── MXCSR flag save/restore ───────────────────────────────────────────────────
function cf_get_flag: DWord;
{$IFDEF CPUX86_64}
var r: DWord;
begin
  asm
    stmxcsr r
  end;
  Result := r;
end;
{$ELSE}
begin
  Result := 0;
end;
{$ENDIF}

procedure cf_set_flag(flag: DWord);
{$IFDEF CPUX86_64}
begin
  asm
    ldmxcsr flag
  end;
end;
{$ELSE}
begin
end;
{$ENDIF}

// ── is_signalingf ─────────────────────────────────────────────────────────────
function cf_is_signalingf(x: Single): Boolean; inline;
var u_sig: Tb32u32;
begin
  u_sig.f := x;
  u_sig.u := u_sig.u xor $00400000;
  Result := (u_sig.u and $7FFFFFFF) > $7FC00000;
end;

// ── isint: returns non-zero if y is an integer ────────────────────────────────
function cf_isint(y: Single): Integer; inline;
var wy_ii: Tb32u32;
    ey_ii, s_ii: Integer;
begin
  wy_ii.f := y;
  ey_ii := Integer((wy_ii.u shr 23) and $FF) - 127;
  s_ii := ey_ii + 9;
  if ey_ii >= 0 then begin
    if s_ii >= 32 then begin Result := 1; Exit; end;
    if (wy_ii.u shl s_ii) = 0 then Result := 1 else Result := 0;
    Exit;
  end;
  if (wy_ii.u shl 1) = 0 then begin Result := 1; Exit; end;
  Result := 0;
end;

// ── fast_two_sum ──────────────────────────────────────────────────────────────
procedure cf_fast_two_sum(var s, t: Double; a, b: Double); inline;
var e_fts: Double;
begin
  s := a + b;
  e_fts := s - a;
  t := b - e_fts;
end;

// ── a_mul: hi+lo = a*b exactly ───────────────────────────────────────────────
procedure cf_a_mul(var hi, lo: Double; a, b: Double); inline;
begin
  hi := a * b;
  lo := pcr_fma(a, b, -hi);
end;

// ── s_mul: (hi+lo) = a*(bh+bl) ────────────────────────────────────────────────
procedure cf_s_mul(var hi, lo: Double; a, bh, bl: Double); inline;
begin
  cf_a_mul(hi, lo, a, bh);
  lo := pcr_fma(a, bl, lo);
end;

// ── d_mul: (hi+lo) = (ah+al)*(bh+bl) ─────────────────────────────────────────
procedure cf_d_mul(var hi, lo: Double; ah, al, bh, bl: Double); inline;
var s_dm, t_dm: Double;
begin
  cf_a_mul(hi, s_dm, ah, bh);
  t_dm := pcr_fma(al, bh, s_dm);
  lo := pcr_fma(ah, bl, t_dm);
end;

// ── p1: approximates log2(1+z) for |z|<=1/64 ─────────────────────────────────
function cf_p1(z: Double): Double;
var z2_p1, z4_p1, c5_p1, c3_p1, c1_p1: Double;
begin
  z2_p1 := z * z;
  c5_p1 := pcr_fma(CF_P1C[6], z, CF_P1C[5]);
  c3_p1 := pcr_fma(CF_P1C[4], z, CF_P1C[3]);
  c1_p1 := pcr_fma(CF_P1C[2], z, CF_P1C[1]);
  z4_p1 := z2_p1 * z2_p1;
  c5_p1 := pcr_fma(CF_P1C[7], z2_p1, c5_p1);
  c1_p1 := pcr_fma(c3_p1, z2_p1, c1_p1);
  c1_p1 := pcr_fma(c5_p1, z4_p1, c1_p1);
  Result := z * c1_p1;
end;

// ── p2: double-double approx of log2(1+zh+zl) ────────────────────────────────
procedure cf_p2(var h, l: Double; zh, zl: Double);
var t_p2: Double;
begin
  h := CF_P2C[17]; // degree 13
  h := pcr_fma(h, zh, CF_P2C[16]);
  h := pcr_fma(h, zh, CF_P2C[15]);
  h := pcr_fma(h, zh, CF_P2C[14]);
  h := pcr_fma(h, zh, CF_P2C[13]);
  h := pcr_fma(h, zh, CF_P2C[12]);
  h := pcr_fma(h, zh, CF_P2C[11]);
  cf_s_mul(h, l, h, zh, zl);
  cf_fast_two_sum(h, t_p2, CF_P2C[10], h);
  l := l + t_p2;
  cf_d_mul(h, l, h, l, zh, zl);
  cf_fast_two_sum(h, t_p2, CF_P2C[8], h);
  l := l + t_p2 + CF_P2C[9];
  cf_d_mul(h, l, h, l, zh, zl);
  cf_fast_two_sum(h, t_p2, CF_P2C[6], h);
  l := l + t_p2 + CF_P2C[7];
  cf_d_mul(h, l, h, l, zh, zl);
  cf_fast_two_sum(h, t_p2, CF_P2C[4], h);
  l := l + t_p2 + CF_P2C[5];
  cf_d_mul(h, l, h, l, zh, zl);
  cf_fast_two_sum(h, t_p2, CF_P2C[2], h);
  l := l + t_p2 + CF_P2C[3];
  cf_d_mul(h, l, h, l, zh, zl);
  cf_fast_two_sum(h, t_p2, CF_P2C[0], h);
  l := l + t_p2 + CF_P2C[1];
  cf_d_mul(h, l, h, l, zh, zl);
end;

// ── q1: approximates 2^z for |z|<=2^-6 ──────────────────────────────────────
function cf_q1(z: Double): Double;
var z2_q1, c3_q1, c0_q1, c2_q1: Double;
begin
  z2_q1 := z * z;
  c3_q1 := pcr_fma(CF_Q1C[4], z, CF_Q1C[3]);
  c0_q1 := pcr_fma(CF_Q1C[1], z, CF_Q1C[0]);
  c2_q1 := pcr_fma(c3_q1, z, CF_Q1C[2]);
  Result := pcr_fma(c2_q1, z2_q1, c0_q1);
end;

// ── q2: double-double approx of 2^(h+l) ──────────────────────────────────────
procedure cf_q2(var qh, ql: Double; h, l: Double);
var h2_q2, c7_q2, c5_q2, t_q2: Double;
begin
  h2_q2 := h * h;
  c7_q2 := pcr_fma(CF_Q2C[11], h, CF_Q2C[10]);
  c5_q2 := pcr_fma(CF_Q2C[9],  h, CF_Q2C[8]);
  c5_q2 := pcr_fma(c7_q2, h2_q2, c5_q2);
  qh := c5_q2 * h;
  cf_fast_two_sum(qh, ql, CF_Q2C[7], qh);
  cf_d_mul(qh, ql, qh, ql, h, l);
  cf_fast_two_sum(qh, t_q2, CF_Q2C[5], qh);
  ql := ql + t_q2 + CF_Q2C[6];
  cf_d_mul(qh, ql, qh, ql, h, l);
  cf_fast_two_sum(qh, t_q2, CF_Q2C[3], qh);
  ql := ql + t_q2 + CF_Q2C[4];
  cf_d_mul(qh, ql, qh, ql, h, l);
  cf_fast_two_sum(qh, t_q2, CF_Q2C[1], qh);
  ql := ql + t_q2 + CF_Q2C[2];
  cf_d_mul(qh, ql, qh, ql, h, l);
  cf_fast_two_sum(qh, t_q2, CF_Q2C[0], qh);
  ql := ql + t_q2;
end;

// ── _log2p1 (fast approximation of log2(1+x)) ────────────────────────────────
function cf_log2p1_fast(x: Double): Double;
var u_lp: Double;
    v_lp: Tb64u64;
    m_lp: UInt64;
    e_lp: Int64;
    i_lp: Integer;
    t_lp, z_lp, p_lp: Double;
begin
  u_lp := 1.0 + x;
  v_lp.f := u_lp;
  m_lp := v_lp.u and $FFFFFFFFFFFFF;
  e_lp := Int64(v_lp.u shr 52) - $3FF;
  if m_lp >= $6A09E667F3BCD then Inc(e_lp);
  v_lp.u := UInt64(Int64(v_lp.u) - e_lp * Int64($10000000000000));
  t_lp := v_lp.f;
  v_lp.f := v_lp.f + 2.0;
  i_lp := Integer(v_lp.u shr 45) - $2002D;
  z_lp := pcr_fma(CF_INV[i_lp], t_lp, -1.0);
  p_lp := cf_p1(z_lp);
  Result := Double(e_lp) + (CF_LOG2INV[i_lp][0] + p_lp);
end;

// ── exp2_1: fast path 2^t → Single (returns -1 if rounding test fails) ────────
function cf_exp2_1(t: Double): Single;
var k_e1: Double;
    r_e1: Double;
    v_e1, err_e1: Tb64u64;
    i_e1: Integer;
    lb_e1, rb_e1: Single;
begin
  k_e1 := pcr_roundeven(t);
  r_e1 := t - k_e1;
  v_e1.f := 3.015625 + r_e1;
  i_e1 := Integer(v_e1.u shr 46) - $10010;
  r_e1 := r_e1 - CF_EXP2T[i_e1];
  v_e1.f := CF_EXP2U[i_e1][0] * cf_q1(r_e1);
  err_e1.f := 5.062616992290714e-13; // 0x1.1dp-41
  v_e1.u := UInt64(Int64(v_e1.u) + Int64(k_e1) * Int64($10000000000000));
  if v_e1.f < 1.175494350822881e-38 then begin // 0x1.00000000008e2p-126
    Result := -1.0;
    Exit;
  end;
  err_e1.u := UInt64(Int64(err_e1.u) + Int64(k_e1) * Int64($10000000000000));
  lb_e1 := Single(v_e1.f) - Single(err_e1.f);
  rb_e1 := Single(v_e1.f) + Single(err_e1.f);
  if lb_e1 <> rb_e1 then begin
    Result := -1.0;
    Exit;
  end;
  Result := lb_e1;
end;

// ── is_exact_or_midpoint ──────────────────────────────────────────────────────
function cf_is_exact_or_midpoint(x, y: Single; var midpoint: Integer): Integer;
var v_iem, w_iem: Tb32u32;
    vd_iem: Tb64u64;
    e_iem: Int32;
    m_iem: UInt64;
    t_iem: Integer;
    my_iem: UInt64;
    y_int_iem: Integer;
    ez_iem: Int32;
    n_iem: UInt32;
    f_iem: Int32;
    t2_iem: Integer;
    n0_iem: UInt32;
    iters_iem: Integer;
    dm_iem: Double;
    s_iem: Double;
begin
  v_iem.f := x;
  w_iem.f := y;
  if ((v_iem.u shl 1) <> 0) and ((w_iem.u shl (32 - 16)) <> 0) then begin
    Result := 0; Exit;
  end;
  if (v_iem.u shl 1) = 0 then begin // x = 0
    Result := 1; Exit;
  end;
  e_iem := Int32(Int32((v_iem.u shl 1) shr 24) - $96);
  if (e_iem < -76) or (30 < e_iem) then begin Result := 0; Exit; end;
  vd_iem.f := 1.0 + Double(x);
  e_iem := Int32(Int32((vd_iem.u shl 1) shr 53) - $433);
  if (y >= 0.0) and (cf_isint(y) <> 0) then begin
    m_iem := vd_iem.u and $FFFFFFFFFFFFF;
    if e_iem >= -1074 then
      m_iem := m_iem or $10000000000000
    else
      Inc(e_iem);
    t_iem := cf_bsf64(m_iem);
    m_iem := m_iem shr t_iem;
    Inc(e_iem, t_iem);
    if (y = 0.0) or (y = 1.0) then begin
      if m_iem > $1000000 then midpoint := 1;
      Result := 1; Exit;
    end;
    if m_iem = 1 then begin
      if (-149 <= Int64(Trunc(y)) * Int64(e_iem))
      and (Int64(Trunc(y)) * Int64(e_iem) < 128) then
        Result := 1
      else
        Result := 0;
      Exit;
    end;
    if (y < 0.0) or (y > 15.0) then begin Result := 0; Exit; end;
    y_int_iem := Integer(Trunc(y));
    if m_iem > CF_XMAX[y_int_iem] then begin Result := 0; Exit; end;
    my_iem := m_iem * m_iem;
    t2_iem := 2;
    while t2_iem < y_int_iem do begin
      my_iem := my_iem * m_iem;
      Inc(t2_iem);
    end;
    t_iem := 1 + cf_bsr64(my_iem);
    ez_iem := e_iem * y_int_iem + t_iem;
    if (ez_iem <= -149) or (128 < ez_iem) then begin Result := 0; Exit; end;
    if my_iem > $1000000 then midpoint := 1;
    if e_iem * y_int_iem >= -149 then Result := 1 else Result := 0;
    Exit;
  end;
  // second branch: y is not a non-negative integer
  n_iem := w_iem.u and $7FFFFF;
  f_iem := Int32(Int32((w_iem.u shl 1) shr 24) - $96);
  if f_iem >= -149 then
    n_iem := n_iem or $800000
  else
    Inc(f_iem);
  t_iem := pcr_bsf32(n_iem);
  n_iem := n_iem shr t_iem;
  Inc(f_iem, t_iem);
  m_iem := vd_iem.u and $FFFFFFFFFFFFF;
  if e_iem >= -1074 then
    m_iem := m_iem or $10000000000000
  else
    Inc(e_iem);
  t_iem := cf_bsf64(m_iem);
  m_iem := m_iem shr t_iem;
  Inc(e_iem, t_iem);
  if y < 0.0 then begin
    if m_iem <> 1 then begin Result := 0; Exit; end;
    t_iem := pcr_bsf32(UInt32(e_iem));
    if -f_iem > t_iem then begin Result := 0; Exit; end;
    if e_iem >= 0 then begin
      if f_iem >= 0 then
        ez_iem := -(e_iem shl f_iem) * Int32(n_iem)
      else
        ez_iem := -sar_i32(e_iem, -f_iem) * Int32(n_iem);
    end else begin
      if f_iem >= 0 then
        ez_iem := ((-e_iem) shl f_iem) * Int32(n_iem)
      else
        ez_iem := sar_i32(-e_iem, -f_iem) * Int32(n_iem);
    end;
    if (-149 <= ez_iem) and (ez_iem < 128) then Result := 1 else Result := 0;
    Exit;
  end;
  // y > 0, non-integer: extract squares from m*2^e
  iters_iem := -f_iem;
  while iters_iem > 0 do begin
    Dec(iters_iem);
    if (e_iem and 1) <> 0 then begin Result := 0; Exit; end;
    e_iem := e_iem div 2;
    dm_iem := Double(m_iem);
    s_iem := Double(Round(pcr_sqrt(dm_iem)));
    if s_iem * s_iem <> dm_iem then begin Result := 0; Exit; end;
    m_iem := UInt64(Trunc(s_iem));
  end;
  if m_iem > 1 then begin
    if n_iem > 15 then begin Result := 0; Exit; end;
    if m_iem > CF_XMAX[n_iem] then begin Result := 0; Exit; end;
  end;
  my_iem := m_iem;
  n0_iem := n_iem;
  while n0_iem > 1 do begin
    my_iem := my_iem * m_iem;
    Dec(n0_iem);
  end;
  t_iem := 1 + pcr_bsr32(UInt32(my_iem));
  if (-149 <= e_iem * Int32(n_iem)) and (e_iem * Int32(n_iem) + t_iem <= 128) then
    Result := 1
  else
    Result := 0;
end;

// ── exp2_2: accurate path for 2^(h+l) → Single ───────────────────────────────
function cf_exp2_2(h, l: Double; x, y: Single; exact: Integer; flag: DWord): Single;
const
  CF_ERR_E22: array[0..1] of Double = (
    8.744365362193872e-26,  // 0x1.b1p-84
    4.861355328424485e-29   // 0x1.edp-95
  );
var k_e22: Double;
    r_e22: Double;
    v_e22, w_e22: Tb64u64;
    i_e22: Integer;
    qh_e22, ql_e22: Double;
    small_e22: Boolean;
    err_e22: Double;
    left_e22, right_e22: Single;
    vtz_e22, wtz_e22: Integer;
begin
  if y = 1.0 then begin
    cf_set_flag(flag);
    Result := 1.0 + x;
    Exit;
  end;
  k_e22 := pcr_roundeven(h);
  // if h+l is tiny, 2^(h+l) rounds to 1
  if (k_e22 = 0.0) and (pcr_fabs(h) <= 4.299566335638736e-08) then begin
    Result := Single(1.0 + h * 0.5);
    Exit;
  end;
  r_e22 := h - k_e22;
  cf_fast_two_sum(h, l, r_e22, l);
  v_e22.f := 3.015625 + h;
  i_e22 := Integer(v_e22.u shr 46) - $10010;
  h := h - CF_EXP2T[i_e22];
  cf_fast_two_sum(h, l, h, l);
  cf_q2(qh_e22, ql_e22, h, l);
  cf_d_mul(qh_e22, ql_e22, CF_EXP2U[i_e22][0], CF_EXP2U[i_e22][1], qh_e22, ql_e22);
  cf_fast_two_sum(qh_e22, ql_e22, qh_e22, ql_e22);
  // rounding test
  w_e22.f := qh_e22;
  if ((w_e22.u + 1) and $FFFFFFF) <= 2 then begin
    small_e22 := (k_e22 = 0.0) and (i_e22 = 16) and (pcr_fabs(h) <= 3.814697265625e-06);
    if small_e22 then err_e22 := CF_ERR_E22[1] else err_e22 := CF_ERR_E22[0];
    v_e22.f := qh_e22 + (ql_e22 - err_e22);
    v_e22.u := UInt64(Int64(v_e22.u) + Int64(k_e22) * Int64($10000000000000));
    w_e22.f := qh_e22 + (ql_e22 + err_e22);
    w_e22.u := UInt64(Int64(w_e22.u) + Int64(k_e22) * Int64($10000000000000));
    if exact <> 0 then begin
      vtz_e22 := cf_bsf64(v_e22.u);
      wtz_e22 := cf_bsf64(w_e22.u);
      cf_set_flag(flag);
      if vtz_e22 >= wtz_e22 then Result := Single(v_e22.f)
      else Result := Single(w_e22.f);
      Exit;
    end;
    left_e22  := Single(v_e22.f);
    right_e22 := Single(w_e22.f);
    if left_e22 <> right_e22 then
      Halt(1);
  end;
  // multiply qh+ql by 2^k
  v_e22.f := qh_e22 + ql_e22;
  w_e22.f := qh_e22;
  if ((w_e22.u shl 36) = 0) and (v_e22.f = qh_e22) and (ql_e22 <> 0.0) then begin
    if ql_e22 > 0.0 then Inc(v_e22.u) else Dec(v_e22.u);
  end;
  v_e22.u := UInt64(Int64(v_e22.u) + Int64(k_e22) * Int64($10000000000000));
  Result := Single(v_e22.f);
end;

// ── log2p1_accurate: double-double approx of log2(1+x) ───────────────────────
procedure cf_log2p1_accurate(var h, l: Double; x: Double);
var v_la: Tb64u64;
    m_la: UInt64;
    e_la: Int64;
    i_la: Integer;
    r_la, zh_la, zl_la, ph_la, pl_la, t_la: Double;
begin
  if 1.0 >= x then begin
    if pcr_fabs(x) >= 1.1102230246251565e-16 then // 2^-53
      cf_fast_two_sum(h, l, 1.0, x)
    else begin
      h := 1.0;
      l := x;
    end;
  end else
    cf_fast_two_sum(h, l, x, 1.0);
  v_la.f := h;
  m_la := v_la.u and $FFFFFFFFFFFFF;
  e_la := Int64(v_la.u shr 52) - $3FF;
  if m_la >= $6A09E667F3BCD then Inc(e_la);
  h := h * CF_SCALE[e_la + 29];
  l := l * CF_SCALE[e_la + 29];
  v_la.f := 2.0 + h;
  i_la := Integer(v_la.u shr 45) - $2002D;
  r_la := CF_INV[i_la];
  zh_la := pcr_fma(r_la, h, -1.0);
  zl_la := r_la * l;
  cf_fast_two_sum(zh_la, zl_la, zh_la, zl_la);
  cf_p2(ph_la, pl_la, zh_la, zl_la);
  cf_fast_two_sum(h, l, Double(e_la), CF_LOG2INV[i_la][0]);
  l := l + CF_LOG2INV[i_la][1];
  cf_fast_two_sum(h, t_la, h, ph_la);
  l := l + t_la + pl_la;
end;

// ── accurate_path ─────────────────────────────────────────────────────────────
function cf_accurate_path(x, y: Single; exact: Integer; flag: DWord): Single;
var h_ap, l_ap: Double;
begin
  cf_log2p1_accurate(h_ap, l_ap, Double(x));
  cf_s_mul(h_ap, l_ap, Double(y), h_ap, l_ap);
  Result := cf_exp2_2(h_ap, l_ap, x, y, exact, flag);
end;

// ── cf_special: handles special/edge cases ────────────────────────────────────
function cf_special(x, y: Single): Single;
var nx_sp, ny_sp: Tb32u32;
    ax_sp, ay_sp: LongWord;
    mone_sp: Tb32u32;
    sy_sp: Integer;
begin
  nx_sp.f := x;
  ny_sp.f := y;
  ax_sp := nx_sp.u shl 1;
  ay_sp := ny_sp.u shl 1;
  mone_sp.f := -1.0;
  if (ax_sp = 0) or (ay_sp = 0) then begin
    if ax_sp = 0 then begin
      if cf_is_signalingf(y) then Result := x + y
      else Result := 1.0;
      Exit;
    end;
    if ay_sp = 0 then begin
      if cf_is_signalingf(x) then begin Result := x + y; Exit; end;
      if x < -1.0 then Result := 0.0 / 0.0
      else Result := 1.0;
      Exit;
    end;
  end;
  if ay_sp >= (LongWord($FF) shl 24) then begin // y = Inf/NaN
    if ax_sp > (LongWord($FF) shl 24) then begin Result := x + y; Exit; end;
    if ay_sp = (LongWord($FF) shl 24) then begin // y = +/-Inf
      if nx_sp.u > mone_sp.u then begin Result := 0.0 / 0.0; Exit; end;
      sy_sp := Integer(ny_sp.u shr 31);
      if nx_sp.u = mone_sp.u then begin
        if sy_sp = 0 then Result := 0.0
        else Result := -y;
        Exit;
      end;
      if x < 0.0 then begin
        if sy_sp = 0 then Result := 0.0 else Result := -y;
        Exit;
      end;
      if x > 0.0 then begin
        if sy_sp <> 0 then Result := 0.0 else Result := y;
        Exit;
      end;
      Result := 1.0;
      Exit;
    end;
    Result := x + y; // y = NaN
    Exit;
  end;
  if nx_sp.u >= (LongWord($FF) shl 23) then begin
    if ax_sp = (LongWord($FF) shl 24) then begin // x = ±Inf
      if (nx_sp.u shr 31) <> 0 then begin Result := 0.0 / 0.0; Exit; end; // -Inf
      if (ny_sp.u shr 31) <> 0 then Result := 1.0 / x
      else Result := x;
      Exit;
    end;
    if ax_sp > (LongWord($FF) shl 24) then begin Result := x + y; Exit; end; // NaN
    if nx_sp.u > mone_sp.u then begin
      Result := 0.0 / 0.0; // x < -1
      Exit;
    end;
    // x = -1
    if (ny_sp.u shr 31) <> 0 then Result := 1.0 / 0.0
    else Result := 0.0;
    Exit;
  end;
  Result := 0.0;
end;

// ── pcr_compoundf: main entry point ──────────────────────────────────────────
function pcr_compoundf(x, n: Single): Single;
var nx_cf, ny_cf: Tb32u32;
    mone_cf: Tb32u32;
    ax_cf, ay_cf: LongWord;
    xd_cf, yd_cf: Double;
    tx_cf, ty_cf, t_cf: Tb64u64;
    flag_cf: DWord;
    l_cf: Double;
    midpoint_cf, exact_cf: Integer;
    res_cf: Single;
begin
  mone_cf.f := -1.0;
  nx_cf.f := x;
  ny_cf.f := n;
  if nx_cf.u >= mone_cf.u then begin
    Result := cf_special(x, n); Exit;
  end;
  ax_cf := nx_cf.u shl 1;
  ay_cf := ny_cf.u shl 1;
  if (ax_cf = 0) or (ax_cf >= (LongWord($FF) shl 24))
  or (ay_cf = 0) or (ay_cf >= (LongWord($FF) shl 24)) then begin
    Result := cf_special(x, n); Exit;
  end;
  xd_cf := Double(x);
  yd_cf := Double(n);
  tx_cf.f := xd_cf;
  ty_cf.f := yd_cf;
  flag_cf := cf_get_flag;
  if ax_cf < $62000000 then begin // |x| < 2^-29
    l_cf := CF_INVLOG2 * (xd_cf - (xd_cf * xd_cf) * 0.5);
  end else begin
    l_cf := cf_log2p1_fast(tx_cf.f);
  end;
  t_cf.f := l_cf * ty_cf.f;
  // detect overflow/underflow
  if (t_cf.u shl 1) >= $80C0000000000000 then begin
    if t_cf.u >= UInt64($C062C00000000000) then begin // t <= -150
      Result := Single(1.1754943508222875e-38) * Single(1.1754943508222875e-38);
      Exit;
    end else if (t_cf.u shr 63) = 0 then begin // t >= 128
      Result := Single(8.507059173023462e+37) * Single(8.507059173023462e+37);
      Exit;
    end;
  end;
  // 2^t rounds to 1
  if (t_cf.u shl 1) <= $7CCE2A8ED5E1A9B2 then begin
    if (t_cf.u shr 63) <> 0 then
      Result := 1.0 - 2.9802322387695312e-08
    else
      Result := 1.0 + 2.9802322387695312e-08;
    Exit;
  end;
  midpoint_cf := 0;
  exact_cf := cf_is_exact_or_midpoint(x, n, midpoint_cf);
  res_cf := cf_exp2_1(t_cf.f);
  if res_cf <> -1.0 then begin
    if (exact_cf <> 0) and (midpoint_cf = 0) then
      cf_set_flag(flag_cf);
    Result := res_cf;
    Exit;
  end;
  cf_set_flag(flag_cf);
  Result := cf_accurate_path(x, n, exact_cf, flag_cf);
end;


// ---------------------------------------------------------------------------
// pcr_sinf / pcr_cosf  (ported from CORE-MATH sinf.c / cosf.c)
// ---------------------------------------------------------------------------

// ---- shared polynomial coefficients (same in sinf.c and cosf.c) -----------
const
  sincos_a: array[0..3] of Double = (
     0.19634954084936204,     // 0x1.921fb54442d17p-3
    -0.0012616486279372187,   // -0x1.4abbce6256a39p-10
     2.432025854080733e-06,   // 0x1.466bc5a518c16p-19
    -2.2318367225754577e-09   // -0x1.32bdc61074ff6p-29
  );
  sincos_b: array[0..3] of Double = (
     0.019276571095877645,    // 0x1.3bd3cc9be45dcp-6
    -6.193103220211784e-05,   // -0x1.03c1f081b0833p-14
     7.958785981094399e-08,   // 0x1.55d3c6fc9ac1fp-24
    -5.4777514393633976e-11   // -0x1.e1d3ff281b40dp-35
  );
  // sinf tb[] = sin(i*pi/16), i=0..31  (one full period)
  sinf_tb: array[0..31] of Double = (
     0.0,                     // sin(0)
     0.19509032201612828,     // sin(pi/16)
     0.3826834323650898,
     0.5555702330196022,
     0.7071067811865476,
     0.8314696123025452,
     0.9238795325112867,
     0.9807852804032304,
     1.0,
     0.9807852804032304,
     0.9238795325112867,
     0.8314696123025452,
     0.7071067811865476,
     0.5555702330196022,
     0.3826834323650898,
     0.19509032201612828,
     0.0,
    -0.19509032201612828,
    -0.3826834323650898,
    -0.5555702330196022,
    -0.7071067811865476,
    -0.8314696123025452,
    -0.9238795325112867,
    -0.9807852804032304,
    -1.0,
    -0.9807852804032304,
    -0.9238795325112867,
    -0.8314696123025452,
    -0.7071067811865476,
    -0.5555702330196022,
    -0.3826834323650898,
    -0.19509032201612828
  );
  // cosf tb[] = cos(i*pi/16), i=0..31
  cosf_tb: array[0..31] of Double = (
     1.0,                     // cos(0)
     0.9807852804032304,      // cos(pi/16)
     0.9238795325112867,
     0.8314696123025452,
     0.7071067811865476,
     0.5555702330196022,
     0.3826834323650898,
     0.19509032201612828,
     0.0,
    -0.19509032201612828,
    -0.3826834323650898,
    -0.5555702330196022,
    -0.7071067811865476,
    -0.8314696123025452,
    -0.9238795325112867,
    -0.9807852804032304,
    -1.0,
    -0.9807852804032304,
    -0.9238795325112867,
    -0.8314696123025452,
    -0.7071067811865476,
    -0.5555702330196022,
    -0.3826834323650898,
    -0.19509032201612828,
     0.0,
     0.19509032201612828,
     0.3826834323650898,
     0.5555702330196022,
     0.7071067811865476,
     0.8314696123025452,
     0.9238795325112867,
     0.9807852804032304
  );
  // 2/pi in fixed-point, little-endian 64-bit limbs
  sincos_ipi: array[0..3] of UInt64 = (
    UInt64($FE5163ABDEBBC562),
    UInt64($DB6295993C439041),
    UInt64($FC2757D1F534DDC0),
    UInt64($A2F9836E4E441529)
  );

// ---- rbig: range-reduction for large arguments (shared sinf/cosf) ---------
// Maps x = 2^(e-127) * m  into  result in [-0.5, 0.5] scaled by pi/2,
// and sets *q to the octant index.
function sincos_rbig(u: LongWord; q: PInteger): Double;
var
  e_exp, k, s, i_val, sgn: Integer;
  m, p3h, p3l, p2l, p1l: UInt64;
  p0, p1, p2, p3: TUInt128;
  a_val: Int64;
  sm: Int64;
begin
  e_exp := Integer((u shr 23) and $FF);
  m := UInt64((u and $007FFFFF) or $800000);  // ~0u>>9 = $007FFFFF, 1<<23 = $800000
  p0 := MulWide(m, sincos_ipi[0]);
  p1 := MulWide(m, sincos_ipi[1]); p1 := p1 + p0.hi;
  p2 := MulWide(m, sincos_ipi[2]); p2 := p2 + p1.hi;
  p3 := MulWide(m, sincos_ipi[3]); p3 := p3 + p2.hi;
  p3h := p3.hi;
  p3l := p3.lo;
  p2l := p2.lo;
  p1l := p1.lo;
  k := e_exp - 124;
  s := k - 23;
  if s < 64 then begin
    i_val := Integer((p3h shl s) or (p3l shr (64 - s)));
    a_val := Int64((p3l shl s) or (p2l shr (64 - s)));
  end else if s = 64 then begin
    i_val := Integer(p3l);
    a_val := Int64(p2l);
  end else begin  // s > 64
    i_val := Integer((p3l shl (s - 64)) or (p2l shr (128 - s)));
    a_val := Int64((p2l shl (s - 64)) or (p1l shr (128 - s)));
  end;
  sgn := sar_i32(Integer(u), 31);   // 0 if positive, -1 if negative
  sm  := sar_i64(a_val, 63);        // sign extension of a
  i_val := i_val - Integer(sm);
  Result := Double(a_val xor Int64(sgn)) * 5.421010862427522e-20;  // *2^-64
  i_val := (i_val xor sgn) - sgn;
  q^ := i_val;
end;

// ---- rltl0: range-reduction for medium arguments (double input) -----------
function sincos_rltl0(x: Double; q: PInteger): Double;
var
  idh, id_d: Double;
  Q_r: Tb64u64;
begin
  idh  := 5.092958178940651 * x;   // 0x1.45f306dc9c883p+2
  id_d := pcr_roundeven(idh);
  Q_r.f := 6755399441055744.0 + id_d;  // 0x1.8p52
  q^ := Integer(LongWord(Q_r.u));
  Result := idh - id_d;
end;

// ---- rltl: range-reduction for medium arguments (float input) -------------
function sincos_rltl(z: Single; q: PInteger): Double;
var
  x, idl, idh, id_d: Double;
  Q_r: Tb64u64;
begin
  x    := Double(z);
  idl  := -3.1558305786379073e-09 * x;  // -0x1.b1bbead603d8bp-29
  idh  :=  5.092958182096481 * x;        // 0x1.45f306ep+2
  id_d := pcr_roundeven(idh);
  Q_r.f := 6755399441055744.0 + id_d;   // 0x1.8p52
  q^ := Integer(LongWord(Q_r.u));
  Result := (idh - id_d) + idl;
end;

// ---- sinf helpers ----------------------------------------------------------

// copysign(1.0f, x) * rh + copysign(1.0f, x) * rl
function sinf_add_sign(x, rh, rl: Single): Single; inline;
var
  t32: Tb32u32;
  sgn: Single;
begin
  t32.f := x;
  if (t32.u shr 31) = 0 then sgn := 1.0 else sgn := -1.0;
  Result := sgn * rh + sgn * rl;
end;

// sinf database lookup for hard cases
function sinf_db(x: Single; r: Double): Single;
const
  n = 4;
  // absolute-value bit patterns of the argument
  db_uarg: array[0..n-1] of LongWord = (
    $46199998,  // |0x1.33333p+13|
    $3F3ADC51,  // |0x1.75b8a2p-1|
    $3FA7832A,  // |0x1.4f0654p+0|
    $4116CBE4   // |0x1.2d97c8p+3|
  );
  // rh bit patterns
  db_urh: array[0..n-1] of LongWord = (
    $BEB1FA5D,  // -0x1.63f4bap-2
    $3F2AB445,  // 0x1.55688ap-1
    $3F7741B6,  // 0x1.ee836cp-1
    $B2CCDE2D   // -0x1.99bc5ap-26
  );
  // rl bit patterns
  db_url: array[0..n-1] of LongWord = (
    $B2000000,  // -0x1p-27
    $B2800000,  // -0x1p-26
    $B2800000,  // -0x1p-26
    $A6000000   // -0x1p-51
  );
var
  t32, trh, trl: Tb32u32;
  ax: LongWord;
  i: Integer;
begin
  t32.f := x;
  ax := t32.u and $7FFFFFFF;
  for i := 0 to n - 1 do begin
    if db_uarg[i] = ax then begin
      trh.u := db_urh[i];
      trl.u := db_url[i];
      Result := sinf_add_sign(x, trh.f, trl.f);
      Exit;
    end;
  end;
  Result := Single(r);
end;

// as_sinf_big: handles |x| > 2^26
function sinf_big(x: Single): Single;
var
  t32: Tb32u32;
  ax: LongWord;
  ia: Integer;
  z_r, z2, z4, aa, bb, s0, c0, r_val: Double;
  t_nan: Tb32u32;
begin
  t32.f := x;
  ax := t32.u shl 1;
  if ax >= $FF000000 then begin    // nan or +-inf
    if (ax shl 8) <> 0 then begin  // nan: propagate
      Result := x + x;
      Exit;
    end;
    // infinity: return NaN
    t_nan.u := $7FC00000;
    Result := t_nan.f;
    Exit;
  end;
  z_r := sincos_rbig(t32.u, @ia);
  z2  := z_r * z_r;
  z4  := z2 * z2;
  aa  := (sincos_a[0] + z2 * sincos_a[1]) + z4 * (sincos_a[2] + z2 * sincos_a[3]);
  bb  := (sincos_b[0] + z2 * sincos_b[1]) + z4 * (sincos_b[2] + z2 * sincos_b[3]);
  s0  := sinf_tb[ia and 31];
  c0  := sinf_tb[(ia + 8) and 31];
  r_val := s0 + z_r * (aa * c0 - bb * (z_r * s0));
  Result := Single(r_val);
end;

// ---- cosf helpers ----------------------------------------------------------

// cosf database lookup for hard cases
function cosf_db(x: Single; r: Double): Single;
const
  n = 5;
  db_uarg: array[0..n-1] of LongWord = (
    $4096CBE4,  // |0x1.2d97c8p+2|
    $5922AA80,  // |0x1.4555p+51|
    $5AA4542C,  // |0x1.48a858p+54|
    $5F18B878,  // |0x1.3170fp+63|
    $6115CB11   // |0x1.2b9622p+67|
  );
  db_urh: array[0..n-1] of LongWord = (
    $324CDE2E,  // 0x1.99bc5cp-27
    $3F08AEBF,  // 0x1.115d7ep-1
    $3EFA40A4,  // 0x1.f48148p-2
    $3F7F14BB,  // 0x1.fe2976p-1
    $3F78142F   // 0x1.f0285ep-1
  );
  db_url: array[0..n-1] of LongWord = (
    $A5800000,  // -0x1p-52
    $B2800000,  // -0x1p-26
    $32000000,  //  0x1p-27
    $32800000,  //  0x1p-26
    $B2800000   // -0x1p-26
  );
var
  t32, trh, trl: Tb32u32;
  ax: LongWord;
  i: Integer;
begin
  t32.f := x;
  ax := t32.u and $7FFFFFFF;
  for i := 0 to n - 1 do begin
    if db_uarg[i] = ax then begin
      trh.u := db_urh[i];
      trl.u := db_url[i];
      Result := trh.f + trl.f;
      Exit;
    end;
  end;
  Result := Single(r);
end;

// as_cosf_big: handles |x| > 2^26
function cosf_big(x: Single): Single;
var
  t32: Tb32u32;
  t64: Tb64u64;
  ax: LongWord;
  ia: Integer;
  z_r, z2, z4, aa, bb, s0, c0, r_val: Double;
  tail: UInt64;
  t_nan: Tb32u32;
begin
  t32.f := x;
  ax := t32.u shl 1;
  if ax >= $FF000000 then begin    // nan or +-inf
    if (ax shl 8) <> 0 then begin  // nan: propagate
      Result := x + x;
      Exit;
    end;
    // infinity: return NaN
    t_nan.u := $7FC00000;
    Result := t_nan.f;
    Exit;
  end;
  z_r := sincos_rbig(t32.u, @ia);
  z2  := z_r * z_r;
  z4  := z2 * z2;
  aa  := (sincos_a[0] + z2 * sincos_a[1]) + z4 * (sincos_a[2] + z2 * sincos_a[3]);
  bb  := (sincos_b[0] + z2 * sincos_b[1]) + z4 * (sincos_b[2] + z2 * sincos_b[3]);
  s0  := cosf_tb[(ia + 8) and 31];   // = -sin(ia*pi/16)
  c0  := cosf_tb[ia and 31];          // = cos(ia*pi/16)
  r_val := c0 + z_r * (aa * s0 - bb * (z_r * c0));
  // tail check for hard cases
  t64.f := r_val;
  tail := (t64.u + UInt64(6)) and ((not UInt64(0)) shr 36);  // 28-bit mask
  if tail <= 12 then begin
    Result := cosf_db(x, r_val);
    Exit;
  end;
  Result := Single(r_val);
end;

// ---- main functions --------------------------------------------------------

function pcr_sinf(x: Single): Single;
var
  t32: Tb32u32;
  ax: LongWord;
  ia: Integer;
  z0_d, z_d: Double;
  z2, z4, aa, bb, s0, c0: Double;
  r_val: Double;
begin
  t32.f := x;
  ax := t32.u shl 1;
  z0_d := Double(x);
  if (ax > $99000000) or (ax < $73000000) then begin
    // |x| > 2^26 or |x| < 2^-12
    if ax < $73000000 then begin
      // |x| < 2^-12
      if ax < $66000000 then begin
        // |x| < 2^-25
        if ax = 0 then begin
          Result := x;
          Exit;
        end;
        // fmaf(-x, |x|, x) via double: for tiny x, sin(x) rounds to x
        Result := Single(Double(-x) * Abs(Double(x)) + Double(x));
        Exit;
      end;
      // 2^-25 <= |x| < 2^-12: use cubic approximation
      // -0x1.555556p-3f = -0.1666666716337204 (f32: $BE2AAAAB)
      Result := (-0.1666666716337204 * x) * (x * x) + x;
      Exit;
    end;
    Result := sinf_big(x);
    Exit;
  end;
  if ax < $822D97C8 then begin
    // medium range: use rltl0 (double-precision reduction)
    if (ax = $7E75B8A2) or (ax = $7F4F0654) then begin
      Result := sinf_db(x, 0.0);
      Exit;
    end;
    z_d := sincos_rltl0(z0_d, @ia);
  end else begin
    // larger medium range: use rltl (float-based reduction)
    if ax = $8C333330 then begin
      Result := sinf_db(x, 0.0);
      Exit;
    end;
    z_d := sincos_rltl(x, @ia);
  end;
  z2   := z_d * z_d;
  z4   := z2 * z2;
  aa   := (sincos_a[0] + z2 * sincos_a[1]) + z4 * (sincos_a[2] + z2 * sincos_a[3]);
  bb   := (sincos_b[0] + z2 * sincos_b[1]) + z4 * (sincos_b[2] + z2 * sincos_b[3]);
  s0   := sinf_tb[ia and 31];
  c0   := sinf_tb[(ia + 8) and 31];
  r_val := s0 + aa * (z_d * c0) - bb * (z2 * s0);
  Result := Single(r_val);
end;

function pcr_cosf(x: Single): Single;
var
  t32: Tb32u32;
  ax: LongWord;
  ia: Integer;
  z0_d, z_d: Double;
  z2, z4, aa, bb, c0, s0: Double;
  r_val: Double;
begin
  t32.f := x;
  ax := t32.u shl 1;
  z0_d := Double(x);
  if (ax > $99000000) or (ax < $73000000) then begin
    // |x| > 2^26 or |x| < 2^-12
    if ax < $73000000 then begin
      // |x| < 2^-12
      if ax < $66000000 then begin
        // |x| < 2^-25
        if ax = 0 then begin
          Result := 1.0;
          Exit;
        end;
        // cos(x) = 1 - 2^-25 for tiny nonzero x (correctly rounded)
        Result := 1.0 - 2.9802322387695312e-08;   // 1 - 0x1p-25f
        Exit;
      end;
      // 2^-25 <= |x| < 2^-12: use quadratic approximation
      Result := -0.5 * x * x + 1.0;
      Exit;
    end;
    Result := cosf_big(x);
    Exit;
  end;
  if ax < $82A41896 then begin
    // medium range: use rltl0
    if ax = $812D97C8 then begin
      Result := cosf_db(x, 0.0);
      Exit;
    end;
    z_d := sincos_rltl0(z0_d, @ia);
  end else begin
    // larger medium range: use rltl
    z_d := sincos_rltl(x, @ia);
  end;
  z2   := z_d * z_d;
  z4   := z2 * z2;
  aa   := (sincos_a[0] + z2 * sincos_a[1]) + z4 * (sincos_a[2] + z2 * sincos_a[3]);
  bb   := (sincos_b[0] + z2 * sincos_b[1]) + z4 * (sincos_b[2] + z2 * sincos_b[3]);
  c0   := cosf_tb[ia and 31];
  s0   := cosf_tb[(ia + 8) and 31];   // = -sin(ia*pi/16)
  r_val := c0 + aa * (z_d * s0) - bb * (z2 * c0);
  Result := Single(r_val);
end;

// ============================================================
//  pcr_sincosf  — port of sincosf.c (cr_sincosf)
//  pcr_tanf     — port of tanf.c    (cr_tanf)
// ============================================================

// ---- tanf helpers (rbig uses k=e-127, different from sincos_rbig) ----------

const
  tanf_cn: array[0..3] of Double = (
    1.5707963267948966,      // 0x1.921fb54442d18p+0
   -0.49720165641032027,     // -0x1.fd226e573289fp-2
    0.02683402276915988,     // 0x1.b7a60c8dac9f6p-6
   -0.00017660096093977045   // -0x1.725beb40f33e5p-13
  );
  tanf_cd: array[0..3] of Double = (
    1.0,                     // 0x1p+0
   -1.1389954387488281,      // -0x1.2395347fb829dp+0
    0.1421268437745497,      // 0x1.2313660f29c36p-3
   -0.0031314039049681057    // -0x1.9a707ab98d1c1p-9
  );

function tanf_rbig(u: LongWord; q: PInteger): Double;
var
  e_exp, k, s_shift, i_val, sgn: Integer;
  m_val, p3h, p3l, p2l, p1l: UInt64;
  p0, p1, p2, p3: TUInt128;
  a_val: Int64;
  sm: Int64;
begin
  e_exp := Integer((u shr 23) and $FF);
  m_val := UInt64((u and $007FFFFF) or $800000);
  p0 := MulWide(m_val, sincos_ipi[0]);
  p1 := MulWide(m_val, sincos_ipi[1]); p1 := p1 + p0.hi;
  p2 := MulWide(m_val, sincos_ipi[2]); p2 := p2 + p1.hi;
  p3 := MulWide(m_val, sincos_ipi[3]); p3 := p3 + p2.hi;
  p3h := p3.hi;
  p3l := p3.lo;
  p2l := p2.lo;
  p1l := p1.lo;
  k       := e_exp - 127;   // tanf uses e-127; sincos_rbig uses e-124
  s_shift := k - 23;
  if s_shift < 64 then begin
    i_val := Integer((p3h shl s_shift) or (p3l shr (64 - s_shift)));
    a_val := Int64((p3l shl s_shift) or (p2l shr (64 - s_shift)));
  end else if s_shift = 64 then begin
    i_val := Integer(p3l);
    a_val := Int64(p2l);
  end else begin  // s_shift > 64
    i_val := Integer((p3l shl (s_shift - 64)) or (p2l shr (128 - s_shift)));
    a_val := Int64((p2l shl (s_shift - 64)) or (p1l shr (128 - s_shift)));
  end;
  sgn   := sar_i32(Integer(u), 31);
  sm    := sar_i64(a_val, 63);
  i_val := i_val - Integer(sm);
  Result := Double(a_val xor Int64(sgn)) * 5.421010862427522e-20;  // *2^-64
  i_val := (i_val xor sgn) - sgn;
  q^ := i_val;
end;

// tanf rltl: multiplies by 2/pi (different from sincos_rltl which uses 16/pi)
function tanf_rltl(z: Single; q: PInteger): Double;
var
  x, idl, idh, id_d: Double;
  Q_r: Tb64u64;
begin
  x    := Double(z);
  idl  := -3.944788223297384e-10 * x;   // -0x1.b1bbead603d8bp-32
  idh  :=  0.6366197727620602 * x;       // 0x1.45f306ep-1
  id_d := pcr_roundeven(idh);
  Q_r.f := 6755399441055744.0 + id_d;   // 0x1.8p52  (fast int extract)
  q^   := Integer(LongWord(Q_r.u));
  Result := (idh - id_d) + idl;
end;

// ---- sincosf helpers -------------------------------------------------------

// Database of 9 hard cases for sincosf
procedure sincosf_database(x: Single; var sout, cout: Single);
const
  sc_db_uarg: array[0..8] of LongWord = (
    $46199998,   // 0x1.33333p+13
    $3F3ADC51,   // 0x1.75b8a2p-1
    $3FA7832A,   // 0x1.4f0654p+0
    $4116CBE4,   // 0x1.2d97c8p+3
    $4096CBE4,   // 0x1.2d97c8p+2
    $5922AA80,   // 0x1.4555p+51
    $5AA4542C,   // 0x1.48a858p+54
    $5F18B878,   // 0x1.3170fp+63
    $6115CB11    // 0x1.2b9622p+67
  );
  sc_db_sh: array[0..8] of Double = (
    -0.34761324524879456,      // -0x1.63f4bap-2
     0.6668131947517395,       //  0x1.55688ap-1
     0.9658464193344116,       //  0x1.ee836cp-1
    -2.384975950064927e-08,    // -0x1.99bc5ap-26
    -1.0,                      // -0x1p+0
    -0.8455373048782349,       // -0x1.b0ea44p-1
     0.8724101781845093,       //  0x1.beac8cp-1
     0.08465760201215744,      //  0x1.5ac1eep-4
    -0.24683333933353424       // -0x1.f983c2p-3
  );
  sc_db_sl: array[0..8] of Double = (
    -7.450580596923828e-09,    // -0x1p-27
    -1.4901161193847656e-08,   // -0x1p-26
    -1.4901161193847656e-08,   // -0x1p-26
    -4.440892098500626e-16,    // -0x1p-51
     2.9802322387695312e-08,   //  0x1p-25
     1.4901161193847656e-08,   //  0x1p-26
     1.4901161193847656e-08,   //  0x1p-26
    -9.313225746154785e-10,    // -0x1p-30
     3.725290298461914e-09     //  0x1p-28
  );
  sc_db_ch: array[0..8] of Double = (
    -0.937637984752655,        // -0x1.e01216p-1
     0.7452248930931091,       //  0x1.7d8e1ep-1
     0.25911521911621094,      //  0x1.09558p-2
    -1.0,                      // -0x1p+0
     1.1924880638503055e-08,   //  0x1.99bc5cp-27
     0.5339164137840271,       //  0x1.115d7ep-1
     0.4887744188308716,       //  0x1.f48148p-2
     0.996410071849823,        //  0x1.fe2976p-1
     0.9690579771995544        //  0x1.f0285ep-1
  );
  sc_db_cl: array[0..8] of Double = (
    -1.4901161193847656e-08,   // -0x1p-26
     1.4901161193847656e-08,   //  0x1p-26
    -7.450580596923828e-09,    // -0x1p-27
     2.9802322387695312e-08,   //  0x1p-25
    -2.220446049250313e-16,    // -0x1p-52
    -1.4901161193847656e-08,   // -0x1p-26
     7.450580596923828e-09,    //  0x1p-27
     1.4901161193847656e-08,   //  0x1p-26
    -1.4901161193847656e-08    // -0x1p-26
  );
var
  t32: Tb32u32;
  ax: LongWord;
  ii: Integer;
begin
  t32.f := x;
  ax := t32.u and $7FFFFFFF;
  for ii := 0 to 8 do begin
    if sc_db_uarg[ii] = ax then begin
      sout := sinf_add_sign(x, Single(sc_db_sh[ii]), Single(sc_db_sl[ii]));
      cout := Single(sc_db_ch[ii] + sc_db_cl[ii]);
      Exit;
    end;
  end;
end;

// Large-argument case for sincosf (|x| >= threshold)
// Uses sincos_rbig (k=e-124 — correct for sincosf)
procedure sincosf_big(x: Single; var sout, cout: Single);
var
  t32: Tb32u32;
  ax: LongWord;
  ia: Integer;
  z_r, z2, z4, aa, bb, s0, c0, s_d, c_d: Double;
  tr: Tb64u64;
  tail: UInt64;
begin
  t32.f := x;
  ax := t32.u shl 1;
  if ax >= $FF000000 then begin
    if (ax shl 8) <> 0 then begin
      sout := x + x;   // NaN: propagate
      cout := x + x;
      Exit;
    end;
    pcr_feraiseexcept_invalid;
    sout := pcr_nanf(nil);
    cout := pcr_nanf(nil);
    Exit;
  end;
  z_r := sincos_rbig(t32.u, @ia);
  z2  := z_r * z_r;
  z4  := z2 * z2;
  aa  := (sincos_a[0] + z2 * sincos_a[1]) + z4 * (sincos_a[2] + z2 * sincos_a[3]);
  bb  := (sincos_b[0] + z2 * sincos_b[1]) + z4 * (sincos_b[2] + z2 * sincos_b[3]);
  bb  := bb * z_r;
  s0  := sinf_tb[ia and 31];
  c0  := sinf_tb[(ia + 8) and 31];
  s_d := s0 + z_r * (aa * c0 - bb * s0);
  c_d := c0 - z_r * (aa * s0 + bb * c0);
  sout := Single(s_d);
  cout := Single(c_d);
  tr.f := c_d;
  tail := (tr.u + 6) and UInt64($0FFFFFFF);   // ~0ull>>36
  if tail <= 12 then
    sincosf_database(x, sout, cout);
end;

// ============================================================
//  pcr_sincosf  — main function
// ============================================================
procedure pcr_sincosf(x: Single; out s, c: Single);
var
  t32: Tb32u32;
  ax: LongWord;
  ia: Integer;
  z0_d, z_d, z2, z4, aa, bb, s0, c0: Double;
begin
  t32.f := x;
  ax := t32.u shl 1;
  z0_d := Double(x);
  if ax < $822D97C8 then begin   // |x| < 0x1.2d97c8p+3
    if ax < $73000000 then begin  // |x| < 0x1p-12
      if ax < $66000000 then begin  // |x| < 0x1p-25
        if ax = 0 then begin
          s := x;
          c := 1.0;
        end else begin
          s := pcr_fmaf(-x, pcr_fabsf(x), x);
          c := 1.0 - 2.9802322387695312e-08;  // 1 - 0x1p-25f
        end;
        Exit;
      end;
      // 2^-25 <= |x| < 2^-12
      s := Single((-0.1666666716337204 * Double(x)) * Double(x * x) + Double(x));
      c := Single((-0.5 * Double(x)) * Double(x) + 1.0);
      Exit;
    end;
    if ax = $812D97C8 then begin
      sincosf_database(x, s, c);
      Exit;
    end;
    z_d := sincos_rltl0(z0_d, @ia);
  end else begin
    if ax > $99000000 then begin
      sincosf_big(x, s, c);
      Exit;
    end;
    if ax = $8C333330 then begin
      sincosf_database(x, s, c);
      Exit;
    end;
    z_d := sincos_rltl(x, @ia);
  end;
  z2 := z_d * z_d;
  z4 := z2 * z2;
  aa := (sincos_a[0] + z2 * sincos_a[1]) + z4 * (sincos_a[2] + z2 * sincos_a[3]);
  bb := (sincos_b[0] + z2 * sincos_b[1]) + z4 * (sincos_b[2] + z2 * sincos_b[3]);
  aa := aa * z_d;
  bb := bb * z2;
  s0 := sinf_tb[ia and 31];
  c0 := sinf_tb[(ia + 8) and 31];
  s  := Single(s0 + (aa * c0 - bb * s0));
  c  := Single(c0 - (aa * s0 + bb * c0));
end;

// ============================================================
//  pcr_tanf  — port of tanf.c (cr_tanf)
// ============================================================
function pcr_tanf(x: Single): Single;
const
  tanf_db_uarg: array[0..7] of LongWord = (
    $3F8A1F62,   // 0x1.143ec4p+0
    $4D56D355,   // 0x1.ada6aap+27
    $57D7B0ED,   // 0x1.af61dap+48
    $5980445E,   // 0x1.0088bcp+52
    $63FC86FE,   // 0x1.f90dfcp+72
    $6A662711,   // 0x1.cc4e22p+85
    $6AD36709,   // 0x1.a6ce12p+86
    $72B505BB    // 0x1.6a0b76p+102
  );
  tanf_db_rh: array[0..7] of Double = (
     1.8670953512191772,      //  0x1.ddf9f6p+0
     0.23828700184822083,     //  0x1.e80304p-3
     0.3445502519607544,      //  0x1.60d1c8p-2
     1.7895326614379883,      //  0x1.ca1edp+0
     0.6748017072677612,      //  0x1.597f9cp-1
    -3.9000706672668457,      // -0x1.f33584p+1
    -0.8855070471763611,      // -0x1.c5612ep-1
    -1.8912676572799683       // -0x1.e42a1ep+0
  );
  tanf_db_rl: array[0..7] of Double = (
    -3.409718953073569e-16,   // -0x1.891d24p-52
     4.358793082672146e-18,   //  0x1.419f46p-58
    -3.268032112391647e-17,   // -0x1.2d6c3ap-55
     2.1771658420191705e-16,  //  0x1.f6053p-53
     1.7449127529366843e-16,  //  0x1.925978p-53
     8.173074377011539e-16,   //  0x1.d7254ap-51
    -1.2783292861530832e-16,  // -0x1.26c33ep-53
    -2.4787918922562627e-16   // -0x1.1dc906p-52
  );
var
  t32: Tb32u32;
  e_exp, i_q, jj: Integer;
  z, z2, z4, n_v, n2, d_v, d2, s0, s1, r1: Double;
  tr: Tb64u64;
  tail: UInt64;
  ax, sgn: LongWord;
  x2: Single;
begin
  t32.f := x;
  e_exp := Integer((t32.u shr 23) and $FF);
  if e_exp < 127 + 28 then begin   // |x| < 2^28
    if e_exp < 115 then begin       // |x| < 2^-13 (exponent < 127-12 = 115)
      if e_exp < 102 then begin     // |x| < 2^-26 (exponent < 127-25 = 102)
        Result := pcr_fmaf(x, pcr_fabsf(x), x);
        Exit;
      end;
      x2 := x * x;
      Result := pcr_fmaf(x, 0.3333333432674408 * x2, x);  // 0x1.555556p-2f
      Exit;
    end;
    z := tanf_rltl(x, @i_q);
  end else if e_exp < $FF then begin
    z := tanf_rbig(t32.u, @i_q);
  end else begin
    if (t32.u shl 9) <> 0 then begin
      Result := x + x;   // NaN
      Exit;
    end;
    pcr_feraiseexcept_invalid;
    Result := pcr_nanf(nil);
    Exit;
  end;
  z2  := z * z;
  z4  := z2 * z2;
  n_v := tanf_cn[0] + z2 * tanf_cn[1];
  n2  := tanf_cn[2] + z2 * tanf_cn[3];
  n_v := n_v + z4 * n2;
  d_v := tanf_cd[0] + z2 * tanf_cd[1];
  d2  := tanf_cd[2] + z2 * tanf_cd[3];
  d_v := d_v + z4 * d2;
  n_v := n_v * z;
  // s[] = {0,1}: s0 = i_q&1, s1 = 1-(i_q&1)
  if (i_q and 1) = 0 then begin
    s0 := 0.0;  s1 := 1.0;
  end else begin
    s0 := 1.0;  s1 := 0.0;
  end;
  r1 := (n_v * s1 - d_v * s0) / (n_v * s0 + d_v * s1);
  tr.f := r1;
  tail := (tr.u + 7) and UInt64($1FFFFFFF);   // ~0ull>>35
  if tail <= 14 then begin
    ax  := t32.u and $7FFFFFFF;
    sgn := t32.u shr 31;
    for jj := 0 to 7 do begin
      if tanf_db_uarg[jj] = ax then begin
        if sgn <> 0 then
          Result := Single(-tanf_db_rh[jj] - tanf_db_rl[jj])
        else
          Result := Single( tanf_db_rh[jj] + tanf_db_rl[jj]);
        Exit;
      end;
    end;
  end;
  Result := Single(r1);
end;

end.
