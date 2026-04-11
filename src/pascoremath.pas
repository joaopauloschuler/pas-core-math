{$I pascoremath.inc}
unit pascoremath;

interface

uses pascoremathtypes, pascoremathhelperfuncs;

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
    // recompute e from au after shift (au is now a 25-bit value for subnormals)
    e := au shr 24;
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
    0.3125
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
    // as_special
    if ux >= $FF000000 then // inf or nan
    begin
      if ux > $FF000000 then begin Result := x + x; Exit; end; // nan
      if t.u shr 31 <> 0 then begin Result := 0.0; Exit; end; // -inf -> 0
      Result := x; Exit; // +inf
    end;
    if t.u >= $C3150000 then // x < -149
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
    0.3125
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
  j_idx := (Int64(Int64(jt.u shl 12) shr 48)); // sign-extend 16-bit field
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

end.
