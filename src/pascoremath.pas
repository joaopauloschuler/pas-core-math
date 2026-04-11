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

end.
