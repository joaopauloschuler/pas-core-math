// ---------------------------------------------------------------------------
// pcr_pow — correctly-rounded binary64 power.
// Ported from core-math/src/binary64/pow/{pow.c, pow.h, dint.h, qint.h} by
// Paul Zimmermann, Tom Hubrecht and Claude Etienne Rohrbacher.
//
// Three-tier Ziv strategy:
//   * fast path: pure double-double, ~2^-67 absolute error
//     (uses pow.h tables: cPowInverse, cPowLogInv*, cPowT1*, cPowT2*, cPowP1, cPowQ1)
//   * slow path: TDInt64, ~2^-122 absolute error
//     (uses dint tables: cPowInverse2*, cPowLogInv2*, cPowT12, cPowT22, cPowP2, cPowQ2)
//   * very-slow path: TQInt64, ~2^-242 absolute error
//     (uses qint tables: cPowInverse3*, cPowLogInv3*, cPowT13, cPowT23, cPowP3, cPowQ3)
//
// Conventions:
//   * TDInt64 .ex shifted by +1 vs pow/dint.h (mantissa in [0.5,1) here vs [1,2) there)
//   * TQInt64 .ex matches pow/qint.h exactly (no offset)
//   * Tables are pre-shifted at generation time; see tools/pow_gen.py.
// ---------------------------------------------------------------------------

{$I pow_const.inc}

const
  // Named constants used by the fast-path log/exp helpers (anti-x87 typecasts).
  cPowLog2H:  Tb64u64 = (u: $3FE62E42FEFA3800);  // 0x1.62e42fefa38p-1
  cPowLog2L:  Tb64u64 = (u: $3D2EF35793C76730);  // 0x1.ef35793c7673p-45
  cPowLog1Threshold: Tb64u64 = (u: $3E70000000000000); // 0x1p-24

  // exp_1 thresholds and constants
  cPowRho0:    Tb64u64 = (u: $C0874910EE4E8A27); // -0x1.74910ee4e8a27p+9
  cPowRho1:    Tb64u64 = (u: $C08483B8CCA421AF); // -0x1.483b8cca421afp+9
  cPowRho2:    Tb64u64 = (u: $40862E42E709A95B); // 0x1.62e42e709a95bp+9
  cPowRho3:    Tb64u64 = (u: $40862E4316EA5DF9); // 0x1.62e4316ea5df9p+9
  cPowInvLog2: Tb64u64 = (u: $40B71547652B82FE); // 0x1.71547652b82fep+12
  cPowExpLog2H:Tb64u64 = (u: $3F262E42FEFA39EF); // 0x1.62e42fefa39efp-13
  cPowExpLog2L:Tb64u64 = (u: $3BBABC9E3B39803F); // 0x1.abc9e3b39803fp-68
  cPowDblMax:  Tb64u64 = (u: $7FEFFFFFFFFFFFFF); // 0x1.fffffffffffffp+1023
  cPowMinSub:  Tb64u64 = (u: $0000000000000001); // 0x1p-1074

// --- 5.02e: q_1 / p_1 — pure-double polynomial helpers --------------------

// Given |z| < 2^-12.905, compute (qh, ql) ≈ exp(z) with relative error < 2^-64.9.
// Mirrors pow.c:120 (q_1).
procedure PowQ1(out qh, ql: Double; z: Double);
var
  q, h0, h1, l1: Double;
begin
  q  := pcr_fma(cPowQ1_4.f, z, cPowQ1_3.f);
  q  := pcr_fma(q,           z, cPowQ1_2.f);
  h0 := pcr_fma(q,           z, cPowQ1_1.f);
  pcr_a_mul(h1, l1, z, h0);                 // exact: h1 + l1 = z * h0
  pcr_fasttwosum(qh, ql, cPowQ1_0.f, h1);   // |Q_1[0]| = 1 dominates
  ql := ql + l1;
end;

// Given |z| < 32*2^-13, compute (ph, pl) ≈ log(1+z) - z + z^2/2 such that
//   |ph + pl - (log(1+z) - z)| < 2^-75.5
// (ph holds -z^2/2's leading part; the polynomial tail is folded into pl).
// Mirrors pow.c:323 (p_1).
procedure PowP1(out ph, pl: Double; z: Double);
var
  wh, wl, t, u, v: Double;
begin
  pcr_a_mul(wh, wl, z, z);                  // exact: wh + wl = z^2
  t := pcr_fma(cPowP1_5.f, z, cPowP1_4.f);
  u := pcr_fma(cPowP1_3.f, z, cPowP1_2.f);
  v := pcr_fma(cPowP1_1.f, z, cPowP1_0.f);
  u := pcr_fma(t, wh, u);
  v := pcr_fma(u, wh, v);
  u := v * wh;
  ph := Double(-0.5) * wh;
  pl := pcr_fma(u, z, Double(-0.5) * wl);
end;

// --- 5.02f: log_1 — fast log(x) returning (h, l, cancel) ------------------
// Mirrors pow.c:581 (log_1). Returns 1 when massive cancellation triggered
// a fast_two_sum normalization (caller may use it as a "rerun in slow path"
// hint), else 0.
function PowLog1(out h, l: Double; x: Double): Integer;
var
  xu: Tb64u64;
  m, mraw: UInt64;
  e: Int64;
  t_u: Tb64u64;
  k, c: UInt32;
  i: UInt64;
  t, ee, r, l1, l2, z, th, tl, ph, pl, hh, ll: Double;
begin
  xu.f := x;
  mraw := xu.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12);
  e    := (xu.u shr 52) and $7FF;

  if e <> 0 then begin
    t_u.u := mraw or (UInt64($3FF) shl 52);
    m := mraw + (UInt64(1) shl 52);
    e := e - $3FF;
  end else begin
    // subnormal
    k := UInt32(pcr_clzll(mraw)) - 11;
    e := -$3FE - Int64(k);
    m := mraw shl k;
    t_u.u := m or (UInt64($3FF) shl 52);
  end;

  t := t_u.f;

  // c = 1 when m >= 0x16a09e667f3bcd (~ sqrt(2)*2^52); else 0
  if m >= UInt64($16A09E667F3BCD) then c := 1 else c := 0;

  e := e + Int64(c);
  ee := Int64(e);  // exact int->double for |e|<=1024

  if c = 1 then begin
    i := m shr 45;
    t := t * Double(0.5);
  end else begin
    i := m shr 44;
  end;

  // i is in [128..181] when c=0, [90..127] when c=1; index into 182-entry table
  r  := cPowInverse  [i - 181].f;
  l1 := cPowLogInvHi [i - 181].f;
  l2 := cPowLogInvLo [i - 181].f;

  z := pcr_fma(r, t, Double(-1.0));

  th := pcr_fma(ee, cPowLog2H.f, l1);
  tl := pcr_fma(ee, cPowLog2L.f, l2);

  // fast_sum(h, l, th, z, tl)
  pcr_fasttwosum(h, l, th, z);
  l := l + tl;

  PowP1(ph, pl, z);

  // fast_sum(h, l, h, ph, l + pl)
  hh := h;
  ll := l + pl;
  pcr_fasttwosum(h, l, hh, ph);
  l := l + ll;

  if (e = 0) and (pcr_fabs(l) > pcr_fabs(h) * cPowLog1Threshold.f) then begin
    pcr_fasttwosum(h, l, h, l);
    PowLog1 := 1;
    exit;
  end;
  PowLog1 := 0;
end;

// --- 5.02g: exp_1 — fast exp(rh+rl)*s with overflow/underflow shortcuts ---
// Mirrors pow.c:954 (exp_1). On NaN input or borderline overflow/underflow
// returns NaN in eh/el so caller can defer to the slow path.
procedure PowExp1(out eh, el: Double; rh, rl, s: Double);
var
  kd, zh, zl, qh, ql: Double;
  kInt, M, i2, i1: Int64;
  d: Tb64u64;
  rh_is_nan: Boolean;
begin
  rh_is_nan := rh <> rh;
  if rh_is_nan or (rh > cPowRho2.f) then begin
    if (not rh_is_nan) and (rh > cPowRho3.f) then begin
      // certain overflow; return s*DBL_MAX so rounding in caller picks +/-Inf or DBL_MAX
      eh := cPowDblMax.f * s;
      el := cPowDblMax.f * s;
    end else begin
      // intermediate region or NaN: defer to slow path
      eh := cNaNDoublePos.f;
      el := cNaNDoublePos.f;
    end;
    exit;
  end;

  if rh < cPowRho1.f then begin
    if rh < cPowRho0.f then begin
      // certain underflow
      eh := Double(0.0) * s;
      el := cPowMinSub.f * (Double(0.5) * s);
    end else begin
      eh := cNaNDoublePos.f;
      el := cNaNDoublePos.f;
    end;
    exit;
  end;

  kd := pcr_roundeven(rh * cPowInvLog2.f);

  zh := pcr_fma(cPowExpLog2H.f, -kd, rh);
  zl := pcr_fma(cPowExpLog2L.f, -kd, rl);

  kInt := Trunc(kd);
  M  := (kInt shr 12) + $3FF;
  i2 := (kInt shr 6) and $3F;
  i1 := kInt and $3F;

  pcr_d_mul(eh, el, cPowT2Hi[i1].f, cPowT2Lo[i1].f, cPowT1Hi[i2].f, cPowT1Lo[i2].f);

  PowQ1(qh, ql, zh + zl);

  pcr_d_mul(eh, el, eh, el, qh, ql);

  d.u := UInt64(M) shl 52;
  d.f := d.f * s;
  eh := eh * d.f;
  el := el * d.f;
end;

// --- 5.02h: dint slow path — q_2 / p_2 / log_2 / exp_2 --------------------

// Arithmetic shift right by 12 for signed Int64 (FPC's `shr` is logical).
// Used in exp_2/log_2 to compute floor(k/2^12). Range: |a| <= 4399162.
function PowAsr12(a: Int64): Int64; inline;
begin
  if a >= 0 then PowAsr12 := a shr 12
  else PowAsr12 := -Int64((QWord(-a-1) shr 12) + 1);
end;

// Mirrors q_2 (pow.c:137). |y| < 2^-12.52, returns r ≈ exp(y) as TDInt64.
procedure PowQ2(out r: TDInt64; const y: TDInt64);
begin
  MulDInt11(r, y, cPowQ2_0);
  AddDInt11(r, cPowQ2_1, r);
  MulDInt11(r, y, r);
  AddDInt11(r, cPowQ2_2, r);
  MulDInt11(r, y, r);
  AddDInt(r, cPowQ2_3, r);
  MulDInt(r, y, r);
  AddDInt(r, cPowQ2_4, r);
  MulDInt(r, y, r);
  AddDInt(r, cPowQ2_5, r);
  MulDInt(r, y, r);
  AddDInt(r, cPowQ2_6, r);
  MulDInt(r, y, r);
  AddDInt(r, cPowQ2_7, r);
end;

// Mirrors p_2 (pow.c:343). |z| <= 2^-13, returns r ≈ log(1+z) as TDInt64.
procedure PowP2(out r: TDInt64; const z: TDInt64);
begin
  MulDInt11(r, z, cPowP2_0);
  AddDInt11(r, cPowP2_1, r);
  MulDInt11(r, z, r);
  AddDInt11(r, cPowP2_2, r);
  MulDInt11(r, z, r);
  AddDInt11(r, cPowP2_3, r);
  MulDInt11(r, z, r);
  AddDInt(r, cPowP2_4, r);
  MulDInt21(r, r, z);
  AddDInt(r, cPowP2_5, r);
  MulDInt21(r, r, z);
  AddDInt(r, cPowP2_6, r);
  MulDInt21(r, r, z);
  AddDInt(r, cPowP2_7, r);
  MulDInt21(r, r, z);
  AddDInt(r, cPowP2_8, r);
  MulDInt21(r, r, z);
end;

// Mirrors log_2 (pow.c:652). Returns r ≈ log(x) with relative error < 2^-122.88.
// xIn is taken by const-ref; caller's x is preserved (we mutate a local copy).
procedure PowLog2(out r: TDInt64; const xIn: TDInt64);
var
  xLoc, z, p: TDInt64;
  E: Int64;
  i, j: UInt32;
begin
  xLoc := xIn;
  // E must be the C-convention binary exponent of x (= Pascal ex - 1) so that
  // E*log(2) and the post-reduction xLoc both end up with the right magnitude.
  E := xLoc.ex - 1;

  if xLoc.hi > UInt64($B504F333F9DE6484) then begin
    Inc(E);
    i := UInt32(xLoc.hi shr (63 + 1 - 7));
  end else
    i := UInt32(xLoc.hi shr (63 - 7));

  xLoc.ex := xLoc.ex - E;

  MulDInt11(z, xLoc, cPowInverse21[i - 90]);  // exact

  // 2nd lookup index. C: j = z.hi >> (63 - 13 - z.ex). Convention compatible
  // since shift is signed by-difference; xLoc.ex was relative, and our +1
  // offset cancels in the subtraction below.
  // Pascal's TDInt64.ex is C ex + 1, so the C shift "63-13-z.ex_C" becomes
  //   63 - 13 - (z.ex_P - 1) = 51 - z.ex_P
  j := UInt32(z.hi shr (51 - z.ex));

  MulDInt11(z, z, cPowInverse22[j - 8128]);

  AddDInt(z, DINT_M_ONE, z);  // z := z - 1; exact (Sterbenz)

  // E * log(2)
  MulDIntInt(r, E, DINT_LOG2);

  PowP2(p, z);
  AddDInt(p, cPowLogInv22[j - 8128], p);
  AddDInt(p, cPowLogInv21[i - 90],  p);
  AddDInt(r, p, r);
end;

// Mirrors exp_2 (pow.c:1061). Returns r ≈ exp(x) with relative error < 2^-121.70.
// In Pascal +1 ex convention, the C threshold "x->ex >= 10" becomes ">= 11";
// the underflow/overflow result ex values shift accordingly (-1076→-1075, 1025→1026).
procedure PowExp2(out r: TDInt64; const xIn: TDInt64);
var
  xLoc, Kdint, y: TDInt64;
  kInt, M, i2, i1: Int64;
begin
  xLoc := xIn;

  if xLoc.ex >= 11 then begin
    CpDInt(r, xLoc);
    if xLoc.sgn = 1 then r.ex := -1075 else r.ex := 1026;
    r.sgn := 0;
    Exit;
  end;

  MulDInt11(Kdint, xLoc, DINT_LOG2_INV);
  kInt := DIntToI(Kdint);

  MulDIntInt(Kdint, kInt, DINT_LOG2);
  Kdint.ex := Kdint.ex - 12;
  Kdint.sgn := Byte(Kdint.sgn = 0);

  AddDInt(y, xLoc, Kdint);

  M  := PowAsr12(kInt);
  i2 := (kInt shr 6) and $3F;   // mask cleans up logical-vs-arithmetic difference
  i1 := kInt and $3F;

  PowQ2(r, y);
  MulDInt(r, cPowT12[i2], r);
  MulDInt(r, cPowT22[i1], r);
  r.ex := r.ex + M;
end;
