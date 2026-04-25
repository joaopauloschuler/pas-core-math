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
  cPowRho0:    Tb64u64 = (u: QWord($C0874910EE4E8A27)); // -0x1.74910ee4e8a27p+9
  cPowRho1:    Tb64u64 = (u: QWord($C08483B8CCA421AF)); // -0x1.483b8cca421afp+9
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

// --- 5.02i: qint very-slow path — q_3 / p_3 / log_3 / exp_3 ---------------

// Mirrors q_3 (pow.c:246). |y| < 2^-12.52, returns r ≈ exp(y) as TQInt64.
// QInt convention here matches C (no ex offset).
procedure PowQ3(out r: TQInt64; const y: TQInt64);
begin
  // Q_3[0]: 64-bit coefficient
  MulQInt11(r, y, cPowQ3_0);
  AddQInt22(r, cPowQ3_1, r);

  // Q_3[2..6]: 128-bit coefficients
  MulQInt22(r, y, r); AddQInt22(r, cPowQ3_2, r);
  MulQInt22(r, y, r); AddQInt22(r, cPowQ3_3, r);
  MulQInt22(r, y, r); AddQInt22(r, cPowQ3_4, r);
  MulQInt22(r, y, r); AddQInt22(r, cPowQ3_5, r);
  MulQInt22(r, y, r); AddQInt22(r, cPowQ3_6, r);

  // Q_3[7..11]: full 256-bit, use mul_qint_33 + add_qint
  MulQInt33(r, y, r); AddQInt(r, cPowQ3_7,  r);
  MulQInt33(r, y, r); AddQInt(r, cPowQ3_8,  r);
  MulQInt33(r, y, r); AddQInt(r, cPowQ3_9,  r);
  MulQInt33(r, y, r); AddQInt(r, cPowQ3_10, r);
  MulQInt33(r, y, r); AddQInt(r, cPowQ3_11, r);

  // Q_3[12..14]: full mul_qint + add_qint
  MulQInt(r, y, r); AddQInt(r, cPowQ3_12, r);
  MulQInt(r, y, r); AddQInt(r, cPowQ3_13, r);
  MulQInt(r, y, r); AddQInt(r, cPowQ3_14, r);
end;

// Mirrors p_3 (pow.c:481). |z| <= 2^-13 with z having only the upper limb set.
// Returns r ≈ log(1+z) (with the leading +z absorbed by the polynomial — C
// note: "the polynomial has zero constant coefficient, we multiply by z").
procedure PowP3(out r: TQInt64; const z: TQInt64);
begin
  // P_3[0..1]: degree 18..17
  MulQInt11(r, cPowP3_0, z);
  AddQInt22(r, cPowP3_1, r);

  // P_3[2..3]: degree 16..15
  MulQInt11(r, r, z); AddQInt22(r, cPowP3_2, r);
  MulQInt11(r, r, z); AddQInt22(r, cPowP3_3, r);

  // P_3[4..7]: degree 14..11
  MulQInt21(r, r, z); AddQInt22(r, cPowP3_4, r);
  MulQInt21(r, r, z); AddQInt22(r, cPowP3_5, r);
  MulQInt21(r, r, z); AddQInt22(r, cPowP3_6, r);
  MulQInt21(r, r, z); AddQInt22(r, cPowP3_7, r);

  // P_3[8..13]: degree 10..5, full 256-bit accuracy
  MulQInt31(r, r, z); AddQInt(r, cPowP3_8,  r);
  MulQInt31(r, r, z); AddQInt(r, cPowP3_9,  r);
  MulQInt31(r, r, z); AddQInt(r, cPowP3_10, r);
  MulQInt31(r, r, z); AddQInt(r, cPowP3_11, r);
  MulQInt31(r, r, z); AddQInt(r, cPowP3_12, r);
  MulQInt31(r, r, z); AddQInt(r, cPowP3_13, r);

  // P_3[14..17]: degree 4..1
  MulQInt41(r, r, z); AddQInt(r, cPowP3_14, r);
  MulQInt41(r, r, z); AddQInt(r, cPowP3_15, r);
  MulQInt41(r, r, z); AddQInt(r, cPowP3_16, r);
  MulQInt41(r, r, z); AddQInt(r, cPowP3_17, r);

  // Final multiplication by z (polynomial has zero constant term).
  MulQInt41(r, r, z);
end;

// Mirrors log_3 (pow.c:802). Returns r ≈ log(x) with relative error < 2^-250.74.
procedure PowLog3(out r: TQInt64; const xIn: TQInt64);
var
  xLoc, z, p: TQInt64;
  E: Int64;
  i, j: UInt32;
begin
  xLoc := xIn;
  E := xLoc.ex;

  if xLoc.r0 > UInt64($B504F333F9DE6484) then begin
    Inc(E);
    i := UInt32(xLoc.r0 shr (63 + 1 - 7));
  end else
    i := UInt32(xLoc.r0 shr (63 - 7));

  xLoc.ex := xLoc.ex - E;

  MulQInt(z, xLoc, cPowInverse31[i - 90]);  // exact

  // Second lookup index: j = z.r0 >> (63 - 13 - z.ex). QInt matches C convention,
  // so the C formula transcribes directly.
  j := UInt32(z.r0 shr (63 - 13 - z.ex));

  MulQInt(z, z, cPowInverse32[j - 8128]);  // exact

  AddQInt(z, QINT_M_ONE, z);  // z := z - 1; exact (Sterbenz)

  // E * log(2)
  MulQIntInt(r, E, QINT_LOG2);

  PowP3(p, z);
  AddQInt(p, cPowLogInv32[j - 8128], p);
  AddQInt(p, cPowLogInv31[i - 90],  p);
  AddQInt(r, p, r);
end;

// Mirrors exp_3 (pow.c:1132). Returns r ≈ exp(x) with relative error < 2^-241.10.
// Caller guarantees |x| < 744.45.
procedure PowExp3(out r: TQInt64; const xIn: TQInt64);
var
  xLoc, Kq, y: TQInt64;
  kInt, M, i2, i1: Int64;
begin
  xLoc := xIn;

  // Round-to-integer multiplier; mul_qint_11 ignores the low limb of x (OK).
  MulQInt11(Kq, xLoc, QINT_LOG2_INV);
  kInt := QIntToI(Kq);

  // Reconstruct -kInt*log(2)/2^12 to full accuracy.
  MulQIntInt(Kq, kInt, QINT_LOG2);
  Kq.ex := Kq.ex - 12;
  Kq.sgn := Byte(Kq.sgn = 0);

  AddQInt(y, xLoc, Kq);  // exact (Sterbenz)

  M  := PowAsr12(kInt);
  i2 := (kInt shr 6) and $3F;
  i1 := kInt and $3F;

  PowQ3(r, y);
  MulQInt(r, cPowT13[i2], r);
  MulQInt(r, cPowT23[i1], r);
  r.ex := r.ex + M;
end;

// --- 5.02j: rounding-boundary detection (is_exact / exact_pow) ------------

// Count trailing zeros (64-bit). Undefined if x = 0.
function pcr_ctzll(x: UInt64): Integer; inline;
begin
  pcr_ctzll := Integer(BsfQWord(x));
end;


// Mirrors C fast_extract (pow.h): e = unbiased exponent (subnormal returns -1023);
// m = significand with implicit bit added for normals.
procedure PowFastExtract(out e: Int64; out m: UInt64; x: Double); inline;
var u: Tb64u64;
begin
  u.f := x;
  e := (u.u shr 52) and $7FF;
  m := (u.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12));
  if e <> 0 then m := m + (UInt64(1) shl 52);
  e := e - $3FF;
end;

// Mirrors C extract (pow.h): returns (e, m) such that x = 2^e * m with m odd.
procedure PowExtract(out e: Int64; out m: UInt64; x: Double); inline;
var u: Tb64u64; t: Integer;
begin
  u.f := x;
  e := (u.u shr 52) and $7FF;
  m := (u.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12));
  if e <> 0 then m := m + (UInt64(1) shl 52);
  if m = 0 then begin t := 0; end else t := pcr_ctzll(m);
  m := m shr t;
  if e <> 0 then
    e := e + Int64(t) - $433
  else
    e := Int64(t) - $432;
end;

// Mirrors is_int(x): true iff x is a finite integer.
function PowIsInt(x: Double): Boolean; inline;
begin
  PowIsInt := (x = pcr_roundeven(x));
end;

// Multiply x by 2^e. Mirrors C pow2(double*, int64_t). e in [-2148, 2046] for
// normal-range results; saturates outside.
procedure PowPow2(var x: Double; e: Int64);
var e2: Tb64u64;
begin
  if (e and 1) <> 0 then x := x * Double(2.0);
  e2.u := UInt64(((e shr 1) + $3FF) and $7FF) shl 52;
  x := (x * e2.f) * e2.f;
end;

// Mirrors round_54: shortcut rounding of a TDInt64 to 54 bits, used by exact_pow.
// Returns G (binary exponent of the LSB of k) and k (54-bit signed mantissa).
// Pascal +1 ex convention: G_C = ex_P - 54.
procedure PowRound54(out Glsb, k: Int64; const x: TDInt64); inline;
begin
  Glsb := x.ex - 54;
  k := Int64((x.hi shr 10) + ((x.hi shr 9) and 1));
end;

// Mirrors C exact_pow (pow.c:1232). Algorithm detectRoundingBoundaryCase from [4].
// On entry z is the dint approximation of x^y from phase 2 (relative err < 2^-117).
// Returns true (and writes r) when the result is exactly representable on 54 bits.
function PowExactPow(out r: Double; x, y: Double; const z: TDInt64;
                     {%H-}exact: Integer): Boolean;
var
  s: Int64;
  m, n, t, kpow, mpow: UInt64;
  E, F, Glsb, k, g: Int64;
  cnt: Integer;
  d: TDInt64;
begin
  PowExactPow := False;
  if z.sgn = 1 then s := -1 else s := 1;

  PowExtract(E, m, x);

  if m = 1 then begin
    // x is a power of 2
    if PowIsInt(Double(E) * y) then begin
      if z.sgn = 1 then r := -1.0 else r := 1.0;
      PowPow2(r, Trunc(Double(E) * y));
      PowExactPow := True;
    end;
    exit;
  end;

  if (y < 0.0) or (y > 34.0) then exit;

  PowExtract(F, n, y);
  if (n > 34) or (F < -5) then exit;

  if F < 0 then begin
    // case (b): check that E is divisible by 2^-F
    if (E and (Int64(UInt64($FFFFFFFFFFFFFFFF) shr (64 + F)))) <> 0 then exit;
    g := Int64((UInt64(E) shr -F)) * Int64(n);
    PowRound54(Glsb, k, z);

    // |2^G*k - z| >= 2^-116*z ?
    if k = 0 then exit;
    cnt := pcr_clzll(UInt64(k));
    d.hi := UInt64(k) shl cnt;
    d.lo := 0;
    d.ex := Glsb + 64 - cnt;       // Pascal +1 convention
    d.sgn := Byte(1 - z.sgn);
    AddDInt(d, z, d);              // exact (Sterbenz)
    d.ex := d.ex + 116;
    if CmpDIntAbs(d, z) >= 0 then exit;

    if Glsb > g then exit;

    // Check whether k is an odd number multiplied by 2^(g-Glsb).
    if (UInt64(k) and (not (UInt64($FFFFFFFFFFFFFFFE) shl (g - Glsb)))) =
       (UInt64(1) shl (g - Glsb)) then begin
      r := Double((Int64(UInt64(k) shr (g - Glsb))) * s);
      PowPow2(r, g);
      PowExactPow := True;
      exit;
    end;
    exit;
  end;

  // case (a): n*2^F integer, 1 <= y <= 34
  t := n shl F;
  kpow := 1;
  mpow := m;
  while t <> 0 do begin
    if (t and 1) <> 0 then begin
      // detect overflow of mpow*kpow
      if (mpow <> 0) and (kpow > UInt64($FFFFFFFFFFFFFFFF) div mpow) then exit;
      kpow := mpow * kpow;
    end;
    t := t shr 1;
    if t <> 0 then begin
      if (mpow <> 0) and (mpow > UInt64($FFFFFFFFFFFFFFFF) div mpow) then exit;
      mpow := mpow * mpow;
    end;
  end;
  if (kpow shr 54) <> 0 then exit;
  r := Double(Int64(kpow) * s);
  PowPow2(r, E * Int64(n shl F));
  PowExactPow := True;
end;

// Returns 1 (or sign-bit-style nonzero) if x^y is exactly representable.
// Mirrors C is_exact (pow.c:1349). 0 = not exact, nonzero = exact.
function PowIsExact(x, y: Double): Integer;
const
  cXMax: array[0..33] of UInt64 = (
    0, UInt64($FFFFFFFFFFFFFFFF),
    94906265, 208063, 9741, 1551, 455, 189, 97,
    59, 39, 27, 21, 15, 13, 11, 9, 7, 7, 5, 5,
    5, 5, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 );
var
  v, w: Tb64u64;
  m, n, my, n0: UInt64;
  e, f, ez: Int64;
  t, y_int, i: Integer;
  dm, sd: Double;
begin
  PowIsExact := 0;
  v.f := x; w.f := y;

  // Early exit: |x|<>1 AND low 42 bits of y nonzero.
  if ((v.u shl 1) <> UInt64($7FE0000000000000)) and ((w.u shl 22) <> 0) then exit;

  if (v.u shl 1) = UInt64($7FE0000000000000) then begin
    PowIsExact := 1; exit;
  end;

  if (y >= 0) and PowIsInt(y) then begin
    m := v.u and UInt64($FFFFFFFFFFFFF);
    e := Int64((v.u shl 1) shr 53) - $433;
    if e >= -1074 then m := m or UInt64($10000000000000)
    else e := e + 1;
    if m = 0 then exit;
    t := pcr_ctzll(m);
    m := m shr t;
    e := e + t;
    if (y = 0.0) or (y = 1.0) then begin PowIsExact := 1; exit; end;
    if m = 1 then begin
      if (-1074 <= y * e) and (y * e < 1024) then PowIsExact := 1;
      exit;
    end;
    if (y < 0) or (33 < y) then exit;
    y_int := Trunc(y);
    if m > cXMax[y_int] then exit;
    my := m * m;
    for i := 2 to y_int - 1 do my := my * m;
    t := 64 - pcr_clzll(m);
    ez := e * y_int + t;
    if (ez <= -1074) or (1024 < ez) then exit;
    if e * y_int >= -1074 then PowIsExact := 1;
    exit;
  end;

  // y = n*2^f (n odd)
  n := w.u and UInt64($FFFFFFFFFFFFF);
  f := Int64((w.u shl 1) shr 53) - $433;
  if f >= -1074 then n := n or UInt64($10000000000000)
  else f := f + 1;
  if n = 0 then exit;
  t := pcr_ctzll(n);
  n := n shr t;
  f := f + t;

  m := v.u and UInt64($FFFFFFFFFFFFF);
  e := Int64((v.u shl 1) shr 53) - $433;
  if e >= -1074 then m := m or UInt64($10000000000000)
  else e := e + 1;
  if m = 0 then exit;
  t := pcr_ctzll(m);
  m := m shr t;
  e := e + t;

  if y < 0 then begin
    if m <> 1 then exit;
    if f >= 0 then begin
      if f < 12 then ez := (-Int64(n) * e) shl f else ez := 1024;
    end else begin
      t := pcr_ctzll(UInt64(e));
      if -f > t then exit;
      ez := (-e shr -f) * Int64(n);
    end;
    if (-1074 <= ez) and (ez < 1024) then PowIsExact := 1;
    exit;
  end;

  // y > 0, non-integer, y = n*2^f with n odd, f < 0.
  // x must be a 2^(-f)-th power; m must be a square repeatedly.
  while True do begin
    f := f + 1; if f > 0 then break;
    if (e and 1) <> 0 then exit;
    e := e div 2;
    dm := Double(m);
    sd := Double(Round(Sqrt(dm)));
    if sd * sd <> dm then exit;
    m := UInt64(Trunc(sd));
  end;

  if m > 1 then begin
    if 33 < n then exit;
    if m > cXMax[n] then exit;
  end;
  my := m; n0 := n;
  while n0 > 1 do begin my := my * m; n0 := n0 - 1; end;
  t := 64 - pcr_clzll(my);
  if (-1074 <= e * Int64(n)) and (e * Int64(n) + t <= 1024) then
    PowIsExact := 1;
end;

// --- 5.02k: main pcr_pow_pas with three Ziv phases ------------------------

// Convert a double to a TQInt64. Mirrors qint_fromd (pow.h:371).
procedure PowQIntFromD(out a: TQInt64; b: Double); inline;
var t: Integer;
begin
  PowFastExtract(a.ex, a.r0, b);
  t := pcr_clzll(a.r0);
  a.sgn := Byte(b < 0.0);
  if t > 11 then a.ex := a.ex - (t - 12);
  a.r0 := a.r0 shl t;
  a.r1 := 0; a.r2 := 0; a.r3 := 0;
end;

// Subnormalize a TQInt64 prior to conversion. Mirrors C subnormalize_qint.
// QInt convention matches C — use C thresholds directly.
procedure PowSubnormalizeQInt(var a: TQInt64);
var
  ex: UInt64;
  hi, md, lo: UInt64;
begin
  if a.ex > -1023 then exit;
  ex := UInt64(-(1011 + a.ex));
  hi := a.r0 shr ex;
  md := (a.r0 shr (ex - 1)) and 1;
  if ((a.r0 and (UInt64($FFFFFFFFFFFFFFFF) shr ex)) or a.r1 or a.r2 or a.r3) <> 0 then
    lo := 1
  else
    lo := 0;
  // round-to-nearest-even (default; we don't query fegetround)
  if lo <> 0 then hi := hi + md
  else hi := hi + (hi and md);
  a.r0 := hi shl ex;
  a.r1 := 0; a.r2 := 0; a.r3 := 0;
  if a.r0 = 0 then begin
    a.ex := a.ex + 1;
    a.r0 := UInt64(1) shl 63;
  end;
end;

// Convert a TQInt64 to a double. Mirrors qint_tod (pow.h:431).
function PowQIntToD(var a: TQInt64): Double;
const
  c2pm53: Tb64u64 = (u: $3CA0000000000000);
  c2pm54: Tb64u64 = (u: $3C90000000000000);
  c1p1023: Tb64u64 = (u: $7FE0000000000000);
var
  r, e: Tb64u64;
  rd: Double;
begin
  PowSubnormalizeQInt(a);
  r.u := (a.r0 shr 11) or (UInt64($3FF) shl 52);
  rd := 0.0;
  if (a.r0 and $400) <> 0 then rd := rd + c2pm53.f;
  if ((a.r0 and $3FF) <> 0) or (a.r1 <> 0) or (a.r2 <> 0) or (a.r3 <> 0) then
    rd := rd + c2pm54.f;
  if a.sgn <> 0 then rd := -rd;
  r.u := r.u or (UInt64(a.sgn) shl 63);
  r.f := r.f + rd;

  if a.ex > -1023 then begin
    if a.ex > 1023 then begin
      if a.ex = 1024 then begin
        r.f := r.f * Double(2.0);
        e.u := c1p1023.u;
      end else begin
        r.f := cPowDblMax.f;
        e.f := cPowDblMax.f;
      end;
    end else
      e.u := UInt64((a.ex + 1023) and $7FF) shl 52;
  end else begin
    if a.ex < -1074 then begin
      if a.ex = -1075 then begin
        r.f := r.f * Double(0.5);
        e.f := cPowMinSub.f;
      end else begin
        r.f := cPowMinSub.f;
        e.f := cPowMinSub.f;
      end;
    end else
      e.u := UInt64(1) shl UInt64(a.ex + 1074);
  end;
  PowQIntToD := r.f * e.f;
end;

// Returns 1 if x is a signaling NaN (assumes x is a NaN).
function PowIsSignaling(x: Double): Boolean; inline;
var u: Tb64u64;
begin
  u.f := x;
  PowIsSignaling := (u.u and (UInt64(1) shl 51)) = 0;
end;

// Main entry point. Correctly-rounded x^y for binary64.
// Mirrors cr_pow (pow.c:1495).
function pcr_pow(x, y: Double): Double;
const
  cInfPos:    Tb64u64 = (u: $7FF0000000000000);
  cInfNeg:    Tb64u64 = (u: QWord($FFF0000000000000));
  cPow1pm100: Tb64u64 = (u: $39B0000000000000); // 0x1p-100
var
  s, x0: Double;
  xu, yu: Tb64u64;
  cancel: Integer;
  lh, ll, rh, rl, res_h, res_l, res_min, res_max, errc: Double;
  ey: Integer;
  cs0, cs1: Double;
  y_parity: Integer;
  ay: Double;
  exact: Integer;
  Xd, Yd, Rd_dint: TDInt64;
  qX, qY, qR, qZ: TQInt64;
  rd: UInt64;
  ex64: UInt64;
  m54: UInt64;
  err_arr: array[0..1] of Tb64u64;
  e_exact: Double;
  r1q, r2q, r3q: UInt64;
begin
  err_arr[0].u := $3BF2700000000000;  // 0x1.27p-64 (~ 2^-63.797)
  err_arr[1].u := $3C55700000000000;  // 0x1.57p-58 (~ 2^-57.579)

  s := 1.0;
  x0 := x;
  xu.f := x; yu.f := y;

  // === NaN / Inf special cases ===
  if (xu.u >= UInt64($7FF0000000000000)) or (yu.u >= UInt64($7FF0000000000000)) then
  begin
    if x <> x then begin   // x is NaN
      if (y = 0.0) and not PowIsSignaling(x) then begin pcr_pow := 1.0; exit; end;
      pcr_pow := x + x; exit;
    end;
    if y <> y then begin   // y is NaN
      if (x = 1.0) and not PowIsSignaling(y) then begin pcr_pow := 1.0; exit; end;
      pcr_pow := y + y; exit;
    end;

    // x = +Inf
    if xu.u = cInfPos.u then begin
      if y = 0.0 then begin pcr_pow := 1.0; exit; end;
      if y < 0.0 then begin pcr_pow := 0.0; exit; end;
      if y > 0.0 then begin pcr_pow := cInfPos.f; exit; end;
    end;
    // x = -Inf
    if xu.u = cInfNeg.u then begin
      if PowIsInt(y) and not PowIsInt(y * 0.5) then begin
        if y < 0.0 then begin pcr_pow := -0.0; exit; end;
        pcr_pow := -cInfPos.f; exit;
      end;
      if y < 0.0 then begin pcr_pow := 0.0; exit; end;
      if y > 0.0 then begin pcr_pow := cInfPos.f; exit; end;
    end;

    // y = +Inf
    if yu.u = cInfPos.u then begin
      if x = 0.0 then begin pcr_pow := 0.0; exit; end;
      if (x = -1.0) or (x = 1.0) then begin pcr_pow := 1.0; exit; end;
      if (-1.0 < x) and (x < 1.0) then begin pcr_pow := 0.0; exit; end;
      if (x < -1.0) or (1.0 < x) then begin pcr_pow := cInfPos.f; exit; end;
    end;
    // y = -Inf
    if yu.u = cInfNeg.u then begin
      if x = 0.0 then begin pcr_pow := cInfPos.f; exit; end;
      if (x = -1.0) or (x = 1.0) then begin pcr_pow := 1.0; exit; end;
      if (-1.0 < x) and (x < 1.0) then begin pcr_pow := cInfPos.f; exit; end;
      if (x < -1.0) or (1.0 < x) then begin pcr_pow := 0.0; exit; end;
    end;
  end;

  // === x <= 0 cascade ===
  if x <= 0.0 then begin
    if y = 0.0 then begin pcr_pow := 1.0; exit; end;

    if xu.u = 0 then begin   // x = +0.0
      if PowIsInt(y) and not PowIsInt(y * 0.5) then begin
        if y < 0.0 then begin pcr_pow := cInfPos.f; exit; end;
        pcr_pow := 0.0; exit;
      end;
      if y > 0.0 then begin pcr_pow := 0.0; exit; end;
      pcr_pow := cInfPos.f; exit;
    end;
    if xu.u = UInt64($8000000000000000) then begin   // x = -0.0
      if PowIsInt(y) and not PowIsInt(y * 0.5) then begin
        if y < 0.0 then begin pcr_pow := -cInfPos.f; exit; end;
        pcr_pow := -0.0; exit;
      end;
      if y > 0.0 then begin pcr_pow := 0.0; exit; end;
      pcr_pow := cInfPos.f; exit;
    end;

    if not PowIsInt(y) then begin pcr_pow := cNaNDoublePos.f; exit; end;

    cs0 := 1.0; cs1 := -1.0;
    ay := y; if ay < 0 then ay := -ay;
    if ay >= Double(9007199254740992.0) then y_parity := 0  // 2^53
    else y_parity := Integer(Trunc(y) and 1);
    if y_parity = 0 then s := cs0 else s := cs1;
    x := -x;
  end;

  // === Phase 1: fast path (double-double) ===
  cancel := PowLog1(lh, ll, x);

  // Avoid spurious underflow/overflow in y*log(x).
  ey := (yu.u shr 52) and $7FF;
  if (ey < $36) or (ey >= $7F5) then begin
    lh := cNaNDoublePos.f; ll := cNaNDoublePos.f;
  end;

  pcr_s_mul(rh, rl, y, lh, ll);

  PowExp1(res_h, res_l, rh, rl, s);

  errc := err_arr[cancel].f;
  res_min := res_h + pcr_fma(errc, -res_h, res_l);
  res_max := res_h + pcr_fma(errc,  res_h, res_l);

  exact := PowIsExact(x, y);

  if res_min = res_max then begin
    pcr_pow := res_max; exit;
  end;

  // Easy cases short-circuit
  if y = 1.0 then begin pcr_pow := s * x; exit; end;
  if y = 2.0 then begin pcr_pow := x * x; exit; end;
  if y = 0.5 then begin pcr_pow := Sqrt(x); exit; end;
  if y = 0.0 then begin pcr_pow := 1.0; exit; end;

  // === Phase 2: dint slow path ===
  DIntFromD(Xd, x);
  DIntFromD(Yd, y);
  Xd.sgn := 0;

  PowLog2(Rd_dint, Xd);
  MulDInt21(Rd_dint, Rd_dint, Yd);
  PowExp2(Rd_dint, Rd_dint);

  // Rounding test (Pascal +1 ex convention: subtract 1 from C thresholds).
  // C R.ex < -1075 ; underflow -> -1074 in Pascal
  if Rd_dint.ex < -1074 then begin
    pcr_pow := 0.5 * (s * cPowMinSub.f); exit;
  end;

  if Rd_dint.ex < -1021 then begin   // C: < -1022
    // subnormal: ex_C = -(1022 + R.ex_C); Pascal R.ex_P = R.ex_C + 1
    ex64 := UInt64(-(1021 + Rd_dint.ex));   // 1 .. 53
    m54 := (Rd_dint.lo shr (10 + ex64)) or (Rd_dint.hi shl (54 - ex64));
    if m54 + 14 > 28 then rd := 1 else rd := 0;
  end else begin
    // normal case: rd = 0 iff lo <= 28 or lo >= 2^64-28
    m54 := (Rd_dint.lo shr 10) or (Rd_dint.hi shl 54);
    if m54 + 28 > 56 then rd := 1 else rd := 0;
  end;

  Rd_dint.sgn := Byte(s = -1.0);

  if rd <> 0 then begin
    pcr_pow := DToD(Rd_dint); exit;
  end;

  // Detect rounding boundary cases
  if PowExactPow(e_exact, x0, y, Rd_dint, exact) then begin
    pcr_pow := e_exact; exit;
  end;

  // === Phase 3: qint very-slow path ===
  PowQIntFromD(qX, x);
  PowQIntFromD(qY, y);
  qX.sgn := 0;

  PowLog3(qR, qX);
  MulQInt41(qR, qR, qY);
  PowExp3(qZ, qR);

  // Extra rounding test
  r1q := (qZ.r0 shl 54) or (qZ.r1 shr 10);
  r2q := (qZ.r1 shl 54) or (qZ.r2 shr 10);
  r3q := (qZ.r2 shl 54) or (qZ.r3 shr 10);
  if ((r1q = 0) and (r2q = 0) and (r3q <= 60)) or
     (((not r1q) = 0) and (((not r2q) = 0)) and (r3q + 120 <= 60)) then
    rd := 0
  else
    rd := 1;

  if rd <> 0 then begin
    qZ.sgn := Byte(s = -1.0);
    qZ.r3 := qZ.r3 and (UInt64($FFFFFFFFFFFFFFFF) shl 10);
    pcr_pow := PowQIntToD(qZ); exit;
  end;

  // qR very small: x^y near 1
  if qR.ex < -56 then begin
    if qR.sgn = 0 then pcr_pow := 1.0 + cPow1pm100.f
    else pcr_pow := 1.0 - cPow1pm100.f;
    exit;
  end;

  // Should not reach here for any valid input.
  pcr_pow := cNaNDoublePos.f;
end;
