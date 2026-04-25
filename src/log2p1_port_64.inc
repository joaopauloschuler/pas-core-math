// ===========================================================================
// pcr_log2p1 — correctly-rounded binary64 log2(1+x).
// Ported from core-math/src/binary64/log2p1/log2p1.c (Paul Zimmermann,
// Tom Hubrecht).
//
// Reuses log_const.inc tables (cLogInverse, cLogInvH/L, cLogInverse2,
// cLogLogInv2, cLogP2, cLogLog2Dint, cLogMOneDint, cLogLog2h/l). LOG2_INV
// (dint, 2^12/log(2)) is local to log2p1 (cLog2p1Log2InvDint).
//
// Algorithm:
//   * fast path (cr_log2p1_fast):
//     - |x| < 2^-5 (e <= -6): degree-11 polynomial p_1a in (x, 0), then ×
//       1/log(2) via local DMul matching log2p1.c d_mul order. Defers to
//       accurate path when e <= -920.
//     - else: form (xh,xl)=1+x, run log10p1.c-style cr_log_fast (degree-7
//       P[7] returning (ph,pl)), correct by log(1 + xl/xh) ≈ xl/xh, then
//       × 1/log(2).
//   * accurate path (cr_log2p1_accurate):
//     - |x| < 2^-105: scaled fma using INVLOG2H/L (with a single explicit
//       subnormal-input exception).
//     - 2^-105 ≤ |x| < 2^-5: Pacc[24] degree-17 polynomial (mixed dd/single)
//       plus 247-entry exception table.
//     - else: dint reduction (LogTwo) + correction by log(1+xl/xh) using
//       2-term Taylor on dint, then mul by LOG2_INV (with Y.ex -= 12).
// ===========================================================================

const
{$I log2p1_const_64.inc}

// Local d_mul matching log2p1.c order (a_mul + fma(ah,bl,lo) + fma(al,bh,lo)).
procedure Log2p1DMul(var hi, lo: Double; ah, al, bh, bl: Double); inline;
var t1, t2: Double;
begin
  pcr_a_mul(hi, t1, ah, bh);
  t2 := pcr_fma(ah, bl, t1);
  lo := pcr_fma(al, bh, t2);
end;

procedure Log2p1SMul(var hi, lo: Double; a, bh, bl: Double); inline;
begin
  pcr_s_mul(hi, lo, a, bh, bl);
end;

// p_1: |z| < 0.00212097..., (h,l) = log(1+z) - z to 2^-78.25 absolute.
procedure Log2p1P1(out h, l: Double; z: Double); inline;
var
  z2h, z2l, p56, p34, ph: Double;
begin
  pcr_a_mul(z2h, z2l, z, z);
  p56 := pcr_fma(cLog2p1P[6].f, z, cLog2p1P[5].f);
  p34 := pcr_fma(cLog2p1P[4].f, z, cLog2p1P[3].f);
  ph  := pcr_fma(p56, z2h, p34);
  ph  := pcr_fma(ph, z, cLog2p1P[2].f);
  ph  := ph * z2h;
  pcr_fasttwosum(h, l, Double(-0.5) * z2h, ph * z);
  l := l + Double(-0.5) * z2l;
end;

// p_1a: |z| < 0.03125, (h,l) = log(1+z) - z to 2^-67.14 absolute.
procedure Log2p1P1a(out h, l: Double; z: Double); inline;
var
  z2h, z2l, z4h: Double;
  p910, p78, p56, p34, p710, p36, ph: Double;
begin
  if pcr_fabs(z) >= Double(2.7105054312137611e-77) then  // 0x1p-255
    pcr_a_mul(z2h, z2l, z, z)
  else begin
    z2h := 0; z2l := 0;
  end;
  z4h := z2h * z2h;
  p910 := pcr_fma(cLog2p1Pa[10].f, z, cLog2p1Pa[9].f);
  p78  := pcr_fma(cLog2p1Pa[8].f,  z, cLog2p1Pa[7].f);
  p56  := pcr_fma(cLog2p1Pa[6].f,  z, cLog2p1Pa[5].f);
  p34  := pcr_fma(cLog2p1Pa[4].f,  z, cLog2p1Pa[3].f);
  p710 := pcr_fma(p910, z2h, p78);
  p36  := pcr_fma(p56,  z2h, p34);
  ph   := pcr_fma(p710, z4h, p36);
  ph   := pcr_fma(ph,   z,   cLog2p1Pa[2].f);
  ph   := ph * z2h;
  pcr_fasttwosum(h, l, Double(-0.5) * z2h, ph * z);
  l := l + Double(-0.5) * z2l;
end;

// log2p1's cr_log_fast: (h,l) ≈ log(2^e * v.f), uses degree-7 P[7] via Log2p1P1.
procedure Log2p1LogFast(var h, l: Double; e: Int32; vf: Double);
var
  v: Tb64u64;
  m: UInt64;
  c, i: Int32;
  y, r, l1, l2, z, ph, pl: Double;
  ee, hh, et, lo: Double;
begin
  v.f := vf;
  m := UInt64($10000000000000) + (v.u and UInt64($FFFFFFFFFFFFF));

  if m >= UInt64($16A09E667F3BCD) then begin
    c := 1; e := e + 1;
    y := vf * Double(0.5);
    i := Int32(m shr 44);
  end else begin
    c := 0;
    y := vf;
    i := Int32(m shr 43);
  end;
  if c = 0 then ;

  r  := cLogInverse[i - cLogOFFSET].f;
  l1 := cLogInvH   [i - cLogOFFSET].f;
  l2 := cLogInvL   [i - cLogOFFSET].f;
  z  := pcr_fma(r, y, Double(-1.0));

  Log2p1P1(ph, pl, z);

  ee := e;
  hh := pcr_fma(ee, cLogLog2h.f, l1);
  // fast_two_sum(h, l, hh, z)
  h  := hh + z;
  et := h - hh;
  l  := z - et;

  // l = ph + (l + (l2 + pl))
  lo := l2 + pl;
  l  := l + lo;
  l  := ph + l;

  l := pcr_fma(ee, cLogLog2l.f, l);
end;

// Accurate-tiny: |x| < 2^-105 ⇒ log2p1(x) ~ x/log(2).
function Log2p1AccTiny(x: Double): Double;
var
  sx, h, l, res, ax: Double;
begin
  ax := pcr_fabs(x);
  // Exception: |x| = 0x0.2c316a14459d8p-1022 (subnormal)
  if ax = cLog2p1AccTinyExc.f then begin
    if x > 0 then
      Result := pcr_fma(Double(2.4099198651028841e-181),
                        Double(2.4099198651028841e-181),
                        cLog2p1AccTinyExcB.f)              // 2^-600 * 2^-600 + b
    else
      Result := pcr_fma(Double(-2.4099198651028841e-181),
                        Double(2.4099198651028841e-181),
                        -cLog2p1AccTinyExcB.f);
    Exit;
  end;
  sx := x * Double(8.1129638414606682e+31);  // 2^106
  Log2p1SMul(h, l, sx, cLog2p1InvLog2H.f, cLog2p1InvLog2L.f);
  res := (h + l) * Double(1.2325951644078309e-32);  // 2^-106
  l := pcr_fma(-res, Double(8.1129638414606682e+31), h) + l;
  Result := pcr_fma(l, Double(1.2325951644078309e-32), res);
end;

// Accurate-small: 2^-105 ≤ |x| < 2^-5, polynomial Pacc[24] + 247-exception table.
function Log2p1AccSmall(x: Double): Double;
var
  a, b, c: Int32;
  h, l, t: Double;
  del: ShortInt;
  i: Int32;
begin
  // Binary search by exceptions[a][0] <= x < exceptions[b][0]
  a := 0; b := 247;
  while a + 1 <> b do begin
    c := (a + b) div 2;
    if cLog2p1ExX[c].f <= x then a := c else b := c;
  end;
  if x = cLog2p1ExX[a].f then begin
    h := cLog2p1ExH[a].f;
    if h > 0 then del :=  cLog2p1ExRnd[a]
    else          del := -cLog2p1ExRnd[a];
    Result := h + h * Double(5.5511151231257827e-17) * Double(del);  // 0x1p-54
    Exit;
  end;

  // Pacc layout: idx 0..13 = 7 dd-pairs (degrees 1..7); idx 14..23 = 10 singles
  // (degrees 8..17). C indexing: Pacc[2*i-2],Pacc[2*i-1] for degree i in 1..7;
  // Pacc[i+6] for degree i in 8..17 (so degrees 8..17 ↦ idx 14..23).
  //
  // Single-double Horner over degrees 16..11
  h := pcr_fma(cLog2p1Pacc[23].f, x, cLog2p1Pacc[22].f);  // degree 16 (idx 22+x*23)
  for i := 15 downto 11 do
    h := pcr_fma(h, x, cLog2p1Pacc[i+6].f);
  l := 0;
  // degrees 10..8 single-double accumulated as dd
  for i := 10 downto 8 do begin
    Log2p1SMul(h, l, x, h, l);
    pcr_fasttwosum(h, t, cLog2p1Pacc[i+6].f, h);
    l := l + t;
  end;
  // degrees 7..1 dd polynomial
  for i := 7 downto 1 do begin
    Log2p1SMul(h, l, x, h, l);
    pcr_fasttwosum(h, t, cLog2p1Pacc[2*i-2].f, h);
    l := l + t + cLog2p1Pacc[2*i-1].f;
  end;
  // final * x
  Log2p1SMul(h, l, x, h, l);
  Result := h + l;
end;

// inv_dint: r ≈ 1/a where a is a Double.
procedure Log2p1InvDInt(out r: TDInt64; ad: Double);
var
  q, am: TDInt64;
begin
  DIntFromD(r, Double(4.0) / ad);
  r.ex := r.ex - 2;
  DIntFromD(am, -ad);
  MulDInt(q, am, r);
  AddDInt(q, DINT_ONE, q);
  MulDInt(q, r, q);
  AddDInt(r, r, q);
end;

procedure Log2p1DivDInt(out r: TDInt64; bd, ad: Double);
var
  bm: TDInt64;
begin
  Log2p1InvDInt(r, ad);
  DIntFromD(bm, bd);
  MulDInt(r, r, bm);
end;

function Log2p1Accurate(x: Double): Double;
var
  X64, Y64, C64, X2: TDInt64;
  ax, xh, xl: Double;
  t: Tb64u64;
  e: Int32;
begin
  ax := pcr_fabs(x);
  if ax < Double(0.03125) then begin  // 2^-5
    if ax < Double(2.4651903288156619e-32) then  // 2^-105
      Result := Log2p1AccTiny(x)
    else
      Result := Log2p1AccSmall(x);
    Exit;
  end;

  // (xh, xl) <- 1 + x
  if x > Double(1.0) then
    pcr_fasttwosum(xh, xl, x, Double(1.0))
  else
    pcr_fasttwosum(xh, xl, Double(1.0), x);

  // log2p1(x) is exact when 1+x = 2^e (xl = 0 ∧ xh power of 2)
  if xl = 0 then begin
    t.f := xh;
    if (t.u shl 12) = 0 then begin
      e := Int32(t.u shr 52) - $3FF;
      Result := Double(e);
      Exit;
    end;
  end;

  // x = 2^e exact: directed-rounding adjustment for e >= 49
  t.f := x;
  if (t.u shl 12) = 0 then begin
    e := Int32(t.u shr 52) - $3FF;
    if e >= 49 then begin
      Result := Double(e) + Double(3.5527136788005009e-15);  // 0x1p-48
      Exit;
    end;
  end;

  DIntFromD(X64, xh);
  LogTwo(Y64, X64);

  Log2p1DivDInt(C64, xl, xh);
  MulDInt(X2, C64, C64);
  X2.ex := X2.ex - 1;
  X2.sgn := 1;
  AddDInt(C64, C64, X2);
  AddDInt(Y64, Y64, C64);

  // multiply by 1/log(2): Y *= LOG2_INV (which approximates 2^12/log(2));
  // then divide by 2^12 via Y.ex -= 12.
  MulDInt(Y64, Y64, cLog2p1Log2InvDint);
  Y64.ex := Y64.ex - 12;

  Result := DToD(Y64);
end;

// Fast path
function Log2p1Fast(out h, l: Double; x: Double; e: Int32): Double;
var
  lo, xh, xl, c: Double;
  vv: Tb64u64;
  ee: Int32;
begin
  if e < -5 then begin
    if e <= -920 then begin
      h := 0; l := 0;
      Result := 1;
      Exit;
    end;
    Log2p1P1a(h, lo, x);
    pcr_fasttwosum(h, l, x, h);
    l := l + lo;
    // multiply by 1/log(2)
    Log2p1DMul(h, l, h, l, cLog2p1InvLog2H.f, cLog2p1InvLog2L.f);
    Result := cLog2p1ErrSmall.f;
    Exit;
  end;

  // (xh, xl) <- 1 + x
  if x > Double(1.0) then begin
    if x < Double(1.7976931348623155e+308) then  // 0x1.fffffffffffffp+1023
      pcr_fasttwosum(xh, xl, x, Double(1.0))
    else begin
      xh := x; xl := Double(1.0);
    end;
  end
  else
    pcr_fasttwosum(xh, xl, Double(1.0), x);

  vv.f := xh;
  ee := Int32(vv.u shr 52) - $3FF;
  vv.u := (UInt64($3FF) shl 52) or (vv.u and UInt64($FFFFFFFFFFFFF));
  Log2p1LogFast(h, l, ee, vv.f);

  // log(xh+xl) = log(xh) + log(1+xl/xh)
  if (xh <= Double(4.4942328371557898e+307))            // 2^1022
     or (pcr_fabs(xl) >= Double(4.0)) then
    c := xl / xh
  else
    c := 0;
  l := l + c;

  // multiply by 1/log(2)
  Log2p1DMul(h, l, h, l, cLog2p1InvLog2H.f, cLog2p1InvLog2L.f);

  Result := cLog2p1ErrMain.f;
  if e = 0 then ;
end;

function pcr_log2p1(x: Double): Double;
var
  v, w, minf: Tb64u64;
  e, k: Int32;
  h, l, err, left, right: Double;
begin
  v.f := x;
  e := Int32((v.u shr 52) and UInt64($7FF)) - $3FF;

  if (e = $400) or (x = 0) or (x <= Double(-1.0)) then begin
    minf.u := UInt64($FFF) shl 52;
    if (e = $400) and (x <> minf.f) then begin
      Result := x + x;  // NaN or +Inf
      Exit;
    end;
    if x <= Double(-1.0) then begin
      if x < Double(-1.0) then
        Result := cNaNDouble
      else
        Result := Double(1.0) / Double(-0.0);
      Exit;
    end;
    Result := x + x;  // ±0
    Exit;
  end;

  // x = 2^n - 1 exact (0 <= n <= 52): log2p1(x) = n+1, no inexact flag.
  if (e >= 0) and (e <= 52) then begin
    if x = cLog2p1T[e].f then begin
      Result := Double(e + 1);
      Exit;
    end;
  end;

  // x = 2^k - 1 exact for -53 <= k <= -1: -1 < x <= -1/2.
  if (e = -1) and (x < 0) then begin
    w.f := Double(1.0) + x;  // 1+x is exact
    if (w.u shl 12) = 0 then begin
      k := Int32(w.u shr 52) - $3FF;
      Result := Double(k);
      Exit;
    end;
  end;

  // Normalize v in [1, 2) — Log2p1Fast doesn't actually use it, but keep
  // shape parallel to log10p1 for clarity (param removed).
  err := Log2p1Fast(h, l, x, e);

  left  := h + (l - err);
  right := h + (l + err);
  if left = right then begin
    Result := left;
    Exit;
  end;

  Result := Log2p1Accurate(x);
end;
