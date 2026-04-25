// ===========================================================================
// pcr_log10p1 — correctly-rounded binary64 log10(1+x).
// Ported from core-math/src/binary64/log10p1/log10p1.c (Paul Zimmermann,
// Tom Hubrecht).
//
// Reuses log_const.inc tables (cLogInverse, cLogInvH/L, cLogInverse2,
// cLogLogInv2, cLogP2, cLogLog2Dint, cLogMOneDint, cLogLog2h/l) — these are
// byte-identical between log.c and log10p1.c. LOG10_INV (dint) is reused
// from log10_port.inc as cLog10OneOverLog10Dint.
//
// Algorithm:
//   * fast path (cr_log10p1_fast):
//     - |x| < 2^-5: degree-7 polynomial p_1 in (x, 0), then × 1/log(10)
//       via local DMul matching log10p1.c d_mul order.
//     - else: form (xh,xl)=1+x, run cr_log_fast (separate from log.c's, uses
//       degree-7 P[7] returning a double-double (ph, pl)), correct by
//       log(1 + xl/xh) ≈ xl/xh, then × 1/log(10).
//   * accurate path (cr_log10p1_accurate):
//     - |x| < 2^-5: Pacc[25] degree-17 polynomial (mixed dd/single) plus
//       222-entry exception table.
//     - else: dint reduction (LogTwo) + correction by log(1+xl/xh) using
//       2-term Taylor on dint, then mul by LOG10_INV. 13-entry exact table
//       short-circuit.
// ===========================================================================

const
{$I log10p1_const.inc}

// ---------- Local d_mul / s_mul matching log10p1.c order --------------------
// log10p1.c d_mul: a_mul(hi,lo,ah,bh); lo = fma(ah,bl,lo); lo = fma(al,bh,lo)
// pcr_d_mul does ah*bl AFTER al*bh — different order. Use local to match C.
procedure Log10p1DMul(var hi, lo: Double; ah, al, bh, bl: Double); inline;
var t1, t2: Double;
begin
  pcr_a_mul(hi, t1, ah, bh);
  t2 := pcr_fma(ah, bl, t1);
  lo := pcr_fma(al, bh, t2);
end;

// log10p1.c s_mul: a_mul(hi,lo,a,bh); lo = fma(a,bl,lo). Same as pcr_s_mul.
procedure Log10p1SMul(var hi, lo: Double; a, bh, bl: Double); inline;
begin
  pcr_s_mul(hi, lo, a, bh, bl);
end;

// ---------- p_1 / p_1a (double-double polynomial in z) ----------------------
// p_1: |z| < 0.00212097..., returns (h,l) = log(1+z) - z to 2^-78.25 absolute.
procedure Log10p1P1(out h, l: Double; z: Double); inline;
var
  z2h, z2l, p56, p34, ph: Double;
begin
  pcr_a_mul(z2h, z2l, z, z);
  p56 := pcr_fma(cLog10p1P[6].f, z, cLog10p1P[5].f);
  p34 := pcr_fma(cLog10p1P[4].f, z, cLog10p1P[3].f);
  ph  := pcr_fma(p56, z2h, p34);
  ph  := pcr_fma(ph, z, cLog10p1P[2].f);
  ph  := ph * z2h;
  pcr_fasttwosum(h, l, Double(-0.5) * z2h, ph * z);
  l := l + Double(-0.5) * z2l;
end;

// p_1a: |z| < 0.03125, returns (h,l) = log(1+z) - z to 2^-67.14 absolute.
procedure Log10p1P1a(out h, l: Double; z: Double); inline;
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
  p910 := pcr_fma(cLog10p1Pa[10].f, z, cLog10p1Pa[9].f);
  p78  := pcr_fma(cLog10p1Pa[8].f,  z, cLog10p1Pa[7].f);
  p56  := pcr_fma(cLog10p1Pa[6].f,  z, cLog10p1Pa[5].f);
  p34  := pcr_fma(cLog10p1Pa[4].f,  z, cLog10p1Pa[3].f);
  p710 := pcr_fma(p910, z2h, p78);
  p36  := pcr_fma(p56,  z2h, p34);
  ph   := pcr_fma(p710, z4h, p36);
  ph   := pcr_fma(ph,   z,   cLog10p1Pa[2].f);
  ph   := ph * z2h;
  pcr_fasttwosum(h, l, Double(-0.5) * z2h, ph * z);
  l := l + Double(-0.5) * z2l;
end;

// ---------- log10p1's cr_log_fast (variant of log.c's cr_log_fast) ----------
// Computes (h,l) ≈ log(2^e * v.f) using degree-7 P[7] returning (ph,pl).
procedure Log10p1LogFast(var h, l: Double; e: Int32; vf: Double);
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
  if c = 0 then ;  // silence FPC

  r  := cLogInverse[i - cLogOFFSET].f;
  l1 := cLogInvH   [i - cLogOFFSET].f;
  l2 := cLogInvL   [i - cLogOFFSET].f;
  z  := pcr_fma(r, y, Double(-1.0));

  Log10p1P1(ph, pl, z);

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

// ---------- Accurate-tiny: |x| < 2^-900 -------------------------------------
function Log10p1AccTiny(x: Double): Double;
var
  sx, h, l, res: Double;
begin
  sx := x * Double(8.1129638414606682e+31);  // 2^106
  Log10p1SMul(h, l, sx, cLog10p1InvLog10H.f, cLog10p1InvLog10L.f);
  res := (h + l) * Double(1.2325951644078309e-32);  // 2^-106
  l := pcr_fma(-res, Double(8.1129638414606682e+31), h) + l;
  Result := pcr_fma(l, Double(1.2325951644078309e-32), res);
end;

// ---------- Accurate-small: 2^-900 <= |x| < 2^-5, polynomial Pacc -----------
function Log10p1AccSmall(x: Double): Double;
var
  a, b, c: Int32;
  h, l, t: Double;
  del: ShortInt;
  i: Int32;
begin
  // Binary search by exceptions[a][0] <= x < exceptions[b][0]
  a := 0; b := 222;
  while a + 1 <> b do begin
    c := (a + b) div 2;
    if cLog10p1ExX[c].f <= x then a := c else b := c;
  end;
  if x = cLog10p1ExX[a].f then begin
    h := cLog10p1ExH[a].f;
    if h > 0 then del :=  cLog10p1ExRnd[a]
    else          del := -cLog10p1ExRnd[a];
    Result := h + h * Double(5.5511151231257827e-17) * Double(del);  // 0x1p-54
    Exit;
  end;

  // degrees 16..10 single-double Horner
  h := pcr_fma(cLog10p1Pacc[24].f, x, cLog10p1Pacc[23].f);  // degree 16 at idx 24
  for i := 15 downto 10 do
    h := pcr_fma(h, x, cLog10p1Pacc[i+7].f);
  // degree 9
  pcr_a_mul(h, l, x, h);
  pcr_fasttwosum(h, t, cLog10p1Pacc[9+7].f, h);
  l := l + t;
  // degrees 8..1: dd polynomial
  for i := 8 downto 1 do begin
    Log10p1SMul(h, l, x, h, l);
    pcr_fasttwosum(h, t, cLog10p1Pacc[2*i-2].f, h);
    l := l + t + cLog10p1Pacc[2*i-1].f;
  end;
  Log10p1SMul(h, l, x, h, l);
  Result := h + l;
end;

// ---------- Accurate path: full dint -----------------------------------------
// inv_dint equivalent: r ≈ 1/a where a is a Double.
procedure Log10p1InvDInt(out r: TDInt64; ad: Double);
var
  q, am: TDInt64;
begin
  // r = 4/a, then divide by 4 (ex -= 2).
  DIntFromD(r, Double(4.0) / ad);
  r.ex := r.ex - 2;
  // Newton: r -> r + r*(1 - a*r)
  DIntFromD(am, -ad);
  MulDInt(q, am, r);           // -a*r
  AddDInt(q, DINT_ONE, q);     // 1 - a*r
  MulDInt(q, r, q);            // r*(1-a*r)
  AddDInt(r, r, q);
end;

// div_dint(r, b, a) = b / a where b and a are doubles.
procedure Log10p1DivDInt(out r: TDInt64; bd, ad: Double);
var
  bm: TDInt64;
begin
  Log10p1InvDInt(r, ad);
  DIntFromD(bm, bd);
  MulDInt(r, r, bm);
end;

function Log10p1Accurate(x: Double): Double;
var
  X64, Y64, C64, X2: TDInt64;
  ax, xh, xl: Double;
  i: Int32;
begin
  ax := pcr_fabs(x);
  if ax < Double(7.4583407312002067e-271) then begin  // 2^-900
    Result := Log10p1AccTiny(x); Exit;
  end;
  if ax < Double(0.03125) then begin  // 2^-5
    Result := Log10p1AccSmall(x); Exit;
  end;

  // 13-entry exact-result table
  for i := 0 to 12 do
    if x = cLog10p1Tacc[3*i].f then begin
      Result := cLog10p1Tacc[3*i+1].f + cLog10p1Tacc[3*i+2].f;
      Exit;
    end;

  // (xh, xl) <- 1 + x
  if x > Double(1.0) then
    pcr_fasttwosum(xh, xl, x, Double(1.0))
  else
    pcr_fasttwosum(xh, xl, Double(1.0), x);

  DIntFromD(X64, xh);
  LogTwo(Y64, X64);

  Log10p1DivDInt(C64, xl, xh);
  MulDInt(X2, C64, C64);
  // multiply X2 by -1/2
  X2.ex := X2.ex - 1;
  X2.sgn := 1;
  // C <- C - C^2/2
  AddDInt(C64, C64, X2);
  AddDInt(Y64, Y64, C64);

  MulDInt(Y64, Y64, cLog10OneOverLog10Dint);

  Result := DToD(Y64);
end;

// ---------- Fast path -------------------------------------------------------
function Log10p1Fast(out h, l: Double; x: Double; e: Int32; v: Tb64u64): Double;
var
  lo, xh, xl, c: Double;
  vv: Tb64u64;
  ee: Int32;
begin
  if e < -5 then begin  // |x| < 2^-5
    if e <= -918 then begin
      h := 0; l := 0;
      Result := 1;
      Exit;
    end;
    Log10p1P1a(h, lo, x);
    pcr_fasttwosum(h, l, x, h);
    l := l + lo;
    // multiply by 1/log(10)
    Log10p1DMul(h, l, h, l, cLog10p1InvLog10H.f, cLog10p1InvLog10L.f);
    Result := cLog10p1ErrSmall.f;
    Exit;
  end;

  // (xh, xl) <- 1 + x
  if x > Double(1.0) then begin
    if x < Double(9007199254740992.0) then  // 2^53
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
  Log10p1LogFast(h, l, ee, vv.f);

  if xh <= Double(4.4942328371557898e+307) then  // 2^1022
    c := xl / xh
  else
    c := 0;
  l := l + c;

  // multiply by 1/log(10)
  Log10p1DMul(h, l, h, l, cLog10p1InvLog10H.f, cLog10p1InvLog10L.f);

  Result := cLog10p1ErrMain.f;
  // Suppress unused-variable warnings (FPC).
  if e = 0 then ;
  if v.u = 0 then ;
end;

// ---------- Main entry ------------------------------------------------------
function pcr_log10p1(x: Double): Double;
var
  v: Tb64u64;
  e: Int32;
  h, l, err, left, right: Double;
  minf: Tb64u64;
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
      if x < Double(-1.0) then begin
        Result := cNaNDouble;
        Exit;
      end else begin
        Result := Double(1.0) / Double(-0.0);
        Exit;
      end;
    end;
    Result := x + x;  // ±0
    Exit;
  end;

  // x = 10^n - 1 exact for 1 <= n <= 15 (no inexact)
  if (e >= 3) and (e <= 49) then begin
    if x = cLog10p1T[e].f then begin
      Result := cLog10p1U[e];
      Exit;
    end;
  end;

  // Normalize v in [1, 2)
  v.u := (UInt64($3FF) shl 52) or (v.u and UInt64($FFFFFFFFFFFFF));
  err := Log10p1Fast(h, l, x, e, v);

  left  := h + (l - err);
  right := h + (l + err);
  if left = right then begin
    Result := left;
    Exit;
  end;

  Result := Log10p1Accurate(x);
end;
