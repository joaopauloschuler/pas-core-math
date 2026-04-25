// pcr_exp2m1 — correctly-rounded binary64 2^x - 1.
// Ported from core-math/src/binary64/exp2m1/exp2m1.c.
// Reuses cExpT0 (= T1, 2^(i/64)) and cExpT1 (= T2, 2^(i/4096)) tables already
// present in pascoremath64.pas for exp/exp2/exp10/expm1.

const
  cExp2m1LN2H:    Tb64u64 = (u:$3FE62E42FEFA39EF);   // log(2) hi
  cExp2m1LN2L:    Tb64u64 = (u:$3C7ABC9E3B39803F);   // log(2) lo
  cExp2m1InvLog2: Tb64u64 = (u:$40B71547652B82FE);   // 2^12 / log(2)
  cExp2m1Log2H:   Tb64u64 = (u:$3F262E42FEFA39EF);   // log(2) * 2^-12 hi
  cExp2m1Log2L:   Tb64u64 = (u:$3BBABC9E3B39803F);   // log(2) * 2^-12 lo
  cExp2m1C2:      Tb64u64 = (u:$3FCEBFBDFF82C58F);   // log(2)^2 / 2

  cExp2m1HugeNeg:    UInt64 = QWord($C04B000000000000);     // x <= -54 or -NaN/-Inf
  cExp2m1HugePos:    UInt64 = $4090000000000000;     // |x| >= 1024
  cExp2m1VTinyAx:    UInt64 = $3CC0527DBD87E24D;     // |x| <= 0x1.0527...p-51
  cExp2m1UltraTiny:  UInt64 = $3970000000000000;     // |x| <= 2^-104
  cExp2m1TinyAx:     UInt64 = $3FC0000000000000;     // |x| <= 0.125
  cExp2m1ScaleP12:   Tb64u64 = (u:$40B0000000000000); // 2^12  = 4096
  cExp2m1ScaleM12:   Tb64u64 = (u:$3F30000000000000); // 2^-12
  cExp2m1ScaleP106:  Tb64u64 = (u:$4690000000000000); // 2^106
  cExp2m1ScaleM106:  Tb64u64 = (u:$3950000000000000); // 2^-106
  cExp2m1Eps54:      Tb64u64 = (u:$3C90000000000000); // 2^-54
  cExp2m1OverflowH:  Tb64u64 = (u:$7FEFFFFFFFFFFFFF); // 0x1.fffffffffffffp+1023
  cExp2m1OverflowL:  Tb64u64 = (u:$7C0FFFFFFFFFFFFF); // 0x1.fffffffffffffp+960
  cExp2m1ErrTiny:    Tb64u64 = (u:$3BD4E00000000000); // 0x1.4ep-66
  cExp2m1ErrMain:    Tb64u64 = (u:$3B84200000000000); // 0x1.42p-71

{$I exp2m1_const.inc}

// --- Local primitives matching exp2m1.c (no extra normalization) ---

procedure E2M1AMul(out h, l: Double; a, b: Double);
var hi: Double;
begin
  hi := a * b;
  l  := pcr_fma(a, b, -hi);
  h  := hi;
end;

procedure E2M1SMul(out h, l: Double; a, bh, bl: Double);
var hi, lo: Double;
begin
  hi := a * bh;
  lo := pcr_fma(a, bh, -hi);
  l  := pcr_fma(a, bl, lo);
  h  := hi;
end;

procedure E2M1DMul(out h, l: Double; ah, al, bh, bl: Double);
var lo, hi: Double;
begin
  hi := ah * bh;
  lo := pcr_fma(ah, bh, -hi);
  lo := pcr_fma(ah, bl, lo);
  l  := pcr_fma(al, bh, lo);
  h  := hi;
end;

procedure E2M1FastSum(var h, l: Double; a, bh, bl: Double); inline;
var s, e: Double;
begin
  s := a + bh;
  e := bh - (s - a);
  h := s;
  l := e + bl;
end;

// q_1: degree-4 fast polynomial for exp(z), |z| < 2^-12.88
procedure E2M1Q1(out hi, lo: Double; zh, zl: Double);
var
  z, q, h, l, qz: Double;
begin
  z := zh + zl;
  q := pcr_fma(cExp2m1Q1[4].f, zh, cExp2m1Q1[3].f);
  q := pcr_fma(q, z, cExp2m1Q1[2].f);
  // fast_two_sum (h, l, Q_1[1], q*z)
  qz := q * z;
  h := cExp2m1Q1[1].f + qz;
  l := qz - (h - cExp2m1Q1[1].f);
  // d_mul (h, l, zh, zl, h, l)
  E2M1DMul(h, l, zh, zl, h, l);
  // fast_sum (h, l, Q_1[0], h, l)
  E2M1FastSum(h, l, cExp2m1Q1[0].f, h, l);
  hi := h; lo := l;
end;

// q_2: degree-7 accurate polynomial for exp(z)
procedure E2M1Q2(out hi, lo: Double; zh, zl: Double);
var
  z, q, h, l, t: Double;
begin
  z := zh + zl;
  q := pcr_fma(cExp2m1Q2[8].f, zh, cExp2m1Q2[7].f);
  q := pcr_fma(q, z, cExp2m1Q2[6].f);
  q := pcr_fma(q, z, cExp2m1Q2[5].f);
  E2M1AMul(h, l, q, z);
  pcr_fasttwosum(h, t, cExp2m1Q2[3].f, h);
  l := l + t + cExp2m1Q2[4].f;
  E2M1DMul(h, l, h, l, zh, zl);
  pcr_fasttwosum(h, t, cExp2m1Q2[2].f, h);
  l := l + t;
  E2M1DMul(h, l, h, l, zh, zl);
  pcr_fasttwosum(h, t, cExp2m1Q2[1].f, h);
  l := l + t;
  E2M1DMul(h, l, h, l, zh, zl);
  pcr_fasttwosum(h, t, cExp2m1Q2[0].f, h);
  l := l + t;
  hi := h; lo := l;
end;

// exp_1: fast exp(xh+xl). |xh| < 1024.
procedure E2M1Exp1(var hi, lo: Double; xh, xl: Double);
var
  kd, kh, kl, yh, yl, qh, ql, t1h, t1l, t2h, t2l, h, l: Double;
  Ki, M, i2, i1: Int64;
  d: Tb64u64;
begin
  kd := pcr_roundeven(xh * cExp2m1InvLog2.f);
  E2M1SMul(kh, kl, kd, cExp2m1Log2H.f, cExp2m1Log2L.f);
  pcr_fasttwosum(yh, yl, xh - kh, xl);
  yl := yl - kl;

  Ki := Trunc(kd);
  M  := SarInt64(Ki, 12) + Int64($3FF);
  i2 := SarInt64(Ki, 6) and $3F;
  i1 := Ki and $3F;

  t1h := cExpT0[i2, 1].f; t1l := cExpT0[i2, 0].f;
  t2h := cExpT1[i1, 1].f; t2l := cExpT1[i1, 0].f;
  E2M1DMul(h, l, t2h, t2l, t1h, t1l);

  E2M1Q1(qh, ql, yh, yl);
  E2M1DMul(h, l, h, l, qh, ql);

  d.u := UInt64(M) shl 52;
  if M = $7FF then
  begin
    h := h * Double(2.0);
    l := l * Double(2.0);
    d.u := UInt64(M - 1) shl 52;
  end;
  hi := h * d.f;
  lo := l * d.f;
end;

// exp_2: accurate 2^x via exp(x*log(2)) with |x| < 745
procedure E2M1Exp2(var hi, lo: Double; x: Double);
var
  kd, yh, yl, qh, ql, t1h, t1l, t2h, t2l, h, l: Double;
  Ki, M, i2, i1: Int64;
  d: Tb64u64;
begin
  kd := pcr_roundeven(x * cExp2m1ScaleP12.f);
  yh := pcr_fma(-kd, cExp2m1ScaleM12.f, x);
  E2M1SMul(yh, yl, yh, cExp2m1LN2H.f, cExp2m1LN2L.f);

  Ki := Trunc(kd);
  M  := SarInt64(Ki, 12) + Int64($3FF);
  i2 := SarInt64(Ki, 6) and $3F;
  i1 := Ki and $3F;
  t1h := cExpT0[i2, 1].f; t1l := cExpT0[i2, 0].f;
  t2h := cExpT1[i1, 1].f; t2l := cExpT1[i1, 0].f;
  E2M1DMul(h, l, t2h, t2l, t1h, t1l);

  E2M1Q2(qh, ql, yh, yl);
  E2M1DMul(h, l, h, l, qh, ql);

  d.u := UInt64(M) shl 52;
  if M = $7FF then
  begin
    h := h * Double(2.0);
    l := l * Double(2.0);
    d.u := UInt64(M - 1) shl 52;
  end;
  hi := h * d.f;
  lo := l * d.f;
end;

// Fast tiny path: |x| <= 0.125, returns absolute error bound.
function Exp2m1FastTiny(out h, l: Double; x: Double): Double;
var
  x2, x4, c8, c6, c4, t, hh, ll, e1, e2, saved: Double;
begin
  x2 := x * x; x4 := x2 * x2;
  c8 := pcr_fma(cExp2m1P[10].f, x, cExp2m1P[9].f);
  c6 := pcr_fma(cExp2m1P[8].f,  x, cExp2m1P[7].f);
  c4 := pcr_fma(cExp2m1P[6].f,  x, cExp2m1P[5].f);
  c8 := pcr_fma(cExp2m1P[11].f, x2, c8);
  c4 := pcr_fma(c6, x2, c4);
  c4 := pcr_fma(c8, x4, c4);

  // a_mul (h, l, c4, x)
  hh := c4 * x;
  ll := pcr_fma(c4, x, -hh);

  // fast_two_sum (h, t, P[4], h)
  saved := hh;
  e1 := cExp2m1P[4].f + saved;
  t := saved - (e1 - cExp2m1P[4].f);
  hh := e1;
  ll := ll + t;

  // s_mul (h, l, x, h, l)
  saved := hh;
  e1 := x * saved;
  e2 := pcr_fma(x, saved, -e1);
  ll := pcr_fma(x, ll, e2);
  hh := e1;

  saved := hh;
  e1 := cExp2m1P[2].f + saved;
  t := saved - (e1 - cExp2m1P[2].f);
  hh := e1;
  ll := ll + t + cExp2m1P[3].f;

  saved := hh;
  e1 := x * saved;
  e2 := pcr_fma(x, saved, -e1);
  ll := pcr_fma(x, ll, e2);
  hh := e1;

  saved := hh;
  e1 := cExp2m1P[0].f + saved;
  t := saved - (e1 - cExp2m1P[0].f);
  hh := e1;
  ll := ll + t + cExp2m1P[1].f;

  saved := hh;
  e1 := x * saved;
  e2 := pcr_fma(x, saved, -e1);
  ll := pcr_fma(x, ll, e2);
  hh := e1;

  h := hh; l := ll;
  Result := cExp2m1ErrTiny.f * hh;
end;

// Fast path. tiny=true iff |x| <= 0.125.
function Exp2m1Fast(out h, l: Double; x: Double; tiny: Boolean): Double;
var
  hh, ll, u: Double;
begin
  if tiny then
  begin
    Result := Exp2m1FastTiny(h, l, x);
    Exit;
  end;
  E2M1SMul(hh, ll, x, cExp2m1LN2H.f, cExp2m1LN2L.f);
  E2M1Exp1(hh, ll, hh, ll);
  if x >= 0 then
    pcr_fasttwosum(hh, u, hh, Double(-1.0))
  else
    pcr_fasttwosum(hh, u, Double(-1.0), hh);
  ll := ll + u;
  h := hh; l := ll;
  Result := cExp2m1ErrMain.f * hh;
end;

function Exp2m1AccurateTiny(x: Double): Double;
var
  a, b, c: Int32;
  h, l, t, x2, x4, c13, c11, exc: Double;
  i: Int32;
begin
  a := 0; b := 59;
  while a + 1 <> b do
  begin
    c := (a + b) div 2;
    if cExp2m1ExAccTiny[c, 0].f <= x then a := c else b := c;
  end;
  if x = cExp2m1ExAccTiny[a, 0].f then
  begin
    Result := cExp2m1ExAccTiny[a, 1].f + cExp2m1ExAccTiny[a, 2].f;
    Exit;
  end;

  x2 := x * x; x4 := x2 * x2;
  c13 := pcr_fma(cExp2m1Q[20].f, x, cExp2m1Q[19].f);
  c11 := pcr_fma(cExp2m1Q[18].f, x, cExp2m1Q[17].f);
  c13 := pcr_fma(cExp2m1Q[21].f, x2, c13);
  exc := cExp2m1Q[16].f * x + c11 * x2 + c13 * x4;
  h := cExp2m1Q[15].f + exc;
  l := exc - (h - cExp2m1Q[15].f);
  E2M1SMul(h, l, x, h, l);
  pcr_fasttwosum(h, t, cExp2m1Q[14].f, h);
  l := l + t;
  for i := 6 downto 0 do
  begin
    E2M1SMul(h, l, x, h, l);
    pcr_fasttwosum(h, t, cExp2m1Q[i*2].f, h);
    l := l + t + cExp2m1Q[i*2 + 1].f;
  end;
  E2M1SMul(h, l, x, h, l);
  Result := h + l;
end;

function Exp2m1Accurate(x: Double): Double;
var
  ix: Tb64u64;
  a, b, c: Int32;
  h, l, u: Double;
begin
  ix.f := x;
  if (ix.u and UInt64($7FFFFFFFFFFFFFFF)) <= UInt64($3FC0000000000000) then
  begin
    Result := Exp2m1AccurateTiny(x);
    Exit;
  end;
  a := 0; b := 93;
  while a + 1 <> b do
  begin
    c := (a + b) div 2;
    if cExp2m1ExAcc[c, 0].f <= x then a := c else b := c;
  end;
  if x = cExp2m1ExAcc[a, 0].f then
  begin
    Result := cExp2m1ExAcc[a, 1].f + cExp2m1ExAcc[a, 2].f;
    Exit;
  end;
  E2M1Exp2(h, l, x);
  if x >= 0 then
    pcr_fasttwosum(h, u, h, Double(-1.0))
  else
    pcr_fasttwosum(h, u, Double(-1.0), h);
  l := l + u;
  Result := h + l;
end;

function pcr_exp2m1(x: Double): Double;
var
  ix, dscale: Tb64u64;
  ux, ax: UInt64;
  h, l, h2, res, err, lb, ub, xs: Double;
  a, b, c: Int32;
  ifloor: Int64;
begin
  ix.f := x; ux := ix.u;
  ax := ux and UInt64($7FFFFFFFFFFFFFFF);

  if ux >= cExp2m1HugeNeg then
  begin
    if (ux shr 52) = $FFF then
    begin
      if ux > UInt64($FFF0000000000000) then Result := x + x
      else Result := Double(-1.0);
      Exit;
    end;
    Result := Double(-1.0) + cExp2m1Eps54.f;
    Exit;
  end
  else if ax >= cExp2m1HugePos then
  begin
    if (ux shr 52) = $7FF then begin Result := x + x; Exit; end;
    Result := cExp2m1OverflowH.f + x * cExp2m1OverflowL.f;
    Exit;
  end
  else if ax <= cExp2m1VTinyAx then
  begin
    if ax <= cExp2m1UltraTiny then
    begin
      if x = 0 then begin Result := x; Exit; end;
      xs := x * cExp2m1ScaleP106.f;
      E2M1AMul(h, l, cExp2m1LN2H.f, xs);
      l  := pcr_fma(cExp2m1LN2L.f, xs, l);
      h2 := h + l;
      h2 := h2 * cExp2m1ScaleM106.f;
      h  := pcr_fma(-h2, cExp2m1ScaleP106.f, h);
      h  := h + l;
      res := pcr_fma(h, cExp2m1ScaleM106.f, h2);
      Result := res;
      Exit;
    end
    else
    begin
      a := 0; b := 56;
      while a + 1 <> b do
      begin
        c := (a + b) div 2;
        if cExp2m1ExVTiny[c, 0].f <= x then a := c else b := c;
      end;
      if x = cExp2m1ExVTiny[a, 0].f then
      begin
        Result := cExp2m1ExVTiny[a, 1].f + cExp2m1ExVTiny[a, 2].f;
        Exit;
      end;
      E2M1AMul(h, l, cExp2m1LN2H.f, x);
      l := pcr_fma(cExp2m1LN2L.f, x, l);
      l := l + cExp2m1C2.f * x * x;
      Result := h + l;
      Exit;
    end;
  end;

  // 2^x - 1 exact for integer x in [-53, 53]
  if (ux shl 17) = 0 then
  begin
    ifloor := Trunc(x);
    if (Double(ifloor) = x) and (ifloor >= -53) and (ifloor <= 53) then
    begin
      if ifloor >= 0 then
      begin
        Result := Double(Int64((Int64(1) shl ifloor) - 1));
        Exit;
      end
      else
      begin
        dscale.u := UInt64(Int64($3FF) + ifloor) shl 52;
        Result := Double(-1.0) + dscale.f;
        Exit;
      end;
    end;
  end;

  err := Exp2m1Fast(h, l, x, ax <= cExp2m1TinyAx);
  lb := h + (l - err);
  ub := h + (l + err);
  if lb = ub then begin Result := lb; Exit; end;
  Result := Exp2m1Accurate(x);
end;
