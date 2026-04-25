// pcr_exp10m1 — correctly-rounded binary64 10^x - 1.
// Ported from core-math/src/binary64/exp10m1/exp10m1.c.
// Reuses E2M1AMul / E2M1SMul / E2M1DMul / E2M1FastSum / E2M1Q1 / E2M1Q2 /
// E2M1Exp1 from exp2m1_port.inc, plus cExpT0 (= T1) and cExpT1 (= T2).
// Q_1 / Q_2 in exp10m1.c are byte-identical to cExp2m1Q1 / cExp2m1Q2.

const
  cExp10m1LN10H:    Tb64u64 = (u:$40026BB1BBB55516);   // log(10) hi
  cExp10m1LN10L:    Tb64u64 = (u:$BCAF48AD494EA3E9);   // log(10) lo
  cExp10m1InvLog210: Tb64u64 = (u:$40CA934F0979A371);  // 2^12 * log(10)/log(2)
  cExp10m1Log210H:  Tb64u64 = (u:$3F134413509F79FF);   // log(2)/log(10)/2^12 hi
  cExp10m1Log210M:  Tb64u64 = (u:$BB89DC1DA9800000);   // mid (30-bit so M*k exact)
  cExp10m1Log210L:  Tb64u64 = (u:$B984FD20DBA1F655);   // lo
  cExp10m1C2:       Tb64u64 = (u:$40053524C73CEA69);   // log(10)^2 / 2

  // Bit-pattern thresholds.
  cExp10m1HugeNeg:    UInt64 = $C03041704C068EF0;      // x = -NaN or x <= -0x1.041704c068efp+4
  cExp10m1HugePosAx:  UInt64 = $40734413509F79FE;      // ax > this -> +NaN or overflow
  cExp10m1VTinyAx:    UInt64 = $3C90000000000000;      // |x| <= 2^-54
  cExp10m1UltraTiny:  UInt64 = $3970000000000000;      // |x| <= 2^-104
  cExp10m1TinyAx:     UInt64 = $3FB0000000000000;      // |x| <= 0.0625
  cExp10m1ScaleP106:  Tb64u64 = (u:$4690000000000000); // 2^106
  cExp10m1ScaleM106:  Tb64u64 = (u:$3950000000000000); // 2^-106
  cExp10m1Eps54:      Tb64u64 = (u:$3C90000000000000); // 2^-54
  cExp10m1OverflowH:  Tb64u64 = (u:$7FEFFFFFFFFFFFFF); // 0x1.fffffffffffffp+1023
  cExp10m1ErrTiny:    Tb64u64 = (u:$3BB0A00000000000); // 0x1.0ap-68
  cExp10m1ErrMain:    Tb64u64 = (u:$3B77A00000000000); // 0x1.7ap-72

  // K=4194304 clamp (max value of round(x*INVLOG2_10) in valid x range).
  cExp10m1Kmax:       Int64 = 4194304;

  // Subnormal exception magic constants (cr_exp10m1, |x| <= 2^-104 path).
  cExp10m1Sub1Mag:    Tb64u64 = (u:$000086C73059343C);   // |x| == 0x0.086c73059343cp-1022
  cExp10m1Sub1A:      Tb64u64 = (u:$1E60010000000000);   //  0x1.001p-537
  cExp10m1Sub1B:      Tb64u64 = (u:$1E50000000000000);   //  0x1p-538
  cExp10m1Sub1C:      Tb64u64 = (u:$000136568740CB56);   //  0x1.36568740cb56p-1026
  cExp10m1Sub2Mag:    Tb64u64 = (u:$00013A7B70D0248C);   // |x| == 0x0.13a7b70d0248cp-1022
  cExp10m1Sub2A:      Tb64u64 = (u:$1E5FFE0000000000);   //  0x0.fffp-537 = 0x1.ffep-538
  cExp10m1Sub2C:      Tb64u64 = (u:$0002D41F3B972FC7);   //  0x1.6a0f9dcb97e38p-1025 (subnormal)

{$I exp10m1_const.inc}

// Integer 10^k for k in 1..15 (exact in Double up to 15).
const
  cExp10m1IntT: array[0..15] of Double = (
    0.0,
    9.0, 99.0, 999.0, 9999.0,
    99999.0, 999999.0, 9999999.0, 99999999.0,
    999999999.0, 9999999999.0, 99999999999.0, 999999999999.0,
    9999999999999.0, 99999999999999.0, 999999999999999.0
  );

// Sign-copy: produce |a|.f with sign of x.
function E10M1CopySign(a, x: Double): Double; inline;
var ua, ux: Tb64u64;
begin
  ua.f := a; ux.f := x;
  ua.u := (ua.u and UInt64($7FFFFFFFFFFFFFFF))
        or (ux.u and UInt64($8000000000000000));
  Result := ua.f;
end;

// exp_2: accurate exp(x*log(10)) using INVLOG2_10 reduction. -16.255 < x < 308.255
procedure E10M1Exp2(out hi, lo: Double; x: Double);
var
  kd, yh, yl, eh, el, qh, ql, t1h, t1l, t2h, t2l, h, l: Double;
  Ki, M, i2, i1: Int64;
  d: Tb64u64;
begin
  kd := pcr_roundeven(x * cExp10m1InvLog210.f);
  if kd = 4194304.0 then kd := 4194303.0;

  yh := pcr_fma(-kd, cExp10m1Log210H.f, x);  // exact
  pcr_fasttwosum(yh, yl, yh, -kd * cExp10m1Log210M.f);
  yl := pcr_fma(-kd, cExp10m1Log210L.f, yl);
  // multiply yh+yl by log(10)
  E2M1SMul(eh, el, yh, cExp10m1LN10H.f, cExp10m1LN10L.f);
  yh := eh;
  yl := pcr_fma(yl, cExp10m1LN10H.f, el);

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

// Fast tiny path: |x| <= 0.0625, |x| > 2^-54. Returns absolute error bound.
function Exp10m1FastTiny(out h, l: Double; x: Double): Double;
var
  x2, x4, c9, c7, c5, t, hh, ll, e1, e2, saved: Double;
begin
  x2 := x * x; x4 := x2 * x2;
  c9 := pcr_fma(cExp10m1P[12].f, x, cExp10m1P[11].f);
  c7 := pcr_fma(cExp10m1P[10].f, x, cExp10m1P[9].f);
  c5 := pcr_fma(cExp10m1P[8].f,  x, cExp10m1P[7].f);
  c9 := pcr_fma(cExp10m1P[13].f, x2, c9);
  c5 := pcr_fma(c7, x2, c5);
  c5 := pcr_fma(c9, x4, c5);

  // a_mul (h, l, c5, x)
  hh := c5 * x;
  ll := pcr_fma(c5, x, -hh);
  // fast_two_sum (h, t, P[6], h)
  saved := hh;
  e1 := cExp10m1P[6].f + saved;
  t  := saved - (e1 - cExp10m1P[6].f);
  hh := e1;
  ll := ll + t;

  // s_mul (h, l, x, h, l) ; fast_two_sum (h, t, P[4], h); l += t + P[5]
  saved := hh;
  e1 := x * saved;
  e2 := pcr_fma(x, saved, -e1);
  ll := pcr_fma(x, ll, e2);
  hh := e1;
  saved := hh;
  e1 := cExp10m1P[4].f + saved;
  t  := saved - (e1 - cExp10m1P[4].f);
  hh := e1;
  ll := ll + t + cExp10m1P[5].f;

  // s_mul ... add P[2]+P[3]
  saved := hh;
  e1 := x * saved;
  e2 := pcr_fma(x, saved, -e1);
  ll := pcr_fma(x, ll, e2);
  hh := e1;
  saved := hh;
  e1 := cExp10m1P[2].f + saved;
  t  := saved - (e1 - cExp10m1P[2].f);
  hh := e1;
  ll := ll + t + cExp10m1P[3].f;

  // s_mul ... add P[0]+P[1]
  saved := hh;
  e1 := x * saved;
  e2 := pcr_fma(x, saved, -e1);
  ll := pcr_fma(x, ll, e2);
  hh := e1;
  saved := hh;
  e1 := cExp10m1P[0].f + saved;
  t  := saved - (e1 - cExp10m1P[0].f);
  hh := e1;
  ll := ll + t + cExp10m1P[1].f;

  // s_mul (h, l, x, h, l) — final multiplication by x
  saved := hh;
  e1 := x * saved;
  e2 := pcr_fma(x, saved, -e1);
  ll := pcr_fma(x, ll, e2);
  hh := e1;

  h := hh; l := ll;
  Result := cExp10m1ErrTiny.f * hh;
end;

// Fast path. tiny=true iff |x| <= 0.0625.
function Exp10m1Fast(out h, l: Double; x: Double; tiny: Boolean): Double;
var
  hh, ll, u: Double;
begin
  if tiny then
  begin
    Result := Exp10m1FastTiny(h, l, x);
    Exit;
  end;
  // Main path: -0x1.041704c068efp+4 < x < -0.0625 or 0.0625 < x <= 0x1.34413509f79fep+8
  E2M1SMul(hh, ll, x, cExp10m1LN10H.f, cExp10m1LN10L.f);
  E2M1Exp1(hh, ll, hh, ll);
  if x >= 0 then
    pcr_fasttwosum(hh, u, hh, Double(-1.0))
  else
    pcr_fasttwosum(hh, u, Double(-1.0), hh);
  ll := ll + u;
  h := hh; l := ll;
  Result := cExp10m1ErrMain.f * hh;
end;

function Exp10m1AccurateTiny(x: Double): Double;
var
  a, b, c: Int32;
  h, l, t, x2, x4, c15, c13, exc: Double;
  i: Int32;
begin
  a := 0; b := 145;
  while a + 1 <> b do
  begin
    c := (a + b) div 2;
    if cExp10m1ExAccTiny[c, 0].f <= x then a := c else b := c;
  end;
  if x = cExp10m1ExAccTiny[a, 0].f then
  begin
    Result := cExp10m1ExAccTiny[a, 1].f + cExp10m1ExAccTiny[a, 2].f;
    Exit;
  end;

  x2 := x * x; x4 := x2 * x2;
  c15 := pcr_fma(cExp10m1Q[23].f, x, cExp10m1Q[22].f);  // degree 15
  c13 := pcr_fma(cExp10m1Q[21].f, x, cExp10m1Q[20].f);  // degree 13
  c15 := pcr_fma(cExp10m1Q[24].f, x2, c15);
  // add Q[19]*x + c13*x^2 + c15*x^4 to Q[18] (degree 11)
  exc := cExp10m1Q[19].f * x + c13 * x2 + c15 * x4;
  h := cExp10m1Q[18].f + exc;
  l := exc - (h - cExp10m1Q[18].f);
  // multiply h+l by x and add Q[17] (degree 10)
  E2M1SMul(h, l, x, h, l);
  pcr_fasttwosum(h, t, cExp10m1Q[17].f, h);
  l := l + t;
  // multiply h+l by x and add Q[16] (degree 9)
  E2M1SMul(h, l, x, h, l);
  pcr_fasttwosum(h, t, cExp10m1Q[16].f, h);
  l := l + t;
  // degrees 8..1 are dd: pairs (Q[14],Q[15]) ... (Q[0],Q[1])
  for i := 7 downto 0 do
  begin
    E2M1SMul(h, l, x, h, l);
    pcr_fasttwosum(h, t, cExp10m1Q[i*2].f, h);
    l := l + t + cExp10m1Q[i*2 + 1].f;
  end;
  // multiply h+l by x
  E2M1SMul(h, l, x, h, l);
  Result := h + l;
end;

function Exp10m1Accurate(x: Double): Double;
var
  ix: Tb64u64;
  a, b, c: Int32;
  h, l, u: Double;
begin
  ix.f := x;
  if (ix.u and UInt64($7FFFFFFFFFFFFFFF)) <= UInt64($3FB0000000000000) then
  begin
    Result := Exp10m1AccurateTiny(x);
    Exit;
  end;
  a := 0; b := 81;
  while a + 1 <> b do
  begin
    c := (a + b) div 2;
    if cExp10m1ExAcc[c, 0].f <= x then a := c else b := c;
  end;
  if x = cExp10m1ExAcc[a, 0].f then
  begin
    Result := cExp10m1ExAcc[a, 1].f + cExp10m1ExAcc[a, 2].f;
    Exit;
  end;
  E10M1Exp2(h, l, x);
  if x >= 0 then
    pcr_fasttwosum(h, u, h, Double(-1.0))
  else
    pcr_fasttwosum(h, u, Double(-1.0), h);
  l := l + u;
  Result := h + l;
end;

function pcr_exp10m1(x: Double): Double;
var
  ix: Tb64u64;
  ux, ax: UInt64;
  h, l, h2, res, err, lb, ub, xs, xa: Double;
  a, b, c: Int32;
  ifloor: Int64;
begin
  ix.f := x; ux := ix.u;
  ax := ux and UInt64($7FFFFFFFFFFFFFFF);

  if ux >= cExp10m1HugeNeg then
  begin
    if (ux shr 52) = $FFF then
    begin
      if ux > UInt64($FFF0000000000000) then Result := x + x
      else Result := Double(-1.0);
      Exit;
    end;
    Result := Double(-1.0) + cExp10m1Eps54.f;
    Exit;
  end
  else if ax > cExp10m1HugePosAx then
  begin
    if (ux shr 52) = $7FF then begin Result := x + x; Exit; end;
    Result := cExp10m1OverflowH.f * x;
    Exit;
  end
  else if ax <= cExp10m1VTinyAx then
  begin
    if ax <= cExp10m1UltraTiny then
    begin
      if x = 0 then begin Result := x; Exit; end;
      // Subnormal exception 1: |x| == 0x0.086c73059343cp-1022
      xa := Abs(x);
      if xa = cExp10m1Sub1Mag.f then
      begin
        Result := pcr_fma(-E10M1CopySign(cExp10m1Sub1A.f, x), cExp10m1Sub1B.f,
                           E10M1CopySign(cExp10m1Sub1C.f, x));
        Exit;
      end;
      // Subnormal exception 2: |x| == 0x0.13a7b70d0248cp-1022
      if xa = cExp10m1Sub2Mag.f then
      begin
        Result := pcr_fma( E10M1CopySign(cExp10m1Sub2A.f, x), cExp10m1Sub1B.f,
                           E10M1CopySign(cExp10m1Sub2C.f, x));
        Exit;
      end;
      xs := x * cExp10m1ScaleP106.f;
      E2M1AMul(h, l, cExp10m1LN10H.f, xs);
      l  := pcr_fma(cExp10m1LN10L.f, xs, l);
      h2 := h + l;
      h2 := h2 * cExp10m1ScaleM106.f;
      h  := pcr_fma(-h2, cExp10m1ScaleP106.f, h);
      h  := h + l;
      res := pcr_fma(h, cExp10m1ScaleM106.f, h2);
      Result := res;
      Exit;
    end
    else
    begin
      a := 0; b := 74;
      while a + 1 <> b do
      begin
        c := (a + b) div 2;
        if cExp10m1ExVTiny[c, 0].f <= x then a := c else b := c;
      end;
      if x = cExp10m1ExVTiny[a, 0].f then
      begin
        Result := cExp10m1ExVTiny[a, 1].f + cExp10m1ExVTiny[a, 2].f;
        Exit;
      end;
      E2M1AMul(h, l, cExp10m1LN10H.f, x);
      l := pcr_fma(cExp10m1LN10L.f, x, l);
      l := pcr_fma(cExp10m1C2.f * x, x, l);
      Result := h + l;
      Exit;
    end;
  end;

  // Integer-x exact path: 10^k - 1 for k in 1..15.
  if (ux shl 15) = 0 then
  begin
    ifloor := Trunc(x);
    if (Double(ifloor) = x) and (ifloor >= 1) and (ifloor <= 15) then
    begin
      Result := cExp10m1IntT[ifloor];
      Exit;
    end;
  end;

  err := Exp10m1Fast(h, l, x, ax <= cExp10m1TinyAx);
  lb := h + (l - err);
  ub := h + (l + err);
  if lb = ub then begin Result := lb; Exit; end;
  Result := Exp10m1Accurate(x);
end;
