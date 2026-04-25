// ===========================================================================
// pcr_erfc — correctly-rounded binary64 erfc(x).
// Ported from core-math/src/binary64/erfc/erfc.c (Sibidanov, Zimmermann,
// Hubrecht, Jeannerod).
//
// Reuse map:
//   - cExpT0 / cExpT1   == erfc.c T1 / T2     (2^(i/64), 2^(i/4096))
//   - cExp2m1Q1         == erfc.c Q_1         (degree-4 fast exp poly)
//   - E2M1Exp1          == erfc.c exp_1       (fast double-double exp)
//   - ErfFast / ErfAccurate (from erf_port.inc) — used for |x| <= ~2.9.
//
// Erfc-specific tables in erfc_const.inc.
// ===========================================================================

const
{$I erfc_const.inc}
  // exp_accurate constants
  cErfcInvLog2Acc: Tb64u64 = (u:$3FF71547652B82FE); // 0x1.71547652b82fep+0
  cErfcLog2Hacc:   Tb64u64 = (u:$3FE62E42FEFA39EF); // 0x1.62e42fefa39efp-1
  cErfcLog2Lacc:   Tb64u64 = (u:$3C7ABC9E3B398000); // 0x1.abc9e3b398p-56
  cErfcLog2Tiny:   Tb64u64 = (u:$3989F97B57A079A1); // 0x1.f97b57a079a19p-103

  // erfc_asympt_fast threshold[6]
  cErfcAFTh0: Tb64u64 = (u:$3FBD500000000000); // 0x1.d5p-4
  cErfcAFTh1: Tb64u64 = (u:$3FC59DA6CA291BA6); // 0x1.59da6ca291ba6p-3
  cErfcAFTh2: Tb64u64 = (u:$3FCBC00000000000); // 0x1.bcp-3
  cErfcAFTh3: Tb64u64 = (u:$3FD0C00000000000); // 0x1.0cp-2
  cErfcAFTh4: Tb64u64 = (u:$3FD3800000000000); // 0x1.38p-2
  cErfcAFTh5: Tb64u64 = (u:$3FD6300000000000); // 0x1.63p-2

  cErfcAFXLim: Tb64u64 = (u:$4039DB1BB14E15CA); // 0x1.9db1bb14e15cap+4
  cErfcAFHLim: Tb64u64 = (u:$044151B9A3FDD5C9); // 0x1.151b9a3fdd5c9p-955
  cErfcAFErrK: Tb64u64 = (u:$3BBD900000000000); // 0x1.d9p-68
  cErfcAFErrSub: Tb64u64 = (u:$0010000000000000); // 0x1p-1022

  // erfc_asympt_accurate threshold[10]
  cErfcAATh0: Tb64u64 = (u:$3FB4500000000000); // 0x1.45p-4
  cErfcAATh1: Tb64u64 = (u:$3FBE000000000000); // 0x1.e0p-4
  cErfcAATh2: Tb64u64 = (u:$3FC3F00000000000); // 0x1.3fp-3
  cErfcAATh3: Tb64u64 = (u:$3FC9500000000000); // 0x1.95p-3
  cErfcAATh4: Tb64u64 = (u:$3FCF500000000000); // 0x1.f5p-3
  cErfcAATh5: Tb64u64 = (u:$3FD3100000000000); // 0x1.31p-2
  cErfcAATh6: Tb64u64 = (u:$3FD7100000000000); // 0x1.71p-2
  cErfcAATh7: Tb64u64 = (u:$3FDBC00000000000); // 0x1.bcp-2
  cErfcAATh8: Tb64u64 = (u:$3FE0B00000000000); // 0x1.0bp-1
  cErfcAATh9: Tb64u64 = (u:$3FE3000000000000); // 0x1.3p-1

  // cr_erfc_fast / cr_erfc thresholds
  cErfcSterbenzTh: Tb64u64 = (u:$3FDE861FBB24C00A); // 0x1.e861fbb24c00ap-2
  cErfcAddErr_104: Tb64u64 = (u:$3990000000000000); // 0x1.4p-104
  cErfcAddErr_102: Tb64u64 = (u:$39A0000000000000); // 0x1.4p-102
  cErfcThr1Fast:   Tb64u64 = (u:$4001713786D9C7C0); // 0x1.713786d9c7c08p+1 (rounded; runtime const used)
  cErfcSatTinyBd:  Tb64u64 = (u:$3C8C5BF891B4EF6A); // 0x1.c5bf891b4ef6ap-54
  cErfcSatTinyPos: Tb64u64 = (u:$3C7C5BF891B4EF6A); // 0x1.c5bf891b4ef6ap-55
  cErfcP54:        Tb64u64 = (u:$3C90000000000000); //  0x1p-54
  cErfcSatTinyA:   Tb64u64 = (u:$0000000000000001); //  0x1p-1074
  cErfcSatTinyB:   Tb64u64 = (u:$3FD0000000000000); //  0.25

  // erfc_asympt_accurate subnormal-output exception
  cErfcSubExpX:    Tb64u64 = (u:$403A8F7BFBD15495); // 0x1.a8f7bfbd15495p+4
  cErfcSubExpA:    Tb64u64 = (u:$0000000000000001); // 0x1p-1074
  cErfcSubExpFmaB: Tb64u64 = (u:$BFD0000000000000); // -0.25
  cErfcSubExpC:    Tb64u64 = (u:$00099EF5883F656C); // 0x1.99ef5883f656cp-1024

  cErfcAccUpper:   Tb64u64 = (u:$400B59FFB450828C); // 0x1.b59ffb450828cp+0

// ----------------------------------------------------------------------------
// Local double-double primitives matching erfc.c semantics exactly
// (no final renormalize, unlike pcr_muldd / pcr_d_mul).
// ----------------------------------------------------------------------------
procedure ErfcAMul(out hi, lo: Double; a, b: Double); inline;
begin
  hi := a * b;
  lo := pcr_fma(a, b, -hi);
end;

procedure ErfcFastTwoSum(out hi, lo: Double; a, b: Double); inline;
var e: Double;
begin
  hi := a + b;
  e  := hi - a;
  lo := b - e;
end;

procedure ErfcTwoSum(out hi, lo: Double; a, b: Double); inline;
var aa, bb, da, db: Double;
begin
  hi := a + b;
  aa := hi - b;
  bb := hi - aa;
  da := a - aa;
  db := b - bb;
  lo := da + db;
end;

procedure ErfcDMul(out hi, lo: Double; ah, al, bh, bl: Double); inline;
var t: Double;
begin
  ErfcAMul(hi, t, ah, bh);
  t  := pcr_fma(ah, bl, t);
  lo := pcr_fma(al, bh, t);
end;

procedure ErfcSMul(out hi, lo: Double; a, bh, bl: Double); inline;
var t: Double;
begin
  ErfcAMul(hi, t, a, bh);
  lo := pcr_fma(a, bl, t);
end;

// Multiply x by 2^e where -2046 <= e <= 2046.  Splits into two factors when
// |e| > 1023 to handle subnormal/overflow gracefully.
function ErfcLdExp(x: Double; e: Int64): Double;
var u: Tb64u64; e1, e2: Int64;
begin
  if (e >= -1022) and (e <= 1023) then
  begin
    u.u := UInt64(e + 1023) shl 52;
    Result := x * u.f;
    Exit;
  end;
  if e > 0 then
  begin
    e1 := 1023; e2 := e - 1023;
  end
  else
  begin
    e1 := -1022; e2 := e + 1022;
  end;
  if e2 > 1023 then e2 := 1023;
  if e2 < -1022 then e2 := -1022;
  u.u := UInt64(e1 + 1023) shl 52;
  Result := x * u.f;
  u.u := UInt64(e2 + 1023) shl 52;
  Result := Result * u.f;
end;

// ----------------------------------------------------------------------------
// exp_accurate: 2^e * (h+l) ~= exp(xh+xl) for -742 <= xh+xl <= -2.92.
// cErfcE2 layout: indices 0..15 are 8 dd pairs (degrees 0..7); indices 16..27
// are 12 single doubles for degrees 8..19.  Erfc.c: E2[2i],E2[2i+1] for dd;
// E2[i+8] for single, where the +8 accounts for 16 dd-half slots.
// ----------------------------------------------------------------------------
procedure ErfcExpAccurate(out h, l: Double; out e: Int32; xh, xl: Double);
var
  k: Int32;
  yh, yl, th, tl: Double;
  i: Int32;
begin
  k := Int32(Trunc(pcr_roundeven(xh * cErfcInvLog2Acc.f)));
  yh := pcr_fma(-Double(k), cErfcLog2Hacc.f, xh);
  ErfcTwoSum(th, tl, -Double(k) * cErfcLog2Lacc.f, xl);
  ErfcFastTwoSum(yh, yl, yh, th);
  yl := pcr_fma(-Double(k), cErfcLog2Tiny.f, yl + tl);

  h := cErfcE2[19 + 8].f;
  for i := 18 downto 16 do
    h := pcr_fma(h, yh, cErfcE2[i + 8].f);
  // degree 15
  ErfcAMul(th, tl, h, yh);
  tl := pcr_fma(h, yl, tl);
  ErfcFastTwoSum(h, l, cErfcE2[15 + 8].f, th);
  l := l + tl;
  for i := 14 downto 8 do
  begin
    ErfcAMul(th, tl, h, yh);
    tl := pcr_fma(h, yl, tl);
    tl := pcr_fma(l, yh, tl);
    ErfcFastTwoSum(h, l, cErfcE2[i + 8].f, th);
    l := l + tl;
  end;
  for i := 7 downto 0 do
  begin
    ErfcAMul(th, tl, h, yh);
    tl := pcr_fma(h, yl, tl);
    tl := pcr_fma(l, yh, tl);
    ErfcFastTwoSum(h, l, cErfcE2[2*i].f, th);
    l := l + (tl + cErfcE2[2*i + 1].f);
  end;
  e := k;
end;

// ----------------------------------------------------------------------------
// erfc_asympt_fast: 0x1.713786d9c7c09p+1 < x < 0x1.b39dc41e48bfdp+4.
// Returns absolute error bound on h+l.
// ----------------------------------------------------------------------------
function ErfcAsymptFast(out h, l: Double; x: Double): Double;
var
  eh, el, uh, ul, yh, yl, zh, zl: Double;
  i, j: Int32;
  thresh: Double;
begin
  if x >= cErfcAFXLim.f then
  begin
    h := Double(0.0); l := Double(0.0);
    Result := Double(1.0); Exit;
  end;

  ErfcAMul(uh, ul, x, x);
  E2M1Exp1(eh, el, -uh, -ul);

  yh := Double(1.0) / x;
  yl := yh * pcr_fma(-x, yh, Double(1.0));

  i := 0;
  while True do
  begin
    case i of
      0: thresh := cErfcAFTh0.f;
      1: thresh := cErfcAFTh1.f;
      2: thresh := cErfcAFTh2.f;
      3: thresh := cErfcAFTh3.f;
      4: thresh := cErfcAFTh4.f;
    else thresh := cErfcAFTh5.f;
    end;
    if not (yh > thresh) then Break;
    Inc(i);
    if i >= 6 then Break;
  end;

  ErfcAMul(uh, ul, yh, yh);
  ul := pcr_fma(Double(2.0) * yh, yl, ul);

  zh := cErfcT[i, 12].f;
  zh := pcr_fma(zh, uh, cErfcT[i, 11].f);
  zh := pcr_fma(zh, uh, cErfcT[i, 10].f);

  ErfcSMul(h, l, zh, uh, ul);
  ErfcFastTwoSum(zh, zl, cErfcT[i, 9].f, h);
  zl := zl + l;

  j := 15;
  while j >= 3 do
  begin
    ErfcDMul(h, l, zh, zl, uh, ul);
    ErfcFastTwoSum(zh, zl, cErfcT[i, (j+1) shr 1].f, h);
    zl := zl + l;
    Dec(j, 2);
  end;
  // degree 1: dd p[0]+p[1]
  ErfcDMul(h, l, zh, zl, uh, ul);
  ErfcFastTwoSum(zh, zl, cErfcT[i, 0].f, h);
  zl := zl + (l + cErfcT[i, 1].f);

  ErfcDMul(uh, ul, zh, zl, yh, yl);
  ErfcDMul(h, l, uh, ul, eh, el);

  if h >= cErfcAFHLim.f then
    Result := cErfcAFErrK.f * h
  else
    Result := cErfcAFErrSub.f;
end;

// ----------------------------------------------------------------------------
// cr_erfc_fast: dispatch fast path. Returns absolute error bound.
// ----------------------------------------------------------------------------
function ErfcFastTop(out h, l: Double; x: Double): Double;
var
  err, t, hh, ll: Double;
const
  THRESHOLD1: Double = 2.9018370415802956; // 0x1.713786d9c7c09p+1
begin
  if x < Double(0.0) then
  begin
    err := ErfFast(hh, ll, -x);
    err := err * hh;
    ErfcFastTwoSum(h, t, Double(1.0), hh);
    l := t + ll;
    Result := err + cErfcAddErr_102.f;
    Exit;
  end;

  if x <= THRESHOLD1 then
  begin
    err := ErfFast(hh, ll, x);
    err := err * hh;
    ErfcFastTwoSum(h, t, Double(1.0), -hh);
    l := t - ll;
    if x >= cErfcSterbenzTh.f then
    begin
      Result := err; Exit;
    end;
    Result := err + cErfcAddErr_104.f;
    Exit;
  end;

  Result := ErfcAsymptFast(h, l, x);
end;

// ----------------------------------------------------------------------------
// erfc_asympt_accurate: 0x1.b59ffb450828cp+0 < x < 0x1.b39dc41e48bfdp+4.
// ----------------------------------------------------------------------------
function ErfcAsymptAccurate(x: Double): Double;
var
  i, j: Int32;
  e: Int32;
  eh, el, uh, ul, yh, yl, zh, zl, h, l, res, corr: Double;
  thresh: Double;
  deg: Int32;
begin
  for i := 0 to 21 do
    if x = cErfcAsymptEx[i, 0].f then
    begin
      Result := cErfcAsymptEx[i, 1].f + cErfcAsymptEx[i, 2].f;
      Exit;
    end;

  if x = cErfcSubExpX.f then
  begin
    Result := pcr_fma(cErfcSubExpA.f, cErfcSubExpFmaB.f, cErfcSubExpC.f);
    Exit;
  end;

  ErfcAMul(uh, ul, x, x);
  ErfcExpAccurate(eh, el, e, -uh, -ul);

  yh := Double(1.0) / x;
  yl := yh * pcr_fma(-x, yh, Double(1.0));

  i := 0;
  while True do
  begin
    case i of
      0: thresh := cErfcAATh0.f;
      1: thresh := cErfcAATh1.f;
      2: thresh := cErfcAATh2.f;
      3: thresh := cErfcAATh3.f;
      4: thresh := cErfcAATh4.f;
      5: thresh := cErfcAATh5.f;
      6: thresh := cErfcAATh6.f;
      7: thresh := cErfcAATh7.f;
      8: thresh := cErfcAATh8.f;
    else thresh := cErfcAATh9.f;
    end;
    if not (yh > thresh) then Break;
    Inc(i);
    if i >= 10 then Break;
  end;

  ErfcAMul(uh, ul, yh, yh);
  ul := pcr_fma(Double(2.0) * yh, yl, ul);

  // polynomial degree = 29 + 2*i, leading coeff at p[20+i]
  zh := cErfcTacc[i, 20 + i].f;
  zl := Double(0.0);
  deg := 27 + 2*i;
  j := deg;
  while j >= 13 do
  begin
    ErfcAMul(h, l, zh, uh);
    l := pcr_fma(zh, ul, l);
    l := pcr_fma(zl, uh, l);
    ErfcTwoSum(zh, zl, cErfcTacc[i, (j-1) div 2 + 6].f, h);
    zl := zl + l;
    Dec(j, 2);
  end;
  j := 11;
  while j >= 1 do
  begin
    ErfcAMul(h, l, zh, uh);
    l := pcr_fma(zh, ul, l);
    l := pcr_fma(zl, uh, l);
    ErfcTwoSum(zh, zl, cErfcTacc[i, j-1].f, h);
    zl := zl + (l + cErfcTacc[i, j].f);
    Dec(j, 2);
  end;

  ErfcAMul(uh, ul, zh, yh);
  ul := pcr_fma(zh, yl, ul);
  ul := pcr_fma(zl, yh, ul);

  ErfcFastTwoSum(uh, ul, uh, ul);
  ErfcAMul(h, l, uh, eh);
  l := pcr_fma(uh, el, l);
  l := pcr_fma(ul, eh, l);

  res := ErfcLdExp(h + l, e);
  if res < Double(2.2250738585072014e-308) then
  begin
    corr := h - ErfcLdExp(res, -e);
    corr := corr + l;
    res := res + ErfcLdExp(corr, e);
  end;
  Result := res;
end;

// ----------------------------------------------------------------------------
// cr_erfc_accurate
// ----------------------------------------------------------------------------
function ErfcAccurateTop(x: Double): Double;
var
  h, l, t: Double;
  i: Int32;
begin
  if x < Double(0.0) then
  begin
    for i := 0 to 16 do
      if x = cErfcAccNegEx[i, 0].f then
      begin
        Result := cErfcAccNegEx[i, 1].f + cErfcAccNegEx[i, 2].f; Exit;
      end;
    ErfAccurate(h, l, -x);
    ErfcFastTwoSum(h, t, Double(1.0), h);
    l := t + l;
    Result := h + l;
    Exit;
  end;
  if x <= cErfcAccUpper.f then
  begin
    for i := 0 to 28 do
      if x = cErfcAccPosEx[i, 0].f then
      begin
        Result := cErfcAccPosEx[i, 1].f + cErfcAccPosEx[i, 2].f; Exit;
      end;
    ErfAccurate(h, l, x);
    ErfcFastTwoSum(h, t, Double(1.0), -h);
    l := t - l;
    Result := h + l;
    Exit;
  end;
  Result := ErfcAsymptAccurate(x);
end;

// ----------------------------------------------------------------------------
// pcr_erfc — top level
// ----------------------------------------------------------------------------
function pcr_erfc(x: Double): Double;
var
  t: Tb64u64;
  at: UInt64;
  h, l, err, left, right: Double;
begin
  t.f := x;
  at := t.u and UInt64($7FFFFFFFFFFFFFFF);

  if t.u >= UInt64($8000000000000000) then
  begin
    if t.u >= UInt64($C017744F8F74E94B) then
    begin
      if t.u >= UInt64($FFF0000000000000) then
      begin
        if t.u = UInt64($FFF0000000000000) then begin Result := Double(2.0); Exit; end;
        Result := x + x; Exit;
      end;
      // Returns 2.0 - 2^-54 (rounds to 2 in RN; runtime expression matches C).
      Result := Double(2.0) - cErfcP54.f; Exit;
    end;
    if x >= -cErfcSatTinyBd.f then
    begin
      Result := pcr_fma(-x, cErfcP54.f, Double(1.0));
      Exit;
    end;
  end
  else
  begin
    if at >= UInt64($403B39DC41E48BFD) then
    begin
      if at >= UInt64($7FF0000000000000) then
      begin
        if at = UInt64($7FF0000000000000) then begin Result := Double(0.0); Exit; end;
        Result := x + x; Exit;
      end;
      Result := cErfcSatTinyA.f * cErfcSatTinyB.f; Exit;
    end;
    if x <= cErfcSatTinyPos.f then
    begin
      Result := pcr_fma(-x, cErfcP54.f, Double(1.0));
      Exit;
    end;
  end;

  err := ErfcFastTop(h, l, x);
  left  := h + (l - err);
  right := h + (l + err);
  if left = right then begin Result := left; Exit; end;

  Result := ErfcAccurateTop(x);
end;
