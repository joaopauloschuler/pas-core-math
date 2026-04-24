// ===========================================================================
// pcr_asinh — correctly-rounded binary64 asinh(x).
// Ported from core-math/src/binary64/asinh/asinh.c by Alexei Sibidanov.
// Reuses main-branch tables (B, R1, R2, L1, L2, L2h, L2l, C) and refine
// tables (T1..T4, LL, RefineCh, RefineCl, L20/L21/L22, VOffset) from
// acosh_port.inc — asinh's C tables are byte-identical to acosh's here.
// Reuses AcoshAddDD (adddd) from acosh_port.inc and AtanhMuldd (muldd_acc
// with the C grouping ahhl + (alhh + ahlh)) from atanh_port.inc.
// ===========================================================================

const
  // as_asinh_zero polynomial ch[12][2]
  cAsinhZeroCh: array[0..11, 0..1] of Tb64u64 = (
    ((u:$BFC5555555555555),(u:$BC65555555555555)),
    ((u:$3FB3333333333333),(u:$3C499999999949DF)),
    ((u:$BFA6DB6DB6DB6DB7),(u:$3C32492496091B0C)),
    ((u:$3F9F1C71C71C71C7),(u:$3C1C71A35CFA0671)),
    ((u:$BF96E8BA2E8BA2E9),(u:$3C317F937248CF81)),
    ((u:$3F91C4EC4EC4EC4F),(u:$BC374E3C1DFD4C3D)),
    ((u:$BF8C999999999977),(u:$BC238E7A467ECC55)),
    ((u:$3F87A87878786C7E),(u:$3C2A83C7BACE55EB)),
    ((u:$BF83FDE50D764083),(u:$BC2D024DF7FA0542)),
    ((u:$3F812EF3CEAE4D12),(u:$BC2BA9C13DEB261F)),
    ((u:$BF7DF3BD104AA267),(u:$BC1546DA9BC5B32A)),
    ((u:$3F7A685FC5DE7A04),(u:$3C140D284A1D67F9)));

  cAsinhZeroCl: array[0..4] of Tb64u64 = (
    (u:$BF77828D553EC800), (u:$3F751712F7BEE368), (u:$BF72E6D98527BCC6),
    (u:$3F70095DA47B392C), (u:$BF63B92D6368192C));

  // Fast-path tiny / small / medium polynomials (for |x| < 0x1.bp-4)
  cAsinhTiny1: Tb64u64 = (u:$BFC5555555555555); // -0x1.5555555555555p-3

  cAsinhTiny2: array[0..1] of Tb64u64 = (
    (u:$BFC5555555555555), (u:$3FB3333327C57C60));

  cAsinhSmall: array[0..3] of Tb64u64 = (
    (u:$BFC5555555555555), (u:$3FB333333332F2FF),
    (u:$BFA6DB6D9A665159), (u:$3F9F186866D775F0));

  cAsinhMedCl: array[0..6] of Tb64u64 = (
    (u:$BFC5555555555555), (u:$3FB3333333333310), (u:$BFA6DB6DB6DA466C),
    (u:$3F9F1C71C2EA7BE4), (u:$BF96E8B651B09D72), (u:$3F91C309FC0E69C2),
    (u:$BF8BAB7833C1E000));

  // Database (35 entries): {|x|, |asinh(|x|)| high, |asinh(|x|)| low}
  cAsinhDb: array[0..34, 0..2] of Tb64u64 = (
    ((u:$3FD00F9476450863),(u:$3FCFCB35067F343C),(u:$3C60000000000000)),
    ((u:$3FD1F0A79315B287),(u:$3FD1B68AAE88FEBA),(u:$3C70000000000000)),
    ((u:$3FD2B9618FF7ACB7),(u:$3FD27781D9AA4E25),(u:$BC70000000000000)),
    ((u:$3FD389EF683F3AA7),(u:$3FD33F52DB6DF1AF),(u:$3C70000000000000)),
    ((u:$3FD3B07E0C779DDA),(u:$3FD364303E1AD8F6),(u:$3C70000000000000)),
    ((u:$3FD48441DF33B6D3),(u:$3FD42E385800F0A4),(u:$3C70000000000000)),
    ((u:$3FD687BD068C1C1E),(u:$3FD616CC75D49226),(u:$BC70000000000000)),
    ((u:$3FD8740C4453A056),(u:$3FD7E4F2AD132A1D),(u:$3C70000000000000)),
    ((u:$3FD891ACDA11167E),(u:$3FD8009D924A3FFD),(u:$3C70000000000000)),
    ((u:$3FDBAFC3479FC9CC),(u:$3FDAE3773250E7D2),(u:$3C70000000000000)),
    ((u:$3FDC59869F17B483),(u:$3FDB7EFA91915C95),(u:$3C70000000000000)),
    ((u:$3FDC8BE879787986),(u:$3FDBAD0485E0FE0A),(u:$BC70000000000000)),
    ((u:$3FDE73B46ABB01E1),(u:$3FDD68039861AB53),(u:$3C70000000000000)),
    ((u:$3FDED6236DA268B0),(u:$3FDDC0CB8F638126),(u:$3C70000000000000)),
    ((u:$3FDF399EBAFC1951),(u:$3FDE1A4F519FAB77),(u:$BC70000000000000)),
    ((u:$3FDF70975AB0D471),(u:$3FDE4BAE8BCD6EA6),(u:$3C70000000000000)),
    ((u:$3FDFBDD4A37760B7),(u:$3FDE90F16EB88C09),(u:$3C70000000000000)),
    ((u:$3FDFEE72EFB4BFDD),(u:$3FDEBC791A88BED8),(u:$3C70000000000000)),
    ((u:$3FE02339D6BDB741),(u:$3FDF0B2264E34555),(u:$3C70000000000000)),
    ((u:$3FE09E7C831B1A23),(u:$3FDFE694C3C89138),(u:$3C70000000000000)),
    ((u:$3FE16D32C862FC3B),(u:$3FE0A9C9334066DB),(u:$BC80000000000000)),
    ((u:$3FE857954132083D),(u:$3FE67425FE575C88),(u:$BC80000000000000)),
    ((u:$3FE8A5C3B60F7E11),(u:$3FE6B23AD4415A17),(u:$BC80000000000000)),
    ((u:$3FE9740EB419DD04),(u:$3FE754AB7535D47D),(u:$3C80000000000000)),
    ((u:$3FEA16D9CC06011A),(u:$3FE7D3755D851062),(u:$BC80000000000000)),
    ((u:$3FEBB635BE2213D1),(u:$3FE91167CAE3CFA9),(u:$3C80000000000000)),
    ((u:$3FED4B21EBF542F0),(u:$3FEA3FC7E4DD47D1),(u:$BC80000000000000)),
    ((u:$3FF7B8516FFD2406),(u:$3FF2F5D3B178914A),(u:$3C90000000000000)),
    ((u:$3FF9295B9116E2E2),(u:$3FF3BFFA88639760),(u:$3C90000000000000)),
    ((u:$3FFFEDC65E327140),(u:$3FF710F91E844F9B),(u:$3C90000000000000)),
    ((u:$40057E377B3F0B4B),(u:$3FFB6E2C73F41415),(u:$3C90000000000000)),
    ((u:$4026056B06A21918),(u:$4008C0A26D055288),(u:$3C90000000000000)),
    ((u:$403843E1B5E5979C),(u:$400F0F978201EB84),(u:$3CA0000000000000)),
    ((u:$409FEE8F69C4CD25),(u:$4020A19AEBB51E90),(u:$BCC0000000000000)),
    ((u:$4170FBC6C02B1C90),(u:$40316369CD53BB69),(u:$3CD0000000000000)));

  // Thresholds
  cAsinhBr_Tinier: UInt64 = $3E57137449123EF7; // 0x1.7137449123ef7p-26
  cAsinhBr_Tiny:   UInt64 = $3E5A000000000000; // 0x1.ap-26
  cAsinhBr_VTiny:  UInt64 = $3F30000000000000; // 0x1p-12
  cAsinhBr_Small:  UInt64 = $3F93000000000000; // 0x1.3p-6
  cAsinhBr_Med:    UInt64 = $3FBB000000000000; // 0x1.bp-4
  cAsinhBr_Half:   UInt64 = $3FD0000000000000; // 0x1p-2
  cAsinhBr_One:    UInt64 = $3FF0000000000000; // 0x1p+0
  cAsinhBr_P26:    UInt64 = $4190000000000000; // 0x1p+26
  cAsinhBr_P52:    UInt64 = $4330000000000000; // 0x1p+52
  cAsinhBr_Inf:    UInt64 = $7FF0000000000000;

  cAsinhFastEps:   Tb64u64 = (u:$3CA7900000000000); // 0x1.79p-53
  cAsinhTinyScale: Tb64u64 = (u:$3C30000000000000); // 0x1p-60
  cAsinhMainEps:   Double  = 1.63e-19;
  cAsinhInvLog2:   Tb64u64 = (u:$3FF71547652B82FE); // 0x1.71547652b82fep+0

// Database lookup (sign applied via |x|)
function AsinhDatabase(x, f: Double): Double;
var
  ax, sgn: Double;
  a, b, m: Int32;
begin
  ax := pcr_fabs(x);
  a := 0; b := High(cAsinhDb);
  m := (a + b) div 2;
  while a <= b do
  begin
    if cAsinhDb[m, 0].f < ax then a := m + 1
    else if cAsinhDb[m, 0].f = ax then
    begin
      sgn := pcr_copysign(Double(1.0), x);
      f := sgn * cAsinhDb[m, 1].f + sgn * cAsinhDb[m, 2].f;
      Break;
    end
    else b := m - 1;
    m := (a + b) div 2;
  end;
  Result := f;
end;

// as_asinh_zero: slow path for |x| < 0x1p-2
// Computes asinh(x) = x + x^3 * P(x^2) via seeded polydd in x^2.
function AsinhZero(x, x2h, x2l: Double): Double;
var
  y0, y1, y2: Double;
  ch, cl, th, tl, hh: Double;
  t, w: Tb64u64;
  i: Int32;
begin
  y2 := x2h * (cAsinhZeroCl[0].f + x2h * (cAsinhZeroCl[1].f
      + x2h * (cAsinhZeroCl[2].f + x2h * (cAsinhZeroCl[3].f
      + x2h * cAsinhZeroCl[4].f))));

  // Seeded polydd: n=12, xh=x2h, xl=x2l, seed *l = y2
  ch := cAsinhZeroCh[11, 0].f + y2;
  cl := ((cAsinhZeroCh[11, 0].f - ch) + y2) + cAsinhZeroCh[11, 1].f;
  for i := 10 downto 0 do
  begin
    ch := AtanhMuldd(x2h, x2l, ch, cl, cl);
    th := ch + cAsinhZeroCh[i, 0].f;
    tl := (cAsinhZeroCh[i, 0].f - th) + ch;
    ch := th;
    cl := cl + tl + cAsinhZeroCh[i, 1].f;
  end;
  y2 := cl;
  y1 := ch;

  // y1 = muldd_acc(y1, y2, x2h, x2l, &y2)
  y1 := AtanhMuldd(y1, y2, x2h, x2l, y2);

  // y1 = mulddd(y1, y2, x, &y2) — no cl term; grouping: ahhl + ahlh
  hh := y1 * x;
  y2 := pcr_fma(x, y1, -hh) + y2 * x;
  y1 := hh;

  pcr_fasttwosum(y0, y1, x, y1);
  pcr_fasttwosum(y1, y2, y1, y2);

  t.f := y1;
  if (t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) = 0 then
  begin
    w.f := y2;
    if ((w.u xor t.u) shr 63) <> 0 then t.u := t.u - 1
    else t.u := t.u + 1;
    y1 := t.f;
  end;
  Result := y0 + y1;
end;

// as_asinh_refine: slow path for 0x1p-2 <= |x|, when main-branch Ziv fails.
// zh/zl is ah/al from caller (~ sqrt(x^2+1)+|x| in double-double).
// a ≈ log2 * |lb| drives the 4-level linear interpolation (i1..i4) into LL[][].
function AsinhRefine(x, zh, zl, a: Double): Double;
var
  t, v, t0, w: Tb64u64;
  ex, e, i1, i2, i3, i4: Int32;
  iv: UInt64;
  ed, el0, el1, el2: Double;
  L0, L1, L2: Double;
  t12, t34, th, tl, dh, dl, shm, slm: Double;
  xh, xl, xh2, sh2, sl2, sh_: Double;
  v0, v1, v2, vtmp: Double;
  er, de: UInt64;
  res, sgn: Double;
begin
  t.f := zh;
  ex := (t.u shr 52) and $7FF;
  if zl = Double(0.0) then
    e := ex - $3FF + 1
  else
    e := ex - $3FF;
  t.u := t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12);
  t.u := t.u or (UInt64($3FF) shl 52);
  ed := e;

  v.f := a - ed + cAcoshVOffset.f;
  iv := (v.u - (UInt64($3FF) shl 52)) shr (52 - 16);
  i1 := Int32((iv shr 12) and $1F);
  i2 := Int32((iv shr  8) and $F);
  i3 := Int32((iv shr  4) and $F);
  i4 := Int32(iv and $F);

  el2 := cAcoshL22.f * ed;
  el1 := cAcoshL21.f * ed;
  el0 := cAcoshL20.f * ed;

  L0 := cAcoshLL[0, i1, 0].f + cAcoshLL[1, i2, 0].f
      + (cAcoshLL[2, i3, 0].f + cAcoshLL[3, i4, 0].f);
  L1 := cAcoshLL[0, i1, 1].f + cAcoshLL[1, i2, 1].f
      + (cAcoshLL[2, i3, 1].f + cAcoshLL[3, i4, 1].f);
  L2 := cAcoshLL[0, i1, 2].f + cAcoshLL[1, i2, 2].f
      + (cAcoshLL[2, i3, 2].f + cAcoshLL[3, i4, 2].f);
  L0 := L0 + el0;

  t12 := cAcoshT1[i1].f * cAcoshT2[i2].f;
  t34 := cAcoshT3[i3].f * cAcoshT4[i4].f;
  th  := t12 * t34;
  tl  := pcr_fma(t12, t34, -th);
  dh  := th * t.f;
  dl  := pcr_fma(th, t.f, -dh);
  shm := tl * t.f;
  slm := pcr_fma(tl, t.f, -shm);

  pcr_fasttwosum(xh, xl, dh - Double(1.0), dl);
  if zl <> Double(0.0) then
  begin
    t.f := zl;
    t.u := UInt64(Int64(t.u) - (Int64(e) shl 52));
    xl := xl + th * t.f;
  end;
  xh := AcoshAddDD(xh, xl, shm, slm, xl);

  sl2 := xh * (cAcoshRefineCl[0].f + xh * (cAcoshRefineCl[1].f + xh * cAcoshRefineCl[2].f));
  xh2 := cAcoshRefineCh[2, 0].f + sl2;
  sh2 := ((cAcoshRefineCh[2, 0].f - xh2) + sl2) + cAcoshRefineCh[2, 1].f;
  for i1 := 1 downto 0 do
  begin
    xh2 := pcr_muldd(xh, xl, xh2, sh2, sh2);
    th  := xh2 + cAcoshRefineCh[i1, 0].f;
    tl  := (cAcoshRefineCh[i1, 0].f - th) + xh2;
    xh2 := th;
    sh2 := sh2 + tl + cAcoshRefineCh[i1, 1].f;
  end;
  sh_ := pcr_muldd(xh, xl, xh2, sh2, sl2);
  sh_ := AcoshAddDD(sh_, sl2, el1, el2, sl2);
  sh_ := AcoshAddDD(sh_, sl2, L1, L2, sl2);

  // v0 = fasttwosum(L0, sh_, &v2); v1 = fasttwosum(v2, sl2, &v2);
  // v0 = fasttwosum(v0, v1, &v1);  v1 = fasttwosum(v1, v2, &v2);
  pcr_fasttwosum(v0, vtmp, L0, sh_);
  pcr_fasttwosum(v1, v2, vtmp, sl2);
  pcr_fasttwosum(v0, v1, v0, v1);
  pcr_fasttwosum(v1, v2, v1, v2);

  sgn := pcr_copysign(Double(2.0), x);
  v0 := v0 * sgn;
  v1 := v1 * sgn;
  v2 := v2 * sgn;

  t.f := v1;
  if (t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) = 0 then
  begin
    w.f := v2;
    if ((w.u xor t.u) shr 63) <> 0 then t.u := t.u - 1
    else t.u := t.u + 1;
    v1 := t.f;
  end;

  t0.f := v0;
  er := ((t.u + 41) and (UInt64($FFFFFFFFFFFFFFFF) shr 12));
  de := ((t0.u shr 52) and $7FF) - ((t.u shr 52) and $7FF);
  res := v0 + v1;
  if (de > 99) or (er < 80) then
    Result := AsinhDatabase(x, res)
  else
    Result := res;
end;

function pcr_asinh(x: Double): Double;
var
  ax: Double;
  ix, t: Tb64u64;
  u: UInt64;
  x2h, x2l, x3h, sl, eps, lb, ub: Double;
  c1, c3, c5, x4: Double;
  ah, al, th, tl, rs: Double;
  ex, e, i1, i2: Int32;
  off: Int32;
  iuint, jshift: UInt64;
  dsigned: Int64;
  ed, r, dx, dx2, f: Double;
  lh, ll, sgn: Double;
begin
  ax := pcr_fabs(x);
  ix.f := ax;
  u := ix.u;

  // |x| < 0x1.bp-4 : small-x polynomial fast path
  if u < cAsinhBr_Med then
  begin
    if u < cAsinhBr_Tinier then
    begin
      // |x| < 0x1.7137449123ef7p-26: asinh(x) rounds to x
      if u = 0 then begin Result := x; Exit; end;
      // C uses fma(-2^-60, x, x); pcr_fma_pascal is not correctly rounded for
      // subnormal outputs, so return x directly (bit-exact for this range).
      Result := pcr_fma(cAsinhTinyScale.f * Double(-1.0), x, x);
      Exit;
    end;
    x2h := x * x;
    x2l := pcr_fma(x, x, -x2h);
    x3h := x2h * x;
    if u < cAsinhBr_Small then
    begin
      if u < cAsinhBr_VTiny then
      begin
        if u < cAsinhBr_Tiny then
          // |x| < 0x1.ap-26
          sl := x3h * cAsinhTiny1.f
        else
          sl := x3h * (cAsinhTiny2[0].f + x2h * cAsinhTiny2[1].f);
      end
      else
        sl := x3h * (cAsinhSmall[0].f + x2h * (cAsinhSmall[1].f
           + x2h * (cAsinhSmall[2].f + x2h * cAsinhSmall[3].f)));
    end
    else
    begin
      // 0x1.3p-6 <= |x| < 0x1.bp-4 : 7-coefficient minimax
      c1 := cAsinhMedCl[1].f + x2h * cAsinhMedCl[2].f;
      c3 := cAsinhMedCl[3].f + x2h * cAsinhMedCl[4].f;
      c5 := cAsinhMedCl[5].f + x2h * cAsinhMedCl[6].f;
      x4 := x2h * x2h;
      sl := x3h * (cAsinhMedCl[0].f + x2h * (c1 + x4 * (c3 + x4 * c5)));
    end;
    eps := cAsinhFastEps.f * x3h;
    lb := x + (sl - eps);
    ub := x + (sl + eps);
    if lb = ub then begin Result := lb; Exit; end;
    Result := AsinhZero(x, x2h, x2l);
    Exit;
  end;

  // |x| >= 0x1.bp-4 : main branch
  x2h := Double(0.0); x2l := Double(0.0);
  off := $3FF;
  if u < cAsinhBr_P26 then
  begin
    // |x| < 2^26
    x2h := x * x;
    x2l := pcr_fma(x, x, -x2h);
    if u < cAsinhBr_One then
      pcr_fasttwosum(th, tl, Double(1.0), x2h)
    else
      pcr_fasttwosum(th, tl, x2h, Double(1.0));
    tl := tl + x2l;
    ah := pcr_sqrt(th);
    rs := Double(0.5) / th;
    al := (tl - pcr_fma(ah, ah, -th)) * (rs * ah);
    pcr_fasttwosum(ah, tl, ah, ax);
    al := al + tl;
  end
  else if u < cAsinhBr_P52 then
  begin
    ah := Double(2.0) * ax;
    al := Double(0.5) / ax;
  end
  else
  begin
    if u >= cAsinhBr_Inf then begin Result := x + x; Exit; end;
    off := $3FE;
    ah := ax;
    al := Double(0.0);
  end;

  t.f := ah;
  ex := (t.u shr 52) and $7FF;
  e  := ex - off;
  t.u := t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12);
  ed := e;
  iuint := t.u shr (52 - 5);
  dsigned := Int64(t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 17));
  jshift := UInt64(
              Int64(t.u)
              + (Int64(cAcoshB[iuint].c0) shl 33)
              + (Int64(cAcoshB[iuint].c1) * (dsigned shr 16))
            ) shr (52 - 10);
  t.u := t.u or (UInt64($3FF) shl 52);
  i1 := Int32(jshift shr 5);
  i2 := Int32(jshift and $1F);
  r := cAcoshR1[i1].f * cAcoshR2[i2].f;
  dx := pcr_fma(r, t.f, -Double(1.0));
  dx2 := dx * dx;
  f := dx2 * ((cAcoshC[0].f + dx * cAcoshC[1].f)
            + dx2 * ((cAcoshC[2].f + dx * cAcoshC[3].f) + dx2 * cAcoshC[4].f));
  lh := cAcoshL2h.f * ed + (cAcoshL1[i1, 1].f + cAcoshL2[i2, 1].f);
  ll := cAcoshL2l.f * ed + cAcoshL1[i1, 0].f + cAcoshL2[i2, 0].f + al / ah + f;
  ll := ll + dx;
  sgn := pcr_copysign(Double(1.0), x);
  lh := lh * sgn;
  ll := ll * sgn;
  eps := cAsinhMainEps;
  lb := lh + (ll - eps);
  ub := lh + (ll + eps);
  if lb = ub then begin Result := lb; Exit; end;
  if ax < Double(0.25) then
  begin
    Result := AsinhZero(x, x2h, x2l);
    Exit;
  end;
  Result := AsinhRefine(x, ah, al, cAsinhInvLog2.f * pcr_fabs(lb));
end;
