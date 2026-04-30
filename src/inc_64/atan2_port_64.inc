// ---------------------------------------------------------------------------
// pcr_atan2 — correctly-rounded binary64 atan2(y, x).
// Ported from core-math/src/binary64/atan2/atan2.c
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
//   (Paul Zimmermann and Alexei Sibidanov, 2024-2025).
//
// Uses TInt64 (192-bit) infrastructure from pascoremathtypes (task 2.07a).
// ---------------------------------------------------------------------------
//
// Notes:
//   * The C reference uses fenv_t/feholdexcept/fetestexcept/feupdateenv to
//     manage FE_OVERFLOW / FE_UNDERFLOW / FE_INEXACT around y/x and the
//     final tint_tod conversion. The pas-core-math port intentionally omits
//     this fenv tracking (CORE_MATH_SUPPORT_ERRNO is out of scope) — the
//     resulting Double value matches the C reference; only the FE flags side
//     effects diverge.
//   * The accurate-path branch when t = y/x is exact and tiny falls back to
//     fma(t, -2^-54, t) just as the C does, so the round-to-nearest result
//     is byte-identical.
//
{$I inc_64/atan2_const_64.inc}

// PI / PI/2 / PI/4 as double-double pairs (matches PI_H/PI_L etc.).
const
  cAtan2PiH:    Tb64u64 = (u: $400921FB54442D18); // 0x1.921fb54442d18p+1
  cAtan2PiL:    Tb64u64 = (u: $3CA1A62633145C07); // 0x1.1a62633145c07p-53
  cAtan2Pi2H:   Tb64u64 = (u: $3FF921FB54442D18); // 0x1.921fb54442d18p+0
  cAtan2Pi2L:   Tb64u64 = (u: $3C91A62633145C07); // 0x1.1a62633145c07p-54
  cAtan2P54Neg: Tb64u64 = (u: QWord($BC90000000000000)); // -0x1p-54
  cAtan2P105:   Tb64u64 = (u: $4680000000000000); //  0x1p+105
  cAtan2PM915:  Tb64u64 = (u: $44C0000000000000); //  0x1p-915
  cAtan2P968:   Tb64u64 = (u: $7C70000000000000); //  0x1p+968
  cAtan2PM968:  Tb64u64 = (u: $0370000000000000); //  0x1p-968 (threshold const)
  cAtan2P1022:  Tb64u64 = (u: $7FD0000000000000); //  0x1p+1022
  cAtan2PM1:    Tb64u64 = (u: $3FE0000000000000); //  0.5
  cAtan2EpsK:   Tb64u64 = (u: $3CC0510000000000); //  0x1.051p-51 fast-eps slope
  cAtan2EpsAdd: Tb64u64 = (u: $3A50000000000000); //  0x1p-90 fast-eps additive
  cAtan2EpsRef: Tb64u64 = (u: $3CD4000000000000); //  0x1.4p-50 refine eps slope
  cAtan2EpsRefA:Tb64u64 = (u: $3CC0000000000000); //  0x1p-51 refine eps additive
  cAtan2TinyT:  Tb64u64 = (u: $3E4D12ED0AF1A27F); //  0x1.d12ed0af1a27fp-27 (exact-t cap)

  // Special-path Inf-Inf finf table:
  //   finf[0] = (0x1p-55, 0x1.921fb54442d18p-1)   for x = -Inf  (kw mask bit)
  //   finf[1] = (0x1p-54, 0x1.2d97c7f3321d2p+1)   for x = +Inf
  cAtan2FInfNeg0: Tb64u64 = (u: $3C80000000000000); //  0x1p-55
  cAtan2FInfNeg1: Tb64u64 = (u: $3FE921FB54442D18); //  0x1.921fb54442d18p-1 (= pi/4)
  cAtan2FInfPos0: Tb64u64 = (u: $3C90000000000000); //  0x1p-54
  cAtan2FInfPos1: Tb64u64 = (u: $4002D97C7F3321D2); //  0x1.2d97c7f3321d2p+1 (= 3pi/4)

  cAtan2_53shl52: UInt64 = UInt64(53) shl 52;       // 53<<52, dxy threshold

// ---------------------------------------------------------------------------
// Local fasttwosum / fastsum (match C's atan2.c definitions exactly).
// pcr_fasttwosum is identical, but kept local for clarity / to avoid var-out
// signature mismatches in tight code.
// ---------------------------------------------------------------------------

procedure A2FastTwoSum(out s, e: Double; x, y: Double); inline;
var ss: Double;
begin
  ss := x + y;
  e  := y - (ss - x);
  s  := ss;
end;

// fastsum: returns s = (xh + yh) high; *e = (xl + yl) + sl  with sl = errfast(xh,yh)
procedure A2FastSum(out sh, e: Double; xh, xl, yh, yl: Double); inline;
var sl: Double;
begin
  A2FastTwoSum(sh, sl, xh, yh);
  e := (xl + yl) + sl;
end;

// ---------------------------------------------------------------------------
// Special-path: NaN / Inf / zero handling.
// Mirrors as_atan2_special() from atan2.c.
// ---------------------------------------------------------------------------

function Atan2Special(y0, x0: Double): Double;
var
  iy, ix: Tb64u64;
  aiy, aix, em, sxbit: UInt64;
  hi_pick, lo_pick: Double;
begin
  iy.f := y0;  ix.f := x0;
  aiy := iy.u shl 1;
  aix := ix.u shl 1;
  em  := UInt64($7FF) shl 53;     // 2 * (0x7ff << 52)
  if (aiy >= em) or (aix >= em) then begin
    // NaN if any > em (since em is exactly 2*Inf magnitude).
    if (aiy > em) or (aix > em) then begin
      Result := y0 + x0;
      Exit;
    end;
    // Now neither is NaN; at least one is +/-Inf.
    if (aiy = em) and (aix = em) then begin
      // both +/-Inf: ix.u>>63 selects finf[]:
      //   x = +Inf -> finf[0] = (0x1p-55, pi/4)   -> +/-pi/4
      //   x = -Inf -> finf[1] = (0x1p-54, 3pi/4)  -> +/-3pi/4
      sxbit := ix.u shr 63;
      if sxbit = 0 then begin
        lo_pick := cAtan2FInfNeg0.f;     // 0x1p-55
        hi_pick := cAtan2FInfNeg1.f;     // pi/4
      end else begin
        lo_pick := cAtan2FInfPos0.f;     // 0x1p-54
        hi_pick := cAtan2FInfPos1.f;     // 3pi/4
      end;
      Result := pcr_copysign(hi_pick, y0) + pcr_copysign(lo_pick, y0);
      Exit;
    end;
    // only one of y, x is +/-Inf
    if aix = em then begin
      if x0 < 0 then
        Result := pcr_copysign(cAtan2PiH.f, y0) + pcr_copysign(cAtan2PiL.f, y0)
      else
        Result := pcr_copysign(Double(0.0), y0);
      Exit;
    end;
    // y = +/-Inf, x finite -> +/-pi/2
    Result := pcr_copysign(cAtan2Pi2H.f, y0) + pcr_copysign(cAtan2Pi2L.f, y0);
    Exit;
  end;

  // Zero handling: aiy == 0 or aix == 0
  if (aiy = 0) and (aix = 0) then begin
    if ix.u = 0 then begin     // atan2(+/-0, +0) = +/-0
      Result := y0;
      Exit;
    end;
    // atan2(+/-0, -0) = +/-pi
    if iy.u = 0 then
      Result := cAtan2PiH.f + cAtan2PiL.f
    else
      Result := -cAtan2PiH.f - cAtan2PiL.f;
    Exit;
  end;
  if aiy = 0 then begin
    if x0 > 0 then begin
      Result := y0 / x0;
      Exit;
    end;
    if (iy.u shr 63) = 0 then
      Result := cAtan2PiH.f + cAtan2PiL.f
    else
      Result := -cAtan2PiH.f - cAtan2PiL.f;
    Exit;
  end;
  // only x is zero — y nonzero, finite
  if y0 > 0 then
    Result := cAtan2Pi2H.f + cAtan2Pi2L.f
  else
    Result := -cAtan2Pi2H.f - cAtan2Pi2L.f;
end;

// ---------------------------------------------------------------------------
// Accurate path — tint_t-based [29,29] rational approximation of atan(z).
// Mirrors atan2_accurate() from atan2.c (lines 200-395).
// ---------------------------------------------------------------------------

function Atan2Accurate(y, x: Double): Double;
var
  t, corr: Double;
  res: Double;
  inv: Boolean;
  z, p, q, tmp, p_tmp, q_tmp: TInt64;
  i: Integer;
  err: UInt64;
  sz: UInt64;
  ay, ax: Double;
begin
  // Compute t = y/x.  We do *not* track FE flags — see header note.
  t := y / x;

  if t = Double(0.0) then begin
    if x > Double(0.0) then begin
      Result := t;            // atan(0) = +/-0 with sign(t)
      Exit;
    end;
    if y > Double(0.0) then
      Result := cAtan2PiH.f + cAtan2PiL.f
    else
      Result := -cAtan2PiH.f - cAtan2PiL.f;
    Exit;
  end;

  corr := pcr_fma(t, x, -y);
  if (corr = Double(0.0)) and (x > Double(0.0)) then begin
    // t exact
    if Abs(t) <= cAtan2TinyT.f then begin
      ay := Abs(y);
      if ay >= cAtan2PM915.f then begin
        Result := pcr_fma(t, cAtan2P54Neg.f, t);
        Exit;
      end;
      // |y| < 2^-915 — scale to avoid underflow in t*x - y
      corr := pcr_fma(t * cAtan2P105.f, x, -y * cAtan2P105.f);
      if corr = Double(0.0) then begin
        Result := pcr_fma(t, cAtan2P54Neg.f, t);
        Exit;
      end;
    end;
  end;

  ax := Abs(x);  ay := Abs(y);
  inv := ay > ax;
  if inv then DivTIntD(z, x, y) else DivTIntD(z, y, x);

  // Tiny-z exact-t shortcut for x>0: |y/x| < 2^-96
  if (not inv) and (x > Double(0.0)) and (z.ex <= -96) then begin
    // z->l -= 2; z->m -= (z->l < 2); z->h -= (z->m < 1);
    // Faithful port — C uses post-decrement comparisons; in practice z->l is
    // a large noise word, so the conditional branches almost never fire.
    tmp := z;
    tmp.l := tmp.l - 2;
    if tmp.l < 2 then tmp.m := tmp.m - 1;
    if tmp.m < 1 then tmp.h := tmp.h - 1;
    Result := TIntToD(tmp, 1, y, x);
    Exit;
  end;

  sz := z.sgn;
  z.sgn := 0;

  CpTInt(p, cAtan2P[29]);
  CpTInt(q, cAtan2Q[29]);
  for i := 28 downto 0 do begin
    MulTInt(p_tmp, p, z);
    MulTInt(q_tmp, q, z);
    AddTInt_NoAlias(p, p_tmp, cAtan2P[i]);
    AddTInt_NoAlias(q, q_tmp, cAtan2Q[i]);
  end;
  // multiply p by z
  MulTInt(p, p, z);
  // divide p by q -> z := p / q
  DivTInt(z, p, q);

  err := 662;
  z.sgn := sz;

  if inv then begin
    if z.sgn = 0 then begin
      z.sgn := 1;
      AddTInt(z, TINT_PI2, z);
    end else begin
      AddTInt(z, TINT_PI2, z);
      z.sgn := 1;
    end;
    err := 524;
  end;

  if x < Double(0.0) then begin
    if z.sgn = 0 then begin
      z.sgn := 1;
      AddTInt(z, TINT_PI, z);
      z.sgn := 1;
    end else begin
      AddTInt(z, TINT_PI, z);
    end;
    err := 266;
  end;

  res := TIntToD(z, err, y, x);
  Result := res;
end;

// ---------------------------------------------------------------------------
// pcr_atan2 — main entry (fast path + Ziv-test dispatch to Atan2Accurate).
// Mirrors cr_atan2 (atan2.c lines 461-586).
// ---------------------------------------------------------------------------

function pcr_atan2(y0, x0: Double): Double;
const
  MASK = UInt64($7FFFFFFFFFFFFFFF);
var
 iy, ix, sgn, jj: Tb64u64;
  aiy, aix, sy, sx, GT, dxy, kw: UInt64;
  jt: Int64;
  ax, ay, x, y, fh, fl, t0, zn, zd, z, z2, dz, eps, rh, rl, lb, ub: Double;
  dh, dl, e, rdh, nh, nl, dt, y1, zh, zl: Double;
  signsel: Double;
begin
  iy.f := y0;  ix.f := x0;
  aiy := iy.u and MASK;
  if (aiy = 0) or (aiy >= UInt64($7FF) shl 52) then begin
    Result := Atan2Special(y0, x0);
    Exit;
  end;
  aix := ix.u and MASK;
  if (aix = 0) or (aix >= UInt64($7FF) shl 52) then begin
    Result := Atan2Special(y0, x0);
    Exit;
  end;

  ax := Abs(x0);
  ay := Abs(y0);
  x  := pcr_fmax(ax, ay);
  y  := pcr_fmin(ax, ay);
  sy := iy.u shr 63;
  sx := ix.u shr 63;
  GT := UInt64(Ord(aix < aiy));
  // dxy = (aix - aiy) xor -GT  — magnitude of difference (off-by-one when GT=1,
  // irrelevant at the 53<<52 threshold).
  dxy := (aix - aiy) xor UInt64(-Int64(GT));
  if dxy >= cAtan2_53shl52 then begin
    Result := Atan2Accurate(y0, x0);
    Exit;
  end;

  // sgn = asgn[GT xor sx xor sy]   asgn[0]=+0, asgn[1]=-0
  if (GT xor sx xor sy) = 0 then sgn.u := 0 else sgn.u := UInt64($8000000000000000);

  kw := (sx shl 2) or (sy shl 1) or GT;

  // jj.f = y/x + (2 + 1/128)
  jj.f := y / x + (Double(2.0) + (Double(1.0)/Double(128.0)));
  jt   := Int64((jj.u shr (52 - 7)) and 127);

  signsel := pcr_copysign(Double(1.0), sgn.f);
  fh := cAtan2F2[2*jt + 1].f * signsel;
  fl := cAtan2F2[2*jt    ].f * signsel;
  fh := fh + cAtan2O[2*kw    ].f;
  fl := fl + cAtan2O[2*kw + 1].f;

  if x < cAtan2PM968.f then begin
    x := x * cAtan2P968.f;
    y := y * cAtan2P968.f;
  end;
  if x > cAtan2P1022.f then begin
    if jt <> 0 then begin
      x := x * cAtan2PM1.f;
      y := y * cAtan2PM1.f;
    end;
  end;

  t0 := cAtan2T2[jt].f;
  {$IFDEF AVX2}
  // zn = (-t0)*x + y ; zd = t0*y + x  (hardware FMA, avoids pcr_fma call/ret).
  asm
    movsd        xmm0, t0
    movsd        xmm1, x
    movsd        xmm2, y
    vfnmadd213sd xmm0, xmm1, xmm2          // xmm0 = -(t0*x) + y
    movsd        zn, xmm0
    movsd        xmm0, t0
    movsd        xmm1, y
    movsd        xmm2, x
    vfmadd213sd  xmm0, xmm1, xmm2          // xmm0 = t0*y + x
    movsd        zd, xmm0
  end ['xmm0', 'xmm1', 'xmm2'];
  {$ELSE}
  zn := pcr_fma(-t0, x, y);
  zd := pcr_fma( t0, y, x);
  {$ENDIF}
  z  := zn / zd;

  z2 := z * z;
  z  := z * signsel;
  dz := (z * z2) * (cAtan2B_0.f + z2 * (cAtan2B_1.f + z2 * cAtan2B_2.f));

  eps := Abs(z) * cAtan2EpsK.f + cAtan2EpsAdd.f;

  A2FastTwoSum(rh, z, fh, z);   // rh = fh + z; z = error
  rl := (fl + dz) + z;
  lb := rh + (rl - eps);
  ub := rh + (rl + eps);

  if lb <> ub then begin
    dh := y * t0;
    dl := pcr_fma(y, t0, -dh);
    A2FastTwoSum(dh, e, x, dh);
    if Abs(dh) <= cAtan2P1022.f then
      rdh := Double(1.0) / dh
    else begin
      // matches "fegetexceptflag/fesetexceptflag(FE_UNDERFLOW)" guard in C —
      // we don't track the flag so just compute 1/dh (correct value).
      rdh := Double(1.0) / dh;
    end;
    dl := dl + e;

    nh := x * t0;
    nl := pcr_fma(x, t0, -nh);
    dt := y - nh;
    y1 := dt + nh;
    if y1 = y then begin
      A2FastTwoSum(nh, nl, dt, -nl);
    end else begin
      A2FastTwoSum(nh, nl, dt, (y - y1) - nl);
    end;
    zh := nh * rdh;
    z2 := zh * zh;
    zl := rdh * (pcr_fma(dh, -zh, nh) + (nl - (nh * rdh) * dl));

    zl := zl + zh * z2 * ((cAtan2B2_0.f + z2 * cAtan2B2_1.f) + (z2 * z2) * (cAtan2B2_2.f + z2 * cAtan2B2_3.f));
    zh := zh * signsel;
    zl := zl * signsel;

    eps := cAtan2EpsRef.f * (Abs(zh) * z2 + cAtan2EpsRefA.f);
    A2FastSum(fh, fl, fh, fl, zh, zl);
    lb := fh + (fl - eps);
    ub := fh + (fl + eps);
    if lb <> ub then begin
      Result := Atan2Accurate(y0, x0);
      Exit;
    end;
  end;
  Result := ub;
end;
