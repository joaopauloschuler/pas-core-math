// ---------------------------------------------------------------------------
// pcr_asinpi — correctly-rounded binary64 asin(x)/pi.
// Ported from core-math/src/binary64/asinpi/asinpi.c (Sibidanov & Zimmermann).
//
// asinpi.c is the *u128-based* asin implementation with a final 1/pi multiply,
// distinct from asin.c's double-double version that pcr_asin ports.
// ---------------------------------------------------------------------------

{$I asinpi_const.inc}

// ---------------------------------------------------------------------------
// Local u128 helpers
// ---------------------------------------------------------------------------

function Asinpi_Muuh(a, b: UInt64): UInt64; inline;
begin
  Result := Mulu64u64(a, b).hi;
end;

// signed-high-64( i64 a * i64 b )
function Asinpi_Mh(a, b: Int64): Int64; inline;
var hi: UInt64;
begin
  hi := Mulu64u64(UInt64(a), UInt64(b)).hi;
  if a < 0 then hi := hi - UInt64(b);
  if b < 0 then hi := hi - UInt64(a);
  Result := Int64(hi);
end;

// (u128)a * b >> 64  (u64 x u128 -> u128 truncated)
function Asinpi_MuU(a: UInt64; const b: TUInt128): TUInt128; inline;
var hilo: UInt64;
begin
  Result := Mulu64u64(a, b.hi);
  hilo := Mulu64u64(a, b.lo).hi;
  Result.lo := Result.lo + hilo;
  if Result.lo < hilo then Inc(Result.hi);
end;

// (u128)a * (u128)b >> 64  (u128 x u128 -> u128 truncated)
function Asinpi_MUU(const a, b: TUInt128): TUInt128; inline;
var add1, add2: UInt64;
begin
  Result := Mulu64u64(a.hi, b.hi);
  add1   := Mulu64u64(a.lo, b.hi).hi;
  Result.lo := Result.lo + add1;
  if Result.lo < add1 then Inc(Result.hi);
  add2   := Mulu64u64(a.hi, b.lo).hi;
  Result.lo := Result.lo + add2;
  if Result.lo < add2 then Inc(Result.hi);
end;

// (u128)a * a >> 64
function Asinpi_SqrU(const a: TUInt128): TUInt128; inline;
var cross, os: TUInt128;
begin
  cross := Mulu64u64(a.lo, a.hi);
  // os = cross >> 63  (full 128 bit shift)
  os.lo := (cross.lo shr 63) or (cross.hi shl 1);
  os.hi := cross.hi shr 63;
  Result := Mulu64u64(a.hi, a.hi);
  AddU128(Result, Result, os);
end;

// (u128) a *= b  truncated to 128 bits  (in-place)
procedure Asinpi_MulU128xU64(var a: TUInt128; b: UInt64); inline;
var lo: TUInt128;
begin
  lo := Mulu64u64(a.lo, b);
  a.hi := a.hi * b + lo.hi;
  a.lo := lo.lo;
end;

// (u128) a += b  (in place)
procedure Asinpi_AddU64(var a: TUInt128; b: UInt64); inline;
begin
  a.lo := a.lo + b;
  if a.lo < b then Inc(a.hi);
end;

// pasin(x) — high 128 bits of asin polynomial
function Asinpi_Pasin(const x: TUInt128): TUInt128;
var
  xh: UInt64;
  t: TUInt128;
  acc, tmp: TUInt128;
begin
  xh := x.hi;
  t := cAsinpi_PasinCh[3];
  t.lo := t.lo + Asinpi_Muuh(xh,
    cAsinpi_PasinB[0] + Asinpi_Muuh(xh,
      cAsinpi_PasinB[1] + Asinpi_Muuh(xh,
        cAsinpi_PasinB[2] + Asinpi_Muuh(xh, cAsinpi_PasinB[3]))));
  // return mUU(x, ch[0] + mUU(x, ch[1] + mUU(x, ch[2] + mUU(x, t))))
  tmp := Asinpi_MUU(x, t);
  AddU128(acc, cAsinpi_PasinCh[2], tmp);
  tmp := Asinpi_MUU(x, acc);
  AddU128(acc, cAsinpi_PasinCh[1], tmp);
  tmp := Asinpi_MUU(x, acc);
  AddU128(acc, cAsinpi_PasinCh[0], tmp);
  Result := Asinpi_MUU(x, acc);
end;

// (i128)dc << ss   ss in [0, 64)
procedure Asinpi_ShlI64ToI128(dc: Int64; ss: Int32; out r: TUInt128); inline;
begin
  if ss = 0 then begin
    r.lo := UInt64(dc);
    if dc < 0 then r.hi := UInt64(-1) else r.hi := 0;
  end else if ss < 64 then begin
    r.lo := UInt64(dc) shl ss;
    if dc < 0 then
      r.hi := (UInt64(dc) shr (64 - ss)) or (UInt64(-1) shl ss)
    else
      r.hi := UInt64(dc) shr (64 - ss);
  end else if ss < 128 then begin
    r.lo := 0;
    r.hi := UInt64(dc) shl (ss - 64);
  end else begin
    r.lo := 0; r.hi := 0;
  end;
end;

// (i128)dc >> ss  arithmetic right shift, ss in [0, 64)
procedure Asinpi_ShrI64ArithToI128(dc: Int64; ss: Int32; out r: TUInt128); inline;
begin
  // (i128)dc has hi = sign-ext, lo = dc (as u64)
  if dc >= 0 then begin
    r.hi := 0;
    if ss = 0 then r.lo := UInt64(dc)
    else r.lo := UInt64(dc) shr ss;
  end else begin
    r.hi := UInt64(-1);
    if ss = 0 then r.lo := UInt64(dc)
    else r.lo := (UInt64(dc) shr ss) or (UInt64(-1) shl (64 - ss));
  end;
end;

// arithmetic right shift of Int64 by ss in (0, 63]
function Asinpi_SarI64(x: Int64; ss: Int32): Int64; inline;
begin
  if x < 0 then
    Result := Int64((UInt64(x) shr ss) or (UInt64(-1) shl (64 - ss)))
  else
    Result := Int64(UInt64(x) shr ss);
end;

// (u128)dc << ss  unsigned, ss in [0, 64)
procedure Asinpi_ShlU64ToU128(dc: UInt64; ss: Int32; out r: TUInt128); inline;
begin
  if ss = 0 then begin
    r.lo := dc; r.hi := 0;
  end else if ss < 64 then begin
    r.lo := dc shl ss;
    r.hi := dc shr (64 - ss);
  end else begin
    r.lo := 0;
    r.hi := dc shl (ss - 64);
  end;
end;

// ---------------------------------------------------------------------------
// asinpi_tiny:  |x| < 2^-53
// ---------------------------------------------------------------------------
function Asinpi_Tiny(x: Double): Double;
const
  kTwo106:   Tb64u64 = (u: $4690000000000000);   //  0x1p106
  kTwoM106:  Tb64u64 = (u: $3950000000000000);   //  0x1p-106
  kTwoM600:  Tb64u64 = (u: $1A70000000000000);   //  0x1p-600
  kH56:      Tb64u64 = (u: $3C7B824198B94A89);   //  0x1.b824198b94a89p-56
  kL110:     Tb64u64 = (u: $391FFFFFFFFFFFFF);   //  0x1.fffffffffffffp-110
  kBcp1022:  Tb64u64 = (u: $001BC03DF34E902C);   //  0x1.bc03df34e902cp-1022
  kBcp1023:  Tb64u64 = (u: $000DE01EF9A74816);   //  0x1.bc03df34e902cp-1023
  kBcp1024:  Tb64u64 = (u: $0006F00F7CD3A40B);   //  0x1.bc03df34e902cp-1024
  kCbp1024:  Tb64u64 = (u: $00072E0BD768E9AD);   //  0x1.cb82f5da3a6b4p-1024
var
  t, t1: Tb64u64;
  au: UInt64;
  e: Int32;
  h, l, sx, sgn, twoM600N: Double;
begin
  twoM600N := -kTwoM600.f;
  t.f := x;
  au := t.u and UInt64($7FFFFFFFFFFFFFFF);
  if x = 0.0 then Exit(x);

  if (au shl 12) = UInt64($59AF9A1194EFE000) then begin
    e := Int32((t.u shr 52) and $7FF);
    h := kH56.f;
    l := kL110.f;
    if x > 0 then sgn := Double(1.0) else sgn := Double(-1.0);
    t1.f := sgn;
    t1.u := t1.u - (UInt64($3C9 - e) shl 52);
    Result := pcr_fma(l, t1.f, h * t1.f);
    Exit;
  end;

  if au = UInt64($35CBA89AF1F855) then begin
    if x > 0 then
      Exit(pcr_fma(kTwoM600.f, twoM600N,    kBcp1022.f))
    else
      Exit(pcr_fma(kTwoM600.f, kTwoM600.f, -kBcp1022.f));
  end;

  if au = UInt64($15CBA89AF1F855) then begin
    if x > 0 then
      Exit(pcr_fma(twoM600N,    kTwoM600.f,  kBcp1024.f))
    else
      Exit(pcr_fma(kTwoM600.f,  kTwoM600.f, -kBcp1024.f));
  end;

  if au = UInt64($168E6482549DB1) then begin
    if x > 0 then
      Exit(pcr_fma(kTwoM600.f,  kTwoM600.f,  kCbp1024.f))
    else
      Exit(pcr_fma(twoM600N,    kTwoM600.f, -kCbp1024.f));
  end;

  if au = UInt64($25CBA89AF1F855) then begin
    if x > 0 then
      Exit(pcr_fma(twoM600N,    kTwoM600.f,  kBcp1023.f))
    else
      Exit(pcr_fma(kTwoM600.f,  kTwoM600.f, -kBcp1023.f));
  end;

  h  := x * cAsinpi_ONE_OVER_PI_H.f;
  sx := x * kTwo106.f;
  l  := pcr_fma(sx, cAsinpi_ONE_OVER_PI_H.f, -h * kTwo106.f);
  l  := pcr_fma(sx, cAsinpi_ONE_OVER_PI_L.f, l);
  Result := pcr_fma(l, kTwoM106.f, h);
end;

// ---------------------------------------------------------------------------
// asinpi_small:  2^-53 <= |x| < 2^-26
// ---------------------------------------------------------------------------
function Asinpi_Small(x: Double): Double;
var
  i: Int32;
  h, l, hh, ll: Double;
begin
  for i := 0 to 17 do begin
    if x = cAsinpi_SmallExcX[i].f then
      Exit(cAsinpi_SmallExcH[i].f + cAsinpi_SmallExcL[i].f);
    if x = -cAsinpi_SmallExcX[i].f then
      Exit(-cAsinpi_SmallExcH[i].f - cAsinpi_SmallExcL[i].f);
  end;
  h := cAsinpi_Small_C1H.f;
  l := pcr_fma(cAsinpi_Small_C3.f, x * x, cAsinpi_Small_C1L.f);
  hh := h * x;
  ll := pcr_fma(h, x, -hh);
  ll := pcr_fma(l, x, ll);
  Result := hh + ll;
end;

// ---------------------------------------------------------------------------
// Multiply fi by INV_PI (1/pi as 128-bit), updating fi in place.
// fi.a = INV_PI_H * fi.hi + (INV_PI_H * fi.lo) >> 64 + (INV_PI_L * fi.hi) >> 64
// ---------------------------------------------------------------------------
procedure Asinpi_MulInvPi(var fi: TUInt128); inline;
var m1, m2: TUInt128;
begin
  m1 := Mulu64u64(cAsinpi_INV_PI_H, fi.lo);
  m2 := Mulu64u64(cAsinpi_INV_PI_L, fi.hi);
  fi := Mulu64u64(cAsinpi_INV_PI_H, fi.hi);
  Asinpi_AddU64(fi, m1.hi);
  Asinpi_AddU64(fi, m2.hi);
end;

// ---------------------------------------------------------------------------
// asinpi_acc — slow / accurate path
// ---------------------------------------------------------------------------
function Asinpi_Acc(x: Double): Double;
var
  rm: TFPURoundingMode;
  t: Tb64u64;
  se, ss, sc, indx, k, nz, ce, ixe: Int32;
  xsign: UInt64;
  ax: Double;
  sm, cm, ixm: UInt64;
  sm2, cm2, fi: TUInt128;
  Sm128, dsm2, dsm3, CmU, DU: TUInt128;
  xx, c0, c2, x2: Double;
  ic, ixx, c: Tb64u64;
  h, dc: Int64;
  vU, v2U, twoV, pp, tmpU: TUInt128;
  vmsk: UInt64;
  rnd: UInt64;
begin
  Result := 0;
  rm := pcr_GetRoundMode;
  t.f := x;
  se := Int32((t.u shr 52) and $7FF) - $3FF;
  xsign := t.u and (UInt64(1) shl 63);
  ax := Abs(x);
  sm := (t.u shl 11) or (UInt64(1) shl 63);
  sm2 := Mulu64u64(sm, sm);

  if ax < Double(0.0131875) then begin
    ss := 2 * se;
    ShrU128(sm2, -14 - ss);
    Sm128.lo := 0;
    Sm128.hi := sm shr 1;
    fi := Asinpi_MuU(sm shr 1, Asinpi_Pasin(sm2));
    AddU128(fi, fi, Sm128);
    se := se + $3FF;
  end else begin
    xx := pcr_fma(x, -x, Double(1.0));
    ixx.f := Double(1.0) / xx;
    c.f   := Sqrt(xx);
    ixx.f := ixx.f * c.f;
    x2 := x * x;
    c0 := cAsinpi_MainCh[0].f + ax * cAsinpi_MainCh[1].f;
    c2 := cAsinpi_MainCh[2].f + ax * cAsinpi_MainCh[3].f;
    c0 := c0 + x2 * c2;
    ic.f := c0 * c.f + Double(64.0);
    indx := Int32((((ic.u and (UInt64(not 0) shr 12)) + (Int64(1) shl (52 - 7))) shr (52 - 6)));
    cm := (c.u shl 11) or (UInt64(1) shl 63);
    ce := Int32((c.u shr 52)) - $3FF;
    cm2 := Mulu64u64(cm, cm);
    ss := 52 + 2 * se;
    ShlU128(sm2, ss);
    sc := 52 + 2 * ce;
    ShlU128(cm2, sc);
    AddU128(sm2, sm2, cm2);
    h := Int64(sm2.hi);
    ixm := (ixx.u and (UInt64(not 0) shr 12)) or (UInt64(1) shl 52);
    ixe := Int32(ixx.u shr 52) - $3FF;
    dc := Asinpi_Mh(h, Int64(ixm));
    // dsm2 = (i128)(dc * (cm>>1)) signed
    dsm2 := Mulu64u64(UInt64(dc), cm shr 1);
    if dc < 0 then dsm2.hi := dsm2.hi - (cm shr 1);
    ShlU128(dsm2, 13);
    SubU128(sm2, sm2, dsm2);
    // dsm3 = (i128)(dc * dc) >= 0
    dsm3 := Mulu64u64(UInt64(dc), UInt64(dc));
    if dc < 0 then begin
      dsm3.hi := dsm3.hi - UInt64(dc);
      dsm3.hi := dsm3.hi - UInt64(dc);
    end;
    sc := 28 - ixe * 2;
    if sc >= 0 then ShrU128(dsm3, sc) else ShlU128(dsm3, -sc);
    AddU128(sm2, sm2, dsm3);
    k := ixe - ce;
    ss := 24 + k;
    CmU.lo := 0;
    CmU.hi := cm;
    Asinpi_ShlI64ToI128(dc, ss, DU);
    SubU128(CmU, CmU, DU);
    // h = sm2 >> 14 (low 64 of result)
    h := Int64((sm2.lo shr 14) or (sm2.hi shl 50));
    dc := Asinpi_Mh(h, Int64(ixm));
    ss := 26 - k;
    if ss >= 0 then begin
      Asinpi_ShrI64ArithToI128(dc, ss, tmpU);
    end else begin
      Asinpi_ShlU64ToU128(UInt64(dc), -ss, tmpU);
    end;
    SubU128(CmU, CmU, tmpU);
    fi.lo := UInt64($D313198A2E037073);
    fi.hi := UInt64($03243F6A8885A308);
    Asinpi_MulU128xU64(fi, UInt64(64 - indx));
    if indx = 0 then begin
      ShrU128(CmU, -ce - 7);
      v2U := Asinpi_SqrU(CmU);
      pp := Asinpi_Pasin(v2U);
      tmpU := Asinpi_MUU(CmU, pp);
      AddU128(CmU, CmU, tmpU);
      ShrU128(CmU, 7);
      SubU128(fi, fi, CmU);
    end else begin
      vU := Asinpi_MuU(sm shr (-se), cAsinpi_AccS[indx - 1]);
      tmpU := Asinpi_MUU(CmU, cAsinpi_AccS[63 - indx]);
      ShrU128(tmpU, -ce);
      SubU128(vU, vU, tmpU);
      if (vU.hi and (UInt64(1) shl 63)) <> 0 then
        vmsk := UInt64(-1)
      else
        vmsk := 0;
      v2U := Asinpi_SqrU(vU);
      AddU128(twoV, vU, vU);
      twoV.lo := twoV.lo and vmsk;
      twoV.hi := twoV.hi and vmsk;
      SubU128(v2U, v2U, twoV);
      ShlU128(v2U, 14);
      pp := Asinpi_Pasin(v2U);
      tmpU := Asinpi_MUU(pp, vU);
      AddU128(vU, vU, tmpU);
      tmpU.lo := pp.lo and vmsk;
      tmpU.hi := pp.hi and vmsk;
      SubU128(vU, vU, tmpU);
      AddU128(fi, fi, vU);
    end;
    se := $3FE;
  end;

  Asinpi_MulInvPi(fi);

  nz := pcr_clzll(fi.hi);
  case rm of
    rmNearest: rnd := (fi.hi shr (10 - nz)) and 1;
    rmDown:    rnd := xsign shr 63;
    rmUp:      rnd := (xsign shr 63) xor 1;
  else
    rnd := 0;
  end;
  t.u := (fi.hi shr (11 - nz)) + ((UInt64(se - nz) shl 52) or xsign or rnd);
  Result := t.f;
end;

// ---------------------------------------------------------------------------
// pcr_asinpi — main entry point
// ---------------------------------------------------------------------------
function pcr_asinpi(x: Double): Double;
var
  rm: TFPURoundingMode;
  t, ic, ixx, c: Tb64u64;
  e, ss, sc, indx, ce, ixe, nz: Int32;
  m: UInt64;
  xsign, sm, cm, ixm, d, v2u, v3u: UInt64;
  fi, m1, m2, fpu, fpd, tmp, V128: TUInt128;
  xx, c0, c2, x2, ax: Double;
  Smh, Cmh, v_i, v2_i, v3_i: Int64;
  Vh, Vl: UInt64;
  rnd: UInt64;
begin
  Result := 0;
  rm := pcr_GetRoundMode;
  t.f := x;
  e := Int32((t.u shr 52) and $7FF) - $3FF;
  xsign := t.u and (UInt64(1) shl 63);
  sm := (t.u shl 11) or (UInt64(1) shl 63);

  if e >= 0 then begin
    m := t.u shl 12;
    if (e = 0) and (m = 0) then Exit(x * Double(0.5));
    if (e = $400) and (m <> 0) then Exit(x + x);
    pcr_feraiseexcept_invalid;
    Result := 0.0 / 0.0;
    Exit;
  end;

  if e < -6 then begin
    if e < -26 then begin
      if e < -53 then Exit(Asinpi_Tiny(x));
      Exit(Asinpi_Small(x));
    end;
    // 2^-26 <= |x| < 2^-6
    v2u := Asinpi_Muuh(sm, sm);
    v3u := Asinpi_Muuh(sm, v2u);
    v2u := v2u shr (-2 * e - 14);
    d := Asinpi_Muuh(v3u,
      cAsinpi_MainB[0] + Asinpi_Muuh(v2u,
        cAsinpi_MainB[1] + Asinpi_Muuh(v2u,
          cAsinpi_MainB[2] + Asinpi_Muuh(v2u,
            cAsinpi_MainB[3] + Asinpi_Muuh(v2u, cAsinpi_MainB[4])))));
    ss := 63 + 2 * e;
    fi.lo := d shl ss;
    fi.hi := (d shr (64 - ss)) + (sm shr 1);

    Asinpi_MulInvPi(fi);

    nz := pcr_clzll(fi.hi);
    if rm = rmNearest then Inc(nz);
    fpu := fi;
    tmp.lo := UInt64(15) shl ss;
    tmp.hi := 0;
    AddU128(fpu, fpu, tmp);
    if (((fi.hi xor fpu.hi) shr (11 - nz)) and 1) <> 0 then
      Exit(Asinpi_Acc(x));
    e := e + $3FF;
  end else begin
    // |x| >= 2^-6
    xx := pcr_fma(x, -x, Double(1.0));
    ixx.f := Double(1.0) / xx;
    c.f := Sqrt(xx);
    ixx.f := ixx.f * c.f;
    ax := Abs(x);
    x2 := x * x;
    c0 := cAsinpi_MainCh[0].f + ax * cAsinpi_MainCh[1].f;
    c2 := cAsinpi_MainCh[2].f + ax * cAsinpi_MainCh[3].f;
    c0 := c0 + x2 * c2;
    c0 := c0 * c.f;
    c0 := c0 + Double(64.0);
    ic.f := c0;
    indx := Int32((((ic.u and (UInt64(not 0) shr 12)) + (Int64(1) shl (52 - 7))) shr (52 - 6)));
    cm := (c.u shl 11) or (UInt64(1) shl 63);
    ce := Int32(c.u shr 52) - $3FF;

    // sm2, cm2 (use locals)
    fi := Mulu64u64(sm, sm); // reuse fi as sm2
    tmp := Mulu64u64(cm, cm); // tmp as cm2
    ss := 52 + 2 * e;
    ShlU128(fi, ss);
    sc := 52 + 2 * ce;
    if sc >= 0 then ShlU128(tmp, sc) else ShrU128(tmp, -sc);
    AddU128(fi, fi, tmp);
    v_i := Int64(fi.hi);  // C: i64 h = sm2.bh
    ixm := (ixx.u and (UInt64(not 0) shr 12)) or (UInt64(1) shl 52);
    ixe := Int32(ixx.u shr 52) - $3FF;

    // Smh = (sm << (6+e)) - sh[64-indx]
    ss := 6 + e;
    Smh := Int64(sm shl ss) - Int64(cAsinpi_MainSh[64 - indx]);
    sc := 6 + ce;
    if sc >= 0 then
      Cmh := Int64(cm shl sc)
    else
      Cmh := Int64(cm shr (-sc));
    Cmh := Cmh - Int64(cAsinpi_MainSh[indx]);
    // Cmh -= mh(h, ixm) >> (34 - ixe)
    Cmh := Cmh - Asinpi_SarI64(Asinpi_Mh(v_i, Int64(ixm)), 34 - ixe);

    v_i  := Asinpi_Mh(Smh, Int64(cAsinpi_MainS[indx])) - Asinpi_Mh(Cmh, Int64(cAsinpi_MainS[64 - indx]));
    v2_i := Asinpi_Mh(v_i, v_i);
    v3_i := Asinpi_Mh(v2_i, v_i);
    // v += mh(v3, a[0] + muuh(v2, a[1] + muuh(v2, a[2] + muuh(v2, a[3]))));
    v_i := v_i + Asinpi_Mh(v3_i,
      Int64(cAsinpi_MainA[0] + Asinpi_Muuh(UInt64(v2_i),
        cAsinpi_MainA[1] + Asinpi_Muuh(UInt64(v2_i),
          cAsinpi_MainA[2] + Asinpi_Muuh(UInt64(v2_i), cAsinpi_MainA[3])))));

    fi.lo := UInt64($D313198A2E037073);
    fi.hi := UInt64($03243F6A8885A308);
    Asinpi_MulU128xU64(fi, UInt64(64 - indx));

    // Vh = v >> 5 (arithmetic), Vl = (u64)v << 59
    Vh := UInt64(Asinpi_SarI64(v_i, 5));
    Vl := UInt64(v_i) shl 59;
    V128.hi := Vh;
    V128.lo := Vl;
    AddU128(fi, fi, V128);

    Asinpi_MulInvPi(fi);

    nz := pcr_clzll(fi.hi);
    if rm = rmNearest then Inc(nz);
    fpu := fi;
    fpd := fi;
    tmp.lo := UInt64(124) shl 55;
    tmp.hi := 0;
    AddU128(fpu, fpu, tmp);
    SubU128(fpd, fpd, tmp);
    if (((fpd.hi xor fpu.hi) shr (11 - nz)) and 1) <> 0 then
      Exit(Asinpi_Acc(x));
    e := $3FE;
  end;

  nz := pcr_clzll(fi.hi);
  case rm of
    rmNearest: rnd := (fi.hi shr (10 - nz)) and 1;
    rmDown:    rnd := xsign shr 63;
    rmUp:      rnd := (xsign shr 63) xor 1;
  else
    rnd := 0;
  end;
  t.u := ((fi.hi shr (11 - nz)) + (UInt64(e - nz) shl 52 or rnd)) or xsign;
  Result := t.f;
end;
