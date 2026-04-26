// ===========================================================================
// pcr_cosh — correctly-rounded binary64 cosh(x).
// Ported from core-math/src/binary64/cosh/cosh.c by Alexei Sibidanov.
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
//
// Shares with other exp-family ports already in pascoremath64.pas:
//   * cExpT0 / cExpT1 lookup tables (byte-identical to cosh's t0[]/t1[])
//   * TanhExpAccurate == cosh's as_exp_accurate (identical ch[3][2],
//     polynomial seed cTanhExpP0/P1/P2, and l2h/l2l/l2ll split)
// Cosh-specific content below:
//   * Fast-path c[5] polynomial for |x| < 0.125
//   * as_cosh_zero ch[4][2] + cl[4] polynomial
//   * as_cosh_database 21-entry table
//   * Main-path ch[4] polynomial
// ===========================================================================

const
  cCoshS:       Tb64u64 = (u: $40B71547652B82FE); //  0x1.71547652b82fep+12
  cCoshMagic:   Tb64u64 = (u: $4198000002000000); //  0x1.8000002p+26
  cCoshMagicH:  Tb64u64 = (u: $4198000000000000); //  0x1.8p+26
  cCoshMask26:  UInt64  = QWord($FFFFFFFFFC000000);

  cCoshAxTiny:  UInt64  = $3FC0000000000000; //  0.125
  cCoshAxP26:   UInt64  = $3E50000000000000; //  2^-26
  cCoshAxHuge:  UInt64  = $408633CE8FB9F87D; //  0x1.633ce8fb9f87dp+9
  cCoshAxMed:   UInt64  = $4014000000000000; //  5.0
  cCoshAxBig:   UInt64  = $40425E4F7B2737FA; //  ~36.736801
  cCoshAxBigBr: UInt64  = $403F666666666666; //  0x1.f666666666666p+4

  // Fast-path tiny error: e = x2 * 0x1.84p-51
  cCoshEps184:  Tb64u64 = (u: $3CC8400000000000); //  0x1.84p-51
  cCoshTinyAdj: Tb64u64 = (u: $3C80000000000000); //  0x1p-55

  // Bracket-epsilon constants for the three main paths
  cCoshErrLrg:  Tb64u64 = (u: $3C5621FB54442D18); //  0.12e-18 (approx)
  cCoshErrMed:  Tb64u64 = (u: $3C02020000000000); //  0x1.202p-63
  cCoshErrSml:  Tb64u64 = (u: $3C1C0A0000000000); //  0x1.c0ap-62

  // Overflow return 0x1p1023 * 2 (computed at runtime to avoid precision traps)
  cCoshOvf:     Tb64u64 = (u: $7FE0000000000000); //  0x1p+1023

  // Fast-path c[5] polynomial (|x| < 0.125 branch)
  cCoshFastC0:  Tb64u64 = (u: $3FE0000000000000); //  0x1p-1
  cCoshFastC1:  Tb64u64 = (u: $3FA5555555555554); //  0x1.5555555555554p-5
  cCoshFastC2:  Tb64u64 = (u: $3F56C16C16C1D0C0); //  0x1.6c16c16c1d0cp-10
  cCoshFastC3:  Tb64u64 = (u: $3EFA01A0075066B4); //  0x1.a01a0075066b4p-16
  cCoshFastC4:  Tb64u64 = (u: $3E927FAFF8DCC1C8); //  0x1.27faff8dcc1c8p-22

  // Main-path ch[4] polynomial (used for all |x| >= 0.125)
  cCoshMainCh0: Tb64u64 = (u: $3FF0000000000000); //  0x1p+0
  cCoshMainCh1: Tb64u64 = (u: $3FE0000000000000); //  0x1p-1
  cCoshMainCh2: Tb64u64 = (u: $3FC5555555AAAAAE); //  0x1.5555555aaaaaep-3
  cCoshMainCh3: Tb64u64 = (u: $3FA55555551C98C0); //  0x1.55555551c98cp-5

  // as_cosh_zero ch[4][2]
  cCoshZeroCh: array[0..3, 0..1] of Tb64u64 = (
    ((u:$3FE0000000000000),(u:QWord($B90C7E8DB669F624))),
    ((u:$3FA5555555555555),(u:$3C45555555556135)),
    ((u:$3F56C16C16C16C17),(u:QWord($BBEF49F4A6E838F2))),
    ((u:$3EFA01A01A01A01A),(u:$3B3A4FFBE15316AA)));

  // as_cosh_zero cl[4]
  cCoshZeroCl: array[0..3] of Tb64u64 = (
    (u:$3E927E4FB7789F5C),
    (u:$3E21EED8EFF9089C),
    (u:$3DA939749CE13DAD),
    (u:$3D2AE9891EFB6691));

  // as_cosh_database 21 × 3 entries
  cCoshDb: array[0..20, 0..2] of Tb64u64 = (
    ((u:$3FB9A5E3CBE1985E),(u:$3FF01492F72F984B),(u:QWord($B940000000000000))),
    ((u:$3FC52A11832E847D),(u:$3FF0381E68CAC923),(u:$3970000000000000)),
    ((u:$3FCBF0305E2C6C37),(u:$3FF061F4C39E16F2),(u:$3940000000000000)),
    ((u:$3FD17326FFC09F68),(u:$3FF099318A43AC80),(u:$3970000000000000)),
    ((u:$3FD3D27BF16D8BDB),(u:$3FF0C6091056E06A),(u:QWord($B940000000000000))),
    ((u:$3FE03923F2B47C07),(u:$3FF219C1989E3373),(u:QWord($BC90000000000000))),
    ((u:$3FEA6031CD5F93BA),(u:$3FF5BFF041B260FE),(u:QWord($B940000000000000))),
    ((u:$3FF104B648F113A1),(u:$3FF9EFDCA62B700A),(u:QWord($B920000000000000))),
    ((u:$3FF1585720F35CD9),(u:$3FFA5BF3ACFDE4B2),(u:$3960000000000000)),
    ((u:$3FFE9CC7ED2E1A7E),(u:$400BB0FF220D8EB5),(u:QWord($BCA0000000000000))),
    ((u:$40043180EA854696),(u:$40191F1122B6B63A),(u:$3990000000000000)),
    ((u:$401725811DCF6782),(u:$40645EA160DDC71F),(u:QWord($B9B0000000000000))),
    ((u:$4025AFD56F7D565B),(u:$40D8FF8E0CCEA7C0),(u:$3A50000000000000)),
    ((u:$402759A2AD4C4D56),(u:$40ECB62EEC26BD78),(u:QWord($BA30000000000000))),
    ((u:$4027FCE95EA5C653),(u:$40F3BF8009648DC0),(u:$3A70000000000000)),
    ((u:$403743D5609348AC),(u:$41F7A87A8BB7FA28),(u:QWord($BE90000000000000))),
    ((u:$404E07E71BFCF06F),(u:$45491EC4412C344F),(u:$3E70000000000000)),
    ((u:$4056474C604CC0D7),(u:$47E7A8F65AD009BD),(u:QWord($C130000000000000))),
    ((u:$40654CD1FEA7663A),(u:$4F3C90810D354618),(u:$4860000000000000)),
    ((u:$4072DA9E5E6AF0B0),(u:$5B127D6FE867D6F6),(u:$5480000000000000)),
    ((u:$407D6479EBA7C971),(u:$6A462A88613629B6),(u:QWord($E370000000000000))));

function CoshDatabase(x, f: Double): Double;
var
  a, b, m: Int32;
  ax: Double;
begin
  a := 0; b := 20;
  ax := Abs(x);
  m := (a + b) div 2;
  while a <= b do
  begin
    if cCoshDb[m, 0].f < ax then a := m + 1
    else if cCoshDb[m, 0].f = ax then
    begin
      f := cCoshDb[m, 1].f + cCoshDb[m, 2].f;
      Break;
    end
    else b := m - 1;
    m := (a + b) div 2;
  end;
  Result := f;
end;

function CoshZero(x: Double): Double;
var
  x2, x2l, y0, y1, y2: Double;
  chp, clp, ch_c, e_tmp: Double;
  i: Int32;
  t_u, w_u: Tb64u64;
begin
  x2  := x * x;
  x2l := pcr_fma(x, x, -x2);

  y2 := x2 * (cCoshZeroCl[0].f + x2 * (cCoshZeroCl[1].f
            + x2 * (cCoshZeroCl[2].f + x2 * cCoshZeroCl[3].f)));

  // Seeded polydd: start with ch = c[n-1][0] + *l; cl = ... + c[n-1][1]
  chp := cCoshZeroCh[3, 0].f + y2;
  clp := ((cCoshZeroCh[3, 0].f - chp) + y2) + cCoshZeroCh[3, 1].f;
  for i := 2 downto 0 do
  begin
    chp := pcr_muldd(x2, x2l, chp, clp, clp);
    ch_c := chp + cCoshZeroCh[i, 0].f;
    e_tmp := (cCoshZeroCh[i, 0].f - ch_c) + chp;
    chp := ch_c;
    clp := (clp + cCoshZeroCh[i, 1].f) + e_tmp;
  end;
  y1 := chp; y2 := clp;

  y1 := pcr_muldd(y1, y2, x2, x2l, y2);
  pcr_fasttwosum(y0, y1, Double(1.0), y1);
  pcr_fasttwosum(y1, y2, y1, y2);

  t_u.f := y1;
  if (t_u.u and UInt64($000FFFFFFFFFFFFF)) = 0 then
  begin
    w_u.f := y2;
    if ((w_u.u xor t_u.u) shr 63) <> 0 then t_u.u := t_u.u - 1
    else t_u.u := t_u.u + 1;
    y1 := t_u.f;
  end;
  if (t_u.u and UInt64($000FFFFFFFFFFFFF)) = UInt64($000FFFFFFFFFFFFF) then
    Result := CoshDatabase(x, y0 + y1)
  else
    Result := y0 + y1;
end;

function pcr_cosh(x: Double): Double;
var
  ix, jt, sp, sm, uh_u, ul_u, v: Tb64u64;
  aix: UInt64;
  ax, v0, t: Double;
  x2, x4, p, eps, lb, ub: Double;
  il, jl, i0_, i1_, j0_, j1_, ie_, je_, eh_, el_, ml_: Int64;
  t0h, t0l, t1h, t1l, th, tl: Double;
  q0h, q0l, q1h, q1l, qh, ql, pm, em, dx, dx2, mx, pp: Double;
  rh, rl: Double;
  fph, fpl, fmh, fml: Double;
begin
  ix.f := x;
  ax := Abs(x);
  v0 := pcr_fma(ax, cCoshS.f, cCoshMagic.f);
  jt.f := v0;
  v.u := jt.u and cCoshMask26;
  t := v.f - cCoshMagicH.f;

  aix := ix.u and UInt64($7FFFFFFFFFFFFFFF);

  // |x| < 0.125
  if aix < cCoshAxTiny then
  begin
    if aix < cCoshAxP26 then
    begin
      Result := pcr_fma(ax, cCoshTinyAdj.f, Double(1.0));
      Exit;
    end;
    x2 := x * x; x4 := x2 * x2;
    p := x2 * ((cCoshFastC0.f + x2 * cCoshFastC1.f)
             + x4 * ((cCoshFastC2.f + x2 * cCoshFastC3.f) + x4 * cCoshFastC4.f));
    eps := x2 * cCoshEps184.f;
    lb := Double(1.0) + (p - eps);
    ub := Double(1.0) + (p + eps);
    if lb = ub then begin Result := lb; Exit; end;
    Result := CoshZero(x);
    Exit;
  end;

  // Large / NaN / Inf / overflow
  if aix > cCoshAxHuge then
  begin
    if aix > UInt64($7FF0000000000000) then begin Result := x + x; Exit; end;
    if aix = UInt64($7FF0000000000000) then begin Result := ax; Exit; end;
    Result := cCoshOvf.f * Double(2.0);
    Exit;
  end;

  // 0.125 <= |x| <= ~711
  il := Int64((jt.u shl 14) shr 40);
  jl := -il;
  i1_ := il and $3F;
  i0_ := (il shr 6) and $3F;
  ie_ := SarInt64(il, 12);
  j1_ := jl and $3F;
  j0_ := SarInt64(jl, 6) and $3F;
  je_ := SarInt64(jl, 12);

  sp.u := UInt64(1022 + ie_) shl 52;
  sm.u := UInt64(1022 + je_) shl 52;

  t0h := cExpT0[i0_, 1].f; t0l := cExpT0[i0_, 0].f;
  t1h := cExpT1[i1_, 1].f; t1l := cExpT1[i1_, 0].f;
  th := t0h * t1h;
  tl := t0h * t1l + t1h * t0l + pcr_fma(t0h, t1h, -th);

  dx  := (ax - cTanhL2hA.f * t) + cTanhL2lA.f * t;
  dx2 := dx * dx;
  mx  := -dx;

  pp := dx * ((cCoshMainCh0.f + dx * cCoshMainCh1.f)
            + dx2 * (cCoshMainCh2.f + dx * cCoshMainCh3.f));

  if aix > cCoshAxMed then
  begin
    if aix > cCoshAxBig then
    begin
      sp.u := UInt64(1021 + ie_) shl 52;
      rh := th;
      rl := tl + th * pp;
      eps := cCoshErrLrg.f * th;
      lb := rh + (rl - eps);
      ub := rh + (rl + eps);
      if lb = ub then begin Result := (lb * sp.f) * Double(2.0); Exit; end;

      th := TanhExpAccurate(ax, t, th, tl, tl);
      pcr_fasttwosum(th, tl, th, tl);
      uh_u.f := th; ul_u.f := tl;
      eh_ := (uh_u.u shr 52) and $7FF;
      el_ := (ul_u.u shr 52) and $7FF;
      ml_ := Int64((ul_u.u + 8) and UInt64($000FFFFFFFFFFFFF));
      th := th + tl;
      th := th * Double(2.0);
      th := th * sp.f;
      if (ml_ <= 16) or ((eh_ - el_) > 103) then
        Result := CoshDatabase(x, th)
      else
        Result := th;
      Exit;
    end;
    // 5 < |x| <= 36.736801
    q0h := cExpT0[j0_, 1].f;
    q1h := cExpT1[j1_, 1].f;
    qh := q0h * q1h;
    th := th * sp.f;
    tl := tl * sp.f;
    qh := qh * sm.f;
    pm := mx * ((cCoshMainCh0.f + mx * cCoshMainCh1.f)
              + dx2 * (cCoshMainCh2.f + mx * cCoshMainCh3.f));
    em := qh + qh * pm;
    rh := th;
    rl := (tl + em) + th * pp;

    eps := cCoshErrMed.f * rh;
    lb := rh + (rl - eps);
    ub := rh + (rl + eps);
    if lb = ub then begin Result := lb; Exit; end;

    th := TanhExpAccurate(ax, t, th, tl, tl);
    if aix > cCoshAxBigBr then
    begin
      rh := th + qh;
      rl := ((th - rh) + qh) + tl;
    end
    else
    begin
      qh := q0h * q1h;
      q0l := cExpT0[j0_, 0].f;
      q1l := cExpT1[j1_, 0].f;
      ql := q0h * q1l + q1h * q0l + pcr_fma(q0h, q1h, -qh);
      qh := qh * sm.f;
      ql := ql * sm.f;
      qh := TanhExpAccurate(-ax, -t, qh, ql, ql);
      rh := th + qh;
      rl := (((th - rh) + qh) + ql) + tl;
    end;
  end
  else
  begin
    // 0.125 <= |x| <= 5
    q0h := cExpT0[j0_, 1].f; q0l := cExpT0[j0_, 0].f;
    q1h := cExpT1[j1_, 1].f; q1l := cExpT1[j1_, 0].f;
    qh := q0h * q1h;
    ql := q0h * q1l + q1h * q0l + pcr_fma(q0h, q1h, -qh);
    th := th * sp.f; tl := tl * sp.f;
    qh := qh * sm.f; ql := ql * sm.f;
    pm := mx * ((cCoshMainCh0.f + mx * cCoshMainCh1.f)
              + dx2 * (cCoshMainCh2.f + mx * cCoshMainCh3.f));
    fph := th; fpl := tl + th * pp;
    fmh := qh; fml := ql + qh * pm;

    rh := fph + fmh;
    rl := ((fph - rh) + fmh) + fml + fpl;
    eps := cCoshErrSml.f * rh;
    lb := rh + (rl - eps);
    ub := rh + (rl + eps);
    if lb = ub then begin Result := lb; Exit; end;

    th := TanhExpAccurate( ax,  t, th, tl, tl);
    qh := TanhExpAccurate(-ax, -t, qh, ql, ql);
    rh := th + qh;
    rl := ((th - rh) + qh) + ql + tl;
  end;

  pcr_fasttwosum(rh, rl, rh, rl);
  uh_u.f := rh; ul_u.f := rl;
  eh_ := (uh_u.u shr 52) and $7FF;
  el_ := (ul_u.u shr 52) and $7FF;
  ml_ := Int64((ul_u.u + 8) and UInt64($000FFFFFFFFFFFFF));
  rh := rh + rl;
  if (ml_ <= 16) or ((eh_ - el_) > 103) then
    Result := CoshDatabase(x, rh)
  else
    Result := rh;
end;
