// ===========================================================================
// pcr_sinh — correctly-rounded binary64 sinh(x).
// Ported from core-math/src/binary64/sinh/sinh.c by Alexei Sibidanov.
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
//
// Shares with previously ported exp-family units:
//   * cExpT0 / cExpT1 tables (byte-identical to sinh's t0[]/t1[])
//   * TanhExpAccurate (identical to sinh's as_exp_accurate)
//   * cTanhL2hA / cTanhL2lA (reduction constants)
//   * cCoshS / cCoshMagic / cCoshMagicH / cCoshMask26 (same magic constants)
//   * cCoshMainCh0..3 (sinh ch[4] is byte-identical to cosh's)
//   * cCoshAxMed / cCoshAxBig / cCoshAxBigBr / cCoshAxHuge / cCoshOvf /
//     cCoshTinyAdj (= 0x1p-55)
//   * cCoshErrMed (= 0x1.202p-63, sinh med eps) and
//     cCoshErrSml (= 0x1.c0ap-62, sinh small eps) are byte-identical reuse
//
// Sinh-specific content below:
//   * Fast-path c[5] polynomial for x0 <= |x| < 0.25
//   * as_sinh_zero ch[5][2] + cl[3] polynomial
//   * as_sinh_database 51-entry table
// ===========================================================================

const
  // Fast-path c[5] polynomial (x0 <= |x| < 0.25 branch)
  cSinhFastC0: Tb64u64 = (u:$3FC5555555555555); // 0x1.5555555555555p-3
  cSinhFastC1: Tb64u64 = (u:$3F8111111111151E); // 0x1.111111111151ep-7
  cSinhFastC2: Tb64u64 = (u:$3F2A01A019D0C767); // 0x1.a01a019d0c767p-13
  cSinhFastC3: Tb64u64 = (u:$3EC71DE444A96E11); // 0x1.71de444a96e11p-19
  cSinhFastC4: Tb64u64 = (u:$3E5AE84653752420); // 0x1.ae8465375242p-26

  // as_sinh_zero ch[5][2] — Pillar B: named scalars (Hi = [i,0], Lo = [i,1]).
  cSinhZeroChH_0: Tb64u64 = (u:$3FC5555555555555);
  cSinhZeroChL_0: Tb64u64 = (u:$3C6555555555552F);
  cSinhZeroChH_1: Tb64u64 = (u:$3F81111111111111);
  cSinhZeroChL_1: Tb64u64 = (u:$3C011111115CF00D);
  cSinhZeroChH_2: Tb64u64 = (u:$3F2A01A01A01A01A);
  cSinhZeroChL_2: Tb64u64 = (u:$3B6A0011C925B85C);
  cSinhZeroChH_3: Tb64u64 = (u:$3EC71DE3A556C734);
  cSinhZeroChL_3: Tb64u64 = (u:QWord($BB6B4E2835532BCD));
  cSinhZeroChH_4: Tb64u64 = (u:$3E5AE64567F54482);
  cSinhZeroChL_4: Tb64u64 = (u:QWord($BAEDEFCF17A6AB79));

  // as_sinh_zero y2 polynomial seed: cl[0..2]
  cSinhZeroCl0: Tb64u64 = (u:$3DE6124613AEF206); // 0x1.6124613aef206p-33
  cSinhZeroCl1: Tb64u64 = (u:$3D6AE7F36BEEA815); // 0x1.ae7f36beea815p-41
  cSinhZeroCl2: Tb64u64 = (u:$3CE95785063CD974); // 0x1.95785063cd974p-49

  // Bracket-epsilon constants
  cSinhFastEps: Tb64u64 = (u:$3CAC000000000000); // 0x1.cp-53 (fast-path)
  cSinhErrLrg:  Tb64u64 = (u:$3C01B60000000000); // 0x1.1b6p-63

  // Input-range thresholds
  cSinhAxTiny:   UInt64 = $3FD0000000000000; // 0.25
  cSinhAxTinier: UInt64 = $3E57137449123EF7; // 0x1.7137449123ef7p-26

  // 51-entry as_sinh_database
  cSinhDb: array[0..50, 0..2] of Tb64u64 = (
    ((u:$3FD364303E1AD8F6),(u:$3FD3B07E0C779DDA),(u:QWord($B95BC00000000000))),
    ((u:$3FD4169F234F23B9),(u:$3FD46B7B3B358F99),(u:QWord($BC70000000000000))),
    ((u:$3FD616CC75D49226),(u:$3FD687BD068C1C1E),(u:$390A000000000000)),
    ((u:$3FDAE3773250E7D2),(u:$3FDBAFC3479FC9CC),(u:QWord($B960000000000000))),
    ((u:$3FDB7EFA91915C95),(u:$3FDC59869F17B483),(u:QWord($B970000000000000))),
    ((u:$3FDD68039861AB53),(u:$3FDE73B46ABB01E1),(u:QWord($B922000000000000))),
    ((u:$3FDE90F16EB88C09),(u:$3FDFBDD4A37760B7),(u:QWord($B93F800000000000))),
    ((u:$3FEA3FC7E4DD47D1),(u:$3FED4B21EBF542F0),(u:$394E000000000000)),
    ((u:$3FEAA3B649A96091),(u:$3FEDD32C5ED1E930),(u:$3958A00000000000)),
    ((u:$3FEC13876341B62E),(u:$3FEFD1D7F1C8170C),(u:$3967200000000000)),
    ((u:$3FF2F5D3B178914A),(u:$3FF7B8516FFD2406),(u:QWord($B972800000000000))),
    ((u:$3FF3FFC12B81CBC2),(u:$3FF9A0FF413A1AF3),(u:$394C000000000000)),
    ((u:$3FF44F65DFF00782),(u:$3FFA38A3C3227609),(u:QWord($B980000000000000))),
    ((u:$3FF7346E3C591A14),(u:$40001E9CFA77B855),(u:$3990000000000000)),
    ((u:$3FFB6E2C73F41415),(u:$40057E377B3F0B4B),(u:QWord($B990000000000000))),
    ((u:$3FFDC5059D4E507D),(u:$4009168C60ED5256),(u:$397C600000000000)),
    ((u:$3FFF737F1E8378C7),(u:$400BFFD3F94F40FB),(u:$397A400000000000)),
    ((u:$4003359640329982),(u:$4015E40DF3F985BE),(u:$3999700000000000)),
    ((u:$40058A4FF5ADAC35),(u:$401D671928665BDD),(u:$3990000000000000)),
    ((u:$4008C0A26D055288),(u:$4026056B06A21918),(u:QWord($B99BE00000000000))),
    ((u:$400BC3C2D0C95F52),(u:$40300FEF7383A978),(u:$39B6100000000000)),
    ((u:$4020A19AEBB51E90),(u:$409FEE8F69C4CD25),(u:$3A04800000000000)),
    ((u:$4023EB8F61734227),(u:$40C4AB1CF45E4E26),(u:$3A50000000000000)),
    ((u:$40243A81752EABE7),(u:$40C81D364845ECFA),(u:QWord($BA50000000000000))),
    ((u:$40316369CD53BB69),(u:$4170FBC6C02B1C90),(u:QWord($BAE9000000000000))),
    ((u:$40320E29EA8B51E2),(u:$41808B8ABBA28ABC),(u:$3B09B00000000000)),
    ((u:$40392A5C27AFBE82),(u:$4223C81F9A247253),(u:$3BC0000000000000)),
    ((u:$403A1E4F11B513D7),(u:$4239A65B6C2E2185),(u:QWord($BB9BC00000000000))),
    ((u:$403C089FCF166171),(u:$4265C452E0E37569),(u:$3BA4000000000000)),
    ((u:$403E42A98B3A0BE5),(u:$429938768CA4F8AA),(u:$3C16D00000000000)),
    ((u:$40404DB52248CBB8),(u:$42D0794072349523),(u:$3C60E80000000000)),
    ((u:$40421BC021EEB97E),(u:$4323065064A170FB),(u:$3CB0880000000000)),
    ((u:$40439FC4D3BB7110),(u:$4368A4E90733B95E),(u:$3CD6E00000000000)),
    ((u:$4043C895D86E96C9),(u:$4370F33837882A60),(u:QWord($BCE2800000000000))),
    ((u:$404E07E71BFCF06F),(u:$45491EC4412C344F),(u:$3E70000000000000)),
    ((u:$404F7216C4B435C9),(u:$458A97E7BE23E65A),(u:QWord($BF00000000000000))),
    ((u:$4056474C604CC0D7),(u:$47E7A8F65AD009BD),(u:QWord($C130800000000000))),
    ((u:$4057A60EE15E3E9D),(u:$48662E4DC3BBF53F),(u:$41CB000000000000)),
    ((u:$4061F0DA93354198),(u:$4CD0BD73B73FC74C),(u:$4655900000000000)),
    ((u:$40654CD1FEA7663A),(u:$4F3C90810D354618),(u:$4862000000000000)),
    ((u:$406556C678D5E976),(u:$4F437E7AC4E7F9B3),(u:$48C0200000000000)),
    ((u:$4067945E34B18A9A),(u:$50E1B0E4936A8C9B),(u:QWord($CA5FA00000000000))),
    ((u:$4072DA9E5E6AF0B0),(u:$5B127D6FE867D6F6),(u:$5480A00000000000)),
    ((u:$40754CEBA01331D5),(u:$5E99A86785B5EF3E),(u:QWord($D812200000000000))),
    ((u:$4079E7B643238A14),(u:$653F5DA7FE652978),(u:$5EC0000000000000)),
    ((u:$407C7206C1B753E4),(u:$68E8670DE0B68CAD),(u:QWord($E237800000000000))),
    ((u:$407D6479EBA7C971),(u:$6A462A88613629B6),(u:QWord($E374000000000000))),
    ((u:$407EB9914D4AC1C8),(u:$6C32B67EFF65DCE8),(u:QWord($E5A0200000000000))),
    ((u:$4080BC04AF1B09F5),(u:$7027B1D97C902985),(u:$6995600000000000)),
    ((u:$40826EE1A46D8C8B),(u:$750FBE20477DF4A7),(u:QWord($EE85500000000000))),
    ((u:$4084A869881F72AC),(u:$7B79EA7540A3D1F9),(u:QWord($F4F2D00000000000)))
  );

function SinhDatabase(x, f: Double): Double;
var
  a, b, m: Int32;
  ax, sgn: Double;
begin
  a := 0; b := 50;
  ax := Abs(x);
  if x >= Double(0.0) then sgn := Double(1.0) else sgn := Double(-1.0);
  m := (a + b) div 2;
  while a <= b do
  begin
    if cSinhDb[m, 0].f < ax then a := m + 1
    else if cSinhDb[m, 0].f = ax then
    begin
      f := sgn * cSinhDb[m, 1].f + sgn * cSinhDb[m, 2].f;
      Break;
    end
    else b := m - 1;
    m := (a + b) div 2;
  end;
  Result := f;
end;

function SinhZero(x: Double): Double;
var
  x2, x2l, y0, y1, y2: Double;
  chp, clp, ch_c, e_tmp: Double;
  t_u, w_u: Tb64u64;
begin
  x2  := x * x;
  x2l := pcr_fma(x, x, -x2);

  y2 := x2 * (cSinhZeroCl0.f + x2 * (cSinhZeroCl1.f + x2 * cSinhZeroCl2.f));

  // Seeded polydd over 5 coefficients — Pillar A unroll of `for i := 3 downto 0`.
  chp := cSinhZeroChH_4.f + y2;
  clp := ((cSinhZeroChH_4.f - chp) + y2) + cSinhZeroChL_4.f;

  chp := pcr_muldd(x2, x2l, chp, clp, clp);
  ch_c := chp + cSinhZeroChH_3.f;
  e_tmp := (cSinhZeroChH_3.f - ch_c) + chp;
  chp := ch_c;
  clp := (clp + cSinhZeroChL_3.f) + e_tmp;

  chp := pcr_muldd(x2, x2l, chp, clp, clp);
  ch_c := chp + cSinhZeroChH_2.f;
  e_tmp := (cSinhZeroChH_2.f - ch_c) + chp;
  chp := ch_c;
  clp := (clp + cSinhZeroChL_2.f) + e_tmp;

  chp := pcr_muldd(x2, x2l, chp, clp, clp);
  ch_c := chp + cSinhZeroChH_1.f;
  e_tmp := (cSinhZeroChH_1.f - ch_c) + chp;
  chp := ch_c;
  clp := (clp + cSinhZeroChL_1.f) + e_tmp;

  chp := pcr_muldd(x2, x2l, chp, clp, clp);
  ch_c := chp + cSinhZeroChH_0.f;
  e_tmp := (cSinhZeroChH_0.f - ch_c) + chp;
  chp := ch_c;
  clp := (clp + cSinhZeroChL_0.f) + e_tmp;
  y1 := chp; y2 := clp;

  y1 := pcr_mulddd_pd(y1, y2, x, y2);
  y1 := pcr_muldd(y1, y2, x2, x2l, y2);
  pcr_fasttwosum(y0, y1, x, y1);
  pcr_fasttwosum(y1, y2, y1, y2);

  t_u.f := y1;
  if (t_u.u and UInt64($000FFFFFFFFFFFFF)) = 0 then
  begin
    w_u.f := y2;
    if ((w_u.u xor t_u.u) shr 63) <> 0 then t_u.u := t_u.u - 1
    else t_u.u := t_u.u + 1;
    y1 := t_u.f;
  end;
  Result := y0 + y1;
end;

function pcr_sinh(x: Double): Double;
var
  ix, jt, sp, sm, uh_u, ul_u, v: Tb64u64;
  aix: UInt64;
  ax, sgn, v0, t: Double;
  x2, x3, x4, p, eps, lb, ub: Double;
  il, jl, i0_, i1_, j0_, j1_, ie_, je_, eh_, el_, ml_: Int64;
  t0h, t0l, t1h, t1l, th, tl: Double;
  q0h, q0l, q1h, q1l, qh, ql, pm, em, dx, dx2, mx, pp: Double;
  rh, rl: Double;
  fph, fpl, fmh, fml: Double;
begin
  ix.f := x;
  aix := ix.u and UInt64($7FFFFFFFFFFFFFFF);
  ax := Abs(x);

  // Fast-path: |x| < 0.25
  if aix < cSinhAxTiny then
  begin
    if aix < cSinhAxTinier then
    begin
      // |x| < x0 = 0x1.7137449123ef7p-26: sinh(x) rounds to x.
      // C uses fma(x, 2^-55, x); emulated pcr_fma may mis-round subnormals,
      // so return x directly (bit-exact in round-to-nearest; loses INEXACT).
      Result := x;
      Exit;
    end;
    x2 := x * x; x3 := x2 * x; x4 := x2 * x2;
    p := x3 * ((cSinhFastC0.f + x2 * cSinhFastC1.f)
             + x4 * ((cSinhFastC2.f + x2 * cSinhFastC3.f) + x4 * cSinhFastC4.f));
    eps := x3 * cSinhFastEps.f;
    lb := x + (p - eps);
    ub := x + (p + eps);
    if lb = ub then begin Result := lb; Exit; end;
    Result := SinhZero(x);
    Exit;
  end;

  // Overflow / NaN / Inf
  if aix > cCoshAxHuge then
  begin
    if aix >= UInt64($7FF0000000000000) then begin Result := x + x; Exit; end;
    if x >= Double(0.0) then Result :=  cCoshOvf.f * Double(2.0)
                        else Result := -cCoshOvf.f * Double(2.0);
    Exit;
  end;

  // Reduction magic (same as cosh exp-family)
  v0 := pcr_fma(ax, cCoshS.f, cCoshMagic.f);
  jt.f := v0;
  v.u := jt.u and cCoshMask26;
  t := v.f - cCoshMagicH.f;

  // 0.25 <= |x| <= ~710.47586
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

  if x >= Double(0.0) then sgn := Double(1.0) else sgn := Double(-1.0);

  if aix > cCoshAxMed then
  begin
    if aix > cCoshAxBig then
    begin
      // |x| > ~36.736801: exp(-|x|) is negligible, compute exp(|x|) only.
      sp.u := UInt64(1021 + ie_) shl 52;
      rh := th;
      rl := tl + th * pp;
      rh := rh * sgn;
      rl := rl * sgn;
      eps := cSinhErrLrg.f * th;
      lb := rh + (rl - eps);
      ub := rh + (rl + eps);
      if lb = ub then begin Result := (lb * sp.f) * Double(2.0); Exit; end;

      th := TanhExpAccurate(ax, t, th, tl, tl);
      pcr_fasttwosum(th, tl, th, tl);
      th := th * sgn;
      tl := tl * sgn;
      uh_u.f := th; ul_u.f := tl;
      eh_ := (uh_u.u shr 52) and $7FF;
      el_ := (ul_u.u shr 52) and $7FF;
      ml_ := Int64((ul_u.u + 8) and UInt64($000FFFFFFFFFFFFF));
      th := th + tl;
      th := th * Double(2.0);
      th := th * sp.f;
      if (ml_ <= 16) or ((eh_ - el_) > 103) then
        Result := SinhDatabase(x, th)
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
    rl := (tl - em) + th * pp;

    rh := rh * sgn;
    rl := rl * sgn;
    eps := cCoshErrMed.f * rh;
    lb := rh + (rl - eps);
    ub := rh + (rl + eps);
    if lb = ub then begin Result := lb; Exit; end;

    th := TanhExpAccurate(ax, t, th, tl, tl);
    if aix > cCoshAxBigBr then
    begin
      // |x| > 31.4: keep the unscaled em approximation for exp(-|x|)
      rh := th - qh;
      rl := ((th - rh) - qh) + tl;
    end
    else
    begin
      // 5 < |x| <= 31.4: refine qh/ql as a double-double
      qh := q0h * q1h;
      q0l := cExpT0[j0_, 0].f;
      q1l := cExpT1[j1_, 0].f;
      ql := q0h * q1l + q1h * q0l + pcr_fma(q0h, q1h, -qh);
      qh := qh * sm.f;
      ql := ql * sm.f;
      qh := TanhExpAccurate(-ax, -t, qh, ql, ql);
      rh := th - qh;
      rl := (((th - rh) - qh) - ql) + tl;
    end;
  end
  else
  begin
    // 0.25 <= |x| <= 5: both exp(±|x|) as double-doubles
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

    rh := fph - fmh;
    rl := ((fph - rh) - fmh) - fml + fpl;
    rh := rh * sgn;
    rl := rl * sgn;
    eps := cCoshErrSml.f * rh;
    lb := rh + (rl - eps);
    ub := rh + (rl + eps);
    if lb = ub then begin Result := lb; Exit; end;

    th := TanhExpAccurate( ax,  t, th, tl, tl);
    qh := TanhExpAccurate(-ax, -t, qh, ql, ql);
    rh := th - qh;
    rl := ((th - rh) - qh) - ql + tl;
  end;

  // Common merge path (both remaining branches)
  pcr_fasttwosum(rh, rl, rh, rl);
  uh_u.f := rh; ul_u.f := rl;
  eh_ := (uh_u.u shr 52) and $7FF;
  el_ := (ul_u.u shr 52) and $7FF;
  ml_ := Int64((ul_u.u + 8) and UInt64($000FFFFFFFFFFFFF));
  rh := rh * sgn;
  rl := rl * sgn;
  rh := rh + rl;
  if (ml_ <= 16) or ((eh_ - el_) > 103) then
    Result := SinhDatabase(x, rh)
  else
    Result := rh;
end;
