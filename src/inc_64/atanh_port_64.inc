// ===========================================================================
// pcr_atanh — correctly-rounded binary64 atanh(x).
// Ported from core-math/src/binary64/atanh/atanh.c by Alexei Sibidanov.
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
// Reuses main-branch tables (B, R1, R2, L1, L2) and refine tables
// (T1..T4, LL, RefineCh, RefineCl, L20/L21/L22, VOffset) from acosh_port.inc.
// ===========================================================================

const
  // Main-branch l1/l2 tables. Values are -log(r1)/2 and -log(r2)/2
  // (half of acosh's tables), so we cannot share them.
  cAtanhL1: array[0..32, 0..1] of Tb64u64 = (
    ((u:$0000000000000000),(u:$0000000000000000)),
    ((u:QWord($BE4532C1269E2038)),(u:$3F862E5000000000)),
    ((u:$3E4CE42D81B54E84),(u:$3F962E3C00000000)),
    ((u:QWord($BE525826F815EC3D)),(u:$3FA0A2AC00000000)),
    ((u:$3E50DB1B1E7CEE11),(u:$3FA62E4A00000000)),
    ((u:QWord($BE51F3A8C6C95003)),(u:$3FABB9DC00000000)),
    ((u:QWord($BE5774CD4FB8C30D)),(u:$3FB0A2B200000000)),
    ((u:$3E2452E56C030A0A),(u:$3FB3687F00000000)),
    ((u:$3E36B63C4966A79A),(u:$3FB62E4100000000)),
    ((u:QWord($BE3B20A21CCB525E)),(u:$3FB8F40A00000000)),
    ((u:$3E54006CFB3D8F85),(u:$3FBBB9D100000000)),
    ((u:QWord($BE5CDB026B310C41)),(u:$3FBE7F9B00000000)),
    ((u:QWord($BE569124FDC0F16D)),(u:$3FC0A2B080000000)),
    ((u:QWord($BE5084656CDC2727)),(u:$3FC2059580000000)),
    ((u:QWord($BE5376FA8B0357FD)),(u:$3FC3687C00000000)),
    ((u:$3E3E56AE55A47B4A),(u:$3FC4CB5E80000000)),
    ((u:$3E5070FF8834EEB4),(u:$3FC62E4400000000)),
    ((u:$3E5623516109F4FE),(u:$3FC7912900000000)),
    ((u:QWord($BE2EC656B95FBDAC)),(u:$3FC8F40B00000000)),
    ((u:$3E3F0CA2E729F510),(u:$3FCA56ED80000000)),
    ((u:QWord($BE57D260A858354A)),(u:$3FCBB9D680000000)),
    ((u:$3E4E7279075503D3),(u:$3FCD1CB900000000)),
    ((u:$3E439E1A0A503873),(u:$3FCE7F9D00000000)),
    ((u:$3E5CD86D7B87C3D6),(u:$3FCFE27D80000000)),
    ((u:$3E5060AB88DE341E),(u:$3FD0A2B240000000)),
    ((u:$3E320A860D3F9390),(u:$3FD1542440000000)),
    ((u:QWord($BE4DACEE95FC2F10)),(u:$3FD2059740000000)),
    ((u:$3E545DE3A86E0ACA),(u:$3FD2B70700000000)),
    ((u:$3E4C164CBFB991AF),(u:$3FD3687B00000000)),
    ((u:$3E5D3F66B24225EF),(u:$3FD419EC40000000)),
    ((u:$3E5FC023EFA144BA),(u:$3FD4CB5F80000000)),
    ((u:$3E3086A8AF6F26C0),(u:$3FD57CD280000000)),
    ((u:QWord($BE105C610CA86C39)),(u:$3FD62E4300000000)));

  cAtanhL2: array[0..32, 0..1] of Tb64u64 = (
    ((u:$0000000000000000),(u:$0000000000000000)),
    ((u:QWord($BE337E152A129E4E)),(u:$3F36320000000000)),
    ((u:QWord($BE53F6C916B8BE9C)),(u:$3F46300000000000)),
    ((u:$3E520505936739D5),(u:$3F50A24000000000)),
    ((u:QWord($BE523E2E8CB541BA)),(u:$3F562DC000000000)),
    ((u:QWord($BDFACB7983AC4F5E)),(u:$3F5BBA0000000000)),
    ((u:$3E36F7C7689C63AE),(u:$3F60A2A000000000)),
    ((u:$3E1F5CA695B4C58B),(u:$3F6368C000000000)),
    ((u:QWord($BE4C6C18BD953226)),(u:$3F662E6000000000)),
    ((u:$3E57A516C34846BD),(u:$3F68F46000000000)),
    ((u:QWord($BE4F3B83DD8B8530)),(u:$3F6BBA0000000000)),
    ((u:QWord($BE0C3459046E4E57)),(u:$3F6E800000000000)),
    ((u:$3D9B5C7E34CB79F6),(u:$3F70A2C000000000)),
    ((u:QWord($BE42487E9AF9A692)),(u:$3F7205C000000000)),
    ((u:$3E5F21BBC4AD79CE),(u:$3F73687000000000)),
    ((u:QWord($BE2550FFC857B731)),(u:$3F74CB7000000000)),
    ((u:$3E487458EC1B7B34),(u:$3F762E2000000000)),
    ((u:$3E5103D4FE83EE81),(u:$3F77911000000000)),
    ((u:$3E4810483D3B398C),(u:$3F78F44000000000)),
    ((u:QWord($BE42085CB340608E)),(u:$3F7A573000000000)),
    ((u:$3E512698A119C42F),(u:$3F7BB9D000000000)),
    ((u:QWord($BE5EDB8C172B4C33)),(u:$3F7D1CC000000000)),
    ((u:QWord($BE58B55B87A5E238)),(u:$3F7E7FE000000000)),
    ((u:$3E5BE5E17763F78A),(u:$3F7FE2B000000000)),
    ((u:QWord($BE1C2D496790073E)),(u:$3F80A2A800000000)),
    ((u:$3E56542F523ABEEC),(u:$3F81541000000000)),
    ((u:QWord($BE5B7FDBE5B193F8)),(u:$3F8205A000000000)),
    ((u:$3E5FA4D42FE30C7C),(u:$3F82B70000000000)),
    ((u:$3E50D46AD04ADC86),(u:$3F83688800000000)),
    ((u:QWord($BE51C22D02D17C4C)),(u:$3F8419F000000000)),
    ((u:$3E1A7D1E330DCCCE),(u:$3F84CB7000000000)),
    ((u:$3E0187025E656BA3),(u:$3F857CD000000000)),
    ((u:QWord($BE4532C1269E2038)),(u:$3F862E5000000000)));

  // Main-branch 5-coefficient polynomial. Differs from acosh only in c[0].
  cAtanhC: array[0..4] of Tb64u64 = (
    (u:QWord($BFF0000000000000)), // -0x1p+0
    (u:$3FF5555555555530), //  0x1.555555555553p+0
    (u:QWord($BFFFFFFFFFFFFFA0)), // -0x1.fffffffffffap+0
    (u:$40099999E33A6366), //  0x1.99999e33a6366p+1
    (u:QWord($C01555559EF9525F)));// -0x1.555559ef9525fp+2

  // as_atanh_zero polynomial ch[13][2]
  cAtanhZeroCh: array[0..12, 0..1] of Tb64u64 = (
    ((u:$3FD5555555555555),(u:$3C75555555555555)),  // 3   :  0x1.5555555555555p-2 ,  0x1.5555555555555p-56
    ((u:$3FC999999999999A),(u:QWord($BC6999999999611C))),  // 5
    ((u:$3FC2492492492492),(u:$3C72492490F76B25)),  // 7
    ((u:$3FBC71C71C71C71C),(u:$3C5C71CD5C38A112)),  // 9
    ((u:$3FB745D1745D1746),(u:QWord($BC67556C4165F4CA))),  // 11
    ((u:$3FB3B13B13B13B14),(u:QWord($BC5B893C3B36052E))),  // 13
    ((u:$3FB1111111111105),(u:$3C54E1AFD723ED1F)),  // 15
    ((u:$3FAE1E1E1E1E2678),(u:QWord($BC5F86EA96FB1435))),  // 17
    ((u:$3FAAF286BC9F90CC),(u:$3C41E51A6E54FDE9)),  // 19
    ((u:$3FA8618618C779B6),(u:QWord($BC3AB913DE95C3BF))),  // 21
    ((u:$3FA642C84AA383EB),(u:$3C5632E747641B12)),  // 23
    ((u:$3FA47AE2D205013C),(u:QWord($BC40C9617E7BCFF2))),  // 25
    ((u:$3FA2F664D60473F9),(u:$3C33ADB3E2B7F35E))); // 27

  // as_atanh_zero polynomial cl[5] (plain Horner, tail)
  cAtanhZeroCl: array[0..4] of Tb64u64 = (
    (u:$3FA1A9A91FD692AF), //  0x1.1a9a91fd692afp-5   // 29
    (u:$3FA06DFBB35E7F44), //  0x1.06dfbb35e7f44p-5   // 31
    (u:$3FA037BED4D7588F), //  0x1.037bed4d7588fp-5   // 33
    (u:$3F95ACA6D6D720D6), //  0x1.5aca6d6d720d6p-6   // 35
    (u:$3FA99EA5700D53A5));//  0x1.99ea5700d53a5p-5   // 37

  // Database (only 2 entries, sign-flipped on negative x)
  cAtanhDb: array[0..1, 0..2] of Tb64u64 = (
    ((u:$3FD2DBB7B1C91363),(u:$3FD36F33D51C264D),(u:$3C70000000000000)),
    ((u:$3FDC493DC899E4A5),(u:$3FDE611AA58AB608),(u:QWord($BC70000000000000))));

  // Main-fast-path constants
  cAtanhFastL2h: Tb64u64 = (u:$3FD62E42FEFA3A00); //  0x1.62e42fefa3ap-2  = log(2)/2
  cAtanhFastL2l: Tb64u64 = (u:QWord($BCD0CA86C3898D00)); // -0x1.0ca86c3898dp-50

  // Thresholds / constants
  cAtanhBr_Small:      UInt64 = $3FD0000000000000; // 0x1p-2 — |x| < 0x1p-2 uses small polynomial path
  cAtanhBr_Tiny:       UInt64 = $3E4D12ED0AF1A27F; // 0x1.d12ed0af1a27fp-27
  cAtanhBr_Abs1:       UInt64 = $3FF0000000000000; // 0x1p+0
  cAtanhBr_Inf:        UInt64 = $7FF0000000000000;
  cAtanhTinyScale:     Tb64u64 = (u:$3C80000000000000); // 0x1p-55
  cAtanhMainEps_Base:  Double  = 38e-24;
  cAtanhMainEps_Dx2:   Tb64u64 = (u:$3CE0000000000000); // 0x1p-49
  cAtanhInvLog2Times2: Tb64u64 = (u:$40071547652B82FE); // 0x1.71547652b82fep+1

function AtanhDatabase(x, f: Double): Double;
var
  ax, sgn: Double;
  a, b, m: Int32;
begin
  ax := pcr_fabs(x);
  sgn := pcr_copysign(Double(1.0), x);
  a := 0; b := High(cAtanhDb);
  m := (a + b) div 2;
  while a <= b do
  begin
    if cAtanhDb[m, 0].f < ax then a := m + 1
    else if cAtanhDb[m, 0].f = ax then
    begin
      f := sgn * cAtanhDb[m, 1].f + sgn * cAtanhDb[m, 2].f;
      Break;
    end
    else b := m - 1;
    m := (a + b) div 2;
  end;
  Result := f;
end;

// Local muldd_acc (matches C operand grouping: ahhl + (alhh + ahlh))
function AtanhMuldd(xh, xl, ch, cl: Double; out l: Double): Double; inline;
var ahlh, alhh, ahhh, ahhl, s: Double;
begin
  ahlh := ch * xl;
  alhh := cl * xh;
  ahhh := ch * xh;
  ahhl := pcr_fma(ch, xh, -ahhh);
  ahhl := ahhl + (alhh + ahlh);
  s := ahhh + ahhl;            // fasttwosum
  l := ahhl - (s - ahhh);
  Result := s;
end;

// as_atanh_zero: accurate path for |x| in [2^-26 .. 0.25]
function AtanhZero(x: Double): Double;
var
  x2, x2l, y0, y1, y2: Double;
  ch, cl, th, tl, hh: Double;
  t, w: Tb64u64;
  i: Int32;
begin
  x2  := x * x;
  x2l := pcr_fma(x, x, -x2);
  y2  := x2 * (cAtanhZeroCl[0].f + x2 * (cAtanhZeroCl[1].f
      + x2 * (cAtanhZeroCl[2].f + x2 * (cAtanhZeroCl[3].f + x2 * cAtanhZeroCl[4].f))));

  // Seeded polydd: n=13, seed *l = y2
  ch := cAtanhZeroCh[12, 0].f + y2;
  cl := ((cAtanhZeroCh[12, 0].f - ch) + y2) + cAtanhZeroCh[12, 1].f;
  for i := 11 downto 0 do
  begin
    ch := AtanhMuldd(x2, x2l, ch, cl, cl);
    th := ch + cAtanhZeroCh[i, 0].f;
    tl := (cAtanhZeroCh[i, 0].f - th) + ch;
    ch := th;
    cl := cl + tl + cAtanhZeroCh[i, 1].f;
  end;
  y2 := cl;
  y1 := ch;

  // Non-normalized mulddd(y1, y2, x, &y2):
  //   hh = y1*x, *y2 = fma(x, y1, -hh) + y2*x; return hh
  hh := y1 * x;
  y2 := pcr_fma(x, y1, -hh) + y2 * x;
  y1 := hh;
  y1 := AtanhMuldd(x2, x2l, y1, y2, y2);
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

// as_atanh_refine: accurate path for 0.25 <= |x| < 1 when fast path fails
function AtanhRefine(x, zh, zl, a: Double): Double;
var
  t, v, t0, w: Tb64u64;
  ex, e, i1, i2, i3, i4: Int32;
  iv: UInt64;
  ed, el0, el1, el2: Double;
  L0, L1, L2: Double;
  t12, t34, th, tl, dh, dl, shm, slm: Double;
  xh, xl: Double;
  xh2, sh2, sl2, sh_: Double;
  v0, v1, v2, vtmp: Double;
  er, de: UInt64;
  res, sgn: Double;
begin
  t.f := zh;
  ex := (t.u shr 52) and $7FF;
  e  := ex - $3FF;
  t.u := t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12);
  t.u := t.u or (UInt64($3FF) shl 52);
  ed  := e;

  v.f := a - ed + cAcoshVOffset.f;
  iv  := (v.u - (UInt64($3FF) shl 52)) shr (52 - 16);
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
  // zl' = zl * 2^-e : subtract e from zl's exponent
  t.f := zl;
  t.u := UInt64(Int64(t.u) - (Int64(e) shl 52));
  xl := xl + th * t.f;
  xh := AcoshAddDD(xh, xl, shm, slm, xl);

  sl2 := xh * (cAcoshRefineCl[0].f + xh * (cAcoshRefineCl[1].f + xh * cAcoshRefineCl[2].f));
  // Seeded polydd: n=3, seed *l = sl2
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
  pcr_fasttwosum(v0, vtmp, L0, sh_);
  pcr_fasttwosum(v1, v2, vtmp, sl2);

  t.f := v1;
  if (t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) = 0 then
  begin
    w.f := v2;
    if ((w.u xor t.u) shr 63) <> 0 then t.u := t.u - 1
    else t.u := t.u + 1;
    v1 := t.f;
  end;

  t0.f := v0;
  er := ((t.u + 1) and (UInt64($FFFFFFFFFFFFFFFF) shr 12));
  de := ((t0.u shr 52) and $7FF) - ((t.u shr 52) and $7FF);
  sgn := pcr_copysign(Double(1.0), x);
  v0 := v0 * sgn;
  v1 := v1 * sgn;
  res := v0 + v1;
  if (de > 104) or (er < 3) then
    Result := AtanhDatabase(x, res)
  else
    Result := res;
end;

function pcr_atanh(x: Double): Double;
var
  ax: Double;
  ix, t: Tb64u64;
  aix: UInt64;
  ph, pl, qh, ql, qzz, iqh, th, tl: Double;
  eps, lb, ub: Double;
  ex, e, i1, i2: Int32;
  iuint, jshift: UInt64;
  dsigned: Int64;
  ed, r, dx, dx2, f: Double;
  lh, ll, sgn: Double;
begin
  ax := pcr_fabs(x);
  ix.f := ax;
  aix := ix.u;

  // |x| >= 1
  if aix >= cAtanhBr_Abs1 then
  begin
    if aix = cAtanhBr_Abs1 then
    begin
      // pole: ±inf
      Result := pcr_copysign(Double(1.0), x) / Double(0.0);
      Exit;
    end;
    if aix > cAtanhBr_Inf then
    begin
      Result := x + x; // NaN propagation
      Exit;
    end;
    Result := pcr_sqrt(Double(-1.0)); // |x| > 1 → NaN
    Exit;
  end;

  // |x| < 0x1p-2
  if aix < cAtanhBr_Small then
  begin
    // Very tiny: atanh(x) rounds to x
    if aix < cAtanhBr_Tiny then
    begin
      Result := pcr_fma(x, cAtanhTinyScale.f, x);
      Exit;
    end;
    // NOTE: the upstream C fast-path here (polynomial + fasttwosum/muldd plus
    // eps-based certification) relies on operation grouping and FP-contract
    // semantics that FPC does not reproduce bit-exactly. AtanhZero is correctly
    // rounded on its own, so we skip the fast path entirely.
    Result := AtanhZero(x);
    Exit;
  end;

  // 0.25 <= |x| < 1 : main branch
  pcr_fasttwosum(ph, pl, Double(1.0), ax);
  qh  := Double(1.0) - ax;
  qzz := Double(1.0) - qh;        // fasttwosub: z = x - s
  ql  := qzz - ax;                 // e = z - y
  iqh := Double(1.0) / qh;
  th  := ph * iqh;
  tl  := pcr_fma(ph, iqh, -th)
       + (pl + ph * (pcr_fma(-qh, iqh, Double(1.0)) - ql * iqh)) * iqh;

  t.f := th;
  ex  := (t.u shr 52) and $7FF;
  e   := ex - $3FF;
  t.u := t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12);
  ed  := e;
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
  r := (Double(0.5) * cAcoshR1[i1].f) * cAcoshR2[i2].f;
  dx := pcr_fma(r, t.f, -Double(0.5));
  dx2 := dx * dx;
  f := dx2 * ((cAtanhC[0].f + dx * cAtanhC[1].f)
            + dx2 * (cAtanhC[2].f + dx * cAtanhC[3].f + dx2 * cAtanhC[4].f));
  // rx/dxl for atanh: rx = r*t.f, dxl = fma(r,t.f,-rx)
  lh := (cAtanhL1[i1, 1].f + cAtanhL2[i2, 1].f) + cAtanhFastL2h.f * ed;
  pcr_fasttwosum(lh, ll, lh, r * t.f - Double(0.5));
  ll := ll + cAtanhFastL2l.f * ed
          + (cAtanhL1[i1, 0].f + cAtanhL2[i2, 0].f)
          + pcr_fma(r, t.f, -(r * t.f))
          + Double(0.5) * tl / th;
  ll := ll + f;
  sgn := pcr_copysign(Double(1.0), x);
  lh := lh * sgn;
  ll := ll * sgn;
  eps := cAtanhMainEps_Base + dx2 * cAtanhMainEps_Dx2.f;
  lb := lh + (ll - eps);
  ub := lh + (ll + eps);
  if lb = ub then begin Result := lb; Exit; end;
  pcr_fasttwosum(th, tl, th, tl);
  Result := AtanhRefine(x, th, tl, cAtanhInvLog2Times2.f * pcr_fabs(lh + ll));
end;
