// ===========================================================================
// pcr_sinpi — correctly-rounded binary64 sin(pi*x).
// Ported from core-math/src/binary64/sinpi/sinpi.c by Alexei Sibidanov.
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
//
// Shares with cospi_port (already in pascoremath64.pas):
//   * Sn/Sm/Cm tables via CospiSincosN / CospiSincosN2
//   * Polynomial tables cSinpiRefSh / cSinpiRefCh and their scalars
//   * sn[]/cn[] coefficients — sinpi and cospi's fast-path sn/cn match
//     (verified against cospi.c:141-142 and sinpi.c:155-156).
// Sinpi-specific (not shareable):
//   * as_sinpi_zero 4+3 coefficient polynomial (cospi's is 2+2)
//   * Exception database (4 entries, different from cospi's 8)
//   * Entry-point logic (different symmetries near integer/half-integer)
// ===========================================================================

const
  // Pillar B (Phase 6.3): cSinpiFastC lifted to named scalars.
  cSinpiFastC_0: Tb64u64 = (u:QWord($C014ABBCE625BE51)); //  -0x1.4abbce625be51p+2
  cSinpiFastC_1: Tb64u64 = (u:$400466BC67754B46);        //   0x1.466bc67754b46p+1
  cSinpiFastC_2: Tb64u64 = (u:QWord($BFE32D2CC12A51F4)); //  -0x1.32d2cc12a51f4p-1
  cSinpiFastC_3: Tb64u64 = (u:$3FB5060540058476);        //   0x1.5060540058476p-4

  // ph/pl = 2*pi double-double (used by |x| <= 2^-5 path)
  cSinpiPh: Tb64u64 = (u: $400921FB54442D18); // 0x1.921fb54442d18p+1
  cSinpiPl: Tb64u64 = (u: $3CA1A62633145C07); // 0x1.1a62633145c07p-53

  // Thresholds
  cSinpiAxT5:    UInt64 = $3FA2000000000000;  // |x| <= 0x1.2p-5
  cSinpiAxPm54:  UInt64 = $3C90000000000000;  // 2^-54

  cSinpi2p106:  Tb64u64 = (u: $4690000000000000); // 0x1p+106
  cSinpi2pm106: Tb64u64 = (u: $3950000000000000); // 0x1p-106

  // Fast-path error and polynomial mixers
  cSinpiFastEps1: Tb64u64 = (u: $3D00000000000000); // 0x1p-47
  cSinpiFastEps2: Tb64u64 = (u: $3990000000000000); // 0x1p-102
  cSinpiFastEr:   Tb64u64 = (u: $3C244A9A66CDB0D6); // 5.5e-19

  // as_sinpi_zero polynomial tables (sinpi.c:89-92). Pillar B (Phase 6.3):
  // arrays lifted to named scalars; the polydd loop (i = 2 downto 0) is
  // unrolled at the call site.
  cSinpiZeroCh_0_0: Tb64u64 = (u:QWord($C014ABBCE625BE53)); // -0x1.4abbce625be53p+2
  cSinpiZeroCh_0_1: Tb64u64 = (u:$3CB05511C68477BE);        //  0x1.05511c68477bep-52
  cSinpiZeroCh_1_0: Tb64u64 = (u:$400466BC6775AAE2);        //  0x1.466bc6775aae2p+1
  cSinpiZeroCh_1_1: Tb64u64 = (u:QWord($BC96DC0CBEFAE1DA)); // -0x1.6dc0cbefae1dap-54
  cSinpiZeroCh_2_0: Tb64u64 = (u:QWord($BFE32D2CCE62BD86)); // -0x1.32d2cce62bd86p-1
  cSinpiZeroCh_2_1: Tb64u64 = (u:$3C8066BD54973829);        //  0x1.066bd54973829p-55
  cSinpiZeroCh_3_0: Tb64u64 = (u:$3FB50783487EE781);        //  0x1.50783487ee781p-4
  cSinpiZeroCh_3_1: Tb64u64 = (u:$3C5832989F39A743);        //  0x1.832989f39a743p-58
  cSinpiZeroCl_0:   Tb64u64 = (u:QWord($BF7E3074FDE861F0)); // -0x1.e3074fde861fp-8
  cSinpiZeroCl_1:   Tb64u64 = (u:$3F3E8F4344534DA6);        //  0x1.e8f4344534da6p-12
  cSinpiZeroCl_2:   Tb64u64 = (u:QWord($BEF6F9CD7B8CB9DB)); // -0x1.6f9cd7b8cb9dbp-16

  // as_sinpi_zero reduction constants (pi split into pi0 + pi1 + pi2, sinpi.c:96)
  cSinpiZeroPi0: Tb64u64 = (u: $4009200000000000); // 0x1.92p+1  (6.3125 exactly)
  cSinpiZeroPi1: Tb64u64 = (u: $3F4FB54442D1846A); // 0x1.fb54442d1846ap-11
  cSinpiZeroPi2: Tb64u64 = (u: QWord($BBED9CCEBA3F91F2)); // -0x1.d9cceba3f91f2p-65

  // Sinpi-specific exception database (sinpi.c:140-145) — 4 entries only.
  // Pillar B (Phase 6.3): arrays lifted to named scalars (Iq stays plain int
  // const since it's not a Tb64u64).
  cSinpiDb4_Iq_0: Int32 = 76;
  cSinpiDb4_Iq_1: Int32 = 108;
  cSinpiDb4_Iq_2: Int32 = 490;
  cSinpiDb4_Iq_3: Int32 = 903;
  cSinpiDb4_X_0: Tb64u64 = (u:QWord($BFD276B3FEF46600)); // -0x1.276b3fef466p-2
  cSinpiDb4_X_1: Tb64u64 = (u:QWord($BFD33CAEA0F24C00)); // -0x1.33caea0f24cp-2
  cSinpiDb4_X_2: Tb64u64 = (u:QWord($BFE8A1E8A3E82C00)); // -0x1.8a1e8a3e82cp-1
  cSinpiDb4_X_3: Tb64u64 = (u:QWord($BFDBDD02D1AD6000)); // -0x1.bdd02d1ad60p-2
  cSinpiDb4_R_0: Tb64u64 = (u:$3FBDB8A79A80C3A0);        // 0x1.db8a79a80c3a0p-4
  cSinpiDb4_R_1: Tb64u64 = (u:$3FC5146C0BC45BCE);        // 0x1.5146c0bc45bcep-3
  cSinpiDb4_R_2: Tb64u64 = (u:$3FE5D6561936B699);        // 0x1.5d6561936b699p-1
  cSinpiDb4_R_3: Tb64u64 = (u:$3FEF72C906962631);        // 0x1.f72c906962631p-1
  cSinpiDb4_D_0: Tb64u64 = (u:$3910000000000000);        //  0x1p-110
  cSinpiDb4_D_1: Tb64u64 = (u:$3920000000000000);        //  0x1p-109
  cSinpiDb4_D_2: Tb64u64 = (u:QWord($BC80000000000000)); // -0x1p-55
  cSinpiDb4_D_3: Tb64u64 = (u:$3C80000000000000);        //  0x1p-55

// ---------------------------------------------------------------------------
// SinpiAsZero — small-|x| (|x| <= 0x1.2p-5) accurate refinement.
// Port of as_sinpi_zero in sinpi.c:86.
// ---------------------------------------------------------------------------
function SinpiAsZero(x: Double): Double;
var
  x2, dx2, x3, dx3: Double;
  fh, fl: Double;
  chp, clp, thp, tlp: Double;
  pi0, pi1, pi2: Double;
  y0, y0l, y1, y2: Double;
  b: Tb64u64;
  t, w: Tb64u64;
  sA, dA, sB, dB: Double;
begin
  x2  := x * x;
  dx2 := pcr_fma(x, x, -x2);
  x3  := x2 * x;
  dx3 := pcr_fma(x2, x, -x3) + dx2 * x;

  // fl = x2*(cl[0] + x2*(cl[1] + x2*cl[2]))
  fl := x2 * (cSinpiZeroCl_0.f + x2 * (cSinpiZeroCl_1.f + x2 * cSinpiZeroCl_2.f));

  // polydd(x2, dx2, 4, ch, &fl): unrolled 4-entry.
  // Pillar A (Phase 6.3): for i := 2 downto 0 unrolled.
  chp := cSinpiZeroCh_3_0.f + fl;
  clp := ((cSinpiZeroCh_3_0.f - chp) + fl) + cSinpiZeroCh_3_1.f;
  // i = 2
  chp := pcr_muldd(x2, dx2, chp, clp, clp);
  thp := chp + cSinpiZeroCh_2_0.f;
  tlp := (cSinpiZeroCh_2_0.f - thp) + chp;
  chp := thp;
  clp := clp + tlp + cSinpiZeroCh_2_1.f;
  // i = 1
  chp := pcr_muldd(x2, dx2, chp, clp, clp);
  thp := chp + cSinpiZeroCh_1_0.f;
  tlp := (cSinpiZeroCh_1_0.f - thp) + chp;
  chp := thp;
  clp := clp + tlp + cSinpiZeroCh_1_1.f;
  // i = 0
  chp := pcr_muldd(x2, dx2, chp, clp, clp);
  thp := chp + cSinpiZeroCh_0_0.f;
  tlp := (cSinpiZeroCh_0_0.f - thp) + chp;
  chp := thp;
  clp := clp + tlp + cSinpiZeroCh_0_1.f;
  fh := chp; fl := clp;

  // fh = muldd_acc(fh, fl, x3, dx3, &fl)
  fh := pcr_muldd(fh, fl, x3, dx3, fl);

  // y0 = pi0 * x; extract "exactly representable" part via magic rounding
  pi0 := cSinpiZeroPi0.f;
  pi1 := cSinpiZeroPi1.f;
  pi2 := cSinpiZeroPi2.f;

  y0 := pi0 * x;
  b.f := y0;
  b.u := b.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12);
  b.u := b.u + (UInt64(85) shl 51);
  y0 := (y0 + b.f) - b.f;
  y0l := pcr_fma(pi0, x, -y0);

  y1 := pi1 * x;
  y2 := pcr_fma(pi1, x, -y1) + pi2 * x;

  // y1 = adddd(y1, y2, y0l, 0, &y2) — inlined (no pcr_adddd helper)
  sA := y1 + y0l; dA := sA - y1;
  y2 := ((y0l - dA) + (y1 + (dA - sA))) + y2;
  y1 := sA;
  // y1 = adddd(y1, y2, fh, fl, &y2) — inlined
  sB := y1 + fh; dB := sB - y1;
  y2 := ((fh - dB) + (y1 + (dB - sB))) + (y2 + fl);
  y1 := sB;

  // y0 = fasttwosum(y0, y1, &y1); y1 = fasttwosum(y1, y2, &y2)
  pcr_fasttwosum(y0, y1, y0, y1);
  pcr_fasttwosum(y1, y2, y1, y2);

  t.f := y1;
  if (t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) = 0 then begin
    w.f := y2;
    if ((w.u xor t.u) shr 63) <> 0 then
      t.u := t.u - 1
    else
      t.u := t.u + 1;
    y1 := t.f;
  end;
  Result := y0 + y1;
end;

// ---------------------------------------------------------------------------
// SinpiRefine — general accurate path (iq, z).
// Math identical to CospiSinpiRefine in pascoremath64.pas, but uses the
// 4-entry sinpi-specific exception database.
// ---------------------------------------------------------------------------
function SinpiRefine(iq: Int32; z: Double): Double;
var
  x, x2, dx2: Double;
  sll, slh, cll, clh: Double;
  chp, clp, thp, tlp: Double;
  sbh, sbl, cbh, cbl: Double;
  csh, csl, sch, scl: Double;
  tsh, tsl, tsl2: Double;
  iqm: Int32;
  sgn, dbx: Double;
  i: Int32;
  t_u: Tb64u64;
begin
  x   := z * cSinpiRefScale.f;  // z * 2^-63
  x2  := x * x;
  dx2 := pcr_fma(x, x, -x2);

  sll := cSinpiRefSllK.f * x2;

  // polydd(x2, dx2, 3, sh, &sll) seeded
  chp := cSinpiRefSh_2_0.f + sll;
  clp := ((cSinpiRefSh_2_0.f - chp) + sll) + cSinpiRefSh_2_1.f;
  chp := pcr_muldd(x2, dx2, chp, clp, clp);
  thp := chp + cSinpiRefSh_1_0.f;
  tlp := (cSinpiRefSh_1_0.f - thp) + chp;
  chp := thp;
  clp := clp + tlp + cSinpiRefSh_1_1.f;
  chp := pcr_muldd(x2, dx2, chp, clp, clp);
  thp := chp + cSinpiRefSh_0_0.f;
  tlp := (cSinpiRefSh_0_0.f - thp) + chp;
  chp := thp;
  clp := clp + tlp + cSinpiRefSh_0_1.f;
  slh := chp; sll := clp;

  slh := pcr_mulddd_pd(slh, sll, x * cSinpiRefXlow.f, sll);

  cll := x2 * (cSinpiRefCll0.f + cSinpiRefCll1.f * x2);

  chp := cSinpiRefCh_1_0.f + cll;
  clp := ((cSinpiRefCh_1_0.f - chp) + cll) + cSinpiRefCh_1_1.f;
  chp := pcr_muldd(x2, dx2, chp, clp, clp);
  thp := chp + cSinpiRefCh_0_0.f;
  tlp := (cSinpiRefCh_0_0.f - thp) + chp;
  chp := thp;
  clp := clp + tlp + cSinpiRefCh_0_1.f;
  clh := chp; cll := clp;

  clh := pcr_muldd(clh, cll, x2, dx2, cll);

  CospiSincosN2(iq, sbh, sbl, cbh, cbl);

  csh := pcr_muldd(clh, cll, sbh, sbl, csl);
  sch := pcr_muldd(slh, sll, cbh, cbl, scl);

  pcr_fasttwosum(tsh, tsl, sch, csh);
  tsl := tsl + csl + scl;
  pcr_fasttwosum(tsh, tsl2, sbh, tsh);
  tsl := sbl + tsl + tsl2;

  // Exception database probe (sinpi-specific 4-entry db)
  t_u.f := tsl;
  if ((t_u.u or (UInt64($FFF) shl 52)) = UInt64($FFFFFFFFFFFFFFFF))
     or ((t_u.u shl 12) = 0) then
  begin
    if iq > 2048 then sgn := -Double(1.0) else sgn := Double(1.0);
    iqm := iq and $7FF;
    // Pillar A (Phase 6.3): for i := 0 to 3 unrolled.
    dbx := cSinpiDb4_X_0.f;
    if ((x = dbx) and (iqm = cSinpiDb4_Iq_0))
       or ((x = -dbx) and (iqm = 2048 - cSinpiDb4_Iq_0)) then
    begin
      Result := sgn * cSinpiDb4_R_0.f + sgn * cSinpiDb4_D_0.f;
      Exit;
    end;
    dbx := cSinpiDb4_X_1.f;
    if ((x = dbx) and (iqm = cSinpiDb4_Iq_1))
       or ((x = -dbx) and (iqm = 2048 - cSinpiDb4_Iq_1)) then
    begin
      Result := sgn * cSinpiDb4_R_1.f + sgn * cSinpiDb4_D_1.f;
      Exit;
    end;
    dbx := cSinpiDb4_X_2.f;
    if ((x = dbx) and (iqm = cSinpiDb4_Iq_2))
       or ((x = -dbx) and (iqm = 2048 - cSinpiDb4_Iq_2)) then
    begin
      Result := sgn * cSinpiDb4_R_2.f + sgn * cSinpiDb4_D_2.f;
      Exit;
    end;
    dbx := cSinpiDb4_X_3.f;
    if ((x = dbx) and (iqm = cSinpiDb4_Iq_3))
       or ((x = -dbx) and (iqm = 2048 - cSinpiDb4_Iq_3)) then
    begin
      Result := sgn * cSinpiDb4_R_3.f + sgn * cSinpiDb4_D_3.f;
      Exit;
    end;
  end;
  Result := tsh + tsl;
end;

// Arithmetic shift right for Int64 — FPC's `shr` is logical, which corrupts
// sign-extension. Sinpi's algorithm signs-extends the mantissa for negative x
// then shifts (cospi doesn't, because it keeps m positive; cos is even).
function SinpiSar(m: Int64; n: Int32): Int64; inline;
begin
  if m < 0 then
    Result := Int64(not ((not UInt64(m)) shr n))
  else
    Result := Int64(UInt64(m) shr n);
end;

// ---------------------------------------------------------------------------
// pcr_sinpi — entry point.
// Port of cr_sinpi in sinpi.c:154.
// ---------------------------------------------------------------------------
function pcr_sinpi(x: Double): Double;
var
  ix: Tb64u64;
  ax: UInt64;
  e, s_shift, si: Int32;
  m0: UInt64;
  sgn_mask: Int64;
  m: Int64;
  iq_u: UInt64;
  iq: Int32;
  ph, pl, zh, zl, t_scale, r, rs, rt: Double;
  x2, x3, x4, eps, lb, ub: Double;
  z, z2, fs, fc, er, rr: Double;
  sh, sl, ch_, cl_: Double;
  k: Int64;
  tg: Int32;
begin
  ix.f := x;
  ax   := ix.u and (UInt64($FFFFFFFFFFFFFFFF) shr 1);

  if ax = 0 then begin Result := x; Exit; end;  // preserves sign of zero

  e  := Int32(ax shr 52);
  m0 := (ix.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) or (UInt64(1) shl 52);
  if (ix.u shr 63) <> 0 then begin
    sgn_mask := -1;
    m := -Int64(m0);
  end else begin
    sgn_mask := 0;
    m := Int64(m0);
  end;
  s_shift := 1063 - e;

  if s_shift < 0 then
  begin
    // |x| >= 2^41
    if e = $7FF then
    begin
      if (ix.u shl 12) = 0 then
      begin
        Result := cNaNDoublePos.f;  // Inf -> qNaN
        Exit;
      end;
      Result := x + x; // NaN
      Exit;
    end;
    s_shift := -s_shift - 1;
    if s_shift > 10 then
    begin
      // Result is +/-0 with sign of x
      if sgn_mask <> 0 then Result := -Double(0.0) else Result := Double(0.0);
      Exit;
    end;
    iq_u := UInt64(m) shl s_shift;
    if (iq_u and 2047) = 0 then
    begin
      if sgn_mask <> 0 then Result := -Double(0.0) else Result := Double(0.0);
      Exit;
    end;
    CospiSincosN(Int32(iq_u), sh, sl, ch_, cl_);
    Result := sh + sl;
    Exit;
  end;

  if ax <= cSinpiAxT5 then
  begin
    // |x| <= 0x1.2p-5
    ph := cSinpiPh.f;
    pl := cSinpiPl.f;
    if ax < cSinpiAxPm54 then
    begin
      // |x| < 2^-54
      if ax < $0350000000000000 then  // |x| < 2^-970
      begin
        // Scale-up path to preserve precision for subnormal results.
        // Uses hardware FMA (vfmadd213sd) because pcr_fma_pascal's Dekker
        // emulation is not correctly rounded for subnormal outputs.
        t_scale := x * cSinpi2p106.f;
        zh := ph * t_scale;
        zl := pcr_fma(ph, t_scale, -zh) + pl * t_scale;
        r  := zh + zl;
        rs := r * cSinpi2pm106.f;
        rt := rs * cSinpi2p106.f;
        Result := pcr_fma((zh - rt) + zl, cSinpi2pm106.f, rs);
        Exit;
      end;
      zh := ph * x;
      zl := pcr_fma(ph, x, -zh) + pl * x;
      Result := zh + zl;
      Exit;
    end;
    // Normal-|x| fast path
    zh := ph * x;
    zl := pcr_fma(ph, x, -zh) + pl * x;
    x2 := x * x;
    x3 := x2 * x;
    x4 := x2 * x2;
    eps := x * (x2 * cSinpiFastEps1.f + cSinpiFastEps2.f);
    zl := zl + x3 * ((cSinpiFastC_0.f + x2 * cSinpiFastC_1.f)
                     + x4 * (cSinpiFastC_2.f + x2 * cSinpiFastC_3.f));
    lb := (zl - eps) + zh;
    ub := (zl + eps) + zh;
    if lb = ub then begin Result := lb; Exit; end;
    Result := SinpiAsZero(x);
    Exit;
  end;

  // Integer / half-integer detection: si = e - 1011, check bits below s1
  si := e - 1011;
  if (si >= 0) and (si < 64) and ((UInt64(m0) shl (si + 1)) = 0) then
  begin
    // x is integer or half-integer
    if (UInt64(m0) shl si) = 0 then
    begin
      // x is an integer -> sin(pi*integer) = +/-0 with sign of x
      if sgn_mask <> 0 then Result := -Double(0.0) else Result := Double(0.0);
      Exit;
    end;
    // half-integer: sign alternates with the integer part parity
    if si = 0 then tg := 0
    else tg := Int32((UInt64(m0) shl (si - 1)) shr 63);
    if tg = 0 then
    begin
      // |x| = 1/2 mod 2: sin(pi*x) = +/-1 with sign of x
      if sgn_mask <> 0 then Result := -Double(1.0) else Result := Double(1.0);
    end
    else
    begin
      // |x| = 3/2 mod 2: sin(pi*x) = -/+1 (flip sign of x)
      if sgn_mask <> 0 then Result := Double(1.0) else Result := -Double(1.0);
    end;
    Exit;
  end;

  // Main fast path
  iq_u := UInt64(SinpiSar(m, s_shift)) and 8191;
  iq_u := (iq_u + 1) shr 1;
  iq := Int32(iq_u);

  k := Int64(UInt64(m) shl (e - 1000));
  z := k;
  z2 := z * z;

  // sn[] and cn[] coefficients (same as cospi's cCospiSn / cCospiCn)
  fs := cCospiSn_0.f + z2 * (cCospiSn_1.f + z2 * cCospiSn_2.f);
  fc := cCospiCn_0.f + z2 * cCospiCn_1.f;

  CospiSincosN(iq, sh, sl, ch_, cl_);

  er := cSinpiFastEr.f; // 5.5e-19
  rr := sl + sh * (z2 * fc) + ch_ * (z * fs);
  lb := (rr - er) + sh;
  ub := (rr + er) + sh;
  if lb = ub then begin Result := lb; Exit; end;

  Result := SinpiRefine(iq, z);
end;
