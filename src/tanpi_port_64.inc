// ===========================================================================
// pcr_tanpi — correctly-rounded binary64 tan(pi*x).
// Ported from core-math/src/binary64/tanpi/tanpi.c by Alexei Sibidanov.
// Tanpi is standalone: its T[][2]/ch/cl tables are distinct from
// cospi/sinpi/tan. No sharing with previously-ported pi-trig functions.
// ===========================================================================

const
{$I tanpi_const_64.inc}

// ---------------------------------------------------------------------------
// TanpiDatabase — port of as_tanpi_database (tanpi.c:88).
// ---------------------------------------------------------------------------
function TanpiDatabase(x, f: Double): Double;
var
  ix, t: Tb64u64;
  aix, p, mask: UInt64;
  ax: Double;
  e: Int32;
  sgn: Double;
  a, b, m: Int32;
  dbx: Double;
begin
  ix.f := x;
  aix := ix.u and (UInt64($FFFFFFFFFFFFFFFF) shr 1);
  p := 0;
  ax := Abs(x);
  e := Int32(aix shr 52) - 1022;
  if e >= 0 then
  begin
    mask := not (UInt64($FFFFFFFFFFFFFFFF) shr (12 + e));
    t.u := aix and mask;
    ax := ax - t.f;
    p := ((aix or (UInt64(1) shl 52)) shr (52 - e)) and 1;
    if p <> 0 then ax := Double(0.5) - ax;
  end;
  if (p xor (ix.u shr 63)) <> 0 then sgn := -Double(1.0) else sgn := Double(1.0);
  a := 0; b := 48; m := 24;
  while a <= b do
  begin
    dbx := cTanpiDbX[m].f;
    if dbx < ax then a := m + 1
    else if dbx = ax then
    begin
      f := sgn * cTanpiDbR[m].f + sgn * cTanpiDbD[m].f;
      Break;
    end
    else b := m - 1;
    m := (a + b) div 2;
  end;
  Result := f;
end;

// ---------------------------------------------------------------------------
// TanpiApply — shared "iq==32 or T[iq] addition formula" post-processing.
// Input: (th, tl) = pi*z + f (small-angle correction), (iq, ms) describe the
// table offset and sign of the 1/2-bit.  Output: (th, tl) = tan(pi*x).
// ---------------------------------------------------------------------------
procedure TanpiApply(iq: Int32; ms: Int64; var th, tl: Double);
var
  nl_in, nh_in: Double;
  nl, nh: Double;
  s2sel: Double;
  ml, mh, mh_old, dm, dn, imh: Double;
  ith: Double;
  th_in, tl_in: Double;
  s_tmp: Double;
begin
  if iq = 32 then
  begin
    ith := -Double(1.0) / th;
    tl := (pcr_fma(ith, th, Double(1.0)) + tl * ith) * ith;
    th := ith;
    Exit;
  end;
  nl_in := cTanpiT[iq, 0].f;
  nh_in := cTanpiT[iq, 1].f;
  if (ms + 1) = 0 then s2sel := -Double(1.0) else s2sel := Double(1.0);
  nh := nh_in * s2sel;
  nl := nl_in * s2sel;
  th_in := th;
  tl_in := tl;
  // (mh, ml) = muldd_acc(th, tl, nh, nl)
  mh := pcr_muldd(th_in, tl_in, nh, nl, ml);
  // fasttwosub(1, mh, &dm): s = 1 - mh_old; dm = (1 - s) - mh_old
  mh_old := mh;
  s_tmp := Double(1.0) - mh_old;
  dm := (Double(1.0) - s_tmp) - mh_old;
  mh := s_tmp;
  ml := dm - ml;
  // fasttwosum(nh, th, &dn): s = nh + th; dn = th - (s - nh)
  s_tmp := nh + th_in;
  dn := th_in - (s_tmp - nh);
  nh := s_tmp;
  nl := nl + dn + tl_in;
  imh := Double(1.0) / mh;
  th := nh * imh;
  tl := pcr_fma(nh, imh, -th) + (nl + nh * (pcr_fma(-mh, imh, Double(1.0)) - ml * imh)) * imh;
end;

// ---------------------------------------------------------------------------
// TanpiHwFma — hardware FMA for the subnormal scale-up branch.
// pcr_fma_pascal's Dekker emulation is not correctly rounded for subnormal
// outputs; use hw FMA directly (same trick as sinpi/atanpi).
// ---------------------------------------------------------------------------
{$IFDEF AVX2}
function TanpiHwFma(a, b, c: Double): Double; assembler; nostackframe;
asm
  vfmadd213sd  xmm0, xmm1, xmm2
end;
{$ELSE}
function TanpiHwFma(a, b, c: Double): Double; inline;
begin
  Result := pcr_fma(a, b, c);
end;
{$ENDIF}

// ---------------------------------------------------------------------------
// pcr_tanpi — entry point. Port of cr_tanpi in tanpi.c:167.
// ---------------------------------------------------------------------------
function pcr_tanpi(x: Double): Double;
var
  ix, ul, uh: Tb64u64;
  ax: UInt64;
  res, th, tl: Double;
  e: Int32;
  s, s1: Int32;
  m: Int64;
  sgn: Int64;
  iq_u: UInt64;
  jq: Int64;
  nl, nh: Double;
  ms, sm: Int64;
  k: Int64;
  z, z2, z3, z4, f, eps: Double;
  lb, ub: Double;
  iq: Int32;
  dx2, dx3: Double;
  chp, clp, thp, tlp: Double;
  i: Int32;
  px0, px1, dv: Double;
  x2: Double;
  er, de: UInt64;
  scu, iscu: Tb64u64;
  sc, isc: Double;
  o, v0b, v0h: Double;
  iq_mask: Int32;
begin
  ix.f := x;
  ax := ix.u and (UInt64($FFFFFFFFFFFFFFFF) shr 1);

  if ax >= (UInt64($3F3) shl 52) then
  begin
    // ---------- |x| >= 2^-12 ----------
    if ax >= (UInt64($42D) shl 52) then
    begin
      // |x| >= 2^46: integer / half-integer structure of tan(pi*x)
      if ax >= (UInt64($7FF) shl 52) then
      begin
        if ax > (UInt64($7FF) shl 52) then begin Result := x + x; Exit; end;
        Result := cNaNDoublePos.f; // +-Inf -> qNaN
        Exit;
      end;
      e := Int32(ax shr 52);
      s := e - 1069;
      if s > 6 then
      begin
        if (ix.u shr 63) <> 0 then Result := -Double(0.0) else Result := Double(0.0);
        Exit;
      end;
      m := Int64(ax);
      if (ix.u shr 63) <> 0 then sgn := -1 else sgn := 0;
      iq_u := (UInt64((m xor sgn) - sgn) shl s) and 127;
      if (iq_u and 31) = 0 then
      begin
        jq := Int64(iq_u shr 5);
        if (jq and 1) <> 0 then
        begin
          if (jq and 2) <> 0 then Result := -Double(1.0) / Double(0.0)
          else                    Result :=  Double(1.0) / Double(0.0);
        end
        else
        begin
          if ((jq xor sgn) and 2) <> 0 then Result := -Double(0.0)
          else                              Result :=  Double(0.0);
        end;
        Exit;
      end;
      if (iq_u and 32) <> 0 then
      begin
        nl := -cTanpiT[32 - (iq_u and 31), 0].f;
        nh := -cTanpiT[32 - (iq_u and 31), 1].f;
      end
      else
      begin
        nl := cTanpiT[iq_u and 31, 0].f;
        nh := cTanpiT[iq_u and 31, 1].f;
      end;
      Result := nh + nl;
      Exit;
    end;

    // ---------- 2^-12 <= |x| < 2^46 ----------
    e  := Int32(ax shr 52);
    s  := 1068 - e;
    s1 := e - 1011;
    m  := Int64((ax and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) or (UInt64(1) shl 52));
    ms := SarInt64(Int64(UInt64(m) shl s1), 63);
    if (ix.u shr 63) <> 0 then sgn := -1 else sgn := 0;
    iq_u := (UInt64(m xor ms) shr s) and 63;
    iq_u := (iq_u + 1) shr 1;
    ms := ms xor sgn;
    sm := (m xor sgn) - sgn;
    k := Int64(UInt64(sm) shl (e - 1005));
    z := k;
    if (UInt64(k) shl 1) = 0 then
    begin
      if k = 0 then
      begin
        if (iq_u and 31) = 0 then
        begin
          jq := SarInt64(sm, s + 6);
          if (jq and 1) <> 0 then
          begin
            if (jq and 2) <> 0 then Result := -Double(1.0) / Double(0.0)
            else                    Result :=  Double(1.0) / Double(0.0);
          end
          else
          begin
            if ((jq xor sgn) and 2) <> 0 then Result := -Double(0.0)
            else                              Result :=  Double(0.0);
          end;
          Exit;
        end;
        iq_mask := Int32((UInt64(m) shl s1) shr 58);
        if iq_mask = $10 then
        begin
          if (ix.u shr 63) <> 0 then Result := -Double(1.0) else Result := Double(1.0);
          Exit;
        end;
        if iq_mask = $30 then
        begin
          if (ix.u shr 63) <> 0 then Result :=  Double(1.0) else Result := -Double(1.0);
          Exit;
        end;
      end;
      if (ix.u shr 63) <> 0 then z := -z;
    end;
    iq := Int32(iq_u);

    // Fast path polynomial f = z^3 * poly(z^2)
    z2 := z * z;
    z4 := z2 * z2;
    z3 := z * z2;
    f := z3 * ((cTanpiC[0].f + z2 * cTanpiC[1].f)
             + z4 * (cTanpiC[2].f + z2 * cTanpiC[3].f));
    if z >= Double(0.0) then
      eps := z3 * cTanpi2pm256.f + cTanpi2pm104.f
    else
      eps := z3 * cTanpi2pm256.f - cTanpi2pm104.f;
    th := pcr_mulddd_pd(cTanpiPh.f, cTanpiPl.f, z, tl);
    // fastsumddd: fasttwosum(th, f) -> (th', tmp); tl += tmp
    pcr_fasttwosum(th, tlp, th, f);
    tl := tl + tlp;

    TanpiApply(iq, ms, th, tl);

    eps := eps + eps * (th * th);
    lb := th + (tl - eps);
    ub := th + (tl + eps);
    if lb = ub then begin Result := lb; Exit; end;

    // --- Refinement path ---
    z := z * cTanpi2pm63.f;
    z2 := z * z;
    dx2 := pcr_fma(z, z, -z2);
    tl := z2 * (cTanpiRefCl[0].f + z2 * (cTanpiRefCl[1].f + z2 * cTanpiRefCl[2].f));
    // polydd(z2, dx2, 5, ch, &tl) — seeded
    chp := cTanpiRefCh[4, 0].f + tl;
    clp := ((cTanpiRefCh[4, 0].f - chp) + tl) + cTanpiRefCh[4, 1].f;
    for i := 3 downto 0 do
    begin
      chp := pcr_muldd(z2, dx2, chp, clp, clp);
      thp := chp + cTanpiRefCh[i, 0].f;
      tlp := (cTanpiRefCh[i, 0].f - thp) + chp;
      chp := thp;
      clp := clp + tlp + cTanpiRefCh[i, 1].f;
    end;
    th := chp; tl := clp;
    th := pcr_mulddd_pd(th, tl, z, tl);
    TanpiApply(iq, ms, th, tl);
    pcr_fasttwosum(th, tl, th, tl);
    res := th;
  end
  else
  begin
    // ---------- |x| < 2^-12 ----------
    if ax = 0 then begin Result := x; Exit; end;
    if ax < (UInt64($3CA) shl 52) then
    begin
      // |x| < 2^-53
      if ax < (UInt64($6B) shl 52) then
      begin
        // |x| < 2^-916: scale-up path
        e := Int32(ax shr 52);
        scu.u  := UInt64(Int64(2045) - e) shl 52;
        iscu.u := UInt64(Int64(1)    + e) shl 52;
        sc  := scu.f;
        isc := iscu.f;
        z := x * sc;
        th := pcr_mulddd_pd(cTanpiPi0.f, cTanpiPi1.f, z, tl);
        res := th * isc;
        if Abs(res) < Double(2.2250738585072014e-308) then
        begin
          // Forced underflow
          if (ix.u shr 63) <> 0 then o := -Double(2.2250738585072014e-308)
          else                       o :=  Double(2.2250738585072014e-308);
          v0b := (o + res) * sc;
          v0h := res * sc;
          tl := tl + (th - v0h);
          v0b := v0b + tl;
          Result := v0b * isc - o;
          Exit;
        end;
      end
      else
      begin
        // 2^-916 <= |x| < 2^-53
        th := pcr_mulddd_pd(cTanpiPi0.f, cTanpiPi1.f, x, tl);
        res := th;
      end;
    end
    else
    begin
      // 2^-53 <= |x| < 2^-12: degree-3 fast path then refine
      x2 := x * x;
      z3 := x * x2;
      f := z3 * (cTanpiSmallC[0].f + x2 * (cTanpiSmallC[1].f + x2 * cTanpiSmallC[2].f));
      px0 := pcr_mulddd_pd(cTanpiPi0.f, cTanpiPi1.f, x, px1);
      pcr_fasttwosum(th, tlp, px0, f);
      tl := px1 + tlp;
      eps := x * (x2 * cTanpiSmEpsA.f + cTanpiSmEpsB.f);
      lb := th + (tl - eps);
      ub := th + (tl + eps);
      if lb = ub then begin Result := lb; Exit; end;

      // Refinement
      dx2 := pcr_fma(x, x, -x2);
      dx3 := pcr_fma(x2, x, -z3) + dx2 * x;
      tl := x2 * (cTanpiSmRefCl[0].f + x2 * cTanpiSmRefCl[1].f);
      // polydd seeded: ch = c[n-1][0] + *l; cl = ((c[n-1][0] - ch) + *l) + c[n-1][1]
      chp := cTanpiSmRefCh[2, 0].f + tl;
      clp := ((cTanpiSmRefCh[2, 0].f - chp) + tl) + cTanpiSmRefCh[2, 1].f;
      for i := 1 downto 0 do
      begin
        chp := pcr_muldd(x2, dx2, chp, clp, clp);
        thp := chp + cTanpiSmRefCh[i, 0].f;
        tlp := (cTanpiSmRefCh[i, 0].f - thp) + chp;
        chp := thp;
        clp := clp + tlp + cTanpiSmRefCh[i, 1].f;
      end;
      th := chp; tl := clp;
      th := pcr_muldd(z3, dx3, th, tl, tl);
      // fasttwosum(px0, th, &dv)
      pcr_fasttwosum(th, dv, px0, th);
      tl := tl + px1 + dv;
      pcr_fasttwosum(th, tl, th, tl);
      res := th;
    end;
  end;

  ul.f := tl;
  uh.f := th;
  er := (ul.u + 6) and (UInt64($FFFFFFFFFFFFFFFF) shr 12);
  de := ((uh.u - ul.u) shr 52) and UInt64($7FF);
  if (er <= 12) or (de > 102) then
  begin
    Result := TanpiDatabase(x, res);
    Exit;
  end;
  Result := res;
end;
