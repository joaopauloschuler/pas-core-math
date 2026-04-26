// ---------------------------------------------------------------------------
// pcr_atan2pi — correctly-rounded binary64 atan2pi(y, x) = atan2(y, x) / pi.
// Ported from core-math/src/binary64/atan2pi/atan2pi.c
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
//   (Paul Zimmermann, 2023-2025), itself derived from atan2.c.
//
// Reuses cAtan2P / cAtan2Q (TInt64 polynomial coefficients) from the atan2
// port — the [29,29] rational approximation of atan(z) is shared. Adds its
// own dedicated 65×10 Pfast / 65-entry Xfast / 65-entry err_fast tables for
// the fast path (per-segment polynomials, distinct from the Sibidanov
// reduction used by atan2's fast path).
//
// fenv / errno / FE flag tracking is intentionally omitted (out of scope per
// Phase 0 design) — the returned Double matches the C reference; only the
// FE side effects diverge.
// ---------------------------------------------------------------------------

{$I inc_64/atan2pi_const_64.inc}

const
  cA2pi_MASK: UInt64 = UInt64($7FFFFFFFFFFFFFFF);
  cA2pi_2pm54: Tb64u64 = (u: $3C90000000000000); // 0x1p-54
  cA2pi_2pm55: Tb64u64 = (u: $3C80000000000000); // 0x1p-55

// ---------------------------------------------------------------------------
// Local helpers — match atan2pi.c's a_mul / s_mul / d_mul / fast_two_sum
// / dd_sum1 byte-for-byte. (pcr_d_mul renormalizes; we cannot reuse it here.)
// ---------------------------------------------------------------------------

procedure A2piFastTwoSum(out hi, lo: Double; a, b: Double); inline;
var e: Double;
begin
  hi := a + b;
  e  := hi - a;
  lo := b - e;
end;

// dd_sum1: assumes |ah| >= |bh|; al implicitly 0.
procedure A2piDDSum1(var h, l: Double; ah, bh, bl: Double); inline;
var hh, ll: Double;
begin
  A2piFastTwoSum(hh, ll, ah, bh);
  h := hh;
  l := ll + bl;
end;

procedure A2piAMul(out hi, lo: Double; a, b: Double); inline;
begin
  hi := a * b;
  lo := pcr_fma(a, b, -hi);
end;

procedure A2piSMul(out hi, lo: Double; a, bh, bl: Double); inline;
begin
  A2piAMul(hi, lo, a, bh);
  lo := pcr_fma(a, bl, lo);
end;

// d_mul: (ah+al)*(bh+bl) -- (al*bl) ignored. Order: a_mul, fma(al,bh,...),
// fma(ah,bl,...). Matches atan2pi.c exactly.
procedure A2piDMul(out hi, lo: Double; ah, al, bh, bl: Double); inline;
var s, t: Double;
begin
  A2piAMul(hi, s, ah, bh);
  t  := pcr_fma(al, bh, s);
  lo := pcr_fma(ah, bl, t);
end;

// ---------------------------------------------------------------------------
// fast_div: h+l ~ bh/ah via Karp-Markstein. Matches atan2pi.c::fast_div.
// ---------------------------------------------------------------------------
procedure A2piFastDiv(out h, l: Double; bh, ah: Double); inline;
var y, eh: Double;
begin
  y  := Double(1.0) / ah;
  h  := bh * y;
  eh := pcr_fma(ah, -h, bh);
  l  := y * eh;
end;

// ---------------------------------------------------------------------------
// Special-path: NaN / Inf / zero handling. Mirrors cr_atan2pi's special
// branches (atan2pi.c lines 796-844).
// ---------------------------------------------------------------------------

function Atan2piSpecial(y, x: Double): Double;
var
  uy, ux: Tb64u64;
  ay, ax: UInt64;
  ey, ex: UInt64;
  isNaNy, isNaNx, isInfy, isInfx: Boolean;
begin
  uy.f := y;  ux.f := x;
  ay := uy.u and cA2pi_MASK;
  ax := ux.u and cA2pi_MASK;
  ey := ay shr 52;
  ex := ax shr 52;

  // NaN/Inf branch
  if (ey = $7FF) or (ex = $7FF) then begin
    isNaNy := (ey = $7FF) and (ay <> (ey shl 52));
    isNaNx := (ex = $7FF) and (ax <> (ex shl 52));
    if isNaNy or isNaNx then begin
      Result := y + x;
      Exit;
    end;
    isInfy := ay = (UInt64($7FF) shl 52);
    isInfx := ax = (UInt64($7FF) shl 52);
    if isInfy and isInfx then begin
      if x < 0 then
        if y > 0 then Result := Double(0.75) else Result := Double(-0.75)
      else
        if y > 0 then Result := Double(0.25) else Result := Double(-0.25);
      Exit;
    end;
    if isInfx then begin
      if x < 0 then begin
        if (uy.u shr 63) <> 0 then Result := Double(-1.0) else Result := Double(1.0);
        Exit;
      end;
      Result := pcr_copysign(Double(0.0), y);
      Exit;
    end;
    // y = +/-Inf, x finite
    if y > 0 then Result := Double(0.5) else Result := Double(-0.5);
    Exit;
  end;

  // Zero branch
  if (ay = 0) and (ax = 0) then begin
    if ux.u = 0 then begin
      Result := y; // atan2pi(+/-0, +0) = +/-0
      Exit;
    end;
    // atan2pi(+/-0, -0) = +/-1
    if uy.u = 0 then Result := Double(1.0) else Result := Double(-1.0);
    Exit;
  end;
  if ay = 0 then begin
    if x > 0 then begin
      Result := y;
      Exit;
    end;
    if uy.u = 0 then Result := Double(1.0) else Result := Double(-1.0);
    Exit;
  end;
  // ax = 0 (only x is zero)
  if y > 0 then Result := Double(0.5) else Result := Double(-0.5);
end;

// ---------------------------------------------------------------------------
// Atan2piAccurate — TInt64-based [29,29] rational approximation of atan,
// then divide by pi. Mirrors atan2pi_accurate (atan2pi.c lines 211-372).
// ---------------------------------------------------------------------------

function Atan2piAccurate(y, x: Double): Double;
var
  absy, absx: Double;
  inv: Boolean;
  z, p, q: TInt64;
  i: Integer;
  err: UInt64;
  sz: UInt64;
  uy, ux: Tb64u64;
  ey, ex: Int32;
begin
  absy := Abs(y);
  absx := Abs(x);
  inv  := absy > absx;

  if absy = absx then begin
    if y > 0 then begin
      if x > 0 then Result := Double(0.25) else Result := Double(0.75);
    end else begin
      if x < 0 then Result := Double(-0.75) else Result := Double(-0.25);
    end;
    Exit;
  end;

  if inv then begin
    uy.f := y;  ux.f := x;
    ex := Int32((ux.u shr 52) and $7FF);
    ey := Int32((uy.u shr 52) and $7FF);
    if ey - ex > 54 then begin
      if x > 0 then begin
        if y > 0 then Result :=  Double(0.5) - cA2pi_2pm55.f
                 else Result := Double(-0.5) + cA2pi_2pm55.f;
      end else begin
        if y > 0 then Result :=  Double(0.5) + cA2pi_2pm54.f
                 else Result := Double(-0.5) - cA2pi_2pm54.f;
      end;
      Exit;
    end;
    DivTIntD(z, x, y);
  end else begin
    if x < 0 then begin
      uy.f := y;  ux.f := x;
      ex := Int32((ux.u shr 52) and $7FF);
      ey := Int32((uy.u shr 52) and $7FF);
      if ey - ex < -54 then begin
        if y > 0 then Result :=  Double(1.0) - cA2pi_2pm54.f
                 else Result := Double(-1.0) + cA2pi_2pm54.f;
        Exit;
      end;
    end;
    DivTIntD(z, y, x);
  end;

  sz := z.sgn;
  z.sgn := 0;

  CpTInt(p, cAtan2P[29]);
  CpTInt(q, cAtan2Q[29]);
  for i := 28 downto 0 do begin
    MulTInt(p, p, z);
    MulTInt(q, q, z);
    AddTInt(p, p, cAtan2P[i]);
    AddTInt(q, q, cAtan2Q[i]);
  end;
  MulTInt(p, p, z);
  DivTInt(z, p, q);

  // multiply by 1/pi
  MulTInt(z, z, TINT_ONE_OVER_PI);
  err := 683;
  z.sgn := sz;

  if inv then begin
    if z.sgn = 0 then begin
      z.sgn := 1;
      AddTInt(z, TINT_ONE_HALF, z);
    end else begin
      AddTInt(z, TINT_ONE_HALF, z);
      z.sgn := 1;
    end;
    err := 346;
  end;

  if x < 0 then begin
    if z.sgn = 0 then begin
      z.sgn := 1;
      AddTInt(z, TINT_ONE, z);
      z.sgn := 1;
    end else begin
      AddTInt(z, TINT_ONE, z);
    end;
    err := 177;
  end;

  Result := TIntToD(z, err, y, x);
end;

// ---------------------------------------------------------------------------
// Atan2piFast — fast path returning (h, l, err). |err| is the relative-error
// bound such that |h+l - atan2pi(y,x)| < err*|h|. Mirrors atan2pi_fast.
// ---------------------------------------------------------------------------

function Atan2piFast(out h, l: Double; y, x: Double): Double;
var
  vy, vx: Tb64u64;
  ay, ax: UInt64;
  inv: Boolean;
  negx, negz: Boolean;
  zh, zl, sz, z, t, hi, lo: Double;
  i, j: Int32;
  pbase: Int32;     // index into cAtan2piPFast = 10*i
  hh, ll: Double;
  tmp: Double;
begin
  vy.f := y;  vx.f := x;
  ay := vy.u shl 1;
  ax := vx.u shl 1;
  inv  := ay > ax;
  negx := (vx.u shr 63) <> 0;
  negz := ((vy.u xor vx.u) shr 63) <> 0;
  if inv then begin
    tmp := y; y := x; x := tmp;
  end;

  // Avoid spurious inexact when |y|=|x|.
  if ax = ay then begin
    if y = x then begin
      if x > 0 then h := Double(0.25) else h := Double(-0.75);
    end else begin
      if x > 0 then h := Double(-0.25) else h := Double(0.75);
    end;
    l := Double(0.0);
    Result := Double(0.0);  // exact
    Exit;
  end;

  A2piFastDiv(zh, zl, y, x);

  sz := Double(1.0);
  if zh < 0 then begin
    zh := -zh;
    zl := -zl;
    sz := Double(-1.0);
  end;

  // 0 <= i <= 64
  i := Int32(Trunc(Double(64.0) * zh));
  if i < 0 then i := 0;
  if i > 64 then i := 64;
  pbase := 10 * i;
  zh := zh - cAtan2piXfast[i].f;
  z  := zh + zl;

  if i = 0 then
    hi := pcr_fma(cAtan2piP9.f, z * z, cAtan2piPFast[pbase + 9].f)
  else
    hi := cAtan2piPFast[pbase + 9].f;

  for j := 6 downto 2 do
    hi := pcr_fma(hi, z, cAtan2piPFast[pbase + 2 + j].f);

  // degree 1: h+l <- h*(zh+zl) + p[2] + p[3]
  A2piSMul(hh, ll, hi, zh, zl);
  A2piFastTwoSum(hi, t, cAtan2piPFast[pbase + 2].f, hh);
  lo := ll + (t + cAtan2piPFast[pbase + 3].f);

  // degree 0: h+l <- (h+l)*(zh+zl) + p[0] + p[1]
  A2piDMul(hh, ll, hi, lo, zh, zl);
  A2piFastTwoSum(hi, t, cAtan2piPFast[pbase + 0].f, hh);
  lo := ll + (t + cAtan2piPFast[pbase + 1].f);
  h := hi; l := lo;

  // multiply by 1/pi
  A2piDMul(hh, ll, h, l, cAtan2piOneOverPiH.f, cAtan2piOneOverPiL.f);
  h := hh; l := ll;

  // restore sign
  h := h * sz;
  l := l * sz;

  // -1/4 < h+l < 1/4
  if inv then begin
    if not negz then
      A2piDDSum1(h, l, Double(0.5), -h, -l)
    else
      A2piDDSum1(h, l, Double(-0.5), -h, -l);
  end;

  if negx then begin
    if not negz then
      A2piDDSum1(h, l, Double(-1.0), h, l)
    else
      A2piDDSum1(h, l, Double(1.0), h, l);
  end;

  Result := cAtan2piErrFast[i].f;
end;

// ---------------------------------------------------------------------------
// pcr_atan2pi — main entry. Fast path + Ziv test, fallback to accurate path.
// Mirrors cr_atan2pi (atan2pi.c lines 790-866).
// ---------------------------------------------------------------------------

function pcr_atan2pi(y, x: Double): Double;
var
  uy, ux: Tb64u64;
  ay, ax: UInt64;
  ey, ex: Int32;
  h, l, err, left, right: Double;
begin
  uy.f := y;  ux.f := x;
  ay := uy.u and cA2pi_MASK;
  ax := ux.u and cA2pi_MASK;
  ey := Int32(ay shr 52);
  ex := Int32(ax shr 52);

  if (ey = $7FF) or (ex = $7FF) or (ay = 0) or (ax = 0) then begin
    Result := Atan2piSpecial(y, x);
    Exit;
  end;

  // Fast-path eligibility: 54 <= ex,ey <= 2044 and |ey-ex| <= 341.
  if (ex >= 54) and (ex <= 2044) and (ey >= 54) and (ey <= 2044)
     and (ey - ex >= -341) and (ey - ex <= 341) then
  begin
    err := Atan2piFast(h, l, y, x);
    {$IFDEF AVX2}
    asm
      movsd        xmm0, h
      movsd        xmm1, err
      movsd        xmm2, l
      vfnmadd213sd xmm0, xmm1, xmm2          // -(h*err) + l
      movsd        left, xmm0
      movsd        xmm0, h
      movsd        xmm1, err
      movsd        xmm2, l
      vfmadd213sd  xmm0, xmm1, xmm2          // h*err + l
      movsd        right, xmm0
    end ['xmm0', 'xmm1', 'xmm2'];
    left  := h + left;
    right := h + right;
    {$ELSE}
    left  := h + pcr_fma(h, -err, l);
    right := h + pcr_fma(h,  err, l);
    {$ENDIF}
    if left = right then begin
      Result := left;
      Exit;
    end;
  end;

  Result := Atan2piAccurate(y, x);
end;
