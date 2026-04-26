// ===========================================================================
// pcr_sincos — correctly-rounded binary64 sin(x)/cos(x) combined.
// Ported from core-math/src/binary64/sincos/sincos.c.
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
// Reuses CosReduce / CosReduceFast / CosEval* primitives from cos_port.inc,
// SinAccurate from sin_port.inc, and CosAccurate from cos_port.inc.
// ===========================================================================

// ---------------------------------------------------------------------------
// sincos_fast — computes h+l ≈ sin(x) and hc+lc ≈ cos(x) in double-double.
// Returns the absolute error bound for both.
// ---------------------------------------------------------------------------
function SinCosFast(out h, l, hc, lc: Double; x: Double): Double;
var
  i: Int32;
  neg, negc: Boolean;
  is_sin: Boolean;
  absx, err1, sh, sl, ch, cl, uh, ul, scale: Double;
  Sh_, Sl_, Ch_, Cl_: Double;
begin
  neg := x < 0;
  negc := False;
  is_sin := True;
  if neg then absx := -x else absx := x;

  i := CosReduceFast(h, l, absx, err1);

  // if i >= 2^10: sin(pi+x) = -sin(x), cos(pi+x) = -cos(x)
  if (i shr 10) <> 0 then begin
    neg := not neg;
    negc := not negc;
  end;
  i := i and $3FF;

  // if i >= 2^9: sin(pi/2+x) = cos(x), cos(pi/2+x) = -sin(x)
  if (i shr 9) <> 0 then begin
    is_sin := not is_sin;
    negc := not negc;
  end;
  i := i and $1FF;

  // if i >= 2^8: sin(pi/2-x) = cos(x), cos(pi/2-x) = sin(x)
  if (i and $100) <> 0 then begin
    is_sin := not is_sin;
    i := $1FF - i;
    h := cCos1pm11.f - h;   // exact
    l := -l;
  end;

  // Now 0 <= i < 256 and 0 <= h+l < 2^-11.
  h := h - cCosSC[i][0].f;
  pcr_a_mul(uh, ul, h, h);
  ul := pcr_fma(h + h, l, ul);
  CosEvalPSfast(sh, sl, h, l, uh, ul);
  CosEvalPCfast(ch, cl, uh, ul);

  if is_sin then begin
    // sin(x) ~ SC[i][2] * (sh+sl) + SC[i][1] * (ch+cl)
    // cos(x) ~ SC[i][2] * (ch+cl) - SC[i][1] * (sh+sl)
    scale := cCosSC[i][2].f; if neg then scale := -scale;
    pcr_s_mul(Sh_, Sl_, scale, sh, sl);
    scale := cCosSC[i][1].f; if neg then scale := -scale;
    pcr_s_mul(Ch_, Cl_, scale, ch, cl);
    pcr_fasttwosum(h, l, Ch_, Sh_);
    l := l + Sl_ + Cl_;

    scale := cCosSC[i][2].f; if negc then scale := -scale;
    pcr_s_mul(Ch_, Cl_, scale, ch, cl);
    scale := cCosSC[i][1].f; if negc then scale := -scale;
    pcr_s_mul(Sh_, Sl_, scale, sh, sl);
    pcr_fasttwosum(hc, lc, Ch_, -Sh_);
    lc := lc + Cl_ - Sl_;
  end else begin
    // sin(x) ~ SC[i][2] * (ch+cl) - SC[i][1] * (sh+sl)
    // cos(x) ~ -SC[i][2] * (sh+sl) - SC[i][1] * (ch+cl)  [when negc, signs flip]
    scale := cCosSC[i][2].f; if neg then scale := -scale;
    pcr_s_mul(Ch_, Cl_, scale, ch, cl);
    scale := cCosSC[i][1].f; if neg then scale := -scale;
    pcr_s_mul(Sh_, Sl_, scale, sh, sl);
    pcr_fasttwosum(h, l, Ch_, -Sh_);
    l := l + Cl_ - Sl_;

    scale := cCosSC[i][2].f; if negc then scale := -scale;
    pcr_s_mul(Sh_, Sl_, scale, sh, sl);
    scale := cCosSC[i][1].f; if negc then scale := -scale;
    pcr_s_mul(Ch_, Cl_, scale, ch, cl);
    pcr_fasttwosum(hc, lc, Ch_, Sh_);
    lc := lc + Sl_ + Cl_;
  end;

  Result := cCosErrCos.f + err1; // 0x1.81p-69 + err1 (bounds both)
end;

// ---------------------------------------------------------------------------
// pcr_sincos — main entry point.
// Ported faithfully from cr_sincos() in sincos.c:2190.
// ---------------------------------------------------------------------------
procedure pcr_sincos(x: Double; out s, c: Double);
var
  t: Tb64u64;
  e: UInt64;
  h, l, hc, lc, err, right, right_c: Double;
begin
  t.f := x;
  e := (t.u shr 52) and $7FF;

  if e = $7FF then begin
    if (t.u shl 1) <> (UInt64($7FF8) shl 49) then begin
      // Inf or signaling NaN — produce qNaN
      s := Double(0.0) / Double(0.0);
      c := s;
      Exit;
    end;
    // Quiet NaN: propagate
    s := x + x;
    c := s;
    Exit;
  end;

  t.u := t.u and $7FFFFFFFFFFFFFFF;
  if t.u <= cCosTinyBits then begin
    // |x| <= 0x1.6a09e667f3bccp-27: sin rounds to x, cos rounds to 1.0.
    // (C computes 1.0 - 2^-54 for nonzero tiny, but that ties-to-even to 1.0.)
    s := x;
    c := Double(1.0);
    Exit;
  end;

  err := SinCosFast(h, l, hc, lc, x);
  s := h + (l - err);
  right := h + (l + err);
  c := hc + (lc - err);
  right_c := hc + (lc + err);
  if (s = right) and (c = right_c) then Exit;

  if c = right_c then begin
    // Fast path succeeded for cos; refine sin only.
    s := SinAccurate(x);
    Exit;
  end;

  // Accurate cos: CosAccurate expects x > 0.
  if x > 0 then c := CosAccurate(x) else c := CosAccurate(-x);

  if s = right then Exit; // fast path succeeded for sin

  s := SinAccurate(x);
end;
