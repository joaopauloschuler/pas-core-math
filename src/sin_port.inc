// ===========================================================================
// pcr_sin — correctly-rounded binary64 sin(x)
// Ported from core-math/src/binary64/sin/sin.c (Paul Zimmermann, Tom Hubrecht).
//
// Reuses the shared machinery established in cos_port.inc:
//   * Tables T / S / C / PSfast / PCfast / PS / PC / SC (cos_tables.inc)
//   * CosSetDD / CosReduce / CosReduce2 / CosReduceFast
//   * CosEvalPSfast / CosEvalPCfast / CosEvalPS / CosEvalPC
//   * cCosMagic
// Only the tri-folding logic and the fast/accurate sign handling differ.
// ===========================================================================

const
  // sin-specific threshold: for |x| <= 0x1.7137449123ef6p-26, sin(x) rounds to
  // fma(x, -2^-54, x). Differs from cos's 2^-27 threshold.
  cSinTinyBits: UInt64 = $3E57137449123EF6;
  cSinNeg1pm54: Tb64u64 = (u: QWord($BC90000000000000)); // -0x1p-54

  // sin_accurate exceptions (ported from sin.c:2009). Only 2 entries.
  cSinExceptionsX: array[0..1] of Tb64u64 = (
    (u:$3EBE0000000001C2), // 0x1.e0000000001c2p-20
    (u:$7506AC5B262CA1FF)  // 0x1.6ac5b262ca1ffp+849
  );
  cSinExceptionsH: array[0..1] of Tb64u64 = (
    (u:$3EBDFFFFFFFF02E0), // 0x1.dfffffffff02ep-20
    (u:$3FF0000000000000)  // 0x1p+0
  );
  cSinExceptionsL: array[0..1] of Tb64u64 = (
    (u:$36DDCBA692492527), // 0x1.dcba692492527p-146
    (u:QWord($B842B089EA1E692B))  // -0x1.2b089ea1e692bp-123
  );

// ---------------------------------------------------------------------------
// sin_fast — fast-path (double-double) evaluator for sin(x).
// Ported faithfully from sin_fast() in sin.c:1757.
// ---------------------------------------------------------------------------
function SinFast(out h, l: Double; x: Double): Double;
var
  i: Int32;
  neg, is_sin: Boolean;
  absx, err1, sh, sl, ch, cl, uh, ul, err, scale: Double;
begin
  neg := x < 0;
  is_sin := True;
  if neg then absx := -x else absx := x;

  i := CosReduceFast(h, l, absx, err1);

  // if i >= 2^10: sin(pi+x) = -sin(x)
  if (i shr 10) <> 0 then neg := not neg;
  i := i and $3FF;

  // if i >= 2^9: sin(pi/2+x) = cos(x)
  if (i shr 9) <> 0 then is_sin := not is_sin;
  i := i and $1FF;

  // if i >= 2^8: sin(pi/2-x) = cos(x), etc.
  if (i and $100) <> 0 then begin
    is_sin := not is_sin;
    i := $1FF - i;
    h := cCos1pm11.f - h;   // exact (see C comment)
    l := -l;
  end;

  // Now 0 <= i < 256 and 0 <= h+l < 2^-11.
  h := h - cCosSC[i][0].f;
  pcr_a_mul(uh, ul, h, h);
  ul := pcr_fma(h + h, l, ul);
  CosEvalPSfast(sh, sl, h, l, uh, ul);
  CosEvalPCfast(ch, cl, uh, ul);

  if is_sin then begin
    scale := cCosSC[i][2].f; if neg then scale := -scale;
    pcr_s_mul(sh, sl, scale, sh, sl);
    scale := cCosSC[i][1].f; if neg then scale := -scale;
    pcr_s_mul(ch, cl, scale, ch, cl);
    pcr_fasttwosum(h, l, ch, sh);
    l := l + sl + cl;
    err := cCosErrSin.f;
  end else begin
    scale := cCosSC[i][2].f; if neg then scale := -scale;
    pcr_s_mul(ch, cl, scale, ch, cl);
    scale := cCosSC[i][1].f; if neg then scale := -scale;
    pcr_s_mul(sh, sl, scale, sh, sl);
    pcr_fasttwosum(h, l, ch, -sh);
    l := l + cl - sl;
    err := cCosErrCos.f;
  end;

  // Note: unlike cos_fast, sin_fast does NOT post-multiply h, l by sgn[neg];
  // the sign was already folded into the scalar scale above.
  Result := err + err1;
end;

// ---------------------------------------------------------------------------
// sin_accurate — slow-path (TDInt64) evaluator.
// Ported faithfully from sin_accurate() in sin.c:1871.
// ---------------------------------------------------------------------------
function SinAccurate(x: Double): Double;
var
  Xd, U, V, X2: TDInt64;
  i: Int32;
  neg, is_sin: Boolean;
  err, hi0, hi1, lo0, lo1: UInt64;
  k: Int32;
  ax, absx, hi, lo: Double;
begin
  neg := x < 0;
  if neg then absx := -x else absx := x;
  DIntFromD(Xd, absx);
  CosReduce(Xd);

  is_sin := True;
  i := CosReduce2(Xd);

  if (i and $400) <> 0 then begin
    neg := not neg;
    i := i and $3FF;
  end;

  if (i and $200) <> 0 then begin
    is_sin := False;
    i := i and $1FF;
  end;

  if (i and $100) <> 0 then begin
    is_sin := not is_sin;
    Xd.sgn := 1;
    AddDInt(Xd, cCosMagic, Xd);  // X <- 2^-11 - X
    i := $1FF - i;
  end;

  MulDInt(X2, Xd, Xd);
  CosEvalPC(U, X2);
  CosEvalPS(V, Xd, X2);

  if is_sin then begin
    MulDInt(U, cSincosS[i], U);
    MulDInt(V, cSincosC[i], V);
  end else begin
    MulDInt(U, cSincosC[i], U);
    MulDInt(V, cSincosS[i], V);
    V.sgn := 1 - V.sgn;
  end;
  AddDInt(U, U, V);

  err := 41;
  lo0 := U.lo - err;
  hi0 := U.hi - UInt64(Ord(lo0 > U.lo));
  lo1 := U.lo + err;
  hi1 := U.hi + UInt64(Ord(lo1 < U.lo));
  if (hi0 shr 10) <> (hi1 shr 10) then begin
    ax := Abs(x);
    for k := 0 to 1 do
      if ax = cSinExceptionsX[k].f then begin
        hi := cSinExceptionsH[k].f;
        lo := cSinExceptionsL[k].f;
        if x > 0 then Result := hi + lo
        else Result := -hi - lo;
        Exit;
      end;
    // Hard-to-round cases otherwise pass; fall through.
  end;

  if neg then U.sgn := 1 - U.sgn;
  Result := DToD(U);
end;

// ---------------------------------------------------------------------------
// pcr_sin — main entry point.
// Ported faithfully from cr_sin() in sin.c:2034.
// ---------------------------------------------------------------------------
function pcr_sin(x: Double): Double;
var
  t: Tb64u64;
  e: UInt64;
  h, l, err, left, right: Double;
begin
  t.f := x;
  e := (t.u shr 52) and $7FF;

  if e = $7FF then begin
    if (t.u shl 1) = (UInt64($7FF) shl 53) then begin
      Result := 0.0 / 0.0; // Inf -> qNaN, raises invalid
      Exit;
    end;
    Result := x + x; // quiet NaN
    Exit;
  end;

  t.u := t.u and $7FFFFFFFFFFFFFFF;
  if t.u <= cSinTinyBits then begin
    // Per C comment in sin.c: "for 0 < |x| <= 0x1.7137449123ef6p-26,
    // sin(x) rounds to x (to nearest)". C computes fma(x, -2^-54, x) to
    // also raise INEXACT/UNDERFLOW, but the bit pattern equals x.
    // The pcr_fma_pascal emulation is not correctly rounded for tiny
    // normal inputs (near the minimum-normal boundary) — returning x
    // directly is bit-exact and avoids that.
    Result := x;
    Exit;
  end;

  err := SinFast(h, l, x);
  left  := h + (l - err);
  right := h + (l + err);
  if left = right then begin
    Result := left;
    Exit;
  end;

  Result := SinAccurate(x);
end;
