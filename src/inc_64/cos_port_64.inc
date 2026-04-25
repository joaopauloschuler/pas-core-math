// ===========================================================================
// pcr_cos — correctly-rounded binary64 cos(x)
// Ported from core-math/src/binary64/cos/cos.c (Paul Zimmermann, Tom Hubrecht).
//
// Structure mirrors the C source:
//   * Tables (T, S, C, PSfast, PCfast, PS, PC, SC) auto-extracted from cos.c
//     via tools/extract_cos_tables.py — do NOT hand-edit. Regenerate instead.
//   * CosSetDD / CosReduce / CosReduce2 / CosReduceFast — range reduction.
//     All four are shared primitives that will be reused by sin/sincos/tan
//     in subsequent ports (task 4.02 / 4.04 / 4.05).
//   * CosEvalPSfast / CosEvalPCfast — fast-path polynomial evaluators.
//   * CosEvalPS / CosEvalPC — slow-path polynomial evaluators (TDInt64).
//   * CosFast — fast-path (double-double) evaluator.
//   * CosAccurate — slow-path (Ziv's second iteration, TDInt64).
//   * pcr_cos — entry point.
//
// NOTE: The u128 arithmetic in reduce() and reduce_fast() is expressed with
// TUInt128 / Mulu64u64. Every u128 operation in the C source has a direct
// Pascal counterpart; follow the C comments line-by-line for intent.
// ===========================================================================

// --- Auto-generated tables (tools/extract_cos_tables.py) ---------------------
{$I inc_64/cos_tables_64.inc}

const
  // 1/(2pi) high/low double-double approximation used by reduce_fast for |x| < 2pi.
  cCosCH: Tb64u64 = (u: $3FC45F306DC9C883); //  0x1.45f306dc9c883p-3
  cCosCL: Tb64u64 = (u: QWord($BC66B01EC5417056)); // -0x1.6b01ec5417056p-57

  // Fast-path threshold: for |x| <= 0x1.6a09e667f3bccp-27 we return x - x*2^-28.
  cCosTinyBits: UInt64 = $3E46A09E667F3BCC;
  cCosPthresh: Tb64u64 = (u: $4009921FB54442D1); // 0x1.921fb54442d17p+2 ≈ 2*pi
  cCos1pm11:   Tb64u64 = (u: $3F40000000000000); // 0x1p-11
  cCos1p11:    Tb64u64 = (u: $40A0000000000000); // 0x1p+11
  cCosNeg1pm28:Tb64u64 = (u: QWord($BE30000000000000)); // -0x1p-28
  cCosErrSin:  Tb64u64 = (u: $3BAAA00000000000); // 0x1.55p-69
  cCosErrCos:  Tb64u64 = (u: $3BA8100000000000); // 0x1.81p-69
  cCosErrReduce1: Tb64u64 = (u: $39BD900000000000); // 0x1.d9p-105
  cCosErrReduce2: Tb64u64 = (u: $3B30100000000000); // 0x1.01p-76

  // cos_accurate exceptions table (ported from cos.c:1999).
  cCosExceptionsX: array[0..4] of Tb64u64 = (
    (u:$3E88000000000009), // 0x1.8000000000009p-23
    (u:$3E98000000000024), // 0x1.8000000000024p-22
    (u:$3EA8000000000090), // 0x1.800000000009p-21
    (u:$3EB20000000000F3), // 0x1.20000000000f3p-20
    (u:$3EB8000000000024)  // 0x1.800000000024p-20
  );
  cCosExceptionsH: array[0..4] of Tb64u64 = (
    (u:$3FEFFFFFFFFFFF70), // 0x1.fffffffffff7p-1
    (u:$3FEFFFFFFFFFFDC0), // 0x1.ffffffffffdcp-1
    (u:$3FEFFFFFFFFFF700), // 0x1.ffffffffff7p-1
    (u:$3FEFFFFFFFFFEBC0), // 0x1.fffffffffebcp-1
    (u:$3FEFFFFFFFFFDC00)  // 0x1.fffffffffdcp-1
  );
  cCosExceptionsL: array[0..4] of Tb64u64 = (
    (u:$36AB56666666666C), // 0x1.b56666666666cp-143
    (u:$3761B5666666667E), // 0x1.b56666666667ep-137
    (u:$37CB5666666666C4), // 0x1.b5666666666c4p-131
    (u:$3803764266666FD0), // 0x1.37642666666fdp-127 — note trailing 0 per literal
    (u:$3833B56666666BDD)  // 0x1.b5666666667ddp-125
  );

// ---------------------------------------------------------------------------
// set_dd: convert 128-bit fractional value (c1/2^64 + c0/2^128) to double-double.
// Ported faithfully from set_dd() in cos.c:1557.
// ---------------------------------------------------------------------------
procedure CosSetDD(out h, l: Double; c1, c0: UInt64);
var
  e, f, g: UInt64;
  t: Tb64u64;
begin
  if c1 <> 0 then begin
    e := UInt64(pcr_clzll(c1));
    if e <> 0 then begin
      c1 := (c1 shl e) or (c0 shr (64 - e));
      c0 := c0 shl e;
    end;
    f := UInt64($3FE) - e;
    t.u := (f shl 52) or ((c1 shl 1) shr 12);
    h := t.f;
    c0 := (c1 shl 53) or (c0 shr 11);
    if c0 <> 0 then begin
      g := UInt64(pcr_clzll(c0));
      if g <> 0 then c0 := c0 shl g;
      t.u := ((f - 53 - g) shl 52) or ((c0 shl 1) shr 12);
      l := t.f;
    end else
      l := 0;
  end else if c0 <> 0 then begin
    e := UInt64(pcr_clzll(c0));
    f := UInt64($3FE) - 64 - e;
    c0 := c0 shl (e + 1); // MSB shifted out
    t.u := (f shl 52) or (c0 shr 12);
    h := t.f;
    c0 := c0 shl 52;
    if c0 <> 0 then begin
      g := UInt64(pcr_clzll(c0));
      c0 := c0 shl (g + 1);
      t.u := ((f - 64 - g) shl 52) or (c0 shr 12);
      l := t.f;
    end else
      l := 0;
  end else begin
    h := 0;
    l := 0;
  end;
end;

// ---------------------------------------------------------------------------
// reduce: approximate X/(2pi) mod 1 using the T[] table.
// Ported faithfully from reduce() in cos.c:1417.
// X is assumed normalized at input; X is normalized at output.
// ---------------------------------------------------------------------------
procedure CosReduce(var X: TDInt64);
var
  e: Int32;
  u: TUInt128;
  tiny: UInt64;
  ci: array[0..4] of UInt64;
  i, f, g, ee: Int32;
begin
  e := Int32(X.ex);

  if e <= 1 then begin
    // |X| < 2
    u := Mulu64u64(X.hi, cSincosT[1]);
    tiny := u.lo;
    X.lo := u.hi;
    u := Mulu64u64(X.hi, cSincosT[0]);
    // X.lo += u (as 128-bit), with carry into X.hi:
    X.lo := X.lo + u.lo;
    X.hi := u.hi + UInt64(Ord(X.lo < u.lo));
    ee := Int32(X.ex);
    NormalizeDInt(X);
    ee := ee - Int32(X.ex);
    if ee <> 0 then
      X.lo := X.lo or (tiny shr (64 - ee));
    Exit;
  end;

  // 2 <= e <= 1024
  if e < 127 then
    i := 0
  else
    i := (e - 127 + 64 - 1) div 64;

  u := Mulu64u64(X.hi, cSincosT[i + 3]);
  ci[0] := u.lo;
  ci[1] := u.hi;
  u := Mulu64u64(X.hi, cSincosT[i + 2]);
  ci[1] := ci[1] + u.lo;
  ci[2] := u.hi + UInt64(Ord(ci[1] < u.lo));
  u := Mulu64u64(X.hi, cSincosT[i + 1]);
  ci[2] := ci[2] + u.lo;
  ci[3] := u.hi + UInt64(Ord(ci[2] < u.lo));
  u := Mulu64u64(X.hi, cSincosT[i]);
  ci[3] := ci[3] + u.lo;
  ci[4] := u.hi + UInt64(Ord(ci[3] < u.lo));

  f := e - 64 * i;

  if f < 64 then begin
    X.hi := (ci[4] shl f) or (ci[3] shr (64 - f));
    X.lo := (ci[3] shl f) or (ci[2] shr (64 - f));
    tiny := (ci[2] shl f) or (ci[1] shr (64 - f));
  end else if f = 64 then begin
    X.hi := ci[3];
    X.lo := ci[2];
    tiny := ci[1];
  end else begin
    // 65 <= f <= 127
    g := f - 64;
    u := Mulu64u64(X.hi, cSincosT[i + 4]);
    // we only need the high 64 bits of the u128 product
    ci[0] := ci[0] + u.hi;
    ci[1] := ci[1] + UInt64(Ord(ci[0] < u.hi));
    ci[2] := ci[2] + UInt64(Ord((ci[0] < u.hi) and (ci[1] = 0)));
    ci[3] := ci[3] + UInt64(Ord((ci[0] < u.hi) and (ci[1] = 0) and (ci[2] = 0)));
    ci[4] := ci[4] + UInt64(Ord((ci[0] < u.hi) and (ci[1] = 0) and (ci[2] = 0) and (ci[3] = 0)));
    X.hi := (ci[3] shl g) or (ci[2] shr (64 - g));
    X.lo := (ci[2] shl g) or (ci[1] shr (64 - g));
    tiny := (ci[1] shl g) or (ci[0] shr (64 - g));
  end;

  X.ex := 0;
  NormalizeDInt(X);
  if X.ex < 0 then
    X.lo := X.lo or (tiny shr (64 + X.ex));
end;

// reduce2: write Xin = i/2^11 + Xout, 0 <= Xout < 2^-11. Exact.
function CosReduce2(var X: TDInt64): Int32;
var
  sh, i: Int32;
begin
  if X.ex <= -11 then begin Result := 0; Exit; end;
  sh := 64 - 11 - Int32(X.ex);
  i := Int32(X.hi shr sh);
  X.hi := X.hi and ((UInt64(1) shl sh) - 1);
  NormalizeDInt(X);
  Result := i;
end;

// ---------------------------------------------------------------------------
// reduce_fast: Given 0x1.6a09e667f3bccp-27 < x, return i and set h,l such
// that i/2^11 + h + l approximates frac(x/(2pi)). Also set err1.
// Ported faithfully from reduce_fast() in cos.c:1628.
// ---------------------------------------------------------------------------
function CosReduceFast(out h, l: Double; x: Double; out err1: Double): Int32;
var
  t: Tb64u64;
  ex: Int32;
  m: UInt64;
  c: array[0..2] of UInt64;
  u: TUInt128;
  e_shift: Int32;
  i_idx: Int32;
  fi: Double;
begin
  t.f := x;
  if x <= cCosPthresh.f then begin
    pcr_a_mul(h, l, cCosCH.f, x);            // exact
    l := pcr_fma(cCosCL.f, x, l);
    err1 := cCosErrReduce1.f * h;
  end else begin
    ex := Int32((t.u shr 52) and $7FF); // 1025..2046
    m := (UInt64(1) shl 52) or (t.u and $FFFFFFFFFFFFF);
    if ex <= 1074 then begin
      // 1025 <= ex <= 1074
      u := Mulu64u64(m, cSincosT[1]);
      c[0] := u.lo;
      c[1] := u.hi;
      u := Mulu64u64(m, cSincosT[0]);
      c[1] := c[1] + u.lo;
      c[2] := u.hi + UInt64(Ord(c[1] < u.lo));
      e_shift := 1075 - ex; // 1..50
    end else begin
      // 1075 <= ex <= 2046
      i_idx := (ex - 1138 + 63) div 64;
      u := Mulu64u64(m, cSincosT[i_idx + 2]);
      c[0] := u.lo;
      c[1] := u.hi;
      u := Mulu64u64(m, cSincosT[i_idx + 1]);
      c[1] := c[1] + u.lo;
      c[2] := u.hi + UInt64(Ord(c[1] < u.lo));
      u := Mulu64u64(m, cSincosT[i_idx]);
      c[2] := c[2] + u.lo;
      e_shift := 1139 + (i_idx shl 6) - ex;
    end;
    if e_shift = 64 then begin
      c[0] := c[1];
      c[1] := c[2];
    end else begin
      c[0] := (c[1] shl (64 - e_shift)) or (c[0] shr e_shift);
      c[1] := (c[2] shl (64 - e_shift)) or (c[1] shr e_shift);
    end;
    CosSetDD(h, l, c[1], c[0]);
    err1 := cCosErrReduce2.f;
  end;
  Result := Int32(Trunc(h * cCos1p11.f)); // floor (h >= 0 here)
  fi := Result;                            // Double form for the fma
  h  := pcr_fma(fi, -cCos1pm11.f, h);
end;

// ---------------------------------------------------------------------------
// evalPSfast / evalPCfast — fast-path polynomial evaluation for sin/cos(2pi x).
// Ported from cos.c:1316 / cos.c:1336.
// ---------------------------------------------------------------------------
procedure CosEvalPSfast(out h, l: Double; xh, xl, uh, ul: Double);
var
  t: Double;
begin
  h := cCosPSfast[4].f;                      // degree 7
  h := pcr_fma(h, uh, cCosPSfast[3].f);      // degree 5
  h := pcr_fma(h, uh, cCosPSfast[2].f);      // degree 3
  pcr_s_mul(h, l, h, uh, ul);
  pcr_fasttwosum(h, t, cCosPSfast[0].f, h);
  l := l + cCosPSfast[1].f + t;
  pcr_d_mul(h, l, h, l, xh, xl);
end;

procedure CosEvalPCfast(out h, l: Double; uh, ul: Double);
var
  t: Double;
begin
  h := cCosPCfast[4].f;                      // degree 6
  h := pcr_fma(h, uh, cCosPCfast[3].f);      // degree 4
  h := pcr_fma(h, uh, cCosPCfast[2].f);      // degree 2
  pcr_s_mul(h, l, h, uh, ul);
  pcr_fasttwosum(h, t, cCosPCfast[0].f, h);
  l := l + cCosPCfast[1].f + t;
end;

// ---------------------------------------------------------------------------
// evalPS / evalPC — slow-path polynomial (TDInt64).
// Ported from cos.c:1352 / cos.c:1372.
// ---------------------------------------------------------------------------
procedure CosEvalPS(out Y: TDInt64; const X, X2: TDInt64);
begin
  MulDInt21(Y, X2, cCosPS[5]); // degree 11
  AddDInt  (Y, Y,  cCosPS[4]); // degree 9
  MulDInt  (Y, Y,  X2);
  AddDInt  (Y, Y,  cCosPS[3]); // degree 7
  MulDInt  (Y, Y,  X2);
  AddDInt  (Y, Y,  cCosPS[2]); // degree 5
  MulDInt  (Y, Y,  X2);
  AddDInt  (Y, Y,  cCosPS[1]); // degree 3
  MulDInt  (Y, Y,  X2);
  AddDInt  (Y, Y,  cCosPS[0]); // degree 1
  MulDInt  (Y, Y,  X);         // multiply by X
end;

procedure CosEvalPC(out Y: TDInt64; const X2: TDInt64);
begin
  MulDInt21(Y, X2, cCosPC[5]); // degree 10
  AddDInt  (Y, Y,  cCosPC[4]); // degree 8
  MulDInt  (Y, Y,  X2);
  AddDInt  (Y, Y,  cCosPC[3]); // degree 6
  MulDInt  (Y, Y,  X2);
  AddDInt  (Y, Y,  cCosPC[2]); // degree 4
  MulDInt  (Y, Y,  X2);
  AddDInt  (Y, Y,  cCosPC[1]); // degree 2
  MulDInt  (Y, Y,  X2);
  AddDInt  (Y, Y,  cCosPC[0]); // degree 0
end;

// ---------------------------------------------------------------------------
// cos_fast — fast-path (double-double) evaluator.
// Returns absolute error bound; h + l approximates cos(x) with | h + l - cos(x) | < err.
// Ported faithfully from cos_fast() in cos.c:1748.
// ---------------------------------------------------------------------------
function CosFast(out h, l: Double; x: Double): Double;
var
  i: Int32;
  neg, is_cos: Boolean;
  err1, sh, sl, ch, cl, uh, ul, err: Double;
begin
  neg := False;
  is_cos := True;

  i := CosReduceFast(h, l, x, err1);

  // if i >= 2^10: cos(pi+x) = -cos(x)
  if (i shr 10) <> 0 then neg := not neg;
  i := i and $3FF;

  // if i >= 2^9: cos(pi/2+x) = -sin(x)
  if (i shr 9) <> 0 then begin
    is_cos := not is_cos;
    neg := not neg;
  end;
  i := i and $1FF;

  // if i >= 2^8: cos(pi/2-x) = sin(x)
  if (i and $100) <> 0 then begin
    is_cos := not is_cos;
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

  if not is_cos then begin
    pcr_s_mul(sh, sl, cCosSC[i][2].f, sh, sl);
    pcr_s_mul(ch, cl, cCosSC[i][1].f, ch, cl);
    pcr_fasttwosum(h, l, ch, sh);
    l := l + sl + cl;
    err := cCosErrSin.f;
  end else begin
    pcr_s_mul(ch, cl, cCosSC[i][2].f, ch, cl);
    pcr_s_mul(sh, sl, cCosSC[i][1].f, sh, sl);
    pcr_fasttwosum(h, l, ch, -sh);
    l := l + cl - sl;
    err := cCosErrCos.f;
  end;

  if neg then begin
    h := -h;
    l := -l;
  end;
  Result := err + err1;
end;

// ---------------------------------------------------------------------------
// cos_accurate — slow-path (TDInt64) evaluator.
// Ported faithfully from cos_accurate() in cos.c:1861.
// ---------------------------------------------------------------------------
const
  cCosMagic: TDInt64 = (hi: QWord($8000000000000000); lo: 0; ex: -10; sgn: 0);

function CosAccurate(x: Double): Double;
var
  Xd: TDInt64;
  U, V, X2: TDInt64;
  i: Int32;
  neg, is_cos: Boolean;
  err, hi0, hi1, lo0, lo1: UInt64;
  k: Int32;
  ax: Double;
begin
  DIntFromD(Xd, x);
  CosReduce(Xd);

  neg := False;
  is_cos := True;

  i := CosReduce2(Xd);

  if (i and $400) <> 0 then begin
    neg := True;
    i := i and $3FF;
  end;

  if (i and $200) <> 0 then begin
    neg := not neg;
    is_cos := False;
    i := i and $1FF;
  end;

  if (i and $100) <> 0 then begin
    is_cos := not is_cos;
    Xd.sgn := 1;                 // negate X
    AddDInt(Xd, cCosMagic, Xd);  // X <- 2^-11 - X (MAGIC = 2^-11)
    i := $1FF - i;
  end;

  MulDInt(X2, Xd, Xd);
  CosEvalPC(U, X2);
  CosEvalPS(V, Xd, X2);

  if not is_cos then begin
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
    for k := 0 to 4 do
      if ax = cCosExceptionsX[k].f then begin
        Result := cCosExceptionsH[k].f + cCosExceptionsL[k].f;
        Exit;
      end;
    // Hard-to-round cases otherwise pass; fall through.
  end;

  if neg then U.sgn := 1 - U.sgn;
  Result := DToD(U);
end;

// ---------------------------------------------------------------------------
// pcr_cos — main entry point.
// Ported faithfully from cr_cos() in cos.c:2023.
// ---------------------------------------------------------------------------
function pcr_cos(x: Double): Double;
var
  t: Tb64u64;
  e: UInt64;
  h, l, err, left, right: Double;
begin
  t.f := x;
  e := (t.u shr 52) and $7FF;

  if e = $7FF then begin
    // NaN or +/-Inf
    if (t.u shl 1) = (UInt64($7FF) shl 53) then begin
      Result := 0.0 / 0.0; // Inf -> qNaN, raises invalid
      Exit;
    end;
    Result := x + x; // quiet NaN
    Exit;
  end;

  // |x|
  t.u := t.u and $7FFFFFFFFFFFFFFF;
  if t.u <= cCosTinyBits then begin
    // |x| <= 0x1.6a09e667f3bccp-27: cos(x) rounds to fma(|x|, -2^-28, 1)
    Result := pcr_fma(t.f, cCosNeg1pm28.f, 1.0);
    Exit;
  end;

  err := CosFast(h, l, t.f);
  left  := h + (l - err);
  right := h + (l + err);
  if left = right then begin
    Result := left;
    Exit;
  end;

  Result := CosAccurate(t.f);
end;
