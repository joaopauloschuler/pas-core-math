// ===========================================================================
// pcr_tan — correctly-rounded binary64 tan(x)
// Ported from core-math/src/binary64/tan/tan.c.
//
// Shares the T, S, C, PSfast, PCfast, PS, PC, SC tables with cos_port.inc
// but needs a *higher-precision* reduce_fast (keeps one extra limb of the
// 2/π product), because it divides by cos2pi(R) which can vanish near π/2.
// Also uses its own fast_div for the fast path and div_dint for the slow.
// ===========================================================================

// ---------------------------------------------------------------------------
// reduce_fast (tan.c:1747) — like CosReduceFast but one extra T[] limb and
// no err1 output. |h + l - frac(x/(2π))| < 2^-104.815.
// ---------------------------------------------------------------------------
function TanReduceFast(out h, l: Double; x: Double): Int32;
var
  t: Tb64u64;
  ex: Int32;
  m: UInt64;
  c: array[0..2] of UInt64;
  u, v: TUInt128;
  e_shift: Int32;
  i_idx: Int32;
  fi: Double;
begin
  t.f := x;
  if x <= cCosPthresh.f then begin
    pcr_a_mul(h, l, cCosCH.f, x);
    l := pcr_fma(cCosCL.f, x, l);
  end else begin
    ex := Int32((t.u shr 52) and $7FF);
    m := (UInt64(1) shl 52) or (t.u and $FFFFFFFFFFFFF);
    if ex <= 1074 then begin
      // 1025 <= ex <= 1074
      v := Mulu64u64(m, cSincosT[2]);
      u := Mulu64u64(m, cSincosT[1]);
      // c[0] = u.lo + v.hi  (low 64 of u+(v>>64))
      c[0] := u.lo + v.hi;
      c[1] := u.hi + UInt64(Ord(c[0] < u.lo));
      u := Mulu64u64(m, cSincosT[0]);
      c[1] := c[1] + u.lo;
      c[2] := u.hi + UInt64(Ord(c[1] < u.lo));
      e_shift := 1075 - ex; // 1..50
    end else begin
      // 1075 <= ex <= 2046
      i_idx := (ex - 1138 + 63) div 64;
      v := Mulu64u64(m, cSincosT[i_idx + 3]);
      u := Mulu64u64(m, cSincosT[i_idx + 2]);
      c[0] := u.lo + v.hi;
      c[1] := u.hi + UInt64(Ord(c[0] < u.lo));
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
  end;
  Result := Int32(Trunc(h * cCos1p11.f));
  fi := Result;
  h := pcr_fma(fi, -cCos1pm11.f, h);
end;

// ---------------------------------------------------------------------------
// fast_div (tan.c:1888): (h+l) = (bh+bl)/(ah+al) with relative error < 2^-96.99.
// ---------------------------------------------------------------------------
procedure TanFastDiv(out h, l: Double; bh, bl, ah, al: Double);
var
  y, eh, el: Double;
begin
  y := 1.0 / ah;
  h := bh * y;
  eh := pcr_fma(ah, -h, bh);
  el := pcr_fma(al, -h, bl);
  l := y * (eh + el);
end;

// ---------------------------------------------------------------------------
// div_dint (from cos.c's dint ops, but tan needs the full div) — we need a
// division of two TDInt64 values. Port from sin/tan's div_dint.
// ---------------------------------------------------------------------------
// NOTE: reuse-checked below.

const
  cTan2pm37:     Tb64u64 = (u: $3DA0000000000000); // 0x1p-37
  cTanErrRelMul: Tb64u64 = (u: $3BD1A00000000000); // 0x1.1ap-66

// ---------------------------------------------------------------------------
// tan_fast — fast-path evaluator for tan(x) in double-double.
// ---------------------------------------------------------------------------
function TanFast(out h, l: Double; x: Double): Double;
var
  i: Int32;
  neg, is_tan: Boolean;
  absx, sh, sl, ch, cl, uh, ul: Double;
  sh0, sl0, ch0, cl0, h1, l1, h2, l2: Double;
begin
  neg := x < 0;
  is_tan := True;
  if neg then absx := -x else absx := x;

  i := TanReduceFast(h, l, absx);

  // tan(pi+x) = tan(x): i >> 10 not used (no sign flip)
  i := i and $3FF;

  // if i >= 2^9: tan(pi/2+x) = -cot(x)
  if (i shr 9) <> 0 then begin
    is_tan := not is_tan;
    neg := not neg;
  end;
  i := i and $1FF;

  // if i >= 2^8: tan(pi/2-x) = cot(x)
  if (i and $100) <> 0 then begin
    is_tan := not is_tan;
    i := $1FF - i;
    h := cCos1pm11.f - h;   // exact
    l := -l;
  end;

  // Defer to slow path: i=0 and h < 2^-37
  if (i = 0) and (h < cTan2pm37.f) then begin
    Result := 1.0; // sentinel: force slow path
    Exit;
  end;

  h := h - cCosSC[i][0].f;
  pcr_fasttwosum(h, l, h, l);
  pcr_a_mul(uh, ul, h, h);
  ul := pcr_fma(h + h, l, ul);
  CosEvalPSfast(sh, sl, h, l, uh, ul);
  CosEvalPCfast(ch, cl, uh, ul);

  pcr_s_mul(sh0, sl0, cCosSC[i][2].f, sh, sl);
  pcr_s_mul(ch0, cl0, cCosSC[i][1].f, ch, cl);
  pcr_fasttwosum(h1, l1, ch0, sh0);
  l1 := l1 + sl0 + cl0;

  pcr_s_mul(ch, cl, cCosSC[i][2].f, ch, cl);
  pcr_s_mul(sh, sl, cCosSC[i][1].f, sh, sl);
  pcr_fasttwosum(h2, l2, ch, -sh);
  l2 := l2 + cl - sl;

  if is_tan then
    TanFastDiv(h, l, h1, l1, h2, l2)
  else
    TanFastDiv(h, l, h2, l2, h1, l1);

  if neg then begin
    h := -h;
    l := -l;
  end;
  // err = h * 2^-65.864 < h * 0x1.1ap-66
  Result := h * cTanErrRelMul.f;
end;

// ---------------------------------------------------------------------------
// Tinv[256]: 64-bit inverse seed table used by InvDInt.
// ---------------------------------------------------------------------------
const
{$I tan_tinv.inc}
  cTanDIntOne: TDInt64 = (hi: $8000000000000000; lo: 0; ex: 1; sgn: 0);

// Subtract u128 "hi128:lo128" = 2^127 from a u64*u64 product, returning 128-bit.
// Used to compute 2^127 - h*t.
procedure U128SubProd(out r: TUInt128; const base: TUInt128; const prod: TUInt128); inline;
begin
  SubU128(r, base, prod);
end;

// r = (u128) t * (some u64)
procedure U128MulU64(out r: TUInt128; t, x: UInt64); inline;
begin
  r := Mulu64u64(t, x);
end;

// ---------------------------------------------------------------------------
// inv_dint (tan.c:333) — compute r ≈ 1/a, relative error < 2^-124.999.
// ---------------------------------------------------------------------------
procedure InvDInt(out r: TDInt64; const a: TDInt64);
var
  h, t: UInt64;
  i: Int32;
  e: TUInt128;
  base: TUInt128;  // 2^127
  ht: TUInt128;
  q: TDInt64;
  a_local: TDInt64;
begin
  a_local := a; // in case r aliases a
  h := a_local.hi; // 2^63 <= h < 2^64

  i := Int32((h shr 55) and $FF);
  t := cTanTinv[i];

  base.hi := UInt64($8000000000000000); base.lo := 0; // 2^127

  // first Newton iteration
  ht := Mulu64u64(h, t);
  SubU128(e, base, ht);                  // e = 2^127 - h*t  (< 2^119)
  // e >> 55: since e < 2^119, (e>>55) < 2^64 — fits in 64-bit.
  // For u128 shifted right by 55 (<64): result_lo = (hi<<9) | (lo>>55), result_hi = hi>>55.
  //   Since e < 2^119, e.hi < 2^55, so hi>>55 == 0. Result fits in 64 bits:
  e := Mulu64u64(t, (e.hi shl 9) or (e.lo shr 55));
  // e is now u128 (t * (prev_e>>55)), t += e>>72.
  // e>>72 = e.hi >> 8  (since 72>=64, we take high 64 of e shifted right (72-64)=8).
  t := t + (e.hi shr 8);

  // second Newton iteration
  ht := Mulu64u64(h, t);
  SubU128(e, base, ht);                  // e < 2^111
  // e>>47 fits in 64 bits since e < 2^111, e>>47 < 2^64.
  // result_lo = (e.hi<<17) | (e.lo>>47); high part dropped.
  e := Mulu64u64(t, (e.hi shl 17) or (e.lo shr 47));
  // t += e>>80 = e.hi >> 16
  t := t + (e.hi shr 16);

  // third Newton iteration
  ht := Mulu64u64(h, t);
  SubU128(e, base, ht);                  // e < 2^95
  // e>>31 fits in 64 bits since e < 2^95, e>>31 < 2^64.
  e := Mulu64u64(t, (e.hi shl 33) or (e.lo shr 31));
  // t += e>>96 = e.hi >> 32
  t := t + (e.hi shr 32);

  r.hi := t;
  r.lo := 0;
  r.ex := 1 - a_local.ex;
  r.sgn := 1;

  // q = -a*r (mul_dint_21 since r.lo=0 -> a_local is b with low zero? actually mul_dint_21
  // expects b.lo=0; here r.lo=0 so call MulDInt21(q, a, r).)
  MulDInt21(q, a_local, r);      // q = a * r, with r.sgn=1 giving q.sgn = a.sgn XOR 1
  r.sgn := 0;                     // restore r's sign to positive
  AddDInt(q, cTanDIntOne, q);    // q = 1 - a*r  (since q was -a*r)
  MulDInt(q, r, q);               // q = r * (1 - a*r)
  AddDInt(r, r, q);               // r = r + r*(1-a*r)
end;

// ---------------------------------------------------------------------------
// div_dint (tan.c:429) — r = b / a, relative error < 2^-123.67.
// ---------------------------------------------------------------------------
procedure DivDInt(out r: TDInt64; const b, a: TDInt64);
begin
  InvDInt(r, a);
  MulDInt(r, r, b);
end;

// ---------------------------------------------------------------------------
// tan_accurate — slow path.
// ---------------------------------------------------------------------------
const
  cTanErr: UInt64 = 86;

  cTanExceptionsX: array[0..1] of Tb64u64 = (
    (u:$3E8DFFFFFFFFFF1F), // 0x1.dffffffffff1fp-22
    (u:$3EADFFFFFFFFFC7C)  // 0x1.dfffffffffc7cp-21
  );
  cTanExceptionsH: array[0..1] of Tb64u64 = (
    (u:$3E8E000000000151), // 0x1.e000000000151p-22
    (u:$3EAE000000000546)  // 0x1.e000000000546p-21
  );
  cTanExceptionsL: array[0..1] of Tb64u64 = (
    (u:$3B3FFFFFFFFFFFFF), // 0x1.fffffffffffffp-76
    (u:$B74658BCEDB6E1D4)  // -0x1.658bcedb6e1d4p-147
  );

function TanAccurate(x: Double): Double;
var
  Xd, U, V, X2, UU, VV, Sin_, Cos_: TDInt64;
  i: Int32;
  neg, is_tan: Boolean;
  err, hi0, hi1, lo0, lo1: UInt64;
  k: Int32;
  ax, absx, hi, lo: Double;
begin
  neg := x < 0;
  if neg then absx := -x else absx := x;
  DIntFromD(Xd, absx);
  CosReduce(Xd);

  is_tan := True;
  i := CosReduce2(Xd);

  i := i and $3FF; // tan(pi+x) = tan(x)

  if (i and $200) <> 0 then begin
    is_tan := False;
    neg := not neg;
    i := i and $1FF;
  end;

  if (i and $100) <> 0 then begin
    is_tan := not is_tan;
    Xd.sgn := 1;
    AddDInt(Xd, cCosMagic, Xd);
    i := $1FF - i;
  end;

  MulDInt(X2, Xd, Xd);
  CosEvalPC(U, X2);
  CosEvalPS(V, Xd, X2);

  // sin2pi(R) via angle-addition
  MulDInt(UU, cSincosS[i], U);
  MulDInt(VV, cSincosC[i], V);
  AddDInt(Sin_, UU, VV);

  // cos2pi(R) via angle-addition
  MulDInt(U, cSincosC[i], U);
  MulDInt(V, cSincosS[i], V);
  V.sgn := 1 - V.sgn;
  AddDInt(Cos_, U, V);

  if is_tan then DivDInt(U, Sin_, Cos_)
  else            DivDInt(U, Cos_, Sin_);

  err := cTanErr;
  lo0 := U.lo - err;
  hi0 := U.hi - UInt64(Ord(lo0 > U.lo));
  lo1 := U.lo + err;
  hi1 := U.hi + UInt64(Ord(lo1 < U.lo));
  if (hi0 shr 10) <> (hi1 shr 10) then begin
    ax := Abs(x);
    for k := 0 to 1 do
      if ax = cTanExceptionsX[k].f then begin
        hi := cTanExceptionsH[k].f;
        lo := cTanExceptionsL[k].f;
        if x > 0 then Result := hi + lo
        else Result := -hi - lo;
        Exit;
      end;
  end;

  if neg then U.sgn := 1 - U.sgn;
  Result := DToD(U);
end;

// ---------------------------------------------------------------------------
// pcr_tan — main entry point.
// ---------------------------------------------------------------------------
function pcr_tan(x: Double): Double;
const
  cTanTinyBits: UInt64 = $3E4D12ED0AF1A27E;
  cTan1pm54:    Tb64u64 = (u: $3C90000000000000); // 0x1p-54
var
  t: Tb64u64;
  e: UInt64;
  h, l, err, left, right: Double;
begin
  t.f := x;
  e := (t.u shr 52) and $7FF;

  if e = $7FF then begin
    if (t.u shl 1) = (UInt64($7FF) shl 53) then begin
      Result := 0.0 / 0.0; // Inf -> qNaN
      Exit;
    end;
    Result := x + x;
    Exit;
  end;

  t.u := t.u and $7FFFFFFFFFFFFFFF;
  if t.u <= cTanTinyBits then begin
    // Tiny: tan(x) ~ x + x^3/3, rounds to fma(x, 2^-54, x).
    // For our non-AVX2 pcr_fma, similar subnormal bug as sin — return x directly.
    Result := x;
    Exit;
  end;

  err := TanFast(h, l, x);
  left := h + (l - err);
  right := h + (l + err);
  if left = right then begin
    Result := left;
    Exit;
  end;

  Result := TanAccurate(x);
end;
