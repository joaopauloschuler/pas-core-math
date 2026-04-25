// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
{$I pascoremath.inc}
unit pascoremathtypes;

interface

uses Math;

type
  TUInt128 = record
    lo, hi: UInt64;
  end;

  Tb32u32 = record
    case boolean of
      false: (f: Single);
      true:  (u: UInt32);
  end;

  Tb64u64 = record
    case boolean of
      false: (f: Double);
      true:  (u: UInt64);
  end;

  // 128-bit significand type used as slow path in every binary64 function.
  // Fields match the little-endian C union in sin.c / dint.h.
  TDInt64 = record
    hi:  UInt64;   // high 64 bits of significand (MSB always 1 when non-zero)
    lo:  UInt64;   // low 64 bits of significand
    ex:  Int64;    // binary exponent (signed)
    sgn: Byte;     // sign: 0 = positive, 1 = negative
  end;

  // 256-bit significand type used exclusively by pow (third Ziv iteration).
  // Arithmetic ported separately in Phase 5.01.
  TQInt64 = record
    r0, r1, r2, r3: UInt64;  // 256-bit significand, r0 = most significant
    ex:             Int64;
    sgn:            Byte;
  end;

  // 192-bit significand type used by atan2 / atan2pi (Phase 2.07 / 2.11).
  // Field order matches the little-endian C union in atan2/tint.h:
  //   m at offset 0 (low 64 of u128 _h), h at offset 8 (high 64 of u128 _h),
  //   l at offset 16 (= u64 _l). Value = (-1)^sgn * (h/2^64 + m/2^128 + l/2^192) * 2^ex.
  // h is the most significant limb; when non-zero, MSB of h is 1 (normalized).
  TInt64 = record
    m:   UInt64;   // middle 64 bits of significand
    h:   UInt64;   // high 64 bits of significand (MSB always 1 when non-zero)
    l:   UInt64;   // low 64 bits of significand
    ex:  Int64;    // binary exponent
    sgn: UInt64;   // sign: 0 = positive, 1 = negative
  end;

operator +(const a: TUInt128; b: UInt64): TUInt128; inline;

procedure AddU128(out r: TUInt128; const a, b: TUInt128); inline;
procedure SubU128(out r: TUInt128; const a, b: TUInt128); inline;
procedure ShlU128(var a: TUInt128; sh: Integer); inline;
procedure ShrU128(var a: TUInt128; sh: Integer); inline;

function Mulu64u64(a, b: UInt64): TUInt128; inline;

// ------- dint64_t arithmetic -------
// All operations ported faithfully from core-math/src/binary64/sin/sin.c.

// Copy
procedure CpDInt(out r: TDInt64; const a: TDInt64); inline;
// True if value is zero (hi = 0)
function DIntZeroP(const a: TDInt64): Boolean; inline;
// Compare absolute values: -1 / 0 / +1
function CmpDIntAbs(const a, b: TDInt64): Integer;
// Add two TDInt64 values (error bounded by 2 ulp_128)
procedure AddDInt(out r: TDInt64; const a, b: TDInt64);
// Multiply two TDInt64 values (error bounded by 6 ulp_128)
procedure MulDInt(out r: TDInt64; const a, b: TDInt64);
// Multiply two TDInt64 values, assuming b.lo = 0 (error bounded by 2 ulp_128)
procedure MulDInt21(out r: TDInt64; const a, b: TDInt64);
// Multiply a TDInt64 by a signed 64-bit integer (ported from mul_dint_2 in
// core-math/src/binary64/log/dint.h). Caller must ensure |b| fits the dint
// range — currently safe for log/log10 where |b| <= 1074.
procedure MulDIntInt(out r: TDInt64; b: Int64; const a: TDInt64);
// Normalize X so that X.hi has its most significant bit set (if X <> 0).
// Used by the 2*pi range-reduction routines (sin/cos/tan/sincos).
procedure NormalizeDInt(var X: TDInt64);
// Convert Double to TDInt64
procedure DIntFromD(out a: TDInt64; b: Double);
// Convert TDInt64 to Double (modifies a via subnormalise; pass a copy if const needed)
function DToD(var a: TDInt64): Double;

// ------- tint64_t (192-bit) arithmetic -------
// All operations ported faithfully from core-math/src/binary64/atan2/tint.h.

procedure CpTInt(out r: TInt64; const a: TInt64); inline;
function TIntZeroP(const a: TInt64): Boolean; inline;
function CmpTIntAbs(const a, b: TInt64): Integer;
// Right shift the significand (h, m, l) by k bits. Does not touch ex/sgn.
procedure RShiftTInt(var r: TInt64; const b: TInt64; k: Integer);
// Left shift the significand (h, m, l) by k bits. Does not touch ex/sgn.
procedure LShiftTInt(var r: TInt64; const b: TInt64; k: Integer);
// r := a + b   (error bounded by 2 ulps in 192-bit)
procedure AddTInt(out r: TInt64; const a, b: TInt64);
// r := a * b   (error bounded by 10 ulps in 192-bit; alias-safe)
procedure MulTInt(out r: TInt64; const a, b: TInt64);
// Convert Double to TInt64 (exact for finite, non-NaN inputs; defined for 0)
procedure TIntFromD(out a: TInt64; b: Double);
// Convert TInt64 to Double with directed rounding driven by err (in ulps of l).
// y, x are pass-through inputs used only for the worst-case panic message.
function TIntToD(const a: TInt64; err: UInt64; y, x: Double): Double;
// r := 1 / A   (relative error < 2^-103.9; A must be non-zero)
procedure InvTInt(out r: TInt64; const A: TInt64);
// r := b / a   (relative error < 2^-185.53)
procedure DivTInt(out r: TInt64; const b, a: TInt64);
// Convenience: r := bd / ad, both Doubles
procedure DivTIntD(out r: TInt64; bd, ad: Double);

const
  cNaNSingle: Single = 0.0/0.0;       // x86 indefinite: 0xFFC00000 (negative quiet NaN)
  cNaNDouble: Double = 0.0/0.0;
  // Positive quiet NaNs matching __builtin_nanf("") and __builtin_nanf("1") in C.
  cNaNSinglePos:  Tb32u32 = (u: $7FC00000); // positive quiet NaN, payload 0
  cNaNSinglePos1: Tb32u32 = (u: $7FC00001); // positive quiet NaN, payload 1
  // Positive quiet NaNs (Double) matching __builtin_nan("") and __builtin_nan("1") in C.
  cNaNDoublePos:  Tb64u64 = (u: $7FF8000000000000); // positive quiet NaN, payload 0
  cNaNDoublePos1: Tb64u64 = (u: $7FF8000000000001); // positive quiet NaN, payload 1
  // 2^(-127): subnormal Single used to trigger IEEE 754 underflow via multiplication
  cUnderflowSingle: Single = 5.877471754111438e-39;

  // Sentinel dint constants (see sin.c ZERO / ONE)
  DINT_ZERO: TDInt64 = (hi: 0; lo: 0; ex: -1076; sgn: 0);
  // hi = $8000000000000000 = 2^63 (MSB set, normalised 1.0 in dint format)
  // ex=1: value = hi/2^64 * 2^ex = (2^63/2^64) * 2^1 = 0.5 * 2 = 1.0
  DINT_ONE:  TDInt64 = (hi: $8000000000000000; lo: 0; ex: 1; sgn: 0);

  // Sentinel TInt64 constants (see atan2/tint.h ZERO / ONE / PI / PI2)
  TINT_ZERO: TInt64 = (m: 0; h: 0; l: 0; ex: -1076; sgn: 0);
  TINT_ONE:  TInt64 = (m: 0; h: $8000000000000000; l: 0; ex: 1; sgn: 0);
  // pi to error < 2^-196.96
  TINT_PI:   TInt64 = (m: $C4C6628B80DC1CD1; h: $C90FDAA22168C234;
                       l: $29024E088A67CC74; ex: 2; sgn: 0);
  // pi/2 to error < 2^-197.96
  TINT_PI2:  TInt64 = (m: $C4C6628B80DC1CD1; h: $C90FDAA22168C234;
                       l: $29024E088A67CC74; ex: 1; sgn: 0);

  // Helpers used inside InvTInt
  cTI_1pm1022: Tb64u64 = (u: $0010000000000000); // 0x1p-1022 (smallest normal)
  cTI_1p53:    Tb64u64 = (u: $4340000000000000); // 0x1p+53

  // Hex-float constants used inside DToD / subnormalise_dint
  cDTD_1pm53:    Tb64u64 = (u: $3CA0000000000000); // 0x1p-53
  cDTD_1pm54:    Tb64u64 = (u: $3C90000000000000); // 0x1p-54
  cDTD_1p1023:   Tb64u64 = (u: $7FE0000000000000); // 0x1p+1023
  cDTD_MaxNorm:  Tb64u64 = (u: $7FEFFFFFFFFFFFFF); // 0x1.fffffffffffffp+1023 = DBL_MAX
  cDTD_MinSub:   Tb64u64 = (u: $0000000000000001); // 0x1p-1074 = smallest subnormal

implementation

// ---------------------------------------------------------------------------
// Internal 64-bit clz helper (no external dependencies, used by dint ops).
// Identical implementation to pcr_clzll in pascoremathhelperfuncs.pas.
// ---------------------------------------------------------------------------
function clzll64(x: UInt64): Integer; inline;
{$IFDEF AVX2}
begin
  Result := 63 - BsrQWord(x);
end;
{$ELSE}
var n: Integer;
begin
  if x = 0 then begin Result := 64; Exit; end;
  n := 0;
  if (x and UInt64($FFFFFFFF00000000)) = 0 then begin n := n + 32; x := x shl 32; end;
  if (x and UInt64($FFFF000000000000)) = 0 then begin n := n + 16; x := x shl 16; end;
  if (x and UInt64($FF00000000000000)) = 0 then begin n := n +  8; x := x shl  8; end;
  if (x and UInt64($F000000000000000)) = 0 then begin n := n +  4; x := x shl  4; end;
  if (x and UInt64($C000000000000000)) = 0 then begin n := n +  2; x := x shl  2; end;
  if (x and UInt64($8000000000000000)) = 0 then n := n + 1;
  Result := n;
end;
{$ENDIF}

// ---------------------------------------------------------------------------
// Internal 128-bit arithmetic helpers (private, used by dint ops)
// ---------------------------------------------------------------------------

// r := a + b  (128-bit addition)
procedure AddU128(out r: TUInt128; const a, b: TUInt128); inline;
var alo: UInt64;
begin
  alo := a.lo;
  r.lo := alo + b.lo;
  r.hi := a.hi + b.hi + UInt64(r.lo < alo);
end;

// r := a - b  (128-bit subtraction)
procedure SubU128(out r: TUInt128; const a, b: TUInt128); inline;
var alo: UInt64;
begin
  alo := a.lo;
  r.lo := alo - b.lo;
  r.hi := a.hi - b.hi - UInt64(alo < b.lo);
end;

// a <<= sh  (0 <= sh; in-place 128-bit left shift)
procedure ShlU128(var a: TUInt128; sh: Integer); inline;
begin
  if sh <= 0 then Exit;
  if sh < 64 then begin
    a.hi := (a.hi shl sh) or (a.lo shr (64 - sh));
    a.lo := a.lo shl sh;
  end else if sh < 128 then begin
    a.hi := a.lo shl (sh - 64);
    a.lo := 0;
  end else begin
    a.hi := 0; a.lo := 0;
  end;
end;

// a >>= sh  (0 <= sh; in-place 128-bit logical right shift)
procedure ShrU128(var a: TUInt128; sh: Integer); inline;
begin
  if sh <= 0 then Exit;
  if sh < 64 then begin
    a.lo := (a.lo shr sh) or (a.hi shl (64 - sh));
    a.hi := a.hi shr sh;
  end else if sh < 128 then begin
    a.lo := a.hi shr (sh - 64);
    a.hi := 0;
  end else begin
    a.lo := 0; a.hi := 0;
  end;
end;

// ---------------------------------------------------------------------------
// TUInt128 operator and Mulu64u64 (existing)
// ---------------------------------------------------------------------------

operator +(const a: TUInt128; b: UInt64): TUInt128; inline;
begin
  Result.lo := a.lo + b;
  Result.hi := a.hi + UInt64(Result.lo < b);
end;

function Mulu64u64(a, b: UInt64): TUInt128; inline;
{$IFDEF AVX2}
var
  rlo, rhi: UInt64;
begin
  asm
    mov  rax, a
    mul  b           // rdx:rax = a * b
    mov  rlo, rax
    mov  rhi, rdx
  end;
  Result.lo := rlo;
  Result.hi := rhi;
end;
{$ELSE}
// Portable fallback: four 32-bit partial products
// done by nanobit in the Lazarus forum: https://forum.lazarus.freepascal.org/index.php/topic,73881.0.html
var
  MulLo, Temp1, Temp2: UInt64;
begin
  MulLo := uint64(uint32(a)) * uint64(uint32(b));
  Temp1 := (a shr 32) * uint64(uint32(b)) + (MulLo shr 32);
  Temp2 := uint64(uint32(a)) * (b shr 32) + uint64(uint32(Temp1));
  Result.lo := (Temp2 shl 32) or (MulLo and $FFFFFFFF);
  Result.hi := (a shr 32) * (b shr 32) + (Temp1 shr 32) + (Temp2 shr 32);
end;
{$ENDIF}

// ---------------------------------------------------------------------------
// dint64_t arithmetic
// ---------------------------------------------------------------------------

procedure CpDInt(out r: TDInt64; const a: TDInt64); inline;
begin
  r := a;
end;

function DIntZeroP(const a: TDInt64): Boolean; inline;
begin
  Result := a.hi = 0;
end;

function CmpDIntAbs(const a, b: TDInt64): Integer;
begin
  if a.hi = 0 then begin
    if b.hi = 0 then Result := 0 else Result := -1;
    Exit;
  end;
  if b.hi = 0 then begin Result := 1; Exit; end;
  if a.ex > b.ex then begin Result := 1; Exit; end;
  if a.ex < b.ex then begin Result := -1; Exit; end;
  // same exponent: compare 128-bit significands as unsigned
  if a.hi > b.hi then Result := 1
  else if a.hi < b.hi then Result := -1
  else if a.lo > b.lo then Result := 1
  else if a.lo < b.lo then Result := -1
  else Result := 0;
end;

// Subnormal rounding used inside DToD — modifies a in place.
// Ported faithfully from subnormalize_dint in sin.c.
procedure SubnormalizeDInt(var a: TDInt64);
var
  ex: UInt64;
  hi, md, lo: UInt64;
  rmode: TFPURoundingMode;
begin
  if a.ex > -1023 then Exit;

  ex := UInt64(-(1011 + a.ex));
  hi := a.hi shr ex;
  md := (a.hi shr (ex - 1)) and 1;
  // logical OR: lo = 1 if any low/residual bits are non-zero
  lo := UInt64(((a.hi and (UInt64($FFFFFFFFFFFFFFFF) shr ex)) <> 0) or (a.lo <> 0));

  rmode := GetRoundMode;
  case rmode of
    rmNearest:
      if lo <> 0 then hi := hi + md
      else hi := hi + (hi and md);
    rmDown:
      if (a.sgn <> 0) and ((md or lo) <> 0) then Inc(hi);
    rmUp:
      if (a.sgn = 0) and ((md or lo) <> 0) then Inc(hi);
    // rmTruncate: truncate towards zero — no correction needed
  end;

  a.hi := hi shl ex;
  a.lo := 0;

  if a.hi = 0 then begin
    Inc(a.ex);
    a.hi := UInt64(1) shl 63;
  end;
end;

// Ported from add_dint in sin.c.
// NOTE: Pascal identifiers are case-insensitive; local 128-bit vars are
// named vA/vB/vC/vD/vE to avoid collision with the parameters a/b.
procedure AddDInt(out r: TDInt64; const a, b: TDInt64);
var
  pa, pb, ptmp: TDInt64;
  vA, vB, vBorig, vC, vD, vE: TUInt128;
  k, ex: UInt64;
  sgn: Byte;
  ch: UInt64;
  cmp: Integer;
  sh: Integer;
begin
  pa := a; pb := b;  // local copies handle aliasing (r may alias a or b)

  // if a is zero (both hi and lo are 0), return b
  if (pa.hi or pa.lo) = 0 then begin
    r := pb;
    Exit;
  end;

  cmp := CmpDIntAbs(pa, pb);

  case cmp of
    0:
      begin
        if (pa.sgn xor pb.sgn) <> 0 then
          r := DINT_ZERO
        else begin
          r := pa;
          Inc(r.ex);
        end;
        Exit;
      end;
    -1:
      begin
        ptmp := pa; pa := pb; pb := ptmp;  // swap so |pa| >= |pb|
      end;
  end;

  // From here |pa| >= |pb|
  vA.hi := pa.hi; vA.lo := pa.lo;
  vB.hi := pb.hi; vB.lo := pb.lo;
  vBorig := vB;  // save original pb.r for Sterbenz case
  k := UInt64(pa.ex - pb.ex);

  if k > 0 then begin
    if k < 128 then
      ShrU128(vB, Integer(k))
    else begin
      vB.hi := 0; vB.lo := 0;
    end;
  end;

  sgn := pa.sgn;
  r.ex := pa.ex;

  if (pa.sgn xor pb.sgn) <> 0 then begin
    // Different signs: vC = vA - vB
    SubU128(vC, vA, vB);

    ch := vC.hi;
    if ch <> 0 then
      ex := UInt64(clzll64(ch))
    else
      ex := 64 + UInt64(clzll64(vC.lo));

    if ex > 0 then begin
      sh := Integer(ex);
      if k = 1 then begin
        // Sterbenz case: vC = (vA << ex) - (vBorig << (ex - 1))
        vD := vA; ShlU128(vD, sh);
        vE := vBorig; ShlU128(vE, sh - 1);
        SubU128(vC, vD, vE);
      end else begin
        vD := vA; ShlU128(vD, sh);
        vE := vB; ShlU128(vE, sh);
        SubU128(vC, vD, vE);
      end;
      Dec(r.ex, Int64(ex));
      ex := UInt64(clzll64(vC.hi));  // now 0 or 1
    end;

    // Final normalization
    ShlU128(vC, Integer(ex));
    Dec(r.ex, Int64(ex));
  end else begin
    // Same signs: vC = vA + vB
    AddU128(vC, vA, vB);

    // Detect 128-bit overflow: vC < vA
    if (vC.hi < vA.hi) or ((vC.hi = vA.hi) and (vC.lo < vA.lo)) then begin
      // vC = (1 << 127) | (vC >> 1)
      vC.lo := (vC.lo shr 1) or (vC.hi shl 63);
      vC.hi := UInt64($8000000000000000) or (vC.hi shr 1);
      Inc(r.ex);
    end;
  end;

  r.sgn := sgn;
  r.hi := vC.hi;
  r.lo := vC.lo;
end;

// Ported from mul_dint in sin.c.
// NOTE: overlap between r and a is allowed (inputs saved to locals first).
procedure MulDInt(out r: TDInt64; const a, b: TDInt64);
var
  m1, m2, rr: TUInt128;
  ah, bh, al, bl: UInt64;
  ex: UInt64;
  rex_a, rex_b: Int64;
  rsgn: Byte;
begin
  // Save inputs before any write to r (r may alias a)
  ah := a.hi; al := a.lo;
  bh := b.hi; bl := b.lo;
  rex_a := a.ex; rex_b := b.ex;
  rsgn := a.sgn xor b.sgn;

  // hi * hi
  rr := Mulu64u64(ah, bh);

  // middle terms: add high 64 bits of (hi*lo) and (lo*hi)
  m1 := Mulu64u64(ah, bl);
  m2 := Mulu64u64(al, bh);

  // rr += (m1 >> 64) + (m2 >> 64)  — no overflow (see C comment)
  rr := rr + m1.hi;
  rr := rr + m2.hi;

  // Normalize: ensure MSB of rr.hi is set
  ex := rr.hi shr 63;
  if ex = 0 then begin
    rr.hi := (rr.hi shl 1) or (rr.lo shr 63);
    rr.lo := rr.lo shl 1;
  end;

  r.hi  := rr.hi;
  r.lo  := rr.lo;
  r.ex  := rex_a + rex_b + Int64(ex) - 1;
  r.sgn := rsgn;
end;

// Ported from mul_dint_21 in core-math/src/binary64/cos/cos.c:
// "Multiply two dint64_t numbers, assuming the low part of b is zero,
//  with error bounded by 2 ulps."
procedure MulDInt21(out r: TDInt64; const a, b: TDInt64);
var
  hi, lo: TUInt128;
  ah, al, bh: UInt64;
  rex_a, rex_b: Int64;
  rsgn: Byte;
  ex: UInt64;
begin
  ah := a.hi; al := a.lo;
  bh := b.hi;
  rex_a := a.ex; rex_b := b.ex;
  rsgn := a.sgn xor b.sgn;

  hi := Mulu64u64(ah, bh);
  lo := Mulu64u64(al, bh);

  // r.r = hi + (lo >> 64)
  hi := hi + lo.hi;

  ex := hi.hi shr 63;
  if ex = 0 then begin
    hi.hi := (hi.hi shl 1) or (hi.lo shr 63);
    hi.lo := hi.lo shl 1;
  end;

  r.hi  := hi.hi;
  r.lo  := hi.lo;
  r.ex  := rex_a + rex_b + Int64(ex) - 1;
  r.sgn := rsgn;
end;

// Ported from mul_dint_2 in core-math/src/binary64/log/dint.h.
// Multiplies a dint by an Int64; result has same convention as input.
procedure MulDIntInt(out r: TDInt64; b: Int64; const a: TDInt64);
var
  c: UInt64;
  t, l, sum: TUInt128;
  m: Integer;
  carry: Boolean;
begin
  if b = 0 then begin
    CpDInt(r, DINT_ZERO);
    Exit;
  end;
  if b < 0 then begin
    c := UInt64(-b);
    r.sgn := Byte(a.sgn xor 1);
  end else begin
    c := UInt64(b);
    r.sgn := a.sgn;
  end;

  // t = a.hi * c (128-bit)
  t := Mulu64u64(a.hi, c);
  if t.hi <> 0 then m := clzll64(t.hi)
  else m := 64;

  // t <<= m
  ShlU128(t, m);

  // l = a.lo * c
  l := Mulu64u64(a.lo, c);
  // l = (l << (m-1)) >> 63  -- round bit alignment.
  // m >= 1 in all reachable cases (a.hi != 0).
  if m >= 1 then ShlU128(l, m - 1);
  ShrU128(l, 63);

  // t = l + t (track carry)
  sum.lo := l.lo + t.lo;
  carry  := sum.lo < l.lo;
  sum.hi := l.hi + t.hi + UInt64(carry);
  carry  := (sum.hi < l.hi) or ((sum.hi = l.hi) and carry);
  t := sum;

  if carry then begin
    // round-half-to-even on the bottom bit, then shift right 1 and set MSB
    if (t.lo and 1) <> 0 then begin
      // t += 1 (with carry-out into bit 128 ignored — already overflowed)
      Inc(t.lo);
      if t.lo = 0 then Inc(t.hi);
    end;
    // t >>= 1
    t.lo := (t.lo shr 1) or (t.hi shl 63);
    t.hi := (t.hi shr 1) or (UInt64(1) shl 63);
    Dec(m);
  end;

  r.hi := t.hi;
  r.lo := t.lo;
  r.ex := a.ex + 64 - Int64(m);
end;

// Ported from normalize() in core-math/src/binary64/cos/cos.c:
// shift left so X.hi has its MSB set (if X <> 0); adjust ex accordingly.
procedure NormalizeDInt(var X: TDInt64);
var
  cnt: Integer;
begin
  if X.hi <> 0 then begin
    cnt := clzll64(X.hi);
    if cnt <> 0 then begin
      X.hi := (X.hi shl cnt) or (X.lo shr (64 - cnt));
      X.lo := X.lo shl cnt;
    end;
    X.ex := X.ex - Int64(cnt);
  end else if X.lo <> 0 then begin
    cnt := clzll64(X.lo);
    X.hi := X.lo shl cnt;
    X.lo := 0;
    X.ex := X.ex - Int64(64 + cnt);
  end;
end;

// Ported from dint_fromd / fast_extract in sin.c.
procedure DIntFromD(out a: TDInt64; b: Double);
var
  xu: Tb64u64;
  e: Int64;
  m: UInt64;
  t: Integer;
begin
  xu.f := b;
  e := Int64((xu.u shr 52) and $7FF);
  m := xu.u and UInt64($000FFFFFFFFFFFFF);
  if e <> 0 then m := m or (UInt64(1) shl 52);
  e := e - $3FE;  // biased_exp - 1022

  t := clzll64(m);
  a.sgn := Byte(b < 0.0);
  a.hi  := m shl t;
  if t > 11 then a.ex := e - Int64(t - 12)
  else a.ex := e;
  a.lo  := 0;
end;

// Ported from dint_tod in sin.c.
// Calls SubnormalizeDInt which may modify a.
function DToD(var a: TDInt64): Double;
var
  ru, eu: Tb64u64;
  rd: Double;
  ex_val: Int64;
begin
  SubnormalizeDInt(a);

  ru.u := (a.hi shr 11) or (UInt64($3FF) shl 52);

  rd := 0.0;
  if ((a.hi shr 10) and 1) <> 0 then
    rd := rd + cDTD_1pm53.f;   // 0x1p-53

  if ((a.hi and $3FF) <> 0) or (a.lo <> 0) then
    rd := rd + cDTD_1pm54.f;   // 0x1p-54

  if a.sgn <> 0 then rd := -rd;

  ru.u := ru.u or (UInt64(a.sgn) shl 63);
  ru.f := ru.f + rd;

  ex_val := a.ex;

  if ex_val > -1022 then begin
    // Normal double result
    if ex_val > 1024 then begin
      if ex_val = 1025 then begin
        ru.f := ru.f * 2.0;
        eu.f := cDTD_1p1023.f;    // 0x1p+1023
      end else begin
        ru.f := cDTD_MaxNorm.f;   // DBL_MAX
        eu.f := cDTD_MaxNorm.f;
      end;
    end else
      eu.u := UInt64((ex_val + 1022) and $7FF) shl 52;
  end else begin
    // Subnormal range
    if ex_val < -1073 then begin
      if ex_val = -1074 then begin
        ru.f := ru.f * 0.5;
        eu.f := cDTD_MinSub.f;    // 0x1p-1074
      end else begin
        ru.f := cDTD_MinSub.f;
        eu.f := cDTD_MinSub.f;
      end;
    end else
      eu.u := UInt64(1) shl UInt64(ex_val + 1073);
  end;

  Result := ru.f * eu.f;
end;

// ---------------------------------------------------------------------------
// TInt64 (192-bit) arithmetic — ported from core-math/src/binary64/atan2/tint.h
// ---------------------------------------------------------------------------

procedure CpTInt(out r: TInt64; const a: TInt64); inline;
begin
  r := a;
end;

function TIntZeroP(const a: TInt64): Boolean; inline;
begin
  Result := a.h = 0;
end;

function CmpTIntAbs(const a, b: TInt64): Integer;
begin
  if a.h = 0 then begin
    if b.h = 0 then Result := 0 else Result := -1;
    Exit;
  end;
  if b.h = 0 then begin Result := 1; Exit; end;
  if a.ex > b.ex then begin Result := 1; Exit; end;
  if a.ex < b.ex then begin Result := -1; Exit; end;
  // same exponent: compare 192-bit significands as unsigned (h:m:l)
  if a.h > b.h then Result := 1
  else if a.h < b.h then Result := -1
  else if a.m > b.m then Result := 1
  else if a.m < b.m then Result := -1
  else if a.l > b.l then Result := 1
  else if a.l < b.l then Result := -1
  else Result := 0;
end;

// Right shift only the (h, m, l) significand. Caller manages ex/sgn.
procedure RShiftTInt(var r: TInt64; const b: TInt64; k: Integer);
var bh, bm, bl: UInt64;
begin
  bh := b.h; bm := b.m; bl := b.l;
  if k = 0 then begin r.h := bh; r.m := bm; r.l := bl; end
  else if k < 64 then begin
    r.h := bh shr k;
    r.m := (bm shr k) or (bh shl (64 - k));
    r.l := (bl shr k) or (bm shl (64 - k));
  end
  else if k = 64 then begin
    r.h := 0;
    r.m := bh;
    r.l := bm;
  end
  else if k < 128 then begin
    r.h := 0;
    r.m := bh shr (k - 64);
    r.l := (bm shr (k - 64)) or (bh shl (128 - k));
  end
  else if k < 192 then begin
    r.h := 0;
    r.m := 0;
    r.l := bh shr (k - 128);
  end
  else begin
    r.h := 0; r.m := 0; r.l := 0;
  end;
end;

// Left shift only the (h, m, l) significand. Caller manages ex/sgn.
procedure LShiftTInt(var r: TInt64; const b: TInt64; k: Integer);
var bh, bm, bl: UInt64;
begin
  bh := b.h; bm := b.m; bl := b.l;
  if k = 0 then begin r.h := bh; r.m := bm; r.l := bl; end
  else if k < 64 then begin
    r.h := (bh shl k) or (bm shr (64 - k));
    r.m := (bm shl k) or (bl shr (64 - k));
    r.l := bl shl k;
  end
  else if k = 64 then begin
    r.h := bm;
    r.m := bl;
    r.l := 0;
  end
  else if k < 128 then begin
    r.h := (bm shl (k - 64)) or (bl shr (128 - k));
    r.m := bl shl (k - 64);
    r.l := 0;
  end
  else if k < 192 then begin
    r.h := bl shl (k - 128);
    r.m := 0;
    r.l := 0;
  end
  else begin
    r.h := 0; r.m := 0; r.l := 0;
  end;
end;

// Internal: 192-bit clz of (h:m:l) treated as unsigned.
// Returns 0..192. Undefined-but-defined: returns 192 for the all-zero input.
function clz192(h, m, l: UInt64): Integer; inline;
begin
  if h <> 0 then Result := clzll64(h)
  else if m <> 0 then Result := 64 + clzll64(m)
  else Result := 128 + clzll64(l);
end;

// Ported from add_tint in tint.h.
procedure AddTInt(out r: TInt64; const a, b: TInt64);
var
  pa, pb, ptmp, t: TInt64;
  sh: UInt64;
  ex, ex1: Integer;
  th, tm, tl: UInt64;
  rl, cl, ch, borrow: UInt64;
  pa_hu, t_hu, r_hu, cl_hu: TUInt128;
  cmp: Integer;
begin
  pa := a; pb := b;  // local copies handle aliasing

  cmp := CmpTIntAbs(pa, pb);
  case cmp of
    0:
      begin
        if (pa.sgn xor pb.sgn) <> 0 then begin
          CpTInt(r, TINT_ZERO);
          Exit;
        end;
        CpTInt(r, pa);
        Inc(r.ex);
        Exit;
      end;
    -1:
      begin
        ptmp := pa; pa := pb; pb := ptmp;  // swap so |pa| >= |pb|
      end;
  end;

  // From here |pa| > |pb|, so pa.ex >= pb.ex
  sh := UInt64(pa.ex - pb.ex);
  // rshift writes only h/m/l; preserve t.ex and t.sgn (unused after).
  t.ex := 0; t.sgn := 0;
  if sh < 192 then RShiftTInt(t, pb, Integer(sh))
  else begin t.h := 0; t.m := 0; t.l := 0; end;

  if (pa.sgn xor pb.sgn) <> 0 then begin
    // Subtract: t := pa - t (192-bit subtraction, no borrow out since |pa| > |pb|)
    tl := pa.l - t.l;
    borrow := UInt64(t.l > pa.l);
    tm := pa.m - t.m - borrow;
    // borrow out of the m subtraction: one if (pa.m < t.m) OR (pa.m == t.m and borrow)
    if pa.m < t.m then borrow := 1
    else if (pa.m = t.m) and (borrow = 1) then borrow := 1
    else borrow := 0;
    th := pa.h - t.h - borrow;
    t.h := th; t.m := tm; t.l := tl;

    ex := clz192(t.h, t.m, t.l);
    if (ex <= 1) or (sh = 0) then begin
      LShiftTInt(r, t, ex);
      r.ex := pa.ex - ex;
    end
    else begin
      // ex >= 2 and sh >= 1: redo with no neglected low bits of pb
      // Shift t := pb << (ex - sh), r := pa << ex, then t := r - t.
      LShiftTInt(t, pb, ex - Integer(sh));
      LShiftTInt(r, pa, ex);
      tl := r.l - t.l;
      borrow := UInt64(t.l > r.l);
      tm := r.m - t.m - borrow;
      if r.m < t.m then borrow := 1
      else if (r.m = t.m) and (borrow = 1) then borrow := 1
      else borrow := 0;
      th := r.h - t.h - borrow;
      t.h := th; t.m := tm; t.l := tl;
      ex1 := clz192(t.h, t.m, t.l);
      LShiftTInt(r, t, ex1);
      r.ex := pa.ex - (ex + ex1);
    end;
  end
  else begin
    // Same signs: 192-bit addition.  In C, _h is the high u128 (m,h pair); we
    // use TUInt128 here for parity with the C's u128 add + overflow detection.
    pa_hu.lo := pa.m; pa_hu.hi := pa.h;
    t_hu.lo  := t.m;  t_hu.hi  := t.h;
    rl := pa.l + t.l;
    cl := UInt64(rl < pa.l);
    AddU128(r_hu, pa_hu, t_hu);
    // ch := (r_hu < pa_hu) as u128
    if (r_hu.hi < pa_hu.hi) or
       ((r_hu.hi = pa_hu.hi) and (r_hu.lo < pa_hu.lo)) then ch := 1
    else ch := 0;
    // r_hu += cl
    cl_hu.lo := cl; cl_hu.hi := 0;
    AddU128(r_hu, r_hu, cl_hu);
    // overflow of r_hu when adding cl: result < cl as u128 ⇒ result.hi = 0 ∧ result.lo < cl
    if (r_hu.hi = 0) and (r_hu.lo < cl) then Inc(ch);
    if ch <> 0 then begin
      // 193-bit overflow: shift result right by 1, insert ch=1 at MSB of h.
      r.l := (r_hu.lo shl 63) or (rl shr 1);
      r.m := (r_hu.hi shl 63) or (r_hu.lo shr 1);
      r.h := (ch shl 63) or (r_hu.hi shr 1);
      r.ex := pa.ex + 1;
    end
    else begin
      r.l := rl;
      r.m := r_hu.lo;
      r.h := r_hu.hi;
      r.ex := pa.ex;
    end;
  end;
  r.sgn := pa.sgn;
end;

// Ported from mul_tint in tint.h.
procedure MulTInt(out r: TInt64; const a, b: TInt64);
var
  ah, am, al, bh, bm, bl: UInt64;
  rh, rm1, rm2, rl1, rl2, rl3: TUInt128;
  rsum: TUInt128;
  rh_v, rm_v, rl_v: UInt64;
  hh, lo, cm: UInt64;
  rex_a, rex_b: Int64;
  rsgn: UInt64;
begin
  rex_a := a.ex; rex_b := b.ex;
  rsgn := a.sgn xor b.sgn;
  ah := a.h; am := a.m; al := a.l;
  bh := b.h; bm := b.m; bl := b.l;

  rh  := Mulu64u64(ah, bh);
  rm1 := Mulu64u64(ah, bm);
  rm2 := Mulu64u64(am, bh);
  rl1 := Mulu64u64(ah, bl);
  rl2 := Mulu64u64(am, bm);
  rl3 := Mulu64u64(al, bh);

  rh_v := rh.hi;
  rm_v := rh.lo;
  rl_v := rm1.lo;

  // Accumulate rm1's high part into rm_v (carry to rh_v)
  hh := rm1.hi;
  rm_v := rm_v + hh;
  if rm_v < hh then Inc(rh_v);

  // Accumulate rm2 (lo into rl_v with carry-out cm; hi into rm_v)
  lo := rm2.lo;
  hh := rm2.hi;
  rl_v := rl_v + lo;
  cm := UInt64(rl_v < lo);
  rm_v := rm_v + hh;
  if rm_v < hh then Inc(rh_v);

  // Accumulate (rl1.hi + rl2.hi + rl3.hi) into (rl_v, cm)
  rsum.lo := rl1.hi; rsum.hi := 0;
  rsum := rsum + rl2.hi;
  rsum := rsum + rl3.hi;
  lo := rsum.lo;
  cm := cm + rsum.hi;
  rl_v := rl_v + lo;
  if rl_v < lo then Inc(cm);

  // Accumulate cm into rm_v (carry to rh_v)
  rm_v := rm_v + cm;
  if rm_v < cm then Inc(rh_v);

  r.ex := rex_a + rex_b;
  r.sgn := rsgn;
  if (rh_v shr 63) = 0 then begin
    // Normalize: shift left 1
    rh_v := (rh_v shl 1) or (rm_v shr 63);
    rm_v := (rm_v shl 1) or (rl_v shr 63);
    rl_v := rl_v shl 1;
    Dec(r.ex);
  end;
  r.h := rh_v; r.m := rm_v; r.l := rl_v;
end;

// Ported from tint_fromd in tint.h. Defined for 0 (yields h=m=l=0).
procedure TIntFromD(out a: TInt64; b: Double);
var
  u: Tb64u64;
  ax: UInt64;
  e: Int64;
  cnt: Integer;
begin
  u.f := b;
  a.sgn := u.u shr 63;
  ax := u.u and UInt64($7FFFFFFFFFFFFFFF);
  e := Int64(ax shr 52);
  if e <> 0 then begin
    a.ex := e - $3FE;
    a.h := (UInt64(1) shl 63) or (ax shl 11);
  end
  else begin
    cnt := clzll64(ax);
    a.ex := -$3F2 - cnt;
    if cnt < 64 then a.h := ax shl cnt else a.h := 0;
  end;
  a.m := 0; a.l := 0;
end;

// Ported from tint_tod in tint.h. Calls Math.Ldexp for the final exponent fold.
function TIntToD(const a: TInt64; err: UInt64; y, x: Double): Double;
const
  S: array[0..1] of Double = (1.0, -1.0);
var
  hh, mm, ll, low, notmm, notll: UInt64;
  ex_val: Int64;
  sh: Integer;
  hf, lf, sf: Double;
  worst: Boolean;
  mid: Boolean;
begin
  // Defined extension over the C: zero significand returns +/-0.0 cleanly.
  if (a.h = 0) and (a.m = 0) and (a.l = 0) then begin
    if a.sgn <> 0 then Result := -0.0 else Result := 0.0;
    Exit;
  end;
  if a.ex >= 1025 then begin
    // overflow
    if a.sgn <> 0 then Result := -1.7976931348623157e+308 - 1.7976931348623157e+308
    else Result := 1.7976931348623157e+308 + 1.7976931348623157e+308;
    Exit;
  end;
  if a.ex <= -1074 then begin
    if a.ex < -1074 then begin
      if a.sgn <> 0 then Result := -5e-324 * 0.5 else Result := 5e-324 * 0.5;
      Exit;
    end;
    mid := (a.h = (UInt64(1) shl 63)) and (a.m = 0) and (a.l = 0);
    if a.sgn <> 0 then begin
      if mid then Result := -5e-324 * 0.5 else Result := -5e-324 * 0.75;
    end
    else begin
      if mid then Result := 5e-324 * 0.5 else Result := 5e-324 * 0.75;
    end;
    Exit;
  end;

  hh := a.h; mm := a.m; ll := a.l;
  ex_val := a.ex;
  low := hh and $7FF;
  notmm := not mm;
  notll := not ll;

  // Worst-case detection — we cannot determine correct rounding.
  if (mm = 0) or (notmm = 0) then begin
    worst :=
      ((mm = 0) and ((low = 0) or (low = $400)) and (ll < err)) or
      ((notmm = 0) and ((low = $3FF) or (low = $7FF)) and (notll < err));
    if worst then begin
      WriteLn('Unexpected worst-case found, please report to core-math@inria.fr:');
      WriteLn('Worst-case of atan2 found: y,x=', y, ',', x);
      Halt(1);
    end;
  end;

  if ex_val <= -1022 then begin
    sh := -1021 - ex_val;  // 1 <= sh <= 52
    ll := (mm shl (64 - sh)) or (ll shr sh) or UInt64(ll > 0);
    mm := (hh shl (64 - sh)) or (mm shr sh);
    hh := hh shr sh;
    low := hh and $7FF;
    ex_val := ex_val + sh;
  end;

  hf := Double(hh shr 11);  // 53-bit significand value
  if err = 0 then lf := 0.0
  else if low < $400 then lf := 0.25
  else if low > $400 then lf := 0.75
  else begin
    if (mm = 0) and (ll = 0) then lf := 0.5
    else lf := 0.75;
  end;

  sf := S[a.sgn];
  // h = fma(l, s, s*h) ; h *= 2^-52 ; result = h * 2^(ex_val-1)
  hf := lf * sf + sf * hf;
  hf := hf * 2.220446049250313e-16;  // 0x1p-52
  Result := hf * Math.Ldexp(1.0, Integer(ex_val - 1));
end;

// Ported from inv_tint in tint.h.
procedure InvTInt(out r: TInt64; const A: TInt64);
var
  q: TInt64;
  ad: Double;
  subnormal: Boolean;
begin
  ad := TIntToD(A, 0, 0.0, 0.0);
  subnormal := Abs(ad) < cTI_1pm1022.f;
  if subnormal then ad := ad * cTI_1p53.f;
  TIntFromD(r, 1.0 / ad);
  if subnormal then Inc(r.ex, 53);
  MulTInt(q, A, r);
  q.sgn := 1 - q.sgn;
  AddTInt(q, TINT_ONE, q);
  MulTInt(q, r, q);
  AddTInt(r, r, q);
end;

// Ported from div_tint in tint.h.
procedure DivTInt(out r: TInt64; const b, a: TInt64);
var
  Y, Z: TInt64;
begin
  InvTInt(Y, a);
  MulTInt(r, Y, b);
  MulTInt(Z, a, r);
  Z.sgn := 1 - Z.sgn;
  AddTInt(Z, b, Z);
  MulTInt(Z, Y, Z);
  AddTInt(r, r, Z);
end;

// Ported from div_tint_d in tint.h.
procedure DivTIntD(out r: TInt64; bd, ad: Double);
var
  A, B: TInt64;
begin
  TIntFromD(A, ad);
  TIntFromD(B, bd);
  DivTInt(r, B, A);
end;

end.
