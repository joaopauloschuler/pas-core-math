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
// Normalize X so that X.hi has its most significant bit set (if X <> 0).
// Used by the 2*pi range-reduction routines (sin/cos/tan/sincos).
procedure NormalizeDInt(var X: TDInt64);
// Convert Double to TDInt64
procedure DIntFromD(out a: TDInt64; b: Double);
// Convert TDInt64 to Double (modifies a via subnormalise; pass a copy if const needed)
function DToD(var a: TDInt64): Double;

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

end.
