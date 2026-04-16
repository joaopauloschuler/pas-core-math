// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//                                                                                                                                                                                                      
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
{$I pascoremath.inc}
unit pascoremathhelperfuncs;

interface

uses
  Math, SysUtils, pascoremathtypes;

// Fused multiply-add — hardware FMA3 on x86-64; 80-bit approximation elsewhere.
// Pure-asm body uses System V AMD64 ABI (params in xmm0/1/2).
function pcr_fmaf(x, y, z: Single): Single; {$IFNDEF AVX2} inline; {$ENDIF}
function pcr_fma(x, y, z: Double): Double; {$IFNDEF AVX2} inline; {$ENDIF}

// Absolute value
function pcr_fabsf(x: Single): Single; inline;
function pcr_fabs(x: Double): Double; inline;

// Copy sign of y to magnitude of x
function pcr_copysignf(x, y: Single): Single; inline;
function pcr_copysign(x, y: Double): Double; inline;

// Square root
function pcr_sqrtf(x: Single): Single; inline;
function pcr_sqrt(x: Double): Double; inline;

// Round to nearest even integer — SSE4.1 ROUNDSD/ROUNDSS on x86-64; bit-manip elsewhere.
// Pure-asm body uses System V AMD64 ABI (param in xmm0, result in xmm0).
function pcr_roundevenf(x: Single): Single; {$IFNDEF AVX2} inline; {$ENDIF}
function pcr_roundeven(x: Double): Double; {$IFNDEF AVX2} inline; {$ENDIF}

// NaN-aware maximum — MAXSS/MAXSD on x86-64; branch fallback elsewhere.
// Pure-asm body uses System V AMD64 ABI (params in xmm0/xmm1).
function pcr_fmaxf(x, y: Single): Single; {$IFNDEF AVX2} inline; {$ENDIF}
function pcr_fmax(x, y: Double): Double; {$IFNDEF AVX2} inline; {$ENDIF}

// NaN-aware minimum — MINSS/MINSD on x86-64; branch fallback elsewhere.
// Pure-asm body uses System V AMD64 ABI (params in xmm0/xmm1).
function pcr_fminf(x, y: Single): Single; {$IFNDEF AVX2} inline; {$ENDIF}
function pcr_fmin(x, y: Double): Double; {$IFNDEF AVX2} inline; {$ENDIF}

// Return a NaN (tagp is ignored, matches C nan()/nanf() signature)
function pcr_nanf(const tagp: PAnsiChar): Single; inline;
function pcr_nan(const tagp: PAnsiChar): Double; inline;

// Raise floating-point exceptions
function pcr_feraiseexcept_invalid():Single; inline;
function pcr_feraiseexcept_divbyzero():Single; inline;

implementation

function pcr_fmaf(x, y, z: Single): Single;
{$IFDEF AVX2}
// Pure-asm: System V AMD64 ABI passes x→xmm0, y→xmm1, z→xmm2; result in xmm0.
// VFMADD213SS: xmm0 = xmm0 * xmm1 + xmm2  (correctly rounded IEEE 754 FMA).
assembler;
asm
  vfmadd213ss xmm0, xmm1, xmm2
end;
{$ELSE}
begin
  // 80-bit fallback: correctly rounded for singles (Extended has enough mantissa bits).
  Result := Single(Extended(x) * Extended(y) + Extended(z));
end;
{$ENDIF}

function pcr_fma(x, y, z: Double): Double;
{$IFDEF AVX2}
// Pure-asm: System V AMD64 ABI passes x→xmm0, y→xmm1, z→xmm2; result in xmm0.
// VFMADD213SD: xmm0 = xmm0 * xmm1 + xmm2  (correctly rounded IEEE 754 FMA).
assembler;
asm
  vfmadd213sd xmm0, xmm1, xmm2
end;
{$ELSE}
begin
  // 80-bit fallback (double-rounding — not true FMA; may lose 1 ULP in rare cases).
  Result := Double(Extended(x) * Extended(y) + Extended(z));
end;
{$ENDIF}

function pcr_fabsf(x: Single): Single;
begin
  Result := Abs(x);
end;

function pcr_fabs(x: Double): Double;
begin
  Result := Abs(x);
end;

function pcr_copysignf(x, y: Single): Single;
var
  vx, vy: Tb32u32;
begin
  vx.f := x;
  vy.f := y;
  vx.u := (vx.u and UInt32($7FFFFFFF)) or (vy.u and UInt32($80000000));
  Result := vx.f;
end;

function pcr_copysign(x, y: Double): Double;
var
  vx, vy: Tb64u64;
begin
  vx.f := x;
  vy.f := y;
  vx.u := (vx.u and UInt64($7FFFFFFFFFFFFFFF)) or (vy.u and UInt64($8000000000000000));
  Result := vx.f;
end;

function pcr_sqrtf(x: Single): Single;
begin
  Result := Sqrt(x);
end;

function pcr_sqrt(x: Double): Double;
begin
  Result := Sqrt(x);
end;

function pcr_roundevenf(x: Single): Single;
{$IFDEF AVX2}
// Pure-asm: System V AMD64 ABI passes x→xmm0; result in xmm0.
// ROUNDSS imm8=12 (0x0C): override MXCSR with round-to-nearest-even, suppress PE.
assembler;
asm
  roundss xmm0, xmm0, 12
end;
{$ELSE}
// Portable fallback: round to nearest even using bit manipulation.
// For |x| >= 2^23 the value is already an integer.
var
  v: Tb32u32;
  e, shift: Int32;
  mask, frac, half: UInt32;
begin
  v.f := x;
  e := Int32((v.u shr 23) and $FF) - 127;  // unbiased exponent
  if e >= 23 then
  begin
    // Already an integer (or inf/nan)
    Result := x;
    Exit;
  end;
  if e < 0 then
  begin
    // |x| < 1: round to 0 or +-1
    if e = -1 then
    begin
      // |x| in [0.5, 1): round to nearest even => 0 if exactly 0.5, else +-1
      // Check if it's exactly +-0.5
      if (v.u and $7FFFFFFF) = $3F000000 then
        Result := 0.0  // exact half => round to even (0)
      else if Abs(x) < 0.5 then
        Result := 0.0
      else
        Result := pcr_copysignf(1.0, x);
    end
    else
      Result := 0.0;
    Exit;
  end;
  // e in [0, 22]: some fractional bits present
  shift := 23 - e;                    // number of fractional bits
  mask  := (1 shl shift) - 1;         // mask for fractional bits
  frac  := v.u and mask;
  half  := 1 shl (shift - 1);         // 0.5 in fractional position

  if frac < half then
  begin
    // Round down: clear fractional bits
    v.u := v.u and (not mask);
  end
  else if frac > half then
  begin
    // Round up
    v.u := (v.u and (not mask)) + (1 shl shift);
  end
  else
  begin
    // Exactly halfway: round to even (check integer bit)
    if (v.u and (1 shl shift)) <> 0 then
      // Int32 part is odd => round up
      v.u := (v.u and (not mask)) + (1 shl shift)
    else
      // Int32 part is even => round down
      v.u := v.u and (not mask);
  end;
  Result := v.f;
end;
{$ENDIF}

function pcr_roundeven(x: Double): Double;
{$IFDEF AVX2}
// Pure-asm: System V AMD64 ABI passes x→xmm0; result in xmm0.
// ROUNDSD imm8=12 (0x0C): override MXCSR with round-to-nearest-even, suppress PE.
assembler;
asm
  roundsd xmm0, xmm0, 12
end;
{$ELSE}
// Portable fallback: round to nearest even using bit manipulation.
var
  v: Tb64u64;
  e, shift: Int32;
  mask, half: UInt64;
  frac: UInt64;
begin
  v.f := x;
  e := Int32((v.u shr 52) and $7FF) - 1023;
  if e >= 52 then begin Result := x; Exit; end;
  if e < 0 then begin
    if e = -1 then begin
      if (v.u and $7FFFFFFFFFFFFFFF) = $3FE0000000000000 then Result := 0.0
      else if Abs(x) < 0.5 then Result := 0.0
      else Result := pcr_copysign(1.0, x);
    end else Result := 0.0;
    Exit;
  end;
  shift := 52 - e;
  mask  := (UInt64(1) shl shift) - 1;
  frac  := v.u and mask;
  half  := UInt64(1) shl (shift - 1);
  if frac < half then v.u := v.u and (not mask)
  else if frac > half then v.u := (v.u and (not mask)) + (UInt64(1) shl shift)
  else begin
    if (v.u and (UInt64(1) shl shift)) <> 0 then
      v.u := (v.u and (not mask)) + (UInt64(1) shl shift)
    else v.u := v.u and (not mask);
  end;
  Result := v.f;
end;
{$ENDIF}

function pcr_fmaxf(x, y: Single): Single;
{$IFDEF AVX2}
// Pure-asm: System V AMD64 ABI passes x→xmm0, y→xmm1; result in xmm0.
// MAXSS returns the larger value; if x (first operand) is NaN, returns y.
assembler;
asm
  maxss xmm0, xmm1
end;
{$ELSE}
begin
  if IsNan(x) then Result := y
  else if IsNan(y) then Result := x
  else if x > y then Result := x
  else Result := y;
end;
{$ENDIF}

function pcr_fmax(x, y: Double): Double;
{$IFDEF AVX2}
// Pure-asm: System V AMD64 ABI passes x→xmm0, y→xmm1; result in xmm0.
// MAXSD returns the larger value; if x (first operand) is NaN, returns y.
assembler;
asm
  maxsd xmm0, xmm1
end;
{$ELSE}
begin
  if IsNan(x) then Result := y
  else if IsNan(y) then Result := x
  else if x > y then Result := x
  else Result := y;
end;
{$ENDIF}

function pcr_fminf(x, y: Single): Single;
{$IFDEF AVX2}
// Pure-asm: System V AMD64 ABI passes x→xmm0, y→xmm1; result in xmm0.
// MINSS returns the smaller value; if x (first operand) is NaN, returns y.
assembler;
asm
  minss xmm0, xmm1
end;
{$ELSE}
begin
  if IsNan(x) then Result := y
  else if IsNan(y) then Result := x
  else if x < y then Result := x
  else Result := y;
end;
{$ENDIF}

function pcr_fmin(x, y: Double): Double;
{$IFDEF AVX2}
// Pure-asm: System V AMD64 ABI passes x→xmm0, y→xmm1; result in xmm0.
// MINSD returns the smaller value; if x (first operand) is NaN, returns y.
assembler;
asm
  minsd xmm0, xmm1
end;
{$ELSE}
begin
  if IsNan(x) then Result := y
  else if IsNan(y) then Result := x
  else if x < y then Result := x
  else Result := y;
end;
{$ENDIF}

function pcr_nanf(const tagp: PAnsiChar): Single;
begin
  Result := Single(NaN);
end;

function pcr_nan(const tagp: PAnsiChar): Double;
begin
  Result := NaN;
end;

function pcr_feraiseexcept_invalid(): Single;
var
  x: Single;
begin
  // Raise FE_INVALID by computing 0/0
  x := 0.0;
  Result := x / x;
end;

function pcr_feraiseexcept_divbyzero(): Single;
var
  x: Single;
begin
  // Raise FE_DIVBYZERO by computing 1/0
  x := 0.0;
  Result := 1.0 / x;
end;

end.
