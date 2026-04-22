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

// Count leading zeros in a 64-bit value (__builtin_clzll equivalent).
// Result is undefined when x = 0.
function pcr_clzll(x: UInt64): Integer; inline;

// Current FPU rounding mode (wraps Math.GetRoundMode for use in ported functions).
function pcr_GetRoundMode: TFPURoundingMode; inline;

// ------- binary64 bit-pattern helpers (task 0.10) -------
// Double equivalents of the cf_*/isint_pf/isodd_pf helpers in pascoremath32.pas.

// Detect signaling NaN via bit-flip at $0008000000000000.
function pcr_is_signaling(x: Double): Boolean; inline;

// Returns true if y is an exact integer (general use).
function pcr_isint_d(y: Double): Boolean; inline;

// pow-variant of integer test (same logic as pcr_isint_d).
function pcr_isint_pd(y: Double): Boolean; inline;

// pow-variant of odd-integer test.
function pcr_isodd_pd(y: Double): Boolean; inline;

// pow-variant of signaling-NaN test.
function pcr_is_signaling_pd(x: Double): Boolean; inline;

// Double x double-double multiply: returns (xh+xl)*ch with error in l.
function pcr_mulddd_pd(xh, xl, ch: Double; out l: Double): Double; inline;

// MXCSR flag save/restore (AVX2: hardware register; otherwise no-op).
function pcr_get_mxcsr: DWord;
procedure pcr_set_mxcsr(flag: DWord);

// ------- double-double and polynomial helpers (task 0.9, promoted from pascoremath32) -------

// Degree-12 polynomial evaluator (used by acosf, asinf and their binary64 analogues).
function pcr_poly12(z: Double; const c: array of Double): Double; inline;

// Double-double × double-double product: returns xh*ch + mixed terms, error in l.
function pcr_muldd(xh, xl, ch, cl: Double; out l: Double): Double; inline;

// Horner evaluation of a flat-array double-double polynomial.
// c is flat: c[k*2] = high part, c[k*2+1] = low part.
function pcr_polydd(xh, xl: Double; n: Int32; const c: array of Double; out l: Double): Double; inline;

// All four primitives below write their var outputs LAST (after all value-param
// reads) so that callers may safely alias value params with var params.

// Error-free sum: s + t = a + b exactly.
procedure pcr_fasttwosum(var s, t: Double; a, b: Double); inline;

// Error-free product: hi + lo = a * b exactly (via FMA).
procedure pcr_a_mul(var hi, lo: Double; a, b: Double); inline;

// Scalar × double-double: (hi + lo) = a * (bh + bl).
procedure pcr_s_mul(var hi, lo: Double; a, bh, bl: Double); inline;

// Double-double × double-double: (hi + lo) = (ah + al) * (bh + bl).
procedure pcr_d_mul(var hi, lo: Double; ah, al, bh, bl: Double); inline;

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

// =============================================================================
// Correctly rounded FMA3 emulation
// Coded by MathMan in the forum: https://forum.lazarus.freepascal.org/index.php/topic,73881.30.html
//
// The implementation is based on https://hal.science/hal-04575249/document
// ============================================================================= 
type
  DW = record
    h: Double;
    l: Double;
  end;

function split( x: Double ):DW;
const
  K: Double = 134217729.0;
var
  splittedx: DW;
  gamma: Double;
begin
  gamma := K * x;
  splittedx.h := ( gamma+( x-gamma ) );
  splittedx.l := x - splittedx.h;
  Result := splittedx;
end;

function DekkerProd( a,b: Double ):DW;
var
  splitteda: DW;
  splittedb: DW;
  product: DW;
begin
  splitteda := split( a );
  splittedb := split( b );
  product.h := a * b;
  product.l := ((( -product.h+splitteda.h*splittedb.h )+( splitteda.h*splittedb.l ))
            + splitteda.l*splittedb.h ) + splitteda.l*splittedb.l;
  Result := product;
end;

function TwoSum( a,b: Double ):DW;
var
  z: DW;
  aprime: Double;
begin
  z.h := a + b;
  aprime := z.h - b;
  z.l := ( a-aprime ) + ( b-( z.h-aprime ) );

  Result := z;
end;

function IsNot1or3timesPowerOf2( x: Double ):Boolean;
const
  P: Double = 2251799813685249.0;
  Q: Double = 2251799813685248.0;
var
  Delta: Double;
begin
  Delta := ( P*x ) - ( Q*x );
  Result := ( Delta<>x );
end;

function pcr_fma_pascal( a,b,c: Double ):Double; inline;
var
  x, s, v: DW;
begin
  x := DekkerProd( a,b );
  s := TwoSum( x.h,c );
  v := TwoSum( x.l,s.l );
  if( ( IsNot1or3timesPowerOf2( v.h ) ) or ( v.l=0 ) ) then
  begin
    Result := s.h + v.h;
  end
  else
  begin
    if( ( UInt8( v.l<0 ) xor UInt8( v.h<0 ) )<>0 ) then
    begin
      Result := s.h + ( 0.875*v.h );
    end
    else
    begin
      Result := s.h + ( 1.125*v.h );
    end;
  end;
end;

// =============================================================================
// END OF Correctly rounded FMA3 emulation
// =============================================================================

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
  // Result := Double(Extended(x) * Extended(y) + Extended(z));
  Result := pcr_fma_pascal(x, y, z);
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

function pcr_clzll(x: UInt64): Integer; inline;
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

function pcr_GetRoundMode: TFPURoundingMode; inline;
begin
  Result := GetRoundMode;
end;

// ---------------------------------------------------------------------------
// binary64 bit-pattern helpers (task 0.10)
// ---------------------------------------------------------------------------

function pcr_is_signaling(x: Double): Boolean; inline;
var u: Tb64u64;
begin
  u.f := x;
  u.u := u.u xor UInt64($0008000000000000);
  Result := (u.u and UInt64($7FFFFFFFFFFFFFFF)) > UInt64($7FF8000000000000);
end;

function pcr_isint_d(y: Double): Boolean; inline;
var wy: Tb64u64;
    ey, s: Int32;
begin
  wy.f := y;
  ey := Int32((wy.u shr 52) and $7FF) - 1023;
  s  := ey + 12;
  if ey >= 0 then begin
    if s >= 64 then begin Result := True; Exit; end;
    Result := (wy.u shl s) = 0;
    Exit;
  end;
  Result := (wy.u shl 1) = 0;
end;

function pcr_isint_pd(y: Double): Boolean; inline;
var wy: Tb64u64;
    ey, s: Int32;
begin
  wy.f := y;
  ey := Int32((wy.u shr 52) and $7FF) - 1023;
  s  := ey + 12;
  if ey >= 0 then begin
    if s >= 64 then begin Result := True; Exit; end;
    Result := (wy.u shl s) = 0;
    Exit;
  end;
  Result := (wy.u shl 1) = 0;
end;

function pcr_isodd_pd(y: Double): Boolean; inline;
var wy: Tb64u64;
    ey, s: Int32;
    oddb: UInt64;
begin
  wy.f := y;
  ey := Int32((wy.u shr 52) and $7FF) - 1023;
  s  := ey + 12;
  oddb := 0;
  if ey >= 0 then begin
    if (s < 64) and ((wy.u shl s) = 0) then
      oddb := (wy.u shr (64 - s)) and 1;
    if s = 64 then
      oddb := wy.u and 1;
  end;
  Result := oddb <> 0;
end;

function pcr_is_signaling_pd(x: Double): Boolean; inline;
var u: Tb64u64;
begin
  u.f := x;
  u.u := u.u xor UInt64($0008000000000000);
  Result := (u.u and UInt64($7FFFFFFFFFFFFFFF)) > UInt64($7FF8000000000000);
end;

function pcr_mulddd_pd(xh, xl, ch: Double; out l: Double): Double; inline;
var ahlh, ahhh, ahhl: Double;
begin
  ahlh := ch * xl;
  ahhh := ch * xh;
  ahhl := pcr_fma(ch, xh, -ahhh);
  ahhl := ahhl + ahlh;
  ch   := ahhh + ahhl;
  l    := (ahhh - ch) + ahhl;
  Result := ch;
end;

function pcr_get_mxcsr: DWord;
{$IFDEF AVX2}
var r: DWord;
begin
  asm
    stmxcsr r
  end;
  Result := r;
end;
{$ELSE}
begin
  Result := 0;
end;
{$ENDIF}

procedure pcr_set_mxcsr(flag: DWord);
{$IFDEF AVX2}
begin
  asm
    ldmxcsr flag
  end;
end;
{$ELSE}
begin
end;
{$ENDIF}

// ---------------------------------------------------------------------------
// Double-double and polynomial helpers (promoted from pascoremath32, task 0.9)
// ---------------------------------------------------------------------------

function pcr_poly12(z: Double; const c: array of Double): Double; inline;
var z2, z4, c0, c2, c4, c6, c8, c10: Double;
begin
  z2 := z * z; z4 := z2 * z2;
  c0  := c[0]  + z * c[1];
  c2  := c[2]  + z * c[3];
  c4  := c[4]  + z * c[5];
  c6  := c[6]  + z * c[7];
  c8  := c[8]  + z * c[9];
  c10 := c[10] + z * c[11];
  c0 := c0 + c2 * z2;
  c4 := c4 + c6 * z2;
  c8 := c8 + z2 * c10;
  c0 := c0 + z4 * (c4 + z4 * c8);
  Result := c0;
end;

function pcr_muldd(xh, xl, ch, cl: Double; out l: Double): Double; inline;
var
  ahlh, alhh, ahhh, ahhl: Double;
begin
  ahlh := ch * xl;
  alhh := cl * xh;
  ahhh := ch * xh;
  ahhl := pcr_fma(ch, xh, -ahhh);
  ahhl := ahhl + alhh + ahlh;
  ch := ahhh + ahhl;
  l := (ahhh - ch) + ahhl;
  Result := ch;
end;

function pcr_polydd(xh, xl: Double; n: Int32; const c: array of Double; out l: Double): Double; inline;
var
  i, i2: Int32;
  ch, cl, th, tl: Double;
begin
  i := n - 1;
  i2 := i * 2;
  ch := c[i2];
  cl := c[i2 + 1];
  Dec(i2,2);
  while i2 >= 8 do begin
    ch := pcr_muldd(xh, xl, ch, cl, cl);
    th := ch + c[i2];
    tl := (c[i2] - th) + ch;
    ch := th;
    cl := cl + tl + c[i2 + 1];
    Dec(i2,2);
    ch := pcr_muldd(xh, xl, ch, cl, cl);
    th := ch + c[i2];
    tl := (c[i2] - th) + ch;
    ch := th;
    cl := cl + tl + c[i2 + 1];
    Dec(i2,2);
    ch := pcr_muldd(xh, xl, ch, cl, cl);
    th := ch + c[i2];
    tl := (c[i2] - th) + ch;
    ch := th;
    cl := cl + tl + c[i2 + 1];
    Dec(i2,2);
    ch := pcr_muldd(xh, xl, ch, cl, cl);
    th := ch + c[i2];
    tl := (c[i2] - th) + ch;
    ch := th;
    cl := cl + tl + c[i2 + 1];
    Dec(i2,2);
  end;
  while i2 >= 0 do begin
    ch := pcr_muldd(xh, xl, ch, cl, cl);
    th := ch + c[i2];
    tl := (c[i2] - th) + ch;
    ch := th;
    cl := cl + tl + c[i2 + 1];
    Dec(i2,2);
  end;
  l := cl;
  Result := ch;
end;

procedure pcr_fasttwosum(var s, t: Double; a, b: Double); inline;
var s_tmp: Double;
begin
  s_tmp := a + b;
  t     := b - (s_tmp - a);
  s     := s_tmp;
end;

procedure pcr_a_mul(var hi, lo: Double; a, b: Double); inline;
var t_am: Double;
begin
  t_am := a * b;
  lo   := pcr_fma(a, b, -t_am);
  hi   := t_am;
end;

procedure pcr_s_mul(var hi, lo: Double; a, bh, bl: Double); inline;
var bl_sm: Double;
begin
  bl_sm := bl;             // save bl before pcr_a_mul may overwrite lo
  pcr_a_mul(hi, lo, a, bh);
  lo := pcr_fma(a, bl_sm, lo);
end;

procedure pcr_d_mul(var hi, lo: Double; ah, al, bh, bl: Double); inline;
var s_dm, t_dm, ah_dm: Double;
begin
  ah_dm := ah;             // save ah before pcr_a_mul may overwrite hi
  pcr_a_mul(hi, s_dm, ah_dm, bh);
  t_dm := pcr_fma(al, bh, s_dm);
  lo   := pcr_fma(ah_dm, bl, t_dm);
end;

end.
