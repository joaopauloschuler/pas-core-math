{$I pascoremath.inc}
unit pascoremathhelperfuncs;

interface

uses
  Math, SysUtils, pascoremathtypes;

// Fused multiply-add (double-rounding approximation via 80-bit Extended)
function pcr_fmaf(x, y, z: Single): Single; inline;
function pcr_fma(x, y, z: Double): Double; inline;

// Absolute value
function pcr_fabsf(x: Single): Single; inline;
function pcr_fabs(x: Double): Double; inline;

// Copy sign of y to magnitude of x
function pcr_copysignf(x, y: Single): Single; inline;
function pcr_copysign(x, y: Double): Double; inline;

// Square root
function pcr_sqrtf(x: Single): Single; inline;
function pcr_sqrt(x: Double): Double; inline;

// Round to nearest even integer (banker's rounding)
function pcr_roundevenf(x: Single): Single; inline;
function pcr_roundeven(x: Double): Double; inline;

// NaN-aware maximum
function pcr_fmaxf(x, y: Single): Single; inline;
function pcr_fmax(x, y: Double): Double; inline;

// NaN-aware minimum
function pcr_fminf(x, y: Single): Single; inline;
function pcr_fmin(x, y: Double): Double; inline;

// Bit scan forward: index of lowest set bit (0-based); undefined for x=0
function pcr_bsf32(x: UInt32): Integer; inline;

// Bit scan reverse: index of highest set bit (0-based) = floor(log2(x)); undefined for x=0
function pcr_bsr32(x: UInt32): Integer; inline;

// Return a NaN (tagp is ignored, matches C nan()/nanf() signature)
function pcr_nanf(const tagp: PAnsiChar): Single; inline;
function pcr_nan(const tagp: PAnsiChar): Double; inline;

// Raise floating-point exceptions
procedure pcr_feraiseexcept_invalid; inline;
procedure pcr_feraiseexcept_divbyzero; inline;

implementation

function pcr_fmaf(x, y, z: Single): Single; inline;
begin
  // Note: double-rounding approximation; not a true IEEE FMA.
  // Uses 80-bit extended precision to reduce rounding error.
  Result := Single(Extended(x) * Extended(y) + Extended(z));
end;

function pcr_fma(x, y, z: Double): Double; inline;
begin
  // Note: double-rounding approximation; not a true IEEE FMA.
  // Uses 80-bit extended precision to reduce rounding error.
  Result := Double(Extended(x) * Extended(y) + Extended(z));
end;

function pcr_fabsf(x: Single): Single; inline;
begin
  Result := Abs(x);
end;

function pcr_fabs(x: Double): Double; inline;
begin
  Result := Abs(x);
end;

function pcr_copysignf(x, y: Single): Single; inline;
var
  vx, vy: Tb32u32;
begin
  vx.f := x;
  vy.f := y;
  vx.u := (vx.u and $7FFFFFFF) or (vy.u and $80000000);
  Result := vx.f;
end;

function pcr_copysign(x, y: Double): Double; inline;
var
  vx, vy: Tb64u64;
begin
  vx.f := x;
  vy.f := y;
  vx.u := (vx.u and $7FFFFFFFFFFFFFFF) or (vy.u and $8000000000000000);
  Result := vx.f;
end;

function pcr_sqrtf(x: Single): Single; inline;
begin
  Result := Sqrt(x);
end;

function pcr_sqrt(x: Double): Double; inline;
begin
  Result := Sqrt(x);
end;

function pcr_roundevenf(x: Single): Single; inline;
// Round to nearest even using bit manipulation.
// For |x| >= 2^23 the value is already an integer.
var
  v: Tb32u32;
  e, shift: Integer;
  mask, frac, half: LongWord;
begin
  v.f := x;
  e := Integer((v.u shr 23) and $FF) - 127;  // unbiased exponent
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
      // Integer part is odd => round up
      v.u := (v.u and (not mask)) + (1 shl shift)
    else
      // Integer part is even => round down
      v.u := v.u and (not mask);
  end;
  Result := v.f;
end;

function pcr_roundeven(x: Double): Double; inline;
var
  v: Tb64u64;
  e: Integer;
  shift: Integer;
  mask, half: UInt64;
  frac: UInt64;
begin
  v.f := x;
  e := Integer((v.u shr 52) and $7FF) - 1023;
  if e >= 52 then
  begin
    Result := x;
    Exit;
  end;
  if e < 0 then
  begin
    if e = -1 then
    begin
      if (v.u and $7FFFFFFFFFFFFFFF) = $3FE0000000000000 then
        Result := 0.0
      else if Abs(x) < 0.5 then
        Result := 0.0
      else
        Result := pcr_copysign(1.0, x);
    end
    else
      Result := 0.0;
    Exit;
  end;
  shift := 52 - e;
  mask  := (UInt64(1) shl shift) - 1;
  frac  := v.u and mask;
  half  := UInt64(1) shl (shift - 1);

  if frac < half then
    v.u := v.u and (not mask)
  else if frac > half then
    v.u := (v.u and (not mask)) + (UInt64(1) shl shift)
  else
  begin
    if (v.u and (UInt64(1) shl shift)) <> 0 then
      v.u := (v.u and (not mask)) + (UInt64(1) shl shift)
    else
      v.u := v.u and (not mask);
  end;
  Result := v.f;
end;

function pcr_fmaxf(x, y: Single): Single; inline;
begin
  if IsNan(x) then Result := y
  else if IsNan(y) then Result := x
  else if x > y then Result := x
  else Result := y;
end;

function pcr_fmax(x, y: Double): Double; inline;
begin
  if IsNan(x) then Result := y
  else if IsNan(y) then Result := x
  else if x > y then Result := x
  else Result := y;
end;

function pcr_fminf(x, y: Single): Single; inline;
begin
  if IsNan(x) then Result := y
  else if IsNan(y) then Result := x
  else if x < y then Result := x
  else Result := y;
end;

function pcr_fmin(x, y: Double): Double; inline;
begin
  if IsNan(x) then Result := y
  else if IsNan(y) then Result := x
  else if x < y then Result := x
  else Result := y;
end;

function pcr_bsf32(x: UInt32): Integer; inline;
{$IFDEF CPUX86_64}
var
  r: LongWord;
begin
  asm
    bsf  eax, x
    mov  r, eax
  end;
  Result := r;
end;
{$ELSE}
var
  i: Integer;
begin
  for i := 0 to 31 do
    if (x and (UInt32(1) shl i)) <> 0 then
    begin
      Result := i;
      Exit;
    end;
  Result := -1;  // undefined for x=0
end;
{$ENDIF}

function pcr_bsr32(x: UInt32): Integer; inline;
{$IFDEF CPUX86_64}
var
  r: LongWord;
begin
  asm
    bsr  eax, x
    mov  r, eax
  end;
  Result := r;
end;
{$ELSE}
var
  i: Integer;
begin
  Result := -1;
  for i := 31 downto 0 do
    if (x and (UInt32(1) shl i)) <> 0 then
    begin
      Result := i;
      Exit;
    end;
end;
{$ENDIF}

function pcr_nanf(const tagp: PAnsiChar): Single; inline;
begin
  Result := Single(NaN);
end;

function pcr_nan(const tagp: PAnsiChar): Double; inline;
begin
  Result := NaN;
end;

procedure pcr_feraiseexcept_invalid; inline;
var
  x: Single;
begin
  // Raise FE_INVALID by computing 0/0
  x := 0.0;
  x := x / x;
end;

procedure pcr_feraiseexcept_divbyzero; inline;
var
  x: Single;
begin
  // Raise FE_DIVBYZERO by computing 1/0
  x := 0.0;
  x := 1.0 / x;
end;

end.
