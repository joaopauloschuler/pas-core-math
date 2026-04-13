// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//                                                                                                                                                                                                      
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
{$I pascoremath.inc}
unit pascoremathtypes;

interface

type
  TUInt128 = record
    lo, hi: UInt64;
  end;

  Tb32u32 = record
    case Boolean of
      False: (f: Single);
      True:  (u: LongWord);
  end;

  Tb64u64 = record
    case Boolean of
      False: (f: Double);
      True:  (u: UInt64);
  end;

operator +(const a: TUInt128; b: UInt64): TUInt128; inline;

function Mulu64u64(a, b: UInt64): TUInt128; inline;

implementation

operator +(const a: TUInt128; b: UInt64): TUInt128; inline;
begin
  Result.lo := a.lo + b;
  Result.hi := a.hi + UInt64(Result.lo < b);  // carry
end;

function Mulu64u64(a, b: UInt64): TUInt128; inline;
{$IFDEF CPUX86_64}
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
var
  a_lo, a_hi, b_lo, b_hi: UInt64;
  p0, p1, p2, p3: UInt64;
  mid: UInt64;
begin
  a_lo := a and $FFFFFFFF;
  a_hi := a shr 32;
  b_lo := b and $FFFFFFFF;
  b_hi := b shr 32;

  p0 := a_lo * b_lo;
  p1 := a_lo * b_hi;
  p2 := a_hi * b_lo;
  p3 := a_hi * b_hi;

  // Accumulate middle terms (bits 32..95)
  mid := (p0 shr 32) + (p1 and $FFFFFFFF) + (p2 and $FFFFFFFF);

  Result.lo := (p0 and $FFFFFFFF) or (mid shl 32);
  Result.hi := p3 + (p1 shr 32) + (p2 shr 32) + (mid shr 32);
end;
{$ENDIF}

end.
