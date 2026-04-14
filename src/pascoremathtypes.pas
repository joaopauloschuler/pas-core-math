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

const
  cNaNSingle: Single = 0.0/0.0;
  cNaNDouble: Double = 0.0/0.0;
  // 2^(-127): subnormal Single used to trigger IEEE 754 underflow via multiplication
  cUnderflowSingle: Single = 5.877471754111438e-39;

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
  MulLo, Temp1, Temp2: UInt64;
begin
  MulLo := uint64(uint32(a)) * uint64(uint32(b));
  Temp1 := (a shr 32) * uint64(uint32(b)) + (MulLo shr 32);
  Temp2 := uint64(uint32(a)) * (b shr 32) + uint64(uint32(Temp1));
  Result.lo := ((Temp2 and $FFFFFFFF) shl 32) or (MulLo and $FFFFFFFF);
  Result.hi := (a shr 32) * (b shr 32) + (Temp1 shr 32) + (Temp2 shr 32);
end;
{$ENDIF}

end.
