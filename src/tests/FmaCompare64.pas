// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// SPDX-License-Identifier: MIT
//
// FmaCompare64: compare hardware pcr_fma (VFMADD213SD) against the pure-Pascal
// pcr_fma_pascal emulation. Build with -dAVX2 so pcr_fma uses the VFMADD path
// and the two implementations are genuinely different code.
//
// Usage:
//   FmaCompare64                 - run bundled FMA test suite
//   FmaCompare64 <a> <b> <c>     - compare one triple (decimal or $hex u64)
{$I ../pascoremath.inc}
program FmaCompare64;

uses
  SysUtils, Math, pascoremathtypes, pascoremathhelperfuncs;

var
  TotalCases: Int64 = 0;
  MismatchCases: Int64 = 0;

function ParseArg(const s: string): Double;
var
  v: Tb64u64;
  code: Integer;
begin
  if (Length(s) > 1) and (s[1] = '$') then
  begin
    Val(s, v.u, code);
    if code <> 0 then
      raise Exception.CreateFmt('Bad hex u64: %s', [s]);
    Result := v.f;
  end
  else
  begin
    Val(s, Result, code);
    if code <> 0 then
      raise Exception.CreateFmt('Bad double: %s', [s]);
  end;
end;

function HexU64(x: Double): string;
var v: Tb64u64;
begin v.f := x; Result := '$' + IntToHex(v.u, 16); end;

function BitsDiff(a, b: Double): Int64;
var va, vb: Tb64u64;
begin
  va.f := a; vb.f := b;
  Result := Abs(Int64(va.u) - Int64(vb.u));
end;

procedure CompareOne(a, b, c: Double; Verbose: Boolean);
var
  hw, sw: Double;
  vh, vs: Tb64u64;
  differ: Boolean;
begin
  Inc(TotalCases);
  hw := pcr_fma(a, b, c);
  sw := pcr_fma_pascal(a, b, c);
  vh.f := hw; vs.f := sw;
  differ := (vh.u <> vs.u) and
            not (IsNaN(hw) and IsNaN(sw));
  if differ then Inc(MismatchCases);
  if Verbose or differ then
  begin
    WriteLn(Format('a=%.17g b=%.17g c=%.17g', [a, b, c]));
    WriteLn('  a=', HexU64(a), '  b=', HexU64(b), '  c=', HexU64(c));
    WriteLn('  hw pcr_fma        = ', Format('%.17g', [hw]), '  ', HexU64(hw));
    WriteLn('  sw pcr_fma_pascal = ', Format('%.17g', [sw]), '  ', HexU64(sw));
    if differ then
      WriteLn('  ** MISMATCH: ulp-distance = ', BitsDiff(hw, sw))
    else
      WriteLn('  MATCH');
  end;
end;

function RandU64: UInt64;
begin
  Result := (UInt64(Random($40000000)) shl 34) xor
            (UInt64(Random($40000000)) shl 4)  xor
             UInt64(Random($10));
end;

// Build a random normal double whose biased exponent lies in [EMinB, EMaxB].
// Keep EMaxB <= 2019 so Veltkamp K*x (K = 2^27+1) cannot overflow.
function RandDouble(EMinB, EMaxB: Integer): Double;
var
  v: Tb64u64;
  eb: UInt64;
begin
  v.u := RandU64;
  eb := UInt64(EMinB + Random(EMaxB - EMinB + 1));
  v.u := (v.u and $800FFFFFFFFFFFFF) or (eb shl 52);
  Result := v.f;
end;

procedure RunSuite;
const
  // Rounding-correctness cases that stay inside the emulation's safe domain
  // (|biased-exp| <= 2019 so Dekker's K*x does not overflow).
  // Covers: trivial exactness, tiebreakers, cancellation, subnormals,
  // v.l=0 branch and exact-halfway (0.875/1.125) branches in pcr_fma_pascal.
  Suite: array[0..12, 0..2] of Double = (
    (1.0, 1.0, 0.0),
    (1.5, 1.5, -2.25),                          // exact -> 0
    (0.1, 0.1, 0.0),
    (1.0 + 1e-16, 1.0 - 1e-16, -1.0),
    (1e150, 1e-150, 1.0),                       // wide-range, no overflow in K*x
    (3.0, 1.0 / 3.0, -1.0),                     // 1/3 not exact -> tiny residual
    (1.0000000000000002, 1.0000000000000002, -1.0),
    (0.5, 0.5, 0.25),                           // 0.5 exactly -> 0.5
    (-1.0, 1.0, 1.0),                           // exact -> 0
    (1e-160, 1e-160, 0.0),                      // toward subnormal tail
    (1.0, 1.0, -0.9999999999999999),
    (0.9999999999999999, 1.0000000000000002, -1.0),
    (3.0, 0.1, -0.3)                            // classic FMA tiebreaker
  );
var
  i, k: Integer;
  a, b, c, p: Double;
  vp, vc: Tb64u64;
begin
  WriteLn('-- Fixed suite --');
  for i := 0 to High(Suite) do
    CompareOne(Suite[i][0], Suite[i][1], Suite[i][2], True);

  // Broad sweep: random normals, exponents bounded so Veltkamp split stays finite.
  WriteLn;
  WriteLn('-- Random sweep (normal range, 200000 triples) --');
  RandSeed := 42;
  for k := 1 to 200000 do
  begin
    a := RandDouble(512, 1534);   // unbiased [-511, +511]
    b := RandDouble(512, 1534);
    c := RandDouble(512, 1534);
    CompareOne(a, b, c, False);
  end;

  // Cancellation sweep: c = -round(a*b) +/- a few ULPs. This is the regime
  // where correct FMA vs naive a*b+c differ the most, and it's what
  // Newton-style refinements (e.g. atanh, rsqrt) actually feed the FMA.
  WriteLn;
  WriteLn('-- Cancellation sweep (c ~= -a*b, 200000 triples) --');
  RandSeed := 12345;
  for k := 1 to 200000 do
  begin
    a := RandDouble(900, 1146);   // unbiased [-123, +123], product exp bounded
    b := RandDouble(900, 1146);
    p := a * b;
    if IsNan(p) or IsInfinite(p) then continue;
    vp.f := -p;
    // Perturb the last few bits of -p to get near-cancellation triples.
    vc.u := vp.u xor UInt64(Random(256));
    c := vc.f;
    CompareOne(a, b, c, False);
  end;

  WriteLn;
  WriteLn(Format('Total: %d   Mismatches: %d', [TotalCases, MismatchCases]));
end;

var
  a, b, c: Double;
begin
  // Mask FP exceptions so pcr_fma_pascal's Dekker product doesn't abort on
  // overflow/underflow inputs (it still returns inf/nan, which is what we want
  // to compare).
  SetExceptionMask([exInvalidOp, exDenormalized, exZeroDivide, exOverflow,
                    exUnderflow, exPrecision]);
  WriteLn('FmaCompare64 - pcr_fma vs pcr_fma_pascal');
  {$IFDEF AVX2}
  WriteLn('Built with -dAVX2: pcr_fma uses VFMADD213SD (hardware).');
  {$ELSE}
  WriteLn('WARNING: built WITHOUT -dAVX2; pcr_fma delegates to pcr_fma_pascal.');
  WriteLn('Rebuild with -dAVX2 for a meaningful comparison.');
  {$ENDIF}

  if ParamCount = 3 then
  begin
    a := ParseArg(ParamStr(1));
    b := ParseArg(ParamStr(2));
    c := ParseArg(ParamStr(3));
    CompareOne(a, b, c, True);
    Halt(Ord(MismatchCases <> 0));
  end
  else if ParamCount = 0 then
  begin
    RunSuite;
    Halt(Ord(MismatchCases <> 0));
  end
  else
  begin
    WriteLn('Usage: FmaCompare64 [ <a> <b> <c> ]');
    WriteLn('  no args   -> run bundled suite + random sweep');
    WriteLn('  3 args    -> compare a single triple (decimal or $hexU64)');
    Halt(2);
  end;
end.
