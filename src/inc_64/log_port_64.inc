// ---------------------------------------------------------------------------
// pcr_log — correctly-rounded binary64 natural logarithm.
// Ported from core-math/src/binary64/log/log.c by Paul Zimmermann and
// Tom Hubrecht. See also dint.h.
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
// Tables are emitted by tools/log_gen.py into log_const.inc.
//
// Algorithm:
//   * fast path (cr_log_fast): degree-6 polynomial reduction with 363-entry
//     inverse / -log(r) tables, ~2^-68 absolute error.
//   * accurate path (cr_log_accurate): TDInt64 reduction + degree-13 dint
//     polynomial.
//
// Convention: pas-core-math TDInt64 .ex is offset by +1 vs log.c's dint.h
// (sin.c convention has significand in [0.5,1); log.c has [1,2)). Constants
// in log_const.inc are pre-shifted, and arithmetic ops (AddDInt/MulDInt/
// MulDIntInt) are convention-agnostic. See log_gen.py for details.
// ---------------------------------------------------------------------------

{$I inc_64/log_const_64.inc}

const
  cLogOFFSET: Int32 = 362;

// --- Fast path: cr_log_fast ----------------------------------------------------
// Returns (h, l) such that h + l ≈ log(2^e * v.f), |error| < 2^-68.22.
// `m` is the ieee mantissa with the implicit bit; `vf` is the v.f double.
procedure LogFast(out h, l: Double; e: Int32; vf: Double);
var
  v: Tb64u64;
  m: UInt64;
  i: Int32;
  y, r, l1, l2, z, ph, z2, p45, p23: Double;
  ee, hh, lo, et: Double;
begin
  v.f := vf;
  m := UInt64($10000000000000) + (v.u and UInt64($FFFFFFFFFFFFF));

  // if m >= 0x16a09e667f3bcd (sqrt(2) * 2^52) divide by 2 to avoid cancellation
  if m >= UInt64($16A09E667F3BCD) then begin
    e := e + 1;
    y := vf * 0.5;
    i := Int32(m shr 44);
  end else begin
    y := vf;
    i := Int32(m shr 43);
  end;

  r  := cLogInverse[i - cLogOFFSET].f;
  l1 := cLogInvH   [i - cLogOFFSET].f;
  l2 := cLogInvL   [i - cLogOFFSET].f;
  z  := pcr_fma(r, y, -1.0);

  z2  := z * z;
  p45 := pcr_fma(cLogP[5].f, z, cLogP[4].f);
  p23 := pcr_fma(cLogP[3].f, z, cLogP[2].f);
  ph  := pcr_fma(p45, z2, p23);
  ph  := pcr_fma(ph, z, cLogP[1].f);
  ph  := ph * z2;

  ee := e;
  // hh = fma(ee, log2_h, l1)
  hh := pcr_fma(ee, cLogLog2h.f, l1);
  // fast_two_sum(h, l, hh, z)
  h := hh + z;
  et := h - hh;          // exact
  l := z - et;           // exact

  // l = ph + (l + l2)
  lo := l + l2;
  l := ph + lo;

  // l = fma(ee, log2_l, l)
  l := pcr_fma(ee, cLogLog2l.f, l);
end;

// --- Accurate path polynomial p_2(z) ------------------------------------------
procedure LogPolyAccurate(out r: TDInt64; const z: TDInt64);
var
  k: Integer;
begin
  CpDInt(r, cLogP2[0]);
  for k := 1 to 12 do begin
    MulDInt(r, z, r);
    AddDInt(r, cLogP2[k], r);
  end;
  MulDInt(r, z, r);
end;

// --- Accurate path: log_2(x) = log(x) ----------------------------------------
procedure LogTwo(out r: TDInt64; var x: TDInt64);
var
  E: Int64;
  i: UInt32;
 z, p: TDInt64;
begin
  // x.ex in pascoremath convention; shift by -1 conceptually, work in log.c idx
  // The C code: E = x->ex; i = x->hi >> 55;
  //   if (x->hi > 0xb504f333f9de6484) { E++; i >>= 1; }
  //   x->ex -= E;
  //
  // In our convention E corresponds to the "log.c ex" = pascoremath_ex - 1.
  // We want to match: after `x.ex -= E` (in log.c terms), x has ex_log = 0,
  // i.e. ex_pas = 1. Compute E using x.ex - 1 (the log.c-equivalent ex).
  E := x.ex - 1;
  i := UInt32(x.hi shr 55);
  if x.hi > UInt64($B504F333F9DE6484) then begin
    E := E + 1;
    i := i shr 1;
  end;
  x.ex := x.ex - E;     // now x has ex_pas = 1, i.e. value in [1/sqrt(2), sqrt(2))

  MulDInt(z, x, cLogInverse2[i - 128]);
  AddDInt(z, cLogMOneDint, z);

  // r = E * LOG2
  MulDIntInt(r, E, cLogLog2Dint);

  LogPolyAccurate(p, z);
  AddDInt(p, cLogLogInv2[i - 128], p);

  AddDInt(r, p, r);
end;

function LogAccurate(x: Double): Double;
var
  X64, Y64: TDInt64;
begin
  if x = 1.0 then begin Result := 0.0; Exit; end;
  DIntFromD(X64, x);
  LogTwo(Y64, X64);
  Result := DToD(Y64);
end;

function pcr_log(x: Double): Double;
var
  v: Tb64u64;
  e: Int32;
  h, l, left, right: Double;
  err: Double;
  minf: Tb64u64;
begin
  v.f := x;
  e := Int32(v.u shr 52) - $3FF;

  // x <= 0 or NaN/Inf or subnormal — handle special cases
  if (e >= $400) or (e = -$3FF) then begin
    minf.u := UInt64($FFF) shl 52;
    if (e = $400) or ((e = $C00) and (x <> minf.f)) then begin
      Result := x + x;  // +Inf or NaN
      Exit;
    end;
    if x <= 0.0 then begin
      if x < 0 then begin
        Result := cNaNDouble;
        Exit;
      end else begin
        // log(+/-0) = -Inf via 1.0 / -0.0
        Result := 1.0 / Double(-0.0);
        Exit;
      end;
    end;
    if e = -$3FF then begin
      v.f := v.f * 4503599627370496.0; // 0x1p52
      e := Int32(v.u shr 52) - $3FF - 52;
    end;
  end;

  // normalize v in [1, 2)
  v.u := (UInt64($3FF) shl 52) or (v.u and UInt64($FFFFFFFFFFFFF));

  // x = 1.0 exactly
  if (v.u = UInt64($3FF0000000000000)) and (e = 0) then begin
    Result := 0.0;
    Exit;
  end;

  LogFast(h, l, e, v.f);

  err := cLogErr.f;
  left  := h + (l - err);
  right := h + (l + err);
  if left = right then begin
    Result := left;
    Exit;
  end;

  Result := LogAccurate(x);
end;
