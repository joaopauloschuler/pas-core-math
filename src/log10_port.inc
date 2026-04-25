// ---------------------------------------------------------------------------
// pcr_log10 — correctly-rounded binary64 base-10 logarithm.
// Ported from core-math/src/binary64/log10/log10.c.
//
// Reuses the log/ tables (cLogInverse, cLogInvH/L, cLogP, cLogInverse2,
// cLogLogInv2, cLogP2, cLogLog2Dint, cLogMOneDint, cLogLog2h/l) declared
// in log_const.inc — log10 is log composed with multiplication by 1/log(10).
//
// Algorithm:
//   * fast path: same polynomial reduction as log, then d_mul by
//     (ONE_OVER_LOG10_H, ONE_OVER_LOG10_L). Absolute error < 2^-68.98.
//   * accurate path: LogTwo (defined in log_port.inc) followed by
//     MulDInt(Y, ONE_OVER_LOG10, Y).
//
// log10(10^n) for 0 <= n <= 22 is detected up-front with a perfect-hash on
// the biased exponent so x=10^n returns n exactly (no spurious inexact).
// ---------------------------------------------------------------------------

const
  cLog10OFFSET: Int32 = 362;
  // 0x1.bcb7b1526e50ep-2 = 1/log(10) high part
  cLog10OneOverLog10H: Tb64u64 = (u:$3FDBCB7B1526E50E);
  // 0x1.95355baaafad3p-57 = 1/log(10) low part
  cLog10OneOverLog10L: Tb64u64 = (u:$3C695355BAAAFAD3);
  // err = 0x1.04p-69
  cLog10Err: Tb64u64 = (u:$3BA0400000000000);

  // ONE_OVER_LOG10 dint — C dint has ex=-2; pascoremath convention is +1, so ex=-1.
  cLog10OneOverLog10Dint: TDInt64 =
    (hi:QWord($DE5BD8A937287195); lo:$355BAAAFAD33DC32; ex:-1; sgn:0);

  // pow10[0..22] = bit patterns of 10^n; entries 23..31 are zero (sentinel).
  cLog10Pow10: array[0..31] of Tb64u64 = (
    (u:$3FF0000000000000), (u:$4024000000000000), (u:$4059000000000000),
    (u:$408F400000000000), (u:$40C3880000000000), (u:$40F86A0000000000),
    (u:$412E848000000000), (u:$416312D000000000), (u:$4197D78400000000),
    (u:$41CDCD6500000000), (u:$4202A05F20000000), (u:$42374876E8000000),
    (u:$426D1A94A2000000), (u:$42A2309CE5400000), (u:$42D6BCC41E900000),
    (u:$430C6BF526340000), (u:$4341C37937E08000), (u:$4376345785D8A000),
    (u:$43ABC16D674EC800), (u:$43E158E460913D00), (u:$4415AF1D78B58C40),
    (u:$444B1AE4D6E2EF50), (u:$4480F0CF064DD592), (u:$0000000000000000),
    (u:$0000000000000000), (u:$0000000000000000), (u:$0000000000000000),
    (u:$0000000000000000), (u:$0000000000000000), (u:$0000000000000000),
    (u:$0000000000000000), (u:$0000000000000000)
  );

// --- Fast path: cr_log10_fast -------------------------------------------------
// Returns (h, l) such that h + l ≈ log10(2^e * v.f), |error| < 2^-68.98.
procedure Log10Fast(var h, l: Double; e: Int32; vf: Double);
var
  v: Tb64u64;
  m: UInt64;
  c: Int32;
  i: Int32;
  y, r, l1, l2, z, ph, z2, p45, p23: Double;
  ee, hh, lo, et: Double;
begin
  v.f := vf;
  m := UInt64($10000000000000) + (v.u and UInt64($FFFFFFFFFFFFF));

  if m >= UInt64($16A09E667F3BCD) then begin
    c := 1; e := e + 1;
    y := vf * 0.5;
    i := Int32(m shr 44);
  end else begin
    c := 0;
    y := vf;
    i := Int32(m shr 43);
  end;

  r  := cLogInverse[i - cLog10OFFSET].f;
  l1 := cLogInvH   [i - cLog10OFFSET].f;
  l2 := cLogInvL   [i - cLog10OFFSET].f;
  z  := pcr_fma(r, y, -1.0);

  z2  := z * z;
  p45 := pcr_fma(cLogP[5].f, z, cLogP[4].f);
  p23 := pcr_fma(cLogP[3].f, z, cLogP[2].f);
  ph  := pcr_fma(p45, z2, p23);
  ph  := pcr_fma(ph, z, cLogP[1].f);
  ph  := ph * z2;

  ee := e;
  hh := pcr_fma(ee, cLogLog2h.f, l1);
  // fast_two_sum(h, l, hh, z)
  h := hh + z;
  et := h - hh;
  l := z - et;

  // l = ph + (l + l2)
  lo := l + l2;
  l := ph + lo;

  // l = fma(ee, log2_l, l)
  l := pcr_fma(ee, cLogLog2l.f, l);

  // Divide by log(10) via d_mul with (ONE_OVER_LOG10_H, ONE_OVER_LOG10_L).
  pcr_d_mul(h, l, h, l, cLog10OneOverLog10H.f, cLog10OneOverLog10L.f);
end;

function Log10Accurate(x: Double): Double;
var
  X64, Y64: TDInt64;
begin
  if x = 1.0 then begin Result := 0.0; Exit; end;
  DIntFromD(X64, x);
  LogTwo(Y64, X64);
  // Multiply by 1/log(10).
  MulDInt(Y64, cLog10OneOverLog10Dint, Y64);
  Result := DToD(Y64);
end;

function pcr_log10(x: Double): Double;
var
  v: Tb64u64;
  e: Int32;
  h, l, left, right, err: Double;
  minf: Tb64u64;
  n: UInt32;
begin
  v.f := x;
  e := Int32(v.u shr 52) - $3FF;

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
        Result := 1.0 / Double(-0.0);
        Exit;
      end;
    end;
    if e = -$3FF then begin
      v.f := v.f * 4503599627370496.0; // 0x1p52
      e := Int32(v.u shr 52) - $3FF - 52;
    end;
  end;

  // Detect x = 10^n exactly for 0 <= n <= 22 (returns n with no inexact).
  n := (UInt32($267CFA2) * UInt32(e + 1)) shr 27;
  if (n <= 31) and (v.u = cLog10Pow10[n].u) then begin
    Result := n;
    Exit;
  end;

  // Normalize v in [1, 2).
  v.u := (UInt64($3FF) shl 52) or (v.u and UInt64($FFFFFFFFFFFFF));

  Log10Fast(h, l, e, v.f);

  err := cLog10Err.f;
  left  := h + (l - err);
  right := h + (l + err);
  if left = right then begin
    Result := left;
    Exit;
  end;

  Result := Log10Accurate(x);
end;
