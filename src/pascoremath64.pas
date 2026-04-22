// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
{$I pascoremath.inc}
unit pascoremath64;

interface

uses Math, pascoremathtypes, pascoremathhelperfuncs, ccoremath64;

// ── Ported functions ──────────────────────────────────────────────────────────
function pcr_rsqrt(x: Double): Double;
function pcr_cbrt(x: Double): Double;

// ── Stub functions (delegate to C reference until ported) ────────────────────
function  pcr_acos(x: Double): Double; inline;
function  pcr_acosh(x: Double): Double; inline;
function  pcr_acospi(x: Double): Double; inline;
function  pcr_asin(x: Double): Double; inline;
function  pcr_asinh(x: Double): Double; inline;
function  pcr_asinpi(x: Double): Double; inline;
function  pcr_atan(x: Double): Double; inline;
function  pcr_atanh(x: Double): Double; inline;
function  pcr_atanpi(x: Double): Double; inline;
// pcr_cbrt declared in ported section above
function  pcr_cos(x: Double): Double; inline;
function  pcr_cosh(x: Double): Double; inline;
function  pcr_cospi(x: Double): Double; inline;
function  pcr_erf(x: Double): Double; inline;
function  pcr_erfc(x: Double): Double; inline;
function  pcr_exp(x: Double): Double; inline;
function  pcr_exp10(x: Double): Double; inline;
function  pcr_exp10m1(x: Double): Double; inline;
function  pcr_exp2(x: Double): Double; inline;
function  pcr_exp2m1(x: Double): Double; inline;
function  pcr_expm1(x: Double): Double; inline;
function  pcr_lgamma(x: Double): Double; inline;
function  pcr_log(x: Double): Double; inline;
function  pcr_log10(x: Double): Double; inline;
function  pcr_log10p1(x: Double): Double; inline;
function  pcr_log1p(x: Double): Double; inline;
function  pcr_log2(x: Double): Double; inline;
function  pcr_log2p1(x: Double): Double; inline;
function  pcr_sin(x: Double): Double; inline;
function  pcr_sinh(x: Double): Double; inline;
function  pcr_sinpi(x: Double): Double; inline;
function  pcr_tan(x: Double): Double; inline;
function  pcr_tanh(x: Double): Double; inline;
function  pcr_tanpi(x: Double): Double; inline;
function  pcr_tgamma(x: Double): Double; inline;
function  pcr_atan2(y, x: Double): Double; inline;
function  pcr_atan2pi(y, x: Double): Double; inline;
function  pcr_hypot(x, y: Double): Double; inline;
function  pcr_pow(x, y: Double): Double; inline;
procedure pcr_sincos(x: Double; out s, c: Double); inline;

implementation

// Constants used in pcr_rsqrt / RsqrtRefine.
// $FFF0000000000000 = -Infinity bit pattern; negative NaN have .u > this.
const
  cRsqrtNegNanThresh: Tb64u64 = (u: $FFF0000000000000);

// ---------------------------------------------------------------------------
// Slow-path refine for pcr_rsqrt.
// Ported from as_rsqrt_refine() in core-math/src/binary64/rsqrt/rsqrt.c.
// Not marked inline (equivalent to C __attribute__((noinline))).
// ---------------------------------------------------------------------------
function RsqrtRefine(rf, a: Double): Double;
var
  ir, ia: Tb64u64;
  nz: Int32;
  e_sub: Int64;
  mode: TFPURoundingMode;
  e: Int32;
  rm, am: UInt64;
  rt, rrt, prrt, rts, tt: TUInt128;
  rth, rtl: UInt64;
  t0, t1: UInt64;
  s_val: Int64;
  dd: Int64;
  am2, am20: UInt64;
  mask: UInt64;
  rt_shifted: TUInt128;
  borrow: UInt64;
  inc_val: UInt64;
begin
  ir.f := rf;
  ia.f := a;

  // Normalize subnormal inputs to fake-normal form for the bit-manipulation below.
  // Caller guarantees ia.u <> 0.
  if ia.u < UInt64(1) shl 52 then begin
    nz    := pcr_clzll(ia.u);
    ia.u  := ia.u shl (nz - 11);
    ia.u  := ia.u and UInt64($000FFFFFFFFFFFFF);  // clear sign+exp
    e_sub := Int64(nz) - 12;                       // fake exponent
    ia.u  := ia.u or UInt64(e_sub shl 52);
  end;

  // If the input is exactly a power of 2 the initial estimate rf is already
  // correctly rounded; skip the expensive refinement.
  if ia.u shl 11 = UInt64(1) shl 63 then begin
    Result := rf;
    Exit;
  end;

  mode := pcr_GetRoundMode;
  e    := Int32((ia.u shr 52) and 1);               // LSB of exponent field
  rm   := (ir.u shl 11 or (UInt64(1) shl 63)) shr 11;
  am   := ((ia.u and UInt64($000FFFFFFFFFFFFF))
            or (UInt64(1) shl 52)) shl (5 - e);

  // rt = rm * am  (128-bit)
  rt   := Mulu64u64(rm, am);
  rth  := rt.hi;
  rtl  := rt.lo;

  // rrt = rtl * rm, then reassemble with corrected high word
  rrt    := Mulu64u64(rtl, rm);
  t0     := rrt.lo;
  t1     := rrt.hi + rth * rm;   // 64×64 low 64 bits (upper bits discarded)
  rrt.hi := t1;
  rrt.lo := t0;

  s_val := Int64(rrt.hi shr 63);   // 0 or 1
  dd    := 1 - 2 * s_val;          // +1 or -1

  // rts = ((rt << 1) ^ (-s)) + s
  // -s as u128: all-zeros when s=0, all-ones when s=1
  rt_shifted := rt;
  ShlU128(rt_shifted, 1);
  mask   := UInt64(0) - UInt64(s_val);    // 0x000...0 or 0xFFF...F
  rts.hi := rt_shifted.hi xor mask;
  rts.lo := rt_shifted.lo xor mask;
  rts.lo := rts.lo + UInt64(s_val);
  rts.hi := rts.hi + UInt64(rts.lo < UInt64(s_val));  // carry

  am2  := am shl 1;
  am20 := UInt64(0) - am;   // two's-complement negation of am

  repeat
    ir.u  := ir.u - UInt64(dd);
    prrt  := rrt;
    am20  := am20 + am2;
    // tt = rts - (u128)am20
    borrow := UInt64(rts.lo < am20);
    tt.lo  := rts.lo - am20;
    tt.hi  := rts.hi - borrow;
    SubU128(rrt, rrt, tt);
  until (prrt.hi xor rrt.hi) and UInt64($8000000000000000) <> 0;

  // Undo the last step if rrt ended up non-negative (bit 127 = 0)
  if rrt.hi and UInt64($8000000000000000) = 0 then begin
    ir.u := ir.u + UInt64(dd);
    rrt  := prrt;
  end;

  if mode = rmNearest then begin
    rm  := (ir.u shl 11 or (UInt64(1) shl 63)) shr 11;
    rt  := Mulu64u64(rm, am);
    rrt := rrt + (am shr 2);   // operator+(TUInt128, UInt64) — handles carry
    AddU128(rrt, rrt, rt);
    inc_val := rrt.hi shr 63;
    ir.u := ir.u + inc_val;
  end else begin
    if mode = rmUp then
      ir.u := ir.u + 1;
  end;

  Result := ir.f;
end;

// ---------------------------------------------------------------------------
// pcr_rsqrt — correctly-rounded reciprocal square root (binary64).
// Ported from cr_rsqrt() in core-math/src/binary64/rsqrt/rsqrt.c.
// ---------------------------------------------------------------------------
function pcr_rsqrt(x: Double): Double;
var
  ix: Tb64u64;
  r, rx, drx, h, dr, rf: Double;
  idr, ir: Tb64u64;
  aidr, mid: UInt64;
begin
  ix.f := x;
  r    := 0.0;  // silence compiler warning

  if ix.u < UInt64(1) shl 52 then begin
    // 0 <= x < 0x1p-1022  (subnormal or +0)
    if ix.u <> 0 then
      r := pcr_sqrt(x) / x   // non-zero subnormal: 1/sqrt(x) = sqrt(x)/x
    else begin
      Result := Infinity;     // x = +0: pole error → +Inf
      Exit;
    end;
  end else if ix.u >= UInt64($7FF0000000000000) then begin
    // NaN, Inf, or negative
    if ix.u shl 1 = 0 then begin
      Result := -Infinity;    // x = -0: pole error → -Inf
      Exit;
    end;
    if ix.u > cRsqrtNegNanThresh.u then begin
      Result := x + x;        // x is a negative NaN: propagate
      Exit;
    end;
    if ix.u shr 63 <> 0 then begin
      // x < 0 (finite negative or -Inf): domain error → positive quiet NaN
      Result := cNaNDoublePos.f;
      Exit;
    end;
    if ix.u shl 12 = 0 then begin
      Result := 0.0;          // x = +Inf: rsqrt(+Inf) = 0
      Exit;
    end;
    Result := x + x;          // x is a positive NaN: propagate
    Exit;
  end else begin
    // Normal: 0x1p-1022 <= x < 2^1024
    if ix.u > UInt64($07FD000000000000) then
      // x > 2^1022: 4/x * (0.25*sqrt(x)) avoids spurious underflow in 1/x
      r := (4.0 / x) * (Double(0.25) * pcr_sqrt(x))
    else
      r := (1.0 / x) * pcr_sqrt(x);
  end;

  // Newton-Raphson refinement step
  rx  := r * x;
  // drx = fma(r, x, -rx) — the exact rounding error of r*x.
  // Non-AVX2: pcr_fma_pascal uses Veltkamp split(x), which overflows (NaN) when
  // K*x > DBL_MAX, i.e. biased-exponent(x) >= 2019 ($7E30...).
  // Fix: for those ~0.66% of inputs, scale: x_s=x*2^{-256}, r_s=r*2^{128},
  // rx_s=rx*2^{-128}; then drx = fma(r_s,x_s,-rx_s)*2^{128}  (exact, pure
  // exponent shift — the same bits, just a different biased exponent).
  {$IFDEF AVX2}
  drx := pcr_fma(r, x, -rx);
  {$ELSE}
  if ix.u >= UInt64($7E30000000000000) then begin
    // Large x: scale to keep K*x_s within double range
    drx := pcr_fma(r  * Double(3.402823669209385e38),    // r * 2^128
                   x  * Double(8.636168555094445e-78),    // x * 2^{-256}
                   -(rx * Double(2.9387358770557188e-39))) // rx * 2^{-128}
           * Double(3.402823669209385e38);                // * 2^128
  end else
    drx := pcr_fma(r, x, -rx);
  {$ENDIF}
  h   := pcr_fma(r, rx, -1.0) + r * drx;
  dr  := (r * Double(0.5)) * h;
  rf  := r - dr;
  dr  := dr - (r - rf);

  idr.f := dr;
  ir.f  := rf;

  // mid=0 and aidr out-of-range identify inputs needing the slow-path refine.
  aidr := (idr.u and UInt64($7FFFFFFFFFFFFFFF))
          - (ir.u and UInt64($7FF0000000000000))
          + UInt64($3FE0000000000000);
  mid  := (aidr - UInt64($3C90000000000000) + 16) shr 5;

  if (mid = 0) or (aidr < UInt64($39B0000000000000))
              or (aidr > UInt64($3C9FFFFFFFFFFF80)) then
    rf := RsqrtRefine(rf, x);

  Result := rf;
end;

// ---------------------------------------------------------------------------
// pcr_cbrt — correctly-rounded cube root (binary64).
// Ported from core-math/src/binary64/cbrt/cbrt.c by Alexei Sibidanov.
// ---------------------------------------------------------------------------
const
  // escale[it] = 2^(it/3), bit patterns
  cCbrtEsc0: Tb64u64 = (u: $3FF0000000000000);  // 1.0
  cCbrtEsc1: Tb64u64 = (u: $3FF428A2F98D728B);  // 2^(1/3)
  cCbrtEsc2: Tb64u64 = (u: $3FF965FEA53D6E3D);  // 2^(2/3)
  // polynomial c[4] approximating z^(1/3) on [1,2]
  cCbrtC0: Double = 0.55282341840164717;   // 0x1.1b0babccfef9cp-1
  cCbrtC1: Double = 0.58711429182669816;   // 0x1.2c9a3e94d1da5p-1
  cCbrtC2: Double = -0.16296967194987905;  // -0x1.4dc30b1a1ddbap-3
  cCbrtC3: Double = 0.023104964110781469;  // 0x1.7a8d3e4ec9b07p-6
  cCbrtU0: Double = 0.33333333333333331;   // 1/3  = 0x1.5555555555555p-2
  cCbrtU1: Double = 0.22222222222222221;   // 2/9  = 0x1.c71c71c71c71cp-3
  // off[0] for nearest half-ULP check
  cCbrtOff:   Double = 1.1102230246251565e-16;  // 2^-53
  cCbrtUlp52: Double = 2.2204460492503131e-16;  // 2^-52
  cCbrtThr75: Double = 2.6469779601696886e-23;  // 2^-75
  cCbrtThr98: Double = 3.1554436208840472e-30;  // 2^-98
  cCbrtThr60: Double = 8.6736173798840355e-19;  // 2^-60
  // hard cases for round-to-nearest
  cCbrtH0In:  Tb64u64 = (u: $4009B78223AA307C);
  cCbrtH0Out: Tb64u64 = (u: $3FF79D15D0E8D59C);
  cCbrtH1In:  Tb64u64 = (u: $401A202BFC89DDFF);
  cCbrtH1Out: Tb64u64 = (u: $3FFDE87AA837820F);
  // directed-rounding hard-case inputs (|zz|)
  cCbrtWIn: array[0..6] of Tb64u64 = (
    (u: $3FF3A9CCD7F022DB), (u: $3FF7845D2FAAC6FE),
    (u: $3FFD1EF81CBBBE71), (u: $4000A2014F62987C),
    (u: $400FE18A044A5501), (u: $401A6BB8C803147B),
    (u: $401AC8538A031CBD));
  // directed-rounding hard-case outputs (positive base)
  cCbrtWOut: array[0..6] of Tb64u64 = (
    (u: $3FF1236160BA9B93), (u: $3FF23115E657E49C),
    (u: $3FF388FB44CDCF5A), (u: $3FF46BCBF47DC1E8),
    (u: $3FF95DECFEC9C904), (u: $3FFE05335A6401DE),
    (u: $3FFE281D87098DE8));

function pcr_cbrt(x: Double): Double;
var
  cvt0, cvt1, cvt2, cvt3, cvt4, cvt5, tmp: Tb64u64;
  hx, mant, sign, ix, isc: UInt64;
  e, et, it: UInt32;
  nz, rm: Int32;
  flag: DWord;
  z, zz, r, rr, z2: Double;
  c0, c2, y, y2, y2l, y3, y3l, h, dy, y1: Double;
  ady, ady0, ady1, azz, off: Double;
  m0, m1: Int64;
  i: Integer;
begin
  cvt0.f := x;
  hx   := cvt0.u;
  mant := hx and UInt64($000FFFFFFFFFFFFF);
  sign := hx shr 63;
  e    := UInt32((hx shr 52) and $7FF);

  if ((e + 1) and $7FF) < 2 then begin
    ix := hx and UInt64($7FFFFFFFFFFFFFFF);
    if (e = $7FF) or (ix = 0) then begin Result := x + x; Exit; end;
    nz   := pcr_clzll(ix) - 11;
    mant := mant shl nz;
    mant := mant and UInt64($000FFFFFFFFFFFFF);
    Dec(e, UInt32(nz) - 1);
  end;

  flag := pcr_get_mxcsr;
  rm   := Int32(Ord(GetRoundMode));  // 0=nearest,1=down,2=up,3=zero

  e      := e + 3072;
  cvt1.u := mant or (UInt64($3FF) shl 52);
  cvt5.u := cvt1.u;
  et     := e div 3;
  it     := e mod 3;
  cvt5.u := cvt5.u + (UInt64(it) shl 52);
  cvt5.u := cvt5.u or (sign shl 63);
  zz     := cvt5.f;

  case it of
    0:    isc := cCbrtEsc0.u;
    1:    isc := cCbrtEsc1.u;
    else  isc := cCbrtEsc2.u;
  end;
  isc    := isc or (sign shl 63);
  cvt2.u := isc;
  z  := cvt1.f;
  r  := Double(1.0) / z;

  // rr = r * rsc[it*2 | sign], rsc = {1,-1, 0.5,-0.5, 0.25,-0.25}
  case it shl 1 or sign of
    0:    rr := r * Double(1.0);
    1:    rr := r * Double(-1.0);
    2:    rr := r * Double(0.5);
    3:    rr := r * Double(-0.5);
    4:    rr := r * Double(0.25);
    else  rr := r * Double(-0.25);
  end;

  // first Newton-Raphson: polynomial initial estimate then cubic correction
  z2 := z * z;
  c0 := cCbrtC0 + z * cCbrtC1;
  c2 := cCbrtC2 + z * cCbrtC3;
  y  := c0 + z2 * c2;
  y2 := y * y;
  h  := y2 * (y * r) - Double(1.0);
  y  := y - (h * y) * (cCbrtU0 - cCbrtU1 * h);
  y  := y * cvt2.f;

  // second Newton-Raphson with double-double error term
  y2  := y * y;
  y2l := pcr_fma(y, y, -y2);
  y3  := y2 * y;
  y3l := pcr_fma(y, y2, -y3) + y * y2l;
  h   := ((y3 - zz) + y3l) * rr;
  dy  := h * (y * cCbrtU0);
  y1  := y - dy;
  dy  := (y - y1) - dy;

  ady := Abs(dy);
  if rm = 0 then off := cCbrtOff else off := Double(0.0);
  ady0 := Abs(ady - off);
  ady1 := Abs(ady - (cCbrtUlp52 + off));

  if (ady0 < cCbrtThr75) or (ady1 < cCbrtThr75) then begin
    // extra refinement pass
    y2  := y1 * y1;
    y2l := pcr_fma(y1, y1, -y2);
    y3  := y2 * y1;
    y3l := pcr_fma(y1, y2, -y3) + y1 * y2l;
    h   := ((y3 - zz) + y3l) * rr;
    dy  := h * (y1 * cCbrtU0);
    y   := y1 - dy;
    dy  := (y1 - y) - dy;
    y1  := y;
    ady  := Abs(dy);
    ady0 := Abs(ady - off);
    ady1 := Abs(ady - (cCbrtUlp52 + off));
    if (ady0 < cCbrtThr98) or (ady1 < cCbrtThr98) then begin
      azz := Abs(zz);
      // round-to-nearest hard cases
      if azz = cCbrtH0In.f then begin
        tmp.u := cCbrtH0Out.u or (sign shl 63);
        y1 := tmp.f;
      end;
      if azz = cCbrtH1In.f then begin
        tmp.u := cCbrtH1Out.u or (sign shl 63);
        y1 := tmp.f;
      end;
      // directed-rounding hard cases
      if rm > 0 then
        for i := 0 to 6 do
          if azz = cCbrtWIn[i].f then begin
            // add 1 ULP before applying sign when rm+sign=2
            if rm + Int32(sign) = 2 then
              tmp.u := cCbrtWOut[i].u + UInt64(1)
            else
              tmp.u := cCbrtWOut[i].u;
            tmp.u := tmp.u or (sign shl 63);
            y1 := tmp.f;
          end;
    end;
  end;

  // scale y1 to the correct exponent: add (et-1365) to biased exponent
  cvt3.f := y1;
  cvt3.u := cvt3.u + (UInt64(et - UInt32(1365)) shl 52);

  // check if we are within 1 ULP of a half-ULP boundary
  m0 := Int64(cvt3.u shl 30);
  m1 := m0 shr 63;  // arithmetic shift: 0 if positive, -1 if negative
  if UInt64(m0 xor m1) <= (UInt64(1) shl 30) then begin
    cvt4.f := y1;
    cvt4.u := (cvt4.u + (UInt64(1) shl 15)) and UInt64($FFFFFFFFFFFF0000);
    if (Abs((cvt4.f - y1) - dy) < cCbrtThr60) or (Abs(zz) = Double(1.0)) then begin
      cvt3.u := (cvt3.u + (UInt64(1) shl 15)) and UInt64($FFFFFFFFFFFF0000);
      pcr_set_mxcsr(flag);
    end;
  end;
  Result := cvt3.f;
end;

// ---------------------------------------------------------------------------
// Stubs — delegate to C reference until each function is ported.
// Replace each stub body with the real Pascal port as phases 1-5 progress.
// ---------------------------------------------------------------------------
function  pcr_acos(x: Double): Double;    begin Result := cr_acos(x);    end;
function  pcr_acosh(x: Double): Double;   begin Result := cr_acosh(x);   end;
function  pcr_acospi(x: Double): Double;  begin Result := cr_acospi(x);  end;
function  pcr_asin(x: Double): Double;    begin Result := cr_asin(x);    end;
function  pcr_asinh(x: Double): Double;   begin Result := cr_asinh(x);   end;
function  pcr_asinpi(x: Double): Double;  begin Result := cr_asinpi(x);  end;
function  pcr_atan(x: Double): Double;    begin Result := cr_atan(x);    end;
function  pcr_atanh(x: Double): Double;   begin Result := cr_atanh(x);   end;
function  pcr_atanpi(x: Double): Double;  begin Result := cr_atanpi(x);  end;
// pcr_cbrt — ported above
function  pcr_cos(x: Double): Double;     begin Result := cr_cos(x);     end;
function  pcr_cosh(x: Double): Double;    begin Result := cr_cosh(x);    end;
function  pcr_cospi(x: Double): Double;   begin Result := cr_cospi(x);   end;
function  pcr_erf(x: Double): Double;     begin Result := cr_erf(x);     end;
function  pcr_erfc(x: Double): Double;    begin Result := cr_erfc(x);    end;
function  pcr_exp(x: Double): Double;     begin Result := cr_exp(x);     end;
function  pcr_exp10(x: Double): Double;   begin Result := cr_exp10(x);   end;
function  pcr_exp10m1(x: Double): Double; begin Result := cr_exp10m1(x); end;
function  pcr_exp2(x: Double): Double;    begin Result := cr_exp2(x);    end;
function  pcr_exp2m1(x: Double): Double;  begin Result := cr_exp2m1(x);  end;
function  pcr_expm1(x: Double): Double;   begin Result := cr_expm1(x);   end;
function  pcr_lgamma(x: Double): Double;  begin Result := cr_lgamma(x);  end;
function  pcr_log(x: Double): Double;     begin Result := cr_log(x);     end;
function  pcr_log10(x: Double): Double;   begin Result := cr_log10(x);   end;
function  pcr_log10p1(x: Double): Double; begin Result := cr_log10p1(x); end;
function  pcr_log1p(x: Double): Double;   begin Result := cr_log1p(x);   end;
function  pcr_log2(x: Double): Double;    begin Result := cr_log2(x);    end;
function  pcr_log2p1(x: Double): Double;  begin Result := cr_log2p1(x);  end;
function  pcr_sin(x: Double): Double;     begin Result := cr_sin(x);     end;
function  pcr_sinh(x: Double): Double;    begin Result := cr_sinh(x);    end;
function  pcr_sinpi(x: Double): Double;   begin Result := cr_sinpi(x);   end;
function  pcr_tan(x: Double): Double;     begin Result := cr_tan(x);     end;
function  pcr_tanh(x: Double): Double;    begin Result := cr_tanh(x);    end;
function  pcr_tanpi(x: Double): Double;   begin Result := cr_tanpi(x);   end;
function  pcr_tgamma(x: Double): Double;  begin Result := cr_tgamma(x);  end;
function  pcr_atan2(y, x: Double): Double;  begin Result := cr_atan2(y, x);   end;
function  pcr_atan2pi(y, x: Double): Double; begin Result := cr_atan2pi(y, x); end;
function  pcr_hypot(x, y: Double): Double;  begin Result := cr_hypot(x, y);   end;
function  pcr_pow(x, y: Double): Double;    begin Result := cr_pow(x, y);     end;
procedure pcr_sincos(x: Double; out s, c: Double);
begin
  cr_sincos(x, @s, @c);
end;

end.
