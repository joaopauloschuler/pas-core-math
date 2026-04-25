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
function pcr_atan(x: Double): Double;
function pcr_log2(x: Double): Double;
function pcr_acos(x: Double): Double;
function pcr_tanh(x: Double): Double;
function pcr_cospi(x: Double): Double;
function pcr_asin(x: Double): Double;
function pcr_exp(x: Double): Double;
function pcr_exp2(x: Double): Double;
function pcr_exp10(x: Double): Double;
function pcr_expm1(x: Double): Double;
function pcr_cos(x: Double): Double;
function pcr_sin(x: Double): Double;
function pcr_sinpi(x: Double): Double;
function pcr_atanpi(x: Double): Double;
function pcr_tan(x: Double): Double;
function pcr_cosh(x: Double): Double;
function pcr_sinh(x: Double): Double;
function pcr_log1p(x: Double): Double;
function pcr_log(x: Double): Double;
function pcr_log10(x: Double): Double;
function pcr_log10p1(x: Double): Double;
function pcr_log2p1(x: Double): Double;
function pcr_erf(x: Double): Double;
function pcr_erfc(x: Double): Double;
function pcr_exp2m1(x: Double): Double;
function pcr_exp10m1(x: Double): Double;
procedure pcr_sincos(x: Double; out s, c: Double);
function pcr_acospi(x: Double): Double;
function pcr_asinpi(x: Double): Double;

// ── Stub functions (delegate to C reference until ported) ────────────────────
// pcr_acos declared in ported section above
function  pcr_acosh(x: Double): Double;
// pcr_acospi declared in ported section above
// pcr_asin declared in ported section above
function  pcr_asinh(x: Double): Double;
// pcr_asinpi declared in ported section above
// pcr_atan declared in ported section above
function  pcr_atanh(x: Double): Double;
// pcr_atanpi declared in ported section above
// pcr_cbrt declared in ported section above
// pcr_cos declared in ported section above
// pcr_cosh declared in ported section above
// pcr_cospi declared in ported section above
// pcr_erf declared in ported section above
// pcr_erfc declared in ported section above
// pcr_exp declared in ported section above
// pcr_exp10 declared in ported section above
// pcr_exp10m1 declared in ported section above
// pcr_exp2 declared in ported section above
// pcr_exp2m1 declared in ported section above
// pcr_expm1 declared in ported section above
function  pcr_lgamma(x: Double): Double; inline;
// pcr_log declared in ported section above
// pcr_log10 declared in ported section above
// pcr_log10p1 declared in ported section above
// pcr_log1p declared in ported section above
// pcr_log2 declared in ported section above
// pcr_log2p1 declared in ported section above
// pcr_sin declared in ported section above
// pcr_sinh declared in ported section above
// pcr_sinpi declared in ported section above
// pcr_tan declared in ported section above
// pcr_tanh declared in ported section above
function  pcr_tanpi(x: Double): Double;
function  pcr_tgamma(x: Double): Double; inline;
function  pcr_atan2(y, x: Double): Double;
function  pcr_atan2pi(y, x: Double): Double;
function  pcr_hypot(x, y: Double): Double; inline;
function  pcr_pow(x, y: Double): Double; inline;
// pcr_sincos declared in ported section above

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
// pcr_atan — correctly-rounded arctangent (binary64).
// Ported from C CORE-MATH cr_atan (Alexei Sibidanov, 2023).
// ---------------------------------------------------------------------------

const
  cAtanAHi: array[0..128] of Tb64u64 = (
    (u:$0000000000000000),(u:$3F89224E047E368E),(u:$3F992346247A91F0),
    (u:$3FA2DBAAE9A05DB0),(u:$3FA927278A3B1162),(u:$3FAF7495EA3F3783),
    (u:$3FB2E239CCFF3831),(u:$3FB60B9F7597FDEC),(u:$3FB936BB8C5B2DA2),
    (u:$3FBC63CE377FC802),(u:$3FBF93183A8DB9E9),(u:$3FC1626D85A91E70),
    (u:$3FC2FCAC73A60640),(u:$3FC4986A74CF4E57),(u:$3FC635C990CE0D36),
    (u:$3FC7D4EC54FB5968),(u:$3FC975F5E0553158),(u:$3FCB1909EFD8B762),
    (u:$3FCCBE4CEB4B4CF2),(u:$3FCE65E3F27C9F2A),(u:$3FD007FA758626AE),
    (u:$3FD0DE53475F3B3C),(u:$3FD1B6103D3597E9),(u:$3FD28F459ECAD74D),
    (u:$3FD36A08355C63DC),(u:$3FD4466D542BAC92),(u:$3FD5248AE1701B17),
    (u:$3FD604775FBB27DF),(u:$3FD6E649F7D78649),(u:$3FD7CA1A832D0F84),
    (u:$3FD8B00196B3D022),(u:$3FD998188E816BF0),(u:$3FDA827999FCEF32),
    (u:$3FDB6F3FC8C61E5B),(u:$3FDC5E87185E67B6),(u:$3FDD506C82A2C800),
    (u:$3FDE450E0D273E7A),(u:$3FDF3C8AD985D9EE),(u:$3FE01B819B5A7CF7),
    (u:$3FE09A4C59BD0D4D),(u:$3FE11AB7190834EC),(u:$3FE19CD3FE8E405D),
    (u:$3FE220B5EF047825),(u:$3FE2A6709A74F289),(u:$3FE32E1889047FFD),
    (u:$3FE3B7C3289ED6F3),(u:$3FE44386DB9CE5DB),(u:$3FE4D17B087B265D),
    (u:$3FE561B82AB7F990),(u:$3FE5F457E4F4812E),(u:$3FE6897514751DB6),
    (u:$3FE7212BE621BE6D),(u:$3FE7BB99ED2990CF),(u:$3FE858DE3B716571),
    (u:$3FE8F9197BF85EEB),(u:$3FE99C6E0F634394),(u:$3FEA43002AE42850),
    (u:$3FEAECF5F9BA35A6),(u:$3FEB9A77C18C1AF2),(u:$3FEC4BB009E77983),
    (u:$3FED00CBC7384D2E),(u:$3FEDB9FA89953FCF),(u:$3FEE776EAFC91706),
    (u:$3FEF395D9F0E3C92),(u:$3FF0000000000000),(u:$3FF065C900AAF2D8),
    (u:$3FF0CE29D0883C99),(u:$3FF139447E6A86EE),(u:$3FF1A73D55278C4B),
    (u:$3FF2183B0C4573FF),(u:$3FF28C66FDAF8F09),(u:$3FF303ED61109E20),
    (u:$3FF37EFD8D87607E),(u:$3FF3FDCA42847507),(u:$3FF48089F8BF42CC),
    (u:$3FF507773C537EAD),(u:$3FF592D11142FA55),(u:$3FF622DB63C8ECC2),
    (u:$3FF6B7DF86265200),(u:$3FF7522CBDD428A8),(u:$3FF7F218E25A7461),
    (u:$3FF89801106CC709),(u:$3FF9444A7462122A),(u:$3FF9F7632FA9E871),
    (u:$3FFAB1C35D8A74EA),(u:$3FFB73EE3C3EF16A),(u:$3FFC3E738086BC0F),
    (u:$3FFD11F0DAE40609),(u:$3FFDEF13B73C1406),(u:$3FFED69B4153A45D),
    (u:$3FFFC95ABAD6CF4A),(u:$4000641E192CEAB3),(u:$4000EA21D716FBF7),
    (u:$40017749711A6679),(u:$40020C36C6A7F38E),(u:$4002A99F50FD4F4F),
    (u:$4003504F333F9DE6),(u:$4004012CE2586A17),(u:$4004BD3D87FE0650),
    (u:$400585AA4E1530FA),(u:$40065BC6CC825147),(u:$40074118E4B6A7C8),
    (u:$400837626D70FDB8),(u:$400940AD30ABC792),(u:$400A5F59E90600DD),
    (u:$400B9633283B6D14),(u:$400CE885653127E7),(u:$400E5A3DE972A377),
    (u:$400FF01305ECD8DC),(u:$4010D7DC7CFF4C9E),(u:$4011D0143E71565F),
    (u:$4012E4FF1626B949),(u:$40141BFEE2424771),(u:$40157BE4EAA5E11B),
    (u:$40170D751908C1B1),(u:$4018DC25C117782B),(u:$401AF73F4CA3310F),
    (u:$401D7398D15E70DB),(u:$4020372FB36B87E2),(u:$402208DBDAE055EF),
    (u:$40244E6C595AFDCC),(u:$4027398C57F3F1AD),(u:$402B1D03C03D2F7F),
    (u:$403046E9FE60A77E),(u:$40345AFFED201B55),(u:$403B267195B1FFAE),
    (u:$40445E2455E4AAA7),(u:$40545EED6854CE99),(u:$0000000000000000));

  cAtanALo: array[0..128] of Tb64u64 = (
    (u:$0000000000000000),(u:$3C1A3CA6C727C59D),(u:$3BF138B0EF96A186),
    (u:$3C436E7F8A3F5E42),(u:$BBFAC986EFB92662),(u:$3C406EC8011EE816),
    (u:$BC5858437D431332),(u:$BC3CEBD13EB7C513),(u:$BC5840CAC0D81DB5),
    (u:$3C5400B0FDAA109E),(u:$3C40E04E06C86E72),(u:$3C4F7AD829163CA7),
    (u:$BC52680735CE2CD8),(u:$BC690559690B42E4),(u:$3C591D29110B41AA),
    (u:$BC4EA90E27182780),(u:$BC2DC82AC14E3E1C),(u:$BC573A10FD13DAAF),
    (u:$BC63A7FFBEABDA0B),(u:$BC6DB6627A24D523),(u:$BC645F97DD3099F6),
    (u:$BC66293F68741816),(u:$BC6AB240D40633E9),(u:$BC2DE34D14E832E0),
    (u:$3C6AF540D9FB4926),(u:$3C6DA60FDBC82AC4),(u:$BC792A601170138A),
    (u:$BC67F1FCA1D5D15B),(u:$BC64E223EA716C7B),(u:$3C7B24C824AC51FC),
    (u:$3C64314CD132BA43),(u:$BC711F1E0817879A),(u:$BC6C3DEA4DBAD538),
    (u:$3C660D1B780EE3EB),(u:$BC4AB5EDB7DFA545),(u:$BC68E1437048B5BD),
    (u:$BC706951C97B050F),(u:$BC414AF9522AB518),(u:$BC7ABA0D7D97D1F2),
    (u:$3C4095BC4EBC2C42),(u:$3C8798826FA27774),(u:$3C8008F6258FC98F),
    (u:$BC5462AF7CEB7DE6),(u:$BC71184DFD78B472),(u:$3C79141876DC40C5),
    (u:$3C8481C20189726C),(u:$3C82E851BD025441),(u:$3C713ADA9B8BC419),
    (u:$BC805B4C3C4CBEE8),(u:$BC85619249BD96F1),(u:$BC6B0A0FBCAFC671),
    (u:$BC819FF2DC66DA45),(u:$3C81320449592D92),(u:$BC81FDDCD2F3DA8E),
    (u:$3C6D44A42E35CC97),(u:$BC7585A178B4A18D),(u:$3C6F95A531B3A970),
    (u:$BC396C2D43CA3392),(u:$BC6A5BED94B05DEF),(u:$3C454509D2BFF511),
    (u:$BC6B4C867CEF300C),(u:$BC1DDFAC663D6BC6),(u:$BC7A510683FF7CB6),
    (u:$3C44FDCD8E4E8710),(u:$0000000000000000),(u:$BC8DEEC7FC9042AD),
    (u:$BC8395AE45E0657D),(u:$3C8332CF301A97F3),(u:$BC86CC8C4B78213B),
    (u:$3C870A90841DA57A),(u:$BC6BA39BAD450EE0),(u:$BC88692946D9F93C),
    (u:$3C63B711BF765B58),(u:$3C7C21387985B081),(u:$BC87DDB19D3D0EFC),
    (u:$BC7F5E354CF971F3),(u:$BC700F0AD675330D),(u:$BC82C93F50AB2C0E),
    (u:$3C7BEC391ADC37D5),(u:$BC69686DDC9FFCF5),(u:$BC78D16529514246),
    (u:$BC8092F51E9C2803),(u:$BC807C06755404C4),(u:$3C802E0D43ABC92B),
    (u:$3C5D0184E48AF6F7),(u:$3C773BE957380BC2),(u:$BC702B6E26C84462),
    (u:$3C525C4F3FFA6E1F),(u:$BC5E302DB3C6823F),(u:$3C73207830326C0E),
    (u:$BC66308CEE7927BF),(u:$BC70147EBF0DF4C5),(u:$BC7168533CC41D8B),
    (u:$BC652A0B0333E9C5),(u:$3C68659EECE35395),(u:$3C820FCAD18CB36F),
    (u:$BC752AFDBD5A8C74),(u:$BC79747A792907D7),(u:$3C790C59393B52C8),
    (u:$3C7AF6934F13A3A8),(u:$BC48534DCAB5AD3E),(u:$BC7555AA8BFCA9A1),
    (u:$BC556B3FEE9CA72B),(u:$3C54B3FDD4FDC06C),(u:$3C6285D367C55DDC),
    (u:$BC48712976F17A16),(u:$BC3ABE8AB65D49FC),(u:$3C5CD9BE81AD764B),
    (u:$3C4742C2922656FA),(u:$BC77C842978BEE09),(u:$3C67BC7DEA7C3C03),
    (u:$3C4AEFBE25B404E9),(u:$BC34BCFAAA95CB2C),(u:$3C50FE741E4EC679),
    (u:$3C5FE74A5B0EC709),(u:$3C50CA1C19F710EF),(u:$3C52867B40BA77D6),
    (u:$3C60FD4E0D4B1547),(u:$3C5C16C9ECC1621D),(u:$3C56B81A36E75E8C),
    (u:$BC57C22045771848),(u:$3C5970503BE105C0),(u:$BC3F299D010AEAD2),
    (u:$3C5D2B61DEFF33EC),(u:$3BF0E84D9567203A),(u:$BBFAD44B44B92653),
    (u:$BC3296D577B5E21D),(u:$3C02DB53886013CA),(u:$0000000000000000));

  cAtanC0: array[0..30] of UInt16 = (
      419,  500,  582,  745,  908, 1234, 1559, 2210,
     2860, 4156, 5444, 7989,10476,15224,19601,27105,
    33036,41266,46469,52375,55587,58906,60612,62325,
    63192,64056,64491,64923,65141,65358,65467);
  cAtanC1: array[0..30] of UInt16 = (
       81,   81,  163,  163,  326,  326,  651,  650,
     1299, 1293, 2569, 2520, 4917, 4576, 8341, 6648,
    10210, 6292, 7926, 4038, 4591, 2172, 2390, 1107,
     1207,  556,  605,  278,  303,  139,  151);
  cAtanC2: array[0..30] of UInt16 = (
        0,    0,    0,    0,    0,    0,    0,    1,
        3,    4,   24,   32,  168,  200,  838,  731,
     1998, 1117, 2048,  849, 1291,  479,  688,  247,
      349,  124,  175,   62,   88,   31,   44);

  // Fast-path polynomial: ch[0..3]
  cAtanCh_0: Double = 1.0;
  cAtanCh_1: Tb64u64 = (u:$BFD555555555552B);  // -0x1.555555555552bp-2
  cAtanCh_2: Tb64u64 = (u:$3FC9999999069C20);  //  0x1.9999999069c2p-3
  cAtanCh_3: Tb64u64 = (u:$BFC248D2C8444AC6);  // -0x1.248d2c8444ac6p-3

  // Small-x polynomial: ch2[0..3]
  cAtanCh2_0: Tb64u64 = (u:$BFD5555555555555);  // -0x1.5555555555555p-2
  cAtanCh2_1: Tb64u64 = (u:$3FC99999999998C1);  //  0x1.99999999998c1p-3
  cAtanCh2_2: Tb64u64 = (u:$BFC249249176AEC0);  // -0x1.249249176aecp-3
  cAtanCh2_3: Tb64u64 = (u:$3FBC711FD121AE80);  //  0x1.c711fd121ae8p-4

  // Fast-path Sterbenz error bounds
  cAtanFmaUb: Tb64u64 = (u:$3CD2000000000000);  //  0x4.8p-52
  cAtanFmaLb: Tb64u64 = (u:$3CC4000000000000);  //  0x2.8p-52
  // fast-path error factor: e = h * 0x3.fp-52
  cAtanEFactor: Tb64u64 = (u:$3CCF800000000000);  //  0x3.fp-52

  // pi/2 high and low parts
  cAtanPiHalfH: Tb64u64 = (u:$3FF921FB54442D18);  //  0x1.921fb54442d18p+0
  cAtanPiHalfL: Tb64u64 = (u:$3C91A62633145C07);  //  0x1.1a62633145c07p-54

  // id scaling: ah = IdHi*id, al = IdLo*id (fast path), IdLo2 (refine), at = IdLo3*id
  cAtanIdHi:   Tb64u64 = (u:$3F8921FB54442D00);  //  0x1.921fb54442dp-7
  cAtanIdLo:   Tb64u64 = (u:$3C88469898CC5170);  //  0x1.8469898cc517p-55
  cAtanIdLo2:  Tb64u64 = (u:$3C88469898CC5180);  //  0x1.8469898cc518p-55 (refine)
  cAtanIdLo3:  Tb64u64 = (u:$B97FC8F8CBB5BF80);  // -0x1.fc8f8cbb5bf8p-104

  // phi scale for refine2 index: |a| * phi_scale + 256.5
  cAtanPhiScale: Tb64u64 = (u:$40545F306DC9C883);  //  0x1.45f306dc9c883p6

  // Refine2 slow-path polynomial: ch[][2] (3 pairs) and cl[] (4 scalars)
  cAtanRefCh0H: Tb64u64 = (u:$BFD5555555555555);  // -0x1.5555555555555p-2
  cAtanRefCh0L: Tb64u64 = (u:$BC75555555555555);  // -0x1.5555555555555p-56
  cAtanRefCh1H: Tb64u64 = (u:$3FC999999999999A);  //  0x1.999999999999ap-3
  cAtanRefCh1L: Tb64u64 = (u:$BC6999999999BCB8);  // -0x1.999999999bcb8p-57
  cAtanRefCh2H: Tb64u64 = (u:$BFC2492492492492);  // -0x1.2492492492492p-3
  cAtanRefCh2L: Tb64u64 = (u:$BC6249242093C016);  // -0x1.249242093c016p-57
  cAtanRefCl0:  Tb64u64 = (u:$3FBC71C71C71C71C);  //  0x1.c71c71c71c71cp-4
  cAtanRefCl1:  Tb64u64 = (u:$BFB745D1745D1265);  // -0x1.745d1745d1265p-4
  cAtanRefCl2:  Tb64u64 = (u:$3FB3B13B115BCBC4);  //  0x1.3b13b115bcbc4p-4
  cAtanRefCl3:  Tb64u64 = (u:$BFB1107C41AD3253);  // -0x1.1107c41ad3253p-4

  // Hard cases db[0..11]: (|x| input, |result| output, correction)
  cAtanDbIn:  array[0..11] of Tb64u64 = (
    (u:$3F80DC89A3B55010),(u:$3F7E3FB41D2D2260),(u:$3FE7BA49F739829F),
    (u:$3FCA933FE176B375),(u:$3F7BB04A79820063),(u:$3F7CD30A9499618B),
    (u:$3F8F44AA37B8E66B),(u:$3F7FD2AC95E57EF9),(u:$3F96419079BBF601),
    (u:$3FE7BA49F739829F),(u:$3FCD768804487B07),(u:$3F7BB04A79820063));
  cAtanDbOut: array[0..11] of Tb64u64 = (
    (u:$3F80DC70AC228717),(u:$3F7E3F9013A852F8),(u:$3FE46AC372243536),
    (u:$3FCA33F32AC5CEB5),(u:$3F7BB02ED5C5E956),(u:$3F7CD2EB65F92A46),
    (u:$3F8F440B04187C87),(u:$3F7FD2829FEBC03A),(u:$3F9640AADE8F5427),
    (u:$3FE46AC372243536),(u:$3FCCF5676F373EC1),(u:$3F7BB02ED5C5E956));
  cAtanDbCor: array[0..11] of Tb64u64 = (
    (u:$3C20000000000000),(u:$3C10000000000000),(u:$3920000000000000),
    (u:$B8F0000000000000),(u:$B8C0000000000000),(u:$B8F0000000000000),
    (u:$B8F0000000000000),(u:$B8F0000000000000),(u:$38D0000000000000),
    (u:$3910000000000000),(u:$B910000000000000),(u:$B8C0000000000000));

function AtanAddDD(xh, xl, ch, cl: Double; out l: Double): Double; inline;
var s, d: Double;
begin
  s := xh + ch;
  d := s - xh;
  l := ((ch - d) + (xh + (d - s))) + (xl + cl);
  Result := s;
end;

function AtanRefine2(x, a: Double): Double;
var
  phi: Tb64u64;
  t_x: Tb64u64;
  ip: Int64;
  ta_u: Tb64u64;
  ta, zta, ztal, zmta, v_r, d_r, ev, r_v, rl, h, hl: Double;
  h2, h2l, h4, h3, h3l: Double;
  fl_r, f_r: Double;
  chp, clp, th, tl: Double;
  ah, al_r, at_r, df: Double;
  id_u: Tb64u64;
  v0, v1, v2: Double;
  ax: Double;
  t0, t1, w_r: Tb64u64;
  j: Integer;
begin
  phi.f := Abs(a) * cAtanPhiScale.f + Double(256.5);
  ip := Int64((phi.u shr 44) and UInt64($FF));

  t_x.f := x;

  if ip = 128 then
  begin
    h  := -Double(1.0) / x;
    hl := pcr_fma(h, x, Double(1.0)) * h;
  end
  else
  begin
    ta_u.u := cAtanAHi[ip].u or (t_x.u and UInt64($8000000000000000));
    ta   := ta_u.f;
    zta  := x * ta;
    ztal := pcr_fma(x, ta, -zta);
    zmta := x - ta;
    v_r  := Double(1.0) + zta;
    d_r  := Double(1.0) - v_r;
    ev   := (d_r + zta) - ((d_r + v_r) - Double(1.0)) + ztal;
    r_v  := Double(1.0) / v_r;
    rl   := (pcr_fma(r_v, -v_r, Double(1.0)) - ev * r_v) * r_v;
    h    := r_v * zmta;
    hl   := pcr_fma(r_v, zmta, -h) + rl * zmta;
  end;

  h2 := pcr_muldd(h, hl, h, hl, h2l);
  h4 := h2 * h2;
  h3 := pcr_muldd(h, hl, h2, h2l, h3l);

  fl_r := h2 * ((cAtanRefCl0.f + h2 * cAtanRefCl1.f) +
                 h4 * (cAtanRefCl2.f + h2 * cAtanRefCl3.f));

  // polydd(h2, h2l, 3, ch_ref, fl_r): inline expansion (n=3, ch pairs reversed)
  chp := cAtanRefCh2H.f + fl_r;
  clp := (cAtanRefCh2H.f - chp) + fl_r + cAtanRefCh2L.f;
  chp := pcr_muldd(h2, h2l, chp, clp, clp);
  th  := chp + cAtanRefCh1H.f;  tl := (cAtanRefCh1H.f - th) + chp;
  chp := th;  clp := clp + tl + cAtanRefCh1L.f;
  chp := pcr_muldd(h2, h2l, chp, clp, clp);
  th  := chp + cAtanRefCh0H.f;  tl := (cAtanRefCh0H.f - th) + chp;
  chp := th;  clp := clp + tl + cAtanRefCh0L.f;
  fl_r := clp;
  f_r  := chp;

  f_r := pcr_muldd(h3, h3l, f_r, fl_r, fl_r);

  if ip = 0 then
  begin
    ah   := h;
    al_r := f_r;
    at_r := fl_r;
  end
  else
  begin
    if ip < 128 then
    begin
      df := cAtanALo[ip].f;
      if x < Double(0.0) then df := -df;
    end
    else
      df := Double(0.0);

    id_u.f := Double(ip);
    id_u.u := id_u.u or (t_x.u and UInt64($8000000000000000));

    ah   := cAtanIdHi.f  * id_u.f;
    al_r := cAtanIdLo2.f * id_u.f;
    at_r := cAtanIdLo3.f * id_u.f;
    al_r := AtanAddDD(al_r, at_r, df,   Double(0.0), at_r);
    al_r := AtanAddDD(al_r, at_r, h,    hl,          at_r);
    al_r := AtanAddDD(al_r, at_r, f_r,  fl_r,        at_r);
  end;

  pcr_fasttwosum(v0, v2, ah, al_r);
  pcr_fasttwosum(v1, v2, v2, at_r);

  ax := Abs(x);
  t0.f := v0;
  t1.f := v1;

  if (((t1.u + UInt64(1)) and UInt64($000FFFFFFFFFFFFF)) <= UInt64(2)) or
     (((t0.u shr 52) and UInt64($7FF)) - ((t1.u shr 52) and UInt64($7FF)) > UInt64(103)) then
  begin
    for j := 0 to 11 do
    begin
      if ax = cAtanDbIn[j].f then
      begin
        if x >= Double(0.0) then
          Result := cAtanDbOut[j].f + cAtanDbCor[j].f
        else
          Result := -(cAtanDbOut[j].f + cAtanDbCor[j].f);
        Exit;
      end;
    end;
    if (t1.u and UInt64($000FFFFFFFFFFFFF)) = UInt64(0) then
    begin
      w_r.f := v2;
      if ((w_r.u xor t1.u) shr 63) <> 0 then
        Dec(t1.u)
      else
        Inc(t1.u);
      v1 := t1.f;
    end;
  end;

  Result := v1 + v0;
end;

function pcr_atan(x: Double): Double;
var
  t: Tb64u64;
  at: UInt64;
  i_row: Int64;
  ii: Int64;
  u, ut, ut2: UInt64;
  ta_u, id_u: Tb64u64;
  h, ah, al_v, f_v: Double;
  h2, h4: Double;
  e, ub, lb: Double;
  x2, x3, x4, fs: Double;
begin
  t.f  := x;
  at   := t.u and UInt64($7FFFFFFFFFFFFFFF);
  i_row := Int64(at shr 51) - Int64(2030);

  if at < UInt64($3F7B21C475E6362A) then
  begin
    if at = UInt64(0) then begin Result := x; Exit; end;
    if at < UInt64($3E40000000000000) then
    begin
      Result := x;
      Exit;
    end;
    x2 := x * x;
    x3 := x * x2;
    x4 := x2 * x2;
    fs := x3 * ((cAtanCh2_0.f + x2 * cAtanCh2_1.f) +
                 x4 * (cAtanCh2_2.f + x2 * cAtanCh2_3.f));
    ub := (fs + fs * cAtanFmaUb.f) + x;
    lb := (fs - fs * cAtanFmaLb.f) + x;
    if ub = lb then begin Result := ub; Exit; end;
    Result := AtanRefine2(x, ub);
    Exit;
  end;

  if at > UInt64($4062DED8E34A9035) then
  begin
    ah   := cAtanPiHalfH.f;
    al_v := cAtanPiHalfL.f;
    if x < Double(0.0) then begin ah := -ah; al_v := -al_v; end;
    if at >= UInt64($434D02967C31CDB5) then
    begin
      if at > UInt64($7FF0000000000000) then begin Result := x + x; Exit; end;
      Result := ah + al_v;
      Exit;
    end;
    h := -Double(1.0) / x;
  end
  else
  begin
    u   := t.u and UInt64($0007FFFFFFFFFFFF);
    ut  := u shr 35;
    ut2 := (ut * ut) shr 16;
    ii  := Int64((UInt64(cAtanC0[Integer(i_row)]) shl 16)
               + ut * UInt64(cAtanC1[Integer(i_row)])
               - ut2 * UInt64(cAtanC2[Integer(i_row)])) shr 25;

    ta_u.u := cAtanAHi[ii].u or (t.u and UInt64($8000000000000000));
    id_u.f := Double(ii);
    id_u.u := id_u.u or (t.u and UInt64($8000000000000000));

    al_v := cAtanALo[ii].f;
    if x < Double(0.0) then al_v := -al_v;
    al_v := al_v + cAtanIdLo.f * id_u.f;

    h    := (x - ta_u.f) / (Double(1.0) + x * ta_u.f);
    ah   := cAtanIdHi.f * id_u.f;
  end;

  h2   := h * h;
  h4   := h2 * h2;
  f_v  := (cAtanCh_0 + h2 * cAtanCh_1.f) + h4 * (cAtanCh_2.f + h2 * cAtanCh_3.f);
  al_v := pcr_fma(h, f_v, al_v);
  e    := h * cAtanEFactor.f;
  ub   := (al_v + e) + ah;
  lb   := (al_v - e) + ah;
  if ub = lb then begin Result := ub; Exit; end;
  Result := AtanRefine2(x, ub);
end;

// ---------------------------------------------------------------------------
// pcr_log2 — correctly-rounded base-2 logarithm (binary64).
// Ported from core-math/src/binary64/log2/log2.c by Alexei Sibidanov.
// ---------------------------------------------------------------------------

const
  // B[] = {c0 (ushort), c1 (short signed)}, 32 entries
  cLog2B_c0: array[0..31] of UInt16 = (
      301,  7189, 13383, 18923, 23845, 28184, 31969, 35231,
    37996, 40288, 42129, 43542, 44546, 45160, 45399, 45281,
    44821, 44032, 42929, 41522, 39825, 37848, 35602, 33097,
    30341, 27345, 24115, 20661, 16989, 13107,  9022,  4740);
  cLog2B_c1: array[0..31] of Int16 = (
    27565, 24786, 22167, 19696, 17361, 15150, 13054, 11064,
     9173,  7372,  5657,  4020,  2457,   962,  -468, -1838,
    -3151, -4412, -5622, -6786, -7905, -8982,-10020,-11020,
   -11985,-12916,-13816,-14685,-15526,-16339,-17126,-17889);

  cLog2R1: array[0..32] of Tb64u64 = (
    (u:$3FF7154800000000),(u:$3FF696AF49200000),
    (u:$3FF61AB3FB680000),(u:$3FF5A18441680000),
    (u:$3FF52ADADB480000),(u:$3FF4B6B7C9080000),(u:$3FF4451B0AA80000),(u:$3FF3D5ED8AE00000),
    (u:$3FF3691834680000),(u:$3FF2FE9B07400000),(u:$3FF2967603680000),(u:$3FF2307AFE500000),
    (u:$3FF1CCA9F7F80000),(u:$3FF16B02F0600000),(u:$3FF10B85E7880000),(u:$3FF0AE04B2E00000),
    (u:$3FF0527F52680000),(u:$3FEFF1EB8C400000),(u:$3FEF42A1F1800000),(u:$3FEE9721D4900000),
    (u:$3FEDEF6B35700000),(u:$3FED4B21BF000000),(u:$3FECAA739BD00000),(u:$3FEC0D32A1500000),
    (u:$3FEB735ECF800000),(u:$3FEADCC9FBD00000),(u:$3FEA497426400000),(u:$3FE9B92F24400000),
    (u:$3FE92C2920600000),(u:$3FE8A205C5800000),(u:$3FE81AF33E300000),(u:$3FE796C35FE00000),
    (u:$3FE7154800000000));

  cLog2R2: array[0..32] of Tb64u64 = (
    (u:$3FF0000000000000),(u:$3FEFFA7000000000),(u:$3FEFF4F000000000),(u:$3FEFEF6000000000),
    (u:$3FEFE9E000000000),(u:$3FEFE45000000000),(u:$3FEFDED000000000),(u:$3FEFD94000000000),
    (u:$3FEFD3C000000000),(u:$3FEFCE4000000000),(u:$3FEFC8C000000000),(u:$3FEFC34000000000),
    (u:$3FEFBDC000000000),(u:$3FEFB84000000000),(u:$3FEFB2C000000000),(u:$3FEFAD4000000000),
    (u:$3FEFA7C000000000),(u:$3FEFA24000000000),(u:$3FEF9CD000000000),(u:$3FEF975000000000),
    (u:$3FEF91E000000000),(u:$3FEF8C6000000000),(u:$3FEF86F000000000),(u:$3FEF817000000000),
    (u:$3FEF7C0000000000),(u:$3FEF769000000000),(u:$3FEF711000000000),(u:$3FEF6BA000000000),
    (u:$3FEF663000000000),(u:$3FEF60C000000000),(u:$3FEF5B5000000000),(u:$3FEF55E000000000),
    (u:$3FEF507000000000));

  cLog2L1Lo: array[0..32] of Tb64u64 = (
    (u:$0000000000000000),(u:$3E51435EDC775B51),(u:$3E5C8F1CBF9E4073),(u:$BE57BF30FA53957B),
    (u:$3E4674D30B6276ED),(u:$3E399FCF0D796ACE),(u:$BE55FFD8B92706D2),(u:$3E5BEF90BC5A116D),
    (u:$3E5FCA73B3D53F0D),(u:$BE441024E560E04E),(u:$3E5484024FAD8461),(u:$3E52D9AB90BA7694),
    (u:$BE52C998EA30BA7B),(u:$BE475FC8682F918E),(u:$BE5F02F268A85FB8),(u:$BE450030EA7FAE4B),
    (u:$BE3BFFB8DA5B849D),(u:$3E3F71993FF95475),(u:$3E288028E67F78FA),(u:$3E5F915F5A0B4E89),
    (u:$3E1C2FE288F968F8),(u:$3E57375A75AE0837),(u:$3E291E48BE920323),(u:$3E2EE7BC0D39A3DB),
    (u:$BE2899E2AC5F778C),(u:$3E5F1F20176130A7),(u:$3E403FC59D34A4F3),(u:$BDF68722010E4653),
    (u:$BE558783D505A6EC),(u:$3E4B212AB9F8D51D),(u:$BE5B3ECD767BE776),(u:$3E4B1AD41F07FC10),
    (u:$0000000000000000));

  cLog2L1Hi: array[0..32] of Tb64u64 = (
    (u:$0000000000000000),(u:$3F9FFE3800000000),(u:$3FB000BC00000000),(u:$3FB7FF9400000000),
    (u:$3FBFFFA600000000),(u:$3FC4000580000000),(u:$3FC7FFE500000000),(u:$3FCBFFC300000000),
    (u:$3FCFFFFD00000000),(u:$3FD2000D00000000),(u:$3FD3FFCC00000000),(u:$3FD5FFE180000000),
    (u:$3FD80013C0000000),(u:$3FDA0024C0000000),(u:$3FDBFFD300000000),(u:$3FDDFFD7C0000000),
    (u:$3FDFFFF980000000),(u:$3FE0FFFD40000000),(u:$3FE20010C0000000),(u:$3FE30018A0000000),
    (u:$3FE3FFF480000000),(u:$3FE5001300000000),(u:$3FE6000EA0000000),(u:$3FE7001120000000),
    (u:$3FE7FFFB40000000),(u:$3FE8FFFAE0000000),(u:$3FE9FFF1E0000000),(u:$3FEB0012A0000000),
    (u:$3FEBFFEB80000000),(u:$3FED0004E0000000),(u:$3FEDFFEB20000000),(u:$3FEEFFD520000000),
    (u:$3FF0000000000000));

  cLog2L2Lo: array[0..32] of Tb64u64 = (
    (u:$0000000000000000),(u:$BE4E2B19F9C7B840),(u:$BE4B1D68137631FE),(u:$BE4AA92227513FC3),
    (u:$BE452918E3AB6F5E),(u:$3E41A0B9B9010A9C),(u:$3E509B87F57867EC),(u:$BE33D6B70C673BE6),
    (u:$3E476D340A6780AB),(u:$3E5E4181E37D9E05),(u:$BE48C181F042B901),(u:$3E49890AE7761D66),
    (u:$BE555977AE613D5F),(u:$3E4C75F49ACF5E56),(u:$BE3E1FC84D0D42BC),(u:$BE55AB8A182ED279),
    (u:$BE5B389A7D8A21D9),(u:$BE57376F1A891FA5),(u:$BE30915581B87A8A),(u:$3E4B7F2EAA894FCA),
    (u:$BE5F8C41BDD38C23),(u:$3E304FFB1F3C8215),(u:$3E5A398CCF2B137D),(u:$3E5373B5BCE07F39),
    (u:$BE4E42040FD9D454),(u:$3E4EA00EB5770526),(u:$BE53E2002BB54BB8),(u:$3E5797EC835DB8D3),
    (u:$3E3BDCD45FEDD285),(u:$BE333B8F09D1A210),(u:$BE458937A8EEA36F),(u:$BE43E76203FF54CD),
    (u:$BE16560F19FC3F41));

  cLog2L2Hi: array[0..32] of Tb64u64 = (
    (u:$0000000000000000),(u:$3F500E4000000000),(u:$3F5FF10000000000),(u:$3F68026000000000),
    (u:$3F6FF68000000000),(u:$3F74019000000000),(u:$3F77FD0000000000),(u:$3F7C04C000000000),
    (u:$3F8000C800000000),(u:$3F81FF8800000000),(u:$3F83FEA800000000),(u:$3F85FE1800000000),
    (u:$3F87FDE800000000),(u:$3F89FE0800000000),(u:$3F8BFE8800000000),(u:$3F8DFF6000000000),
    (u:$3F90004800000000),(u:$3F91010C00000000),(u:$3F91FF1000000000),(u:$3F93002C00000000),
    (u:$3F93FE8C00000000),(u:$3F95000000000000),(u:$3F95FEB400000000),(u:$3F97008400000000),
    (u:$3F97FF9400000000),(u:$3F98FECC00000000),(u:$3F9A012400000000),(u:$3F9B00B400000000),
    (u:$3F9C007400000000),(u:$3F9D006000000000),(u:$3F9E007800000000),(u:$3F9F00BC00000000),
    (u:$3FA0009600000000));

  cLog2C0: Tb64u64 = (u:$BFD62E41D56C6400);
  cLog2C1: Tb64u64 = (u:$3FC47FD2632D2D32);
  cLog2C2: Tb64u64 = (u:$BFB5504497831BA7);
  cLog2C3: Tb64u64 = (u:$3FA7A3314C5BEF3C);
  cLog2L2H: Tb64u64 = (u:$3FF7154800000000);
  cLog2L2L: Tb64u64 = (u:$BE9AD47A2F472159);

  // Refine-path tables
  cLog2T_t1: array[0..16] of Tb64u64 = (
    (u:$3FF0000000000000), (u:$3FEEA4AFA0000000), (u:$3FED5818E0000000), (u:$3FEC199BE0000000),
    (u:$3FEAE89F98000000), (u:$3FE9C49180000000), (u:$3FE8ACE540000000), (u:$3FE7A11470000000),
    (u:$3FE6A09E68000000), (u:$3FE5AB07E0000000), (u:$3FE4BFDAD8000000), (u:$3FE3DEA650000000),
    (u:$3FE306FE08000000), (u:$3FE2387A70000000), (u:$3FE172B840000000), (u:$3FE0B55870000000),
    (u:$3FE0000000000000));
  cLog2T_t2: array[0..15] of Tb64u64 = (
    (u:$3FF0000000000000), (u:$3FEFE9D968000000), (u:$3FEFD3C228000000), (u:$3FEFBDBA38000000),
    (u:$3FEFA7C180000000), (u:$3FEF91D800000000), (u:$3FEF7BFDB0000000), (u:$3FEF663278000000),
    (u:$3FEF507658000000), (u:$3FEF3AC948000000), (u:$3FEF252B38000000), (u:$3FEF0F9C20000000),
    (u:$3FEEFA1BF0000000), (u:$3FEEE4AAA0000000), (u:$3FEECF4830000000), (u:$3FEEB9F488000000));
  cLog2T_t3: array[0..15] of Tb64u64 = (
    (u:$3FF0000000000000), (u:$3FEFFE9D20000000), (u:$3FEFFD3A58000000), (u:$3FEFFBD798000000),
    (u:$3FEFFA74E8000000), (u:$3FEFF91248000000), (u:$3FEFF7AFB8000000), (u:$3FEFF64D38000000),
    (u:$3FEFF4EAC8000000), (u:$3FEFF38868000000), (u:$3FEFF22618000000), (u:$3FEFF0C3D0000000),
    (u:$3FEFEF61A0000000), (u:$3FEFEDFF78000000), (u:$3FEFEC9D68000000), (u:$3FEFEB3B60000000));
  cLog2T_t4: array[0..15] of Tb64u64 = (
    (u:$3FF0000000000000), (u:$3FEFFFE9D0000000), (u:$3FEFFFD3A0000000), (u:$3FEFFFBD78000000),
    (u:$3FEFFFA748000000), (u:$3FEFFF9118000000), (u:$3FEFFF7AE8000000), (u:$3FEFFF64C0000000),
    (u:$3FEFFF4E90000000), (u:$3FEFFF3860000000), (u:$3FEFFF2238000000), (u:$3FEFFF0C08000000),
    (u:$3FEFFEF5D8000000), (u:$3FEFFEDFA8000000), (u:$3FEFFEC980000000), (u:$3FEFFEB350000000));

  cLog2LL0: array[0..16, 0..2] of Tb64u64 = (
    ((u:$0000000000000000), (u:$0000000000000000), (u:$0000000000000000)),
    ((u:$3FB000001FDA0000), (u:$3D55ED58A7FF2C40), (u:$39F512ACBB219717)),
    ((u:$3FBFFFFFDA070000), (u:$BD0B38AE1C2D5400), (u:$BA210F11908C5C8D)),
    ((u:$3FC7FFFFEFB50000), (u:$3D4785D91BB08320), (u:$3A2F5DED30CA48C0)),
    ((u:$3FD0000004A60000), (u:$3D5380ABF8FE7EB0), (u:$BA2EBD3F567DF886)),
    ((u:$3FD4000009760000), (u:$3D4D321385BF0A10), (u:$BA1C95A3D59B6C33)),
    ((u:$3FD80000081B4000), (u:$3D5B6E8F25FF2610), (u:$BA24072D39D44270)),
    ((u:$3FDC00000F4F4000), (u:$3D4A5FED334580A0), (u:$B9F6F8FDABADD77A)),
    ((u:$3FDFFFFFF9DE0000), (u:$BD4F85A54FB4A600), (u:$B9F50083897EE638)),
    ((u:$3FE1FFFFFA35C000), (u:$3D5795DB9D185140), (u:$BA1594F07D4D693F)),
    ((u:$3FE3FFFFF9CBC000), (u:$3D32FB09A99477A0), (u:$3A1B2C6248382152)),
    ((u:$3FE5FFFFF6DF2000), (u:$BD53ECF24077C3A0), (u:$BA2382C6B4E2DAA6)),
    ((u:$3FE800000552E000), (u:$3D47044806F72380), (u:$BA26DA3EB126C999)),
    ((u:$3FE9FFFFFC182000), (u:$3D45CF031D5D4CF0), (u:$3A2302C44DA79F4B)),
    ((u:$3FEBFFFFF6B62000), (u:$3D548D11D13372A8), (u:$3A2F045032543C0D)),
    ((u:$3FEDFFFFF7A42000), (u:$BD40CBF3E62FA860), (u:$3A21CFAF84227211)),
    ((u:$3FF0000000000000), (u:$0000000000000000), (u:$0000000000000000)));

  cLog2LL1: array[0..16, 0..2] of Tb64u64 = (
    ((u:$0000000000000000), (u:$0000000000000000), (u:$0000000000000000)),
    ((u:$3F7000024A000000), (u:$BD5799935D030650), (u:$BA0CACAD0B61C964)),
    ((u:$3F8000014A880000), (u:$BD586DD926D35EE8), (u:$BA2277AB54C00E64)),
    ((u:$3F87FFFF7B380000), (u:$3D15FFC7FF00E600), (u:$3A25870383C1225D)),
    ((u:$3F9000004B940000), (u:$BD44AA9DA7EA0120), (u:$BA247DECFE981022)),
    ((u:$3F94000064380000), (u:$BD5015D79FD24C70), (u:$3A2EC26C231DDB4C)),
    ((u:$3F97FFFF8FF80000), (u:$BD5A62A2B755F048), (u:$BA2EEBE936E9A3F9)),
    ((u:$3F9C000048180000), (u:$3D1D1A1BBCE5F580), (u:$BA24DD25E141AD16)),
    ((u:$3FA0000050EE0000), (u:$BD5B32450688D5E0), (u:$3A2577A1C80A65B0)),
    ((u:$3FA2000014760000), (u:$BD3BAC8F6C103D40), (u:$BA1157A2472C9532)),
    ((u:$3FA3FFFFF2440000), (u:$BD4FBAB0C9D05AA0), (u:$BA13680FFABF1A5B)),
    ((u:$3FA5FFFFB1CE0000), (u:$3D4EF9EC2BA49140), (u:$BA224850720B73AB)),
    ((u:$3FA7FFFFD9600000), (u:$3D514E3A91ABEFF0), (u:$BA140E3565C114DD)),
    ((u:$3FAA0000321C0000), (u:$3D3E33C6A0732340), (u:$BA27E8C937A981D8)),
    ((u:$3FABFFFFC56A0000), (u:$BD52076781FF0460), (u:$BA296E054EC8F31E)),
    ((u:$3FADFFFFDBA40000), (u:$BD4EA8B4476EDDE0), (u:$B9E07BFB9FDE770F)),
    ((u:$3FB000001FDA0000), (u:$3D55ED58A7FF2C40), (u:$39F512ACBB219717)));

  cLog2LL2: array[0..16, 0..2] of Tb64u64 = (
    ((u:$0000000000000000), (u:$0000000000000000), (u:$0000000000000000)),
    ((u:$3F30002866000000), (u:$BD5329E31412B688), (u:$3A2279DDC3585863)),
    ((u:$3F3FFFED31000000), (u:$3D5546DCB40B0518), (u:$BA26D048F41DEBEC)),
    ((u:$3F48000388000000), (u:$3D5FCE6789678C88), (u:$3A24959DC2FEA030)),
    ((u:$3F50000668400000), (u:$3D5C5CAEC0AA6620), (u:$BA289D58A2E0D287)),
    ((u:$3F54000937000000), (u:$BD5DF91BF5DB20F8), (u:$3A2D4FC35005498A)),
    ((u:$3F58000A2D800000), (u:$3D3E7464F400A700), (u:$BA26365A147CE4A8)),
    ((u:$3F5C00094A400000), (u:$BD5C55C526E65A68), (u:$BA2A7EA2CF16D2CE)),
    ((u:$3F60000345400000), (u:$3D48E815121136D0), (u:$BA168F8F7A924FD6)),
    ((u:$3F620000F6400000), (u:$3D57FB19A44DA298), (u:$BA2A05FC6ED4124D)),
    ((u:$3F63FFFDB7000000), (u:$3D5BFAC3B1ABE2E8), (u:$BA09AB1F45A1EFB1)),
    ((u:$3F66000516A00000), (u:$BD402544D9FAF830), (u:$3A2286C39ADFFA8C)),
    ((u:$3F67FFFFF4200000), (u:$BD52878B893DA790), (u:$3A2B55C5E09022B1)),
    ((u:$3F6A00056F400000), (u:$BD59D7CCF4019608), (u:$3A2F4B05553FE340)),
    ((u:$3F6BFFFE65200000), (u:$BD5EE85F47C2B5A0), (u:$3A253FA913D4E9A6)),
    ((u:$3F6E0001F7600000), (u:$3D42E594FB2AF5E0), (u:$B9CAFAEC77CA8D86)),
    ((u:$3F7000024A000000), (u:$BD5799935D030650), (u:$BA0CACAD0B61C964)));

  cLog2LL3: array[0..16, 0..2] of Tb64u64 = (
    ((u:$0000000000000000), (u:$0000000000000000), (u:$0000000000000000)),
    ((u:$3EF0014690000000), (u:$BD4743D8EA40A4D0), (u:$B9E75086A21C12AB)),
    ((u:$3F00014C18000000), (u:$3D48E7FA3CF59080), (u:$3A2A11BDCDB1ACE3)),
    ((u:$3F07FF17C8000000), (u:$3D5118DBBCCE8470), (u:$3A009806ECD19760)),
    ((u:$3F0FFFCBB8000000), (u:$BD5EAA70B2051E38), (u:$39E0D287C5AA12B4)),
    ((u:$3F14004294000000), (u:$3D5602B04D79B0C8), (u:$39F8D2C8DDF31648)),
    ((u:$3F1800A218000000), (u:$BD50A23112178630), (u:$BA2C4652435118F1)),
    ((u:$3F1BFF9304000000), (u:$BD4385DC475C16B0), (u:$B9FE7288F74BBBDA)),
    ((u:$3F1FFFF810000000), (u:$3D4355EDCFF4E190), (u:$3A0F8DAAE6C1DB11)),
    ((u:$3F22002FF2000000), (u:$3D23FBB7CD366140), (u:$3A293C0E5B4DCA92)),
    ((u:$3F23FFAC90000000), (u:$BD39BDE416DE25C0), (u:$BA2F4DF1D74776B3)),
    ((u:$3F25FFE340000000), (u:$BD521735D0955DC0), (u:$BA21B0E10765842F)),
    ((u:$3F28001B52000000), (u:$3D40C0949E968810), (u:$BA1D3771C0B220F3)),
    ((u:$3F2A0054C8000000), (u:$3D426936197B1050), (u:$3A2E59CAA334554A)),
    ((u:$3F2BFFD6F0000000), (u:$3D427097B0B27D40), (u:$3A278A37896A0D1F)),
    ((u:$3F2E00132C000000), (u:$3D3DF1E084120280), (u:$3A1C61F449A4B5D5)),
    ((u:$3F30002866000000), (u:$BD5329E31412B688), (u:$3A2279DDC3585863)));

  cLog2Cy0H: Tb64u64 = (u:$3FF71547652B82FE);
  cLog2Cy0L: Tb64u64 = (u:$3C7777D0FFDA0D24);
  cLog2Cy1H: Tb64u64 = (u:$BFE71547652B82FE);
  cLog2Cy1L: Tb64u64 = (u:$BC6777D0FFDA0D24);
  cLog2Cy2H: Tb64u64 = (u:$3FDEC709DC3A03FD);
  cLog2Cy2L: Tb64u64 = (u:$3C7D27E96BE541E5);
  cLog2Cl0:  Tb64u64 = (u:$BFD71547652B82FE);
  cLog2Cl1:  Tb64u64 = (u:$3FD2776C50F1FF14);
  cLog2Cl2:  Tb64u64 = (u:$BFCEC709DC3ECA5D);

function Log2AddDD(xh, xl, ch, cl: Double; out l: Double): Double; inline;
var s, d: Double;
begin
  s := xh + ch;
  d := s - xh;
  l := ((ch - d) + (xh + (d - s))) + (xl + cl);
  Result := s;
end;

function Log2Refine(x, a: Double): Double;
var
  t, v, w: Tb64u64;
  ex, e, i1, i2, i3, i4: Int32;
  k: Int32;
  i: UInt64;
  ed: Double;
  L0, L1, L2: Double;
  t12, t34, th, tl, dh, dl, sh, sl: Double;
  xh, xl, chp, clp, tth, ttl: Double;
  v0, v1, v2: Double;
begin
  t.f := x;
  ex := Int32(t.u shr 52);
  e  := ex - $3FF;
  if ex = 0 then
  begin
    k := pcr_clzll(t.u);
    e := e - (k - 12);
    t.u := t.u shl (k - 11);
  end;
  t.u := t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12);
  t.u := t.u or (UInt64($3FF) shl 52);
  ed := e;
  v.f := a - ed + Double(1.0000305175781250);  // 0x1.00008p+0
  i := (v.u - (UInt64($3FF) shl 52)) shr (52 - 16);
  i1 := Int32(i shr 12);
  i2 := Int32((i shr 8) and $F);
  i3 := Int32((i shr 4) and $F);
  i4 := Int32(i and $F);

  L0 := cLog2LL0[i1,0].f + cLog2LL1[i2,0].f + (cLog2LL2[i3,0].f + cLog2LL3[i4,0].f) + ed;
  L1 := cLog2LL0[i1,1].f + cLog2LL1[i2,1].f + (cLog2LL2[i3,1].f + cLog2LL3[i4,1].f);
  L2 := cLog2LL0[i1,2].f + cLog2LL1[i2,2].f + (cLog2LL2[i3,2].f + cLog2LL3[i4,2].f);

  t12 := cLog2T_t1[i1].f * cLog2T_t2[i2].f;
  t34 := cLog2T_t3[i3].f * cLog2T_t4[i4].f;
  th  := t12 * t34;
  tl  := pcr_fma(t12, t34, -th);
  dh  := th * t.f;
  dl  := pcr_fma(th, t.f, -dh);
  sh  := tl * t.f;
  sl  := pcr_fma(tl, t.f, -sh);

  pcr_fasttwosum(xh, xl, dh - Double(1.0), dl);
  xh := Log2AddDD(xh, xl, sh, sl, xl);

  // polynomial: sl = xh*(cl0 + xh*(cl1 + xh*cl2))
  sl := xh * (cLog2Cl0.f + xh * (cLog2Cl1.f + xh * cLog2Cl2.f));

  // polydd(xh, xl, 3, cy, &sl) — inline with initial *l = sl
  // i=2: ch = cy[2].h + sl, cl = ((cy[2].h - ch) + sl) + cy[2].l
  chp := cLog2Cy2H.f + sl;
  clp := ((cLog2Cy2H.f - chp) + sl) + cLog2Cy2L.f;
  // i=1
  chp := pcr_muldd(xh, xl, chp, clp, clp);
  tth := chp + cLog2Cy1H.f;  ttl := (cLog2Cy1H.f - tth) + chp;
  chp := tth;  clp := clp + ttl + cLog2Cy1L.f;
  // i=0
  chp := pcr_muldd(xh, xl, chp, clp, clp);
  tth := chp + cLog2Cy0H.f;  ttl := (cLog2Cy0H.f - tth) + chp;
  chp := tth;  clp := clp + ttl + cLog2Cy0L.f;
  sh := chp;  sl := clp;

  sh := pcr_muldd(xh, xl, sh, sl, sl);
  sh := Log2AddDD(sh, sl, L1, L2, sl);

  pcr_fasttwosum(v0, v2, L0, sh);
  pcr_fasttwosum(v1, v2, v2, sl);

  t.f := v1;
  if (t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) = 0 then
  begin
    w.f := v2;
    if ((w.u xor t.u) shr 63) <> 0 then
      Dec(t.u)
    else
      Inc(t.u);
    v1 := t.f;
  end;
  Result := v1 + v0;
end;

function pcr_log2(x: Double): Double;
var
  t: Tb64u64;
  ex, e, i1, i2, k: Int32;
  ir: UInt64;
  d: Int64;
  j, sum: UInt64;
  r, o, dxl, dxh, dx, dx2, f_poly, lt, lh, ll, ed: Double;
  eps, lb, ub: Double;
begin
  t.f := x;
  ex  := Int32(t.u shr 52);
  e   := ex - $3FF;

  if ex = 0 then
  begin
    if t.u = 0 then begin Result := -Infinity; Exit; end;
    k := pcr_clzll(t.u);
    e := e - (k - 12);
    t.u := t.u shl (k - 11);
    ex := 0;
  end;

  if ex >= $7FF then
  begin
    if (t.u shl 1) = 0 then begin Result := -Infinity; Exit; end;  // -0
    if (t.u shl 1) > (UInt64($7FF) shl 53) then
    begin Result := x + x; Exit; end;  // NaN
    if (t.u shr 63) <> 0 then
    begin Result := cNaNDouble; Exit; end;  // x < 0
    Result := x;  // +Inf
    Exit;
  end;

  t.u := t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12);
  if t.u = 0 then begin Result := e; Exit; end;

  ed  := e;
  ir  := t.u shr (52 - 5);
  d   := Int64(t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 17));

  // sum = t.u + (B.c0<<33) + B.c1*(d>>16), signed 64-bit safe
  sum := t.u + (UInt64(cLog2B_c0[Integer(ir)]) shl 33)
             + UInt64(Int64(cLog2B_c1[Integer(ir)]) * (d shr 16));
  j   := sum shr (52 - 10);
  t.u := t.u or (UInt64($3FF) shl 52);
  i1  := Int32(j shr 5);
  i2  := Int32(j and $1F);

  r    := cLog2R1[i1].f * cLog2R2[i2].f;
  o    := r * t.f;
  dxl  := pcr_fma(r, t.f, -o);
  dxh  := o - cLog2L2H.f;
  dx   := dxh + dxl;
  dx2  := dx * dx;
  f_poly := dx2 * ((cLog2C0.f + dx * cLog2C1.f) + dx2 * (cLog2C2.f + dx * cLog2C3.f));

  lt := (cLog2L1Hi[i1].f + cLog2L2Hi[i2].f) + ed;
  lh := lt + dxh;
  ll := (lt - lh) + dxh;
  ll := ll + ((cLog2L1Lo[i1].f + cLog2L2Lo[i2].f) + dxl) + dxh * cLog2L2L.f;
  ll := ll + f_poly;

  eps := Double(1.6e-22);
  lb  := lh + (ll - eps);
  ub  := lh + (ll + eps);
  if lb = ub then begin Result := lb; Exit; end;

  Result := Log2Refine(x, ub);
end;

// ---------------------------------------------------------------------------
// pcr_acos — correctly-rounded binary64 arc cosine.
// Ported from core-math/src/binary64/acos/acos.c by Alexei Sibidanov.
// ---------------------------------------------------------------------------

const
  // cc[33][8]: asin approximation coefficients; cc[j][0..1] = (fh, fl) at anchor j,
  // cc[j][2..7] = polynomial coefficients in t = x^2 - j/128.
  cAcosCC: array[0..32, 0..7] of Tb64u64 = (
    ((u:$3FF0000000000000),(u:$0000000000000000),(u:$3FC5555555555555),(u:$3FB33333333333E4),(u:$3FA6DB6DB6D31F82),(u:$3F9F1C71F6889397),(u:$3F96E874B7045B46),(u:$3F91F753132271E2)),
    ((u:$3FF0055A27E0D033),(u:$BC9D9BA10494C062),(u:$3FC57C00CB5D6C4D),(u:$3FB37881F5649A75),(u:$3FA759AF49D494DD),(u:$3FA002E1864DDA2E),(u:$3F97C2D5D468CDD9),(u:$3F9292834C025357)),
    ((u:$3FF00ABE0C129E1E),(u:$3C67CEB0EE49D42A),(u:$3FC5A3385D5C7BA5),(u:$3FB3BF51056F6637),(u:$3FA7DBA76B124B37),(u:$3FA07BE4B02E94C4),(u:$3F98A6BB92513F01),(u:$3F936AFD4C615AEC)),
    ((u:$3FF0102BCFFD6ACD),(u:$BC8C2294C65D2E86),(u:$3FC5CAFF17351901),(u:$3FB407ABBC04FEB3),(u:$3FA86179B8005949),(u:$3FA0F97520BD4E72),(u:$3F99950C5C89F3DF),(u:$3F944F2344E7B664)),
    ((u:$3FF015A397CF0F1C),(u:$BC8EEBD6CCFE3EE3),(u:$3FC5F3581BE7B08B),(u:$3FB4519DDF1AE531),(u:$3FA8EB4B6EE35E92),(u:$3FA17BC85414CD46),(u:$3F9A8E5895E3FCF9),(u:$3F953FAFDC629400)),
    ((u:$3FF01B2588811EEB),(u:$3C47193E5D0A915F),(u:$3FC61C46A67205D2),(u:$3FB49D33A6EEAE0B),(u:$3FA979438563C014),(u:$3FA20316AE977F05),(u:$3F9B9339AFB53AA4),(u:$3F963D6B02C42D0A)),
    ((u:$3FF020B1C7DF0575),(u:$BC8DD547E329C1E5),(u:$3FC645CE0AB901BB),(u:$3FB4EA79C34FC7A6),(u:$3FAA0B8AC08940EC),(u:$3FA28F9BABD0629B),(u:$3F9CA452CF90A55E),(u:$3F97492B016730EF)),
    ((u:$3FF026487C8C5D71),(u:$BC95FD9B68DC3B6E),(u:$3FC66FF1B67D5D70),(u:$3FB5397D613373EB),(u:$3FAAA24BCE3AEB4A),(u:$3FA3219610B150AC),(u:$3F9DC251825103F1),(u:$3F9863D5A3932532)),
    ((u:$3FF02BE9CE0B87CD),(u:$3C7E5D09DA2E0F04),(u:$3FC69AB5325BC359),(u:$3FB58A4C3097AAB3),(u:$3FAB3DB3605F46F2),(u:$3FA3B94821742CAB),(u:$3F9EEDEE7DA72A15),(u:$3F998E6179A3E9A0)),
    ((u:$3FF03195E4C483F1),(u:$BC95DB10AD66EACB),(u:$3FC6C61C22D908F0),(u:$3FB5DCF46AB9F2CB),(u:$3FABDDF049BB1F4D),(u:$3FA456F7DB6AC768),(u:$3FA013F738BD7BB3),(u:$3F9AC9D739783D21)),
    ((u:$3FF0374CEA0C0C9F),(u:$BC9917BFF5241C76),(u:$3FC6F22A497B2EC0),(u:$3FB63184D8A79DB5),(u:$3FAC83339CAF946E),(u:$3FA4FAEF331019D4),(u:$3FA0B8917547D678),(u:$3F9C17533F147E1C)),
    ((u:$3FF03D0F082AFCC8),(u:$BC9018BBCDDB49EB),(u:$3FC71EE385EFDF06),(u:$3FB6880CDA2D3885),(u:$3FAD2DB0CBFAE54D),(u:$3FA5A57C56B50C5E),(u:$3FA16535A40098B2),(u:$3F9D780730B8EBB8)),
    ((u:$3FF042DC6A65FFBF),(u:$BC8C7EA28DCE95D1),(u:$3FC74C4BD7412F9E),(u:$3FB6E09C6D2B72BC),(u:$3FADDD9DCDA253DE),(u:$3FA656F1F62B5001),(u:$3FA21A5AE2AC77EE),(u:$3F9EED3BCA067F0E)),
    ((u:$3FF048B53D05907B),(u:$3C9634FFFED6E2A6),(u:$3FC77A675D1978BE),(u:$3FB73B4435583415),(u:$3FAE9333402EBBF3),(u:$3FA70FA78FD9F73F),(u:$3FA2D8804D934FE1),(u:$3FA03C29691A281C)),
    ((u:$3FF04E99AD5E4BCD),(u:$BC9E97A72FE827E0),(u:$3FC7A93A5917200C),(u:$3FB7981584731C05),(u:$3FAF4EAC9268FAE2),(u:$3FA7CFF9C3B19721),(u:$3FA3A02D9C1E0145),(u:$3FA10D64A0E56953)),
    ((u:$3FF05489E9D99995),(u:$3C8D177637EC6A2B),(u:$3FC7D8C930314681),(u:$3FB7F72262F532E4),(u:$3FB0082416E39013),(u:$3FA8984AAC80DDF4),(u:$3FA471F3CAF18EB8),(u:$3FA1EB1CCE6DD570)),
    ((u:$3FF05A8621FEB16B),(u:$BC7E5B33B1407C5F),(u:$3FC809186C2E57DD),(u:$3FB8587D99442DC8),(u:$3FB06C23D1DFCB7F),(u:$3FA969024036DD22),(u:$3FA54E6DD4D2AF33),(u:$3FA2D62F439F2A31)),
    ((u:$3FF0608E867BFF30),(u:$3C8CBEF5D8580027),(u:$3FC83A2CBD2D8BA2),(u:$3FB8BC3AB9724C6E),(u:$3FB0D377EF1E0C39),(u:$3FAA428EB7ADDF84),(u:$3FA636417BC01FF2),(u:$3FA3CF8ACC7EB2A0)),
    ((u:$3FF066A34930EC8D),(u:$BC9480F445FEDAD1),(u:$3FC86C0AFB447A74),(u:$3FB9226E29948D9C),(u:$3FB13E44A9B5A3A6),(u:$3FAB2564FEA8B3FE),(u:$3FA72A2023D92458),(u:$3FA4D8313CEC3485)),
    ((u:$3FF06CC49D38146C),(u:$BC8B55394F4FC07B),(u:$3FC89EB82831FEED),(u:$3FB98B2D2EB9BB23),(u:$3FB1ACB01E9AB414),(u:$3FAC12012CBD00C6),(u:$3FA82AC7C1D15C38),(u:$3FA5F13925C6EDCA)),
    ((u:$3FF072F2B6F1E601),(u:$BC92DCBB05419970),(u:$3FC8D2397127AEBB),(u:$3FB9F68DF88DA51D),(u:$3FB21EE26A4F62A1),(u:$3FAD08E707F7AE6F),(u:$3FA93903DEE3FEB0),(u:$3FA71BCFB5C57B59)),
    ((u:$3FF0792DCC0FBD20),(u:$BC75BF23EE4F9D54),(u:$3FC9069430AB508A),(u:$3FBA64A7ADB4CD85),(u:$3FB29505C8B48349),(u:$3FAE0AA2921CFA60),(u:$3FAA55AEB46F4322),(u:$3FA8593ACAD3BECE)),
    ((u:$3FF07F76139F761D),(u:$3C9FA1046481BB82),(u:$3FC93BCDF091CCA6),(u:$3FBAD59278EDC42F),(u:$3FB30F46B7261652),(u:$3FAF17C8A17C843E),(u:$3FAB81B2619E15B5),(u:$3FA9AADB395F5AE4)),
    ((u:$3FF085CBC61783C1),(u:$3C90A6E9EFA20176),(u:$3FC971EC6C1531E4),(u:$3FBB496797068912),(u:$3FB38DD419140184),(u:$3FB0187BC3357BBB),(u:$3FACBE0A3DCAFE26),(u:$3FAB122F4FA499D0)),
    ((u:$3FF08C2F1D638E4C),(u:$3C7B47C159534A3D),(u:$3FC9A8F592078624),(u:$3FBBC04165B57AB2),(u:$3FB410DF5F4BED1D),(u:$3FB0AB6BDF478C71),(u:$3FAE0BC44A945C64),(u:$3FAC90D59BCD5701)),
    ((u:$3FF092A054F1A2FC),(u:$BC92F657224E9830),(u:$3FC9E0EF87243A2C),(u:$3FBC3A3B7366A278),(u:$3FB4989CB22E2175),(u:$3FB1450E5BA7AD39),(u:$3FAF6C02C8F0EF93),(u:$3FAE288FFC8D182C)),
    ((u:$3FF0991FA9BFFBF4),(u:$BC5CA1140A1ABBF4),(u:$3FCA19E0A8823B80),(u:$3FBCB772900F9C24),(u:$3FB525431F0CBB2E),(u:$3FB1E5C2D06804E1),(u:$3FB06FFEFA7AA6B8),(u:$3FAFDB4704DCA347)),
    ((u:$3FF09FAD5A6B68F9),(u:$3C7AA1F06E92964E),(u:$3FCA53CF8E28C50E),(u:$3FBD3804DF1DE350),(u:$3FB5B70CC8FA98DC),(u:$3FB28DEF298C979B),(u:$3FB13482F6347EEB),(u:$3FB0D586DE48358C)),
    ((u:$3FF0A649A73E61F2),(u:$3C874AC0D817E9C7),(u:$3FCA8EC30DC93891),(u:$3FBDBC11EA950625),(u:$3FB64E371D5616D3),(u:$3FB33E0023936249),(u:$3FB204426263066A),(u:$3FB1CD12E4629723)),
    ((u:$3FF0ACF4D240CCC4),(u:$3C9DA890F3B40BD3),(u:$3FCACAC23DA07797),(u:$3FBE43BAB7741A98),(u:$3FB6EB030C631819),(u:$3FB3F669D2EB516E),(u:$3FB2E0006AE505AE),(u:$3FB2D58204457C82)),
    ((u:$3FF0B3AF1F4880BB),(u:$3C7F450FB78D32BA),(u:$3FCB07D4778263AF),(u:$3FBECF21DB7BE0EF),(u:$3FB78DB5465013E4),(u:$3FB4B7A8376F0996),(u:$3FB3C88F9F2EF221),(u:$3FB3F02AD9EB9753)),
    ((u:$3FF0BA78D40A9260),(u:$BC957B07A441E242),(u:$3FCB46015C126262),(u:$3FBF5E6B94713F3D),(u:$3FB836967D0AFECF),(u:$3FB5823FDD1707B9),(u:$3FB4BED355269DC2),(u:$3FB51E83065121CF)),
    ((u:$3FF0C152382D7366),(u:$BC9EE6913347C2A6),(u:$3FCB8550D62BFB6E),(u:$3FBFF1BDE0FA3CAD),(u:$3FB8E5F3AB550989),(u:$3FB656BE8B38EBAF),(u:$3FB5C3C13008A099),(u:$3FB662225A1B4F77)));

  // sin(pi/64 * j), double-double form, 0 <= j <= 32.
  cAcosSHi: array[0..32] of Tb64u64 = (
    (u:$0000000000000000),(u:$3FA91F65F10DD814),(u:$3FB917A6BC29B42C),(u:$3FC2C8106E8E613A),
    (u:$3FC8F8B83C69A60B),(u:$3FCF19F97B215F1B),(u:$3FD294062ED59F06),(u:$3FD58F9A75AB1FDD),
    (u:$3FD87DE2A6AEA963),(u:$3FDB5D1009E15CC0),(u:$3FDE2B5D3806F63B),(u:$3FE073879922FFEE),
    (u:$3FE1C73B39AE68C8),(u:$3FE30FF7FCE17035),(u:$3FE44CF325091DD6),(u:$3FE57D69348CECA0),
    (u:$3FE6A09E667F3BCD),(u:$3FE7B5DF226AAFAF),(u:$3FE8BC806B151741),(u:$3FE9B3E047F38741),
    (u:$3FEA9B66290EA1A3),(u:$3FEB728345196E3E),(u:$3FEC38B2F180BDB1),(u:$3FECED7AF43CC773),
    (u:$3FED906BCF328D46),(u:$3FEE212104F686E5),(u:$3FEE9F4156C62DDA),(u:$3FEF0A7EFB9230D7),
    (u:$3FEF6297CFF75CB0),(u:$3FEFA7557F08A517),(u:$3FEFD88DA3D12526),(u:$3FEFF621E3796D7E),
    (u:$3FF0000000000000));
  cAcosSLo: array[0..32] of Tb64u64 = (
    (u:$0000000000000000),(u:$BC2912BD0D569A90),(u:$BC3E2718D26ED688),(u:$3C513000A89A11E0),
    (u:$BC626D19B9FF8D82),(u:$BC642DEEF11DA2C4),(u:$BC75D28DA2C4612D),(u:$BC1EFDC0D58CF620),
    (u:$BC672CEDD3D5A610),(u:$3C65B362CB974183),(u:$3C5E0D891D3C6841),(u:$BC8A5A014347406C),
    (u:$3C8B25DD267F6600),(u:$BC6EFCC626F74A6F),(u:$3C68076A2CFDC6B3),(u:$BC875720992BFBB2),
    (u:$BC8BDD3413B26456),(u:$BC70F537ACDF0AD7),(u:$BC82C5E12ED1336D),(u:$BC830EE286712474),
    (u:$3C39F630E8B6DAC8),(u:$BC8BC69F324E6D61),(u:$BC76E0B1757C8D07),(u:$BC5E7B6BB5AB58AE),
    (u:$3C7457E610231AC2),(u:$BC8014C76C126527),(u:$3C8760B1E2E3F81E),(u:$3C752C7ADC6B4989),
    (u:$3C7562172A361FD3),(u:$BC87A0A8CA13571F),(u:$BC887DF6378811C7),(u:$BC6C57BC2E24AA15),
    (u:$0000000000000000));

  // c[5][2]: inner asin polynomial (after substitution) for the accurate path.
  cAcosCHi: array[0..4] of Tb64u64 = (
    (u:$3FF0000000000000),(u:$3FC5555555555555),(u:$3FB3333333333333),
    (u:$3FA6DB6DB6DB6DB7),(u:$3F9F1C71C71C6D5B));
  cAcosCLo: array[0..4] of Tb64u64 = (
    (u:$B93FC2C76456515B),(u:$3C65555555623513),(u:$3C49997E3427441B),
    (u:$BC1CB95FF08658E6),(u:$3C3B125BCCDCC89E));
  // ct[3]: outer polynomial tail (degree 5..7) used as seed for polydd.
  cAcosCt: array[0..2] of Tb64u64 = (
    (u:$3F96E8BA2EC8CB69),(u:$3F91C4EA7A15C997),(u:$3F8CA8355D39BB67));

  cAcosOffH: array[0..1] of Tb64u64 = ((u:$0000000000000000),(u:$400921FB54442D18));
  cAcosOffL: array[0..1] of Tb64u64 = ((u:$0000000000000000),(u:$3CA1A62633145C07));

  cAcosPiHalfH:  Tb64u64 = (u:$3FF921FB54442D18);  //  0x1.921fb54442d18p+0
  cAcosPiHalfL:  Tb64u64 = (u:$3C91A62633145C07);  //  0x1.1a62633145c07p-54
  cAcosRefScale: Tb64u64 = (u:$40345F306DC9C883);  //  0x1.45f306dc9c883p+4 (=64/pi)
  cAcosPi64H:    Tb64u64 = (u:$3FA921FB54442D00);  //  pi/64 hi
  cAcosPi64M:    Tb64u64 = (u:$3CA8469898CC5180);  //  pi/64 mid
  cAcosPi64L:    Tb64u64 = (u:$B99FC8F8CBB5BF6C);  //  pi/64 lo
  cAcosC2fK:     Tb64u64 = (u:$3FB8000000000000);  //  0x1.8p-4 (anti-cancellation add)
  cAcosEps1:     Tb64u64 = (u:$3CB8C00000000000);  //  0x1.8cp-52
  cAcosEps2:     Tb64u64 = (u:$3960000000000000);  //  0x1p-105
  cAcosP5:       Tb64u64 = (u:$4040000000000000);  //  32.0 = 0x1p5
  cAcosN7:       Tb64u64 = (u:$3F80000000000000);  //  0x1p-7
  cAcosP7:       Tb64u64 = (u:$4060000000000000);  //  128.0 = 0x1p7

// Rare worst-case x inputs in AcosRefine where the generic ±1-ULP bump below
// would disagree with the correctly-rounded result. Patterns from acos.c.
const
  cAcosWcIn: array[0..6] of Tb64u64 = (
    (u:$BFC771164BFD1F84),(u:$BFE4510EE8EB4E67),(u:$BFD011C543F23A17),
    (u:$3FEFFFFFFFFFFDC0),(u:$3FB53EA6C7255E88),(u:$3F4FD737BE914578),
    (u:$3FEFFFFFFFFFFF70));
  cAcosWcHi: array[0..6] of Tb64u64 = (
    (u:$3FFC14601DAAF657),(u:$400211C0E2C2559E),(u:$3FFD318C90D9E8B7),
    (u:$3E98000000000024),(u:$3FF7CDACB6BBE707),(u:$3FF91E006D41D8D8),
    (u:$3E88000000000009));
  cAcosWcLo: array[0..6] of Tb64u64 = (
    (u:$BC90000000000000),(u:$BCA0000000000000),(u:$BC90000000000000),
    (u:$3B30000000000000),(u:$3C90000000000000),(u:$3CA8000000000000),
    (u:$3B20000000000000));

function AcosRefine(x, phi: Double): Double;
var
  s2, dx2, c2h, c2l, c2f, ch, cl: Double;
  jf: Int64;
  Ch_v, Cl_v, Sh_v, Sl_v, ax_r: Double;
  dsh, dsl, dch, dcl: Double;
  Sc, dSc, Cs, dCs, v, dv: Double;
  sgn, jtd: Double;
  v2, dv2: Double;
  fh, fl: Double;
  chp, clp, th, tl: Double;
  ph, pl, ps: Double;
  sh_tmp, d_tmp, sl_tmp: Double;
  t_u, w_u, xu: Tb64u64;
  e: Int64;
  m, ebit: UInt64;
  k: Int32;
begin
  s2  := x * x;
  dx2 := pcr_fma(x, x, -s2);
  // c2h, c2l = fasttwosub(1, s2, &c2l) == fasttwosum(1, -s2)
  pcr_fasttwosum(c2h, c2l, Double(1.0), -s2);
  c2l := c2l - dx2;
  pcr_fasttwosum(c2h, c2l, c2h, c2l);

  c2f := pcr_fma(x, -x, Double(1.0));
  ch  := Sqrt(c2f);
  cl  := (c2l - pcr_fma(ch, ch, -c2f)) * ((Double(0.5) / c2f) * ch);

  jf := Trunc(pcr_roundeven(Abs(phi - cAcosPiHalfH.f) * cAcosRefScale.f));

  Ch_v := cAcosSHi[32 - jf].f;  Cl_v := cAcosSLo[32 - jf].f;
  Sh_v := cAcosSHi[jf].f;       Sl_v := cAcosSLo[jf].f;

  ax_r := Abs(x);
  dsh := ax_r - Sh_v;  dsl := -Sl_v;
  dch := ch   - Ch_v;  dcl := cl - Cl_v;

  Sc  := pcr_fma(Sh_v, dch, cAcosC2fK.f) - cAcosC2fK.f;
  dSc := pcr_fma(Sh_v, dch, -Sc);

  Cs  := pcr_fma(Ch_v, dsh, cAcosC2fK.f) - cAcosC2fK.f;
  dCs := pcr_fma(Ch_v, dsh, -Cs);

  v  := Cs - Sc;
  dv := (Ch_v * dsl + Cl_v * dsh) - (Sh_v * dcl + Sl_v * dch) - (dSc - dCs);
  pcr_fasttwosum(v, dv, v, dv);

  if x >= Double(0.0) then sgn := Double(1.0) else sgn := -Double(1.0);
  jtd := Double(32) - jf * sgn;           // 0..64

  v2 := pcr_muldd(v, dv, v, dv, dv2);
  v  := v  * (-sgn);
  dv := dv * (-sgn);

  fl := v2 * (cAcosCt[0].f + v2 * (cAcosCt[1].f + v2 * cAcosCt[2].f));
  // fh = polydd(v2, dv2, 5, c, &fl) with incoming *l = fl:
  //   i=4: ch = fasttwosum(c[4][0], fl_in, &fl_out); cl = c[4][1] + fl_out
  pcr_fasttwosum(chp, fl, cAcosCHi[4].f, fl);
  clp := cAcosCLo[4].f + fl;
  // i=3..0: ch = muldd(v2, dv2, ch, cl, &cl); ch = fastsum(c[i][0], c[i][1], ch, cl, &cl)
  chp := pcr_muldd(v2, dv2, chp, clp, clp);
  pcr_fasttwosum(th, tl, cAcosCHi[3].f, chp);
  chp := th;  clp := (cAcosCLo[3].f + clp) + tl;
  chp := pcr_muldd(v2, dv2, chp, clp, clp);
  pcr_fasttwosum(th, tl, cAcosCHi[2].f, chp);
  chp := th;  clp := (cAcosCLo[2].f + clp) + tl;
  chp := pcr_muldd(v2, dv2, chp, clp, clp);
  pcr_fasttwosum(th, tl, cAcosCHi[1].f, chp);
  chp := th;  clp := (cAcosCLo[1].f + clp) + tl;
  chp := pcr_muldd(v2, dv2, chp, clp, clp);
  pcr_fasttwosum(th, tl, cAcosCHi[0].f, chp);
  chp := th;  clp := (cAcosCLo[0].f + clp) + tl;
  fh := chp;  fl := clp;

  fh := pcr_muldd(v, dv, fh, fl, fl);

  ph := jtd * cAcosPi64H.f;
  pl := cAcosPi64M.f * jtd;
  ps := cAcosPi64L.f * jtd;
  // pl = sum(fh, fl, pl, ps, &ps): uses (non-fast) twosum
  sh_tmp := fh + pl;
  d_tmp  := sh_tmp - fh;
  sl_tmp := (pl - d_tmp) + (fh + (d_tmp - sh_tmp));
  ps     := (fl + ps) + sl_tmp;
  pl     := sh_tmp;
  pcr_fasttwosum(ph, pl, ph, pl);
  pcr_fasttwosum(pl, ps, pl, ps);
  pcr_fasttwosum(ph, pl, ph, pl);
  pcr_fasttwosum(pl, ps, pl, ps);

  t_u.f := pl;
  e := Int64((t_u.u shr 52) and UInt64($7FF)) - 1023;
  e := 52 - (107 + e);
  if e < 0 then e := 0;
  if e > 52 then e := 52;
  m := (UInt64(1) shl 52) - (UInt64(1) shl e);
  if e = 0 then e := 64;
  ebit := UInt64(1) shl (e - 1);
  if ((t_u.u + ebit) and m) = 0 then
  begin
    xu.f := x;
    for k := 0 to 6 do
      if xu.u = cAcosWcIn[k].u then
      begin
        Result := cAcosWcHi[k].f + cAcosWcLo[k].f;
        Exit;
      end;
    w_u.f := ps;
    if ((w_u.u xor t_u.u) shr 63) <> 0 then
      Dec(t_u.u)
    else
      Inc(t_u.u);
    pl := t_u.f;
  end;

  Result := ph + pl;
end;

function pcr_acos(x: Double): Double;
var
  ix: Tb64u64;
  ax: UInt64;
  k: Int64;
  j: Int64;
  f0h, f0l, t, z, zl, jd: Double;
  t2, d_poly, fh, fl, eps, lb, ub, sum_sh, fastsum_sl: Double;
begin
  ix.f := x;
  ax := ix.u shl 1;

  if ax > UInt64($7FC0000000000000) then
  begin
    // |x| > 0.5 branch: k = sign bit of x
    k := ix.u shr 63;
    f0h := cAcosOffH[k].f;
    f0l := cAcosOffL[k].f;
    if ax >= UInt64($7FE0000000000000) then
    begin
      if ax = UInt64($7FE0000000000000) then
      begin Result := f0h + f0l; Exit; end;                   // |x| = 1
      if ax > UInt64($FFE0000000000000) then
      begin Result := x + x; Exit; end;                       // NaN
      Result := Double(0.0) / Double(0.0); Exit;              // |x| > 1: NaN (negative, like C)
    end;
    t  := Double(2.0) - Double(2.0) * Abs(x);
    jd := pcr_roundeven(t * cAcosP5.f);                       // t * 32
    if x >= Double(0.0) then z := Sqrt(t) else z := -Sqrt(t);
    zl := pcr_fma(z, z, -t) * ((-Double(0.5) / t) * z);
    t  := Double(0.25) * t - jd * cAcosN7.f;                  // 0.25*t - jd/128
  end
  else
  begin
    f0h := cAcosPiHalfH.f;
    f0l := cAcosPiHalfL.f;
    // For |x| <= 0x1.cb3b399d747f2p-55 acos rounds to pi/2 in RN
    if ax <= UInt64($7919676733AE8FE4) then
    begin Result := f0h + f0l; Exit; end;
    t  := x * x;
    jd := pcr_roundeven(t * cAcosP7.f);                       // t * 128
    t  := pcr_fma(x, x, -cAcosN7.f * jd);
    z  := -x;
    zl := Double(0.0);
  end;

  j := Trunc(jd);
  t2 := t * t;
  d_poly := t * ((cAcosCC[j,2].f + t * cAcosCC[j,3].f) +
                 t2 * ((cAcosCC[j,4].f + t * cAcosCC[j,5].f) +
                       t2 * (cAcosCC[j,6].f + t * cAcosCC[j,7].f)));
  fh := cAcosCC[j,0].f;
  fl := cAcosCC[j,1].f + d_poly;
  fh := pcr_muldd(z, zl, fh, fl, fl);
  // fastsum(f0h, f0l, fh, fl, &fl_out): sh = fasttwosum(f0h, fh, &sl); fl_out = (f0l + fl_in) + sl
  fastsum_sl := Double(0.0);
  pcr_fasttwosum(sum_sh, fastsum_sl, f0h, fh);
  fl := (f0l + fl) + fastsum_sl;
  fh := sum_sh;

  eps := Abs(z * t) * cAcosEps1.f + cAcosEps2.f;
  lb  := fh + (fl - eps);
  ub  := fh + (fl + eps);
  if lb <> ub then
    Result := AcosRefine(x, lb)
  else
    Result := lb;
end;

// ---------------------------------------------------------------------------
// pcr_tanh — correctly-rounded binary64 hyperbolic tangent.
// Ported from core-math/src/binary64/tanh/tanh.c by Alexei Sibidanov.
// ---------------------------------------------------------------------------

const
  // t0[64][2]: (lo, hi) pair; hi ~ 2^(i0/64) on the coarse grid.
  // cExpT0: shared by exp/exp2/exp10/expm1/tanh/sinh families.
  // Each entry is a double-double approximation of 2^(i/64) (lo, hi).
  cExpT0: array[0..63, 0..1] of Tb64u64 = (
    ((u:$0000000000000000),(u:$3FF0000000000000)),
    ((u:$BC719083535B085E),(u:$3FF02C9A3E778061)),
    ((u:$3C8D73E2A475B466),(u:$3FF059B0D3158574)),
    ((u:$3C6186BE4BB28500),(u:$3FF0874518759BC8)),
    ((u:$3C98A62E4ADC610A),(u:$3FF0B5586CF9890F)),
    ((u:$3C403A1727C57B52),(u:$3FF0E3EC32D3D1A2)),
    ((u:$BC96C51039449B3A),(u:$3FF11301D0125B51)),
    ((u:$BC932FBF9AF1369E),(u:$3FF1429AAEA92DE0)),
    ((u:$BC819041B9D78A76),(u:$3FF172B83C7D517B)),
    ((u:$3C8E5B4C7B4968E4),(u:$3FF1A35BEB6FCB75)),
    ((u:$3C9E016E00A2643C),(u:$3FF1D4873168B9AA)),
    ((u:$3C8DC775814A8494),(u:$3FF2063B88628CD6)),
    ((u:$3C99B07EB6C70572),(u:$3FF2387A6E756238)),
    ((u:$3C82BD339940E9DA),(u:$3FF26B4565E27CDD)),
    ((u:$3C8612E8AFAD1256),(u:$3FF29E9DF51FDEE1)),
    ((u:$3C90024754DB41D4),(u:$3FF2D285A6E4030B)),
    ((u:$3C86F46AD23182E4),(u:$3FF306FE0A31B715)),
    ((u:$3C932721843659A6),(u:$3FF33C08B26416FF)),
    ((u:$BC963AEABF42EAE2),(u:$3FF371A7373AA9CB)),
    ((u:$BC75E436D661F5E2),(u:$3FF3A7DB34E59FF7)),
    ((u:$3C8ADA0911F09EBC),(u:$3FF3DEA64C123422)),
    ((u:$BC5EF3691C309278),(u:$3FF4160A21F72E2A)),
    ((u:$3C489B7A04EF80D0),(u:$3FF44E086061892D)),
    ((u:$3C73C1A3B69062F0),(u:$3FF486A2B5C13CD0)),
    ((u:$3C7D4397AFEC42E2),(u:$3FF4BFDAD5362A27)),
    ((u:$BC94B309D25957E4),(u:$3FF4F9B2769D2CA7)),
    ((u:$BC807ABE1DB13CAC),(u:$3FF5342B569D4F82)),
    ((u:$3C99BB2C011D93AC),(u:$3FF56F4736B527DA)),
    ((u:$3C96324C054647AC),(u:$3FF5AB07DD485429)),
    ((u:$3C9BA6F93080E65E),(u:$3FF5E76F15AD2148)),
    ((u:$BC9383C17E40B496),(u:$3FF6247EB03A5585)),
    ((u:$BC9BB60987591C34),(u:$3FF6623882552225)),
    ((u:$BC9BDD3413B26456),(u:$3FF6A09E667F3BCD)),
    ((u:$BC6BBE3A683C88AA),(u:$3FF6DFB23C651A2F)),
    ((u:$BC816E4786887A9A),(u:$3FF71F75E8EC5F74)),
    ((u:$BC90245957316DD4),(u:$3FF75FEB564267C9)),
    ((u:$BC841577EE049930),(u:$3FF7A11473EB0187)),
    ((u:$3C705D02BA15797E),(u:$3FF7E2F336CF4E62)),
    ((u:$BC9D4C1DD41532D8),(u:$3FF82589994CCE13)),
    ((u:$BC9FC6F89BD4F6BA),(u:$3FF868D99B4492ED)),
    ((u:$3C96E9F156864B26),(u:$3FF8ACE5422AA0DB)),
    ((u:$3C85CC13A2E3976C),(u:$3FF8F1AE99157736)),
    ((u:$BC675FC781B57EBC),(u:$3FF93737B0CDC5E5)),
    ((u:$BC9D185B7C1B85D0),(u:$3FF97D829FDE4E50)),
    ((u:$3C7C7C46B071F2BE),(u:$3FF9C49182A3F090)),
    ((u:$BC9359495D1CD532),(u:$3FFA0C667B5DE565)),
    ((u:$BC9D2F6EDB8D41E2),(u:$3FFA5503B23E255D)),
    ((u:$3C90FAC90EF7FD32),(u:$3FFA9E6B5579FDBF)),
    ((u:$3C97A1CD345DCC82),(u:$3FFAE89F995AD3AD)),
    ((u:$BC62805E3084D708),(u:$3FFB33A2B84F15FB)),
    ((u:$BC75584F7E54AC3A),(u:$3FFB7F76F2FB5E47)),
    ((u:$3C823DD07A2D9E84),(u:$3FFBCC1E904BC1D2)),
    ((u:$3C811065895048DE),(u:$3FFC199BDD85529C)),
    ((u:$3C92884DFF483CAC),(u:$3FFC67F12E57D14B)),
    ((u:$3C7503CBD1E949DC),(u:$3FFCB720DCEF9069)),
    ((u:$BC9CBC3743797A9C),(u:$3FFD072D4A07897C)),
    ((u:$3C82ED02D75B3706),(u:$3FFD5818DCFBA487)),
    ((u:$3C9C2300696DB532),(u:$3FFDA9E603DB3285)),
    ((u:$BC91A5CD4F184B5C),(u:$3FFDFC97337B9B5F)),
    ((u:$3C839E8980A9CC90),(u:$3FFE502EE78B3FF6)),
    ((u:$BC9E9C23179C2894),(u:$3FFEA4AFA2A490DA)),
    ((u:$3C9DC7F486A4B6B0),(u:$3FFEFA1BEE615A27)),
    ((u:$3C99D3E12DD8A18A),(u:$3FFF50765B6E4540)),
    ((u:$3C874853F3A5931E),(u:$3FFFA7C1819E90D8)));

  // t1[64][2]: fine grid, (lo, hi) pair.
  // cExpT1: shared by exp/exp2/tanh families. Double-double approx of 2^(i/4096).
  cExpT1: array[0..63, 0..1] of Tb64u64 = (
    ((u:$0000000000000000),(u:$3FF0000000000000)),
    ((u:$3C9AE8E38C59C72A),(u:$3FF000B175EFFDC7)),
    ((u:$BC57B5D0D58EA8F4),(u:$3FF00162F3904052)),
    ((u:$3C94115CB6B16A8E),(u:$3FF0021478E11CE6)),
    ((u:$BC8D7C96F201BB2E),(u:$3FF002C605E2E8CF)),
    ((u:$3C984711D4C35EA0),(u:$3FF003779A95F959)),
    ((u:$BC80484245243778),(u:$3FF0042936FAA3D8)),
    ((u:$BC94B237DA2025FA),(u:$3FF004DADB113DA0)),
    ((u:$BC75E00E62D6B30E),(u:$3FF0058C86DA1C0A)),
    ((u:$3C9A1D6CEDBB9480),(u:$3FF0063E3A559473)),
    ((u:$BC94ACF197A00142),(u:$3FF006EFF583FC3D)),
    ((u:$BC6EAF2EA42391A6),(u:$3FF007A1B865A8CA)),
    ((u:$3C7DA93F90835F76),(u:$3FF0085382FAEF83)),
    ((u:$BC86A79084AB093C),(u:$3FF00905554425D4)),
    ((u:$3C986364F8FBE8F8),(u:$3FF009B72F41A12B)),
    ((u:$BC882E8E14E3110E),(u:$3FF00A6910F3B6FD)),
    ((u:$BC84F6B2A7609F72),(u:$3FF00B1AFA5ABCBF)),
    ((u:$BC7E1A258EA8F71A),(u:$3FF00BCCEB7707EC)),
    ((u:$3C74362CA5BC26F2),(u:$3FF00C7EE448EE02)),
    ((u:$3C9095A56C919D02),(u:$3FF00D30E4D0C483)),
    ((u:$BC6406AC4E81A646),(u:$3FF00DE2ED0EE0F5)),
    ((u:$3C9B5A6902767E08),(u:$3FF00E94FD0398E0)),
    ((u:$BC991B2060859320),(u:$3FF00F4714AF41D3)),
    ((u:$3C8427068AB22306),(u:$3FF00FF93412315C)),
    ((u:$3C9C1D0660524E08),(u:$3FF010AB5B2CBD11)),
    ((u:$BC9E7BDFB3204BE8),(u:$3FF0115D89FF3A8B)),
    ((u:$3C8843AA8B9CBBC6),(u:$3FF0120FC089FF63)),
    ((u:$BC734104EE7EDAE8),(u:$3FF012C1FECD613B)),
    ((u:$BC72B6AEB6176892),(u:$3FF0137444C9B5B5)),
    ((u:$3C7A8CD33B8A1BB2),(u:$3FF01426927F5278)),
    ((u:$3C72EDC08E5DA99A),(u:$3FF014D8E7EE8D2F)),
    ((u:$3C857BA2DC7E0C72),(u:$3FF0158B4517BB88)),
    ((u:$3C9B61299AB8CDB8),(u:$3FF0163DA9FB3335)),
    ((u:$BC990565902C5F44),(u:$3FF016F0169949ED)),
    ((u:$3C870FC41C5C2D54),(u:$3FF017A28AF25567)),
    ((u:$3C94B9A6E145D76C),(u:$3FF018550706AB62)),
    ((u:$BC7008EFF5142BFA),(u:$3FF019078AD6A19F)),
    ((u:$BC977669F033C7DE),(u:$3FF019BA16628DE2)),
    ((u:$BC909BB78EEEAD0A),(u:$3FF01A6CA9AAC5F3)),
    ((u:$3C9371231477ECE6),(u:$3FF01B1F44AF9F9E)),
    ((u:$3C75E7626621EB5A),(u:$3FF01BD1E77170B4)),
    ((u:$BC9BC72B100828A4),(u:$3FF01C8491F08F08)),
    ((u:$BC6CE39CBBAB8BBE),(u:$3FF01D37442D5070)),
    ((u:$3C816996709DA2E2),(u:$3FF01DE9FE280AC8)),
    ((u:$BC8C11F5239BF536),(u:$3FF01E9CBFE113EF)),
    ((u:$3C8E1D4EB5EDC6B4),(u:$3FF01F4F8958C1C6)),
    ((u:$BC9AFB99946EE3F0),(u:$3FF020025A8F6A35)),
    ((u:$BC98F06D8A148A32),(u:$3FF020B533856324)),
    ((u:$BC82BF310FC54EB6),(u:$3FF02168143B0281)),
    ((u:$BC9C95A035EB4176),(u:$3FF0221AFCB09E3E)),
    ((u:$BC9491793E46834C),(u:$3FF022CDECE68C4F)),
    ((u:$BC73E8D0D9C49090),(u:$3FF02380E4DD22AD)),
    ((u:$BC9314AA16278AA4),(u:$3FF02433E494B755)),
    ((u:$3C848DAF888E9650),(u:$3FF024E6EC0DA046)),
    ((u:$3C856DC8046821F4),(u:$3FF02599FB483385)),
    ((u:$3C945B42356B9D46),(u:$3FF0264D1244C719)),
    ((u:$BC7082EF51B61D7E),(u:$3FF027003103B10E)),
    ((u:$3C72106ED0920A34),(u:$3FF027B357854772)),
    ((u:$BC9FD4CF26EA5D0E),(u:$3FF0286685C9E059)),
    ((u:$BC909F8775E78084),(u:$3FF02919BBD1D1D8)),
    ((u:$3C564CBBA902CA28),(u:$3FF029CCF99D720A)),
    ((u:$3C94383EF231D206),(u:$3FF02A803F2D170D)),
    ((u:$3C94A47A505B3A46),(u:$3FF02B338C811703)),
    ((u:$3C9E471202234680),(u:$3FF02BE6E199C811)));

  // as_tanh_database: 12-entry sorted-by-|x| table (x, f, dlo).
  cTanhDb: array[0..11, 0..2] of Tb64u64 = (
    ((u:$3FCAC343B179FEC4),(u:$3FCA612499C53078),(u:$3C60000000000000)),
    ((u:$3FD00764A988BF73),(u:$3FCF676484C0703B),(u:$B970000000000000)),
    ((u:$3FD17D1E8A63711F),(u:$3FD110E96A6C2D96),(u:$B970000000000000)),
    ((u:$3FD291C601A05276),(u:$3FD210B7D0C03743),(u:$3C70000000000000)),
    ((u:$3FD36F33D51C264D),(u:$3FD2DBB7B1C91363),(u:$B950000000000000)),
    ((u:$3FD43EAEA23649C3),(u:$3FD39877ED028641),(u:$BC90000000000000)),
    ((u:$3FDD88D7550B2826),(u:$3FDB9A3637366AFD),(u:$3C70000000000000)),
    ((u:$3FDE611AA58AB608),(u:$3FDC493DC899E4A6),(u:$BC90000000000000)),
    ((u:$3FE01EFE7AC8C15D),(u:$3FDDC3FE1B524821),(u:$B970000000000000)),
    ((u:$3FE1005EC0BCCABB),(u:$3FDF20B1C8557DED),(u:$BC90000000000000)),
    ((u:$3FE33DFEB0FA4BFE),(u:$3FE1372F9EE76E99),(u:$3C90000000000000)),
    ((u:$3FE49F24AC5CAC35),(u:$3FE22C495FF06104),(u:$B970000000000000)));

  // as_exp_accurate: ch[3][2]
  cTanhExpCh: array[0..2, 0..1] of Tb64u64 = (
    ((u:$3FF0000000000000),(u:$3A16C16BD194535D)),
    ((u:$3FE0000000000000),(u:$BA28259D904FD34F)),
    ((u:$3FC5555555555555),(u:$3C653E93E9F26E62)));
  // Inner exp polynomial seed coefficients (0x1.555...p-5 etc.)
  cTanhExpP0: Tb64u64 = (u:$3FA5555555555555);
  cTanhExpP1: Tb64u64 = (u:$3F811111113E93E9);
  cTanhExpP2: Tb64u64 = (u:$3F56C16C169400A7);
  // log(2) split for accurate path: l2h + l2l + l2ll
  cTanhL2hA:  Tb64u64 = (u:$3F262E42FF000000);   //  0x1.62e42ffp-13
  cTanhL2lA:  Tb64u64 = (u:$3D0718432A1B0E26);   //  0x1.718432a1b0e26p-47
  cTanhL2llA: Tb64u64 = (u:$3999FF0342542FC3);   //  0x1.9ff0342542fc3p-102

  // as_tanh_zero: ch[10][2]
  cTanhZeroCh: array[0..9, 0..1] of Tb64u64 = (
    ((u:$BFD5555555555555),(u:$BC75555555555555)),
    ((u:$3FC1111111111111),(u:$3C41111111110916)),
    ((u:$BFABA1BA1BA1BA1C),(u:$3C47917917A46F2C)),
    ((u:$3F9664F4882C10FA),(u:$BC09A52A06F1E599)),
    ((u:$BF8226E355E6C23D),(u:$3C2C297394C24E38)),
    ((u:$3F6D6D3D0E157DE0),(u:$BC0311087E5B1526)),
    ((u:$BF57DA36452B75E1),(u:$BBE2868CDE54EA0C)),
    ((u:$3F4355824803667B),(u:$3BD2CD8FC406C3F7)),
    ((u:$BF2F57D7734C821D),(u:$3B9DA22861B4CA80)),
    ((u:$3F1967E18AD3FACF),(u:$BBB0831108273A74)));

  // as_tanh_zero: cl[6]
  cTanhZeroCl: array[0..5] of Tb64u64 = (
    (u:$BF0497D8E6462927),(u:$3EF0B1318C243BD7),(u:$BEDB0F2935E9A120),
    (u:$3EC5E9444536E654),(u:$BEB174FF2A31908C),(u:$3E9749698C8D338D));

  // Medium-path polynomial in x2 (after x^3 factor) for |x| in [2^-30, 0.25)
  cTanhMedC: array[0..7] of Tb64u64 = (
    (u:$BFD5555555555554),(u:$3FC1111111110D61),(u:$BFABA1BA1B983D8B),(u:$3F9664F4820E99F0),
    (u:$BF8226E11E4AC7CF),(u:$3F6D6C4AB70668B6),(u:$BF57BBECB57CE996),(u:$3F41451443697DD8));

  // Large/medium exp polynomial ch[4]
  cTanhChOuter: array[0..3] of Tb64u64 = (
    (u:$4000000000000000),(u:$4000000000000000),(u:$3FF55555557E54FF),(u:$3FE55555553A12F4));

  cTanhSBig:   Tb64u64 = (u:$C0C71547652B82FE);   // -0x1.71547652b82fep+13
  cTanhMagic:  Tb64u64 = (u:$4188000004000000);   //  0x1.8000004p+25
  cTanhP25:    Tb64u64 = (u:$4188000000000000);   //  0x1.8p+25  (= 50331648)
  cTanhL2hM:   Tb64u64 = (u:$BF162E42FF000000);   // -0x1.62e42ffp-14
  cTanhL2lM:   Tb64u64 = (u:$BCF718432A1B0E26);   // -0x1.718432a1b0e26p-48
  cTanhL2L:    Tb64u64 = (u:$BF162E42FEFA39EF);   // -0x1.62e42fefa39efp-14 (large path)
  cTanhP55:    Tb64u64 = (u:$3C80000000000000);   //  0x1p-55
  cTanh1a52:   Tb64u64 = (u:$3CBA000000000000);   //  0x1.ap-52
  cTanh1p62:   Tb64u64 = (u:$3C10000000000000);   //  0x1p-62
  cTanh11p49:  Tb64u64 = (u:$3CE1000000000000);   //  0x1.1p-49
  cTanh1p3_55: Tb64u64 = (u:$3C83000000000000);   //  0x1.3p-55
  cTanhMask27: UInt64  = $FFFFFFFFF8000000;
  cTanhAixLge: UInt64  = $40330FC1931F09CA;
  cTanhAixMed: UInt64  = $400D76C8B4395810;
  cTanhAixSml: UInt64  = $3FD0000000000000;
  cTanhAixT30: UInt64  = $3E10000000000000;
  cTanhAixT32: UInt64  = $3DF0000000000000;

function TanhDatabase(x, f: Double): Double;
var
  a, b, m: Int32;
  ax, sgn: Double;
begin
  a := 0; b := 11;
  ax := Abs(x);
  m := (a + b) div 2;
  while a <= b do
  begin
    if cTanhDb[m, 0].f < ax then a := m + 1
    else if cTanhDb[m, 0].f = ax then
    begin
      if x >= Double(0.0) then sgn := Double(1.0) else sgn := -Double(1.0);
      f := sgn * cTanhDb[m, 1].f + sgn * cTanhDb[m, 2].f;
      Break;
    end
    else b := m - 1;
    m := (a + b) div 2;
  end;
  Result := f;
end;

function TanhExpAccurate(x, t, th, tl: Double; out l: Double): Double;
var
  dx, dxl, dxll, dxh, fl: Double;
  chp, clp, thp, tlp: Double;
  fh, zh, zl, uh, ul, vh, vl: Double;
begin
  dx   := x - cTanhL2hA.f * t;
  dxl  := cTanhL2lA.f * t;
  dxll := cTanhL2llA.f * t + pcr_fma(cTanhL2lA.f, t, -dxl);
  dxh  := dx + dxl;
  dxl  := ((dx - dxh) + dxl) + dxll;

  // Seed fl = dxh*(p0 + dxh*(p1 + dxh*p2))
  fl := dxh * (cTanhExpP0.f + dxh * (cTanhExpP1.f + dxh * cTanhExpP2.f));

  // polydd(dxh, dxl, 3, ch, &fl) — seeded
  chp := cTanhExpCh[2,0].f + fl;
  clp := ((cTanhExpCh[2,0].f - chp) + fl) + cTanhExpCh[2,1].f;
  // i = 1
  chp := pcr_muldd(dxh, dxl, chp, clp, clp);
  thp := chp + cTanhExpCh[1,0].f; tlp := (cTanhExpCh[1,0].f - thp) + chp;
  chp := thp; clp := clp + tlp + cTanhExpCh[1,1].f;
  // i = 0
  chp := pcr_muldd(dxh, dxl, chp, clp, clp);
  thp := chp + cTanhExpCh[0,0].f; tlp := (cTanhExpCh[0,0].f - thp) + chp;
  chp := thp; clp := clp + tlp + cTanhExpCh[0,1].f;

  fh := chp; fl := clp;
  fh := pcr_muldd(dxh, dxl, fh, fl, fl);
  fh := pcr_muldd(th,  tl,  fh, fl, fl);

  zh := th + fh; zl := (th - zh) + fh;
  uh := zh + tl; ul := ((zh - uh) + tl) + zl;
  vh := uh + fl; vl := ((uh - vh) + fl) + ul;
  l := vl;
  Result := vh;
end;

function TanhZero(x: Double): Double;
var
  x2, x2l, y0, y1, y2: Double;
  chp, clp, thp, tlp: Double;
  s_tmp, z_tmp: Double;
  i: Int32;
  t_u, w_u: Tb64u64;
begin
  x2  := x * x;
  x2l := pcr_fma(x, x, -x2);

  y2 := x2 * (cTanhZeroCl[0].f + x2 * (cTanhZeroCl[1].f + x2 * (cTanhZeroCl[2].f
        + x2 * (cTanhZeroCl[3].f + x2 * (cTanhZeroCl[4].f + x2 * cTanhZeroCl[5].f)))));

  // polydd (n=10) seeded with y2
  chp := cTanhZeroCh[9,0].f + y2;
  clp := ((cTanhZeroCh[9,0].f - chp) + y2) + cTanhZeroCh[9,1].f;
  for i := 8 downto 0 do
  begin
    chp := pcr_muldd(x2, x2l, chp, clp, clp);
    thp := chp + cTanhZeroCh[i,0].f;
    tlp := (cTanhZeroCh[i,0].f - thp) + chp;
    chp := thp;
    clp := clp + tlp + cTanhZeroCh[i,1].f;
  end;
  y1 := chp; y2 := clp;

  // y1 = mulddd(y1, y2, x, &y2)
  y1 := pcr_mulddd_pd(y1, y2, x, y2);
  // y1 = muldd_acc(y1, y2, x2, x2l, &y2)
  y1 := pcr_muldd(y1, y2, x2, x2l, y2);

  // y0 = fasttwosum(x, y1, &y1)
  pcr_fasttwosum(y0, y1, x, y1);
  // y1 = fasttwosum(y1, y2, &y2)
  s_tmp := y1 + y2;
  z_tmp := s_tmp - y1;
  y2 := y2 - z_tmp;
  y1 := s_tmp;

  t_u.f := y1;
  if (t_u.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) = 0 then
  begin
    w_u.f := y2;
    if ((w_u.u xor t_u.u) shr 63) <> 0 then
      Dec(t_u.u)
    else
      Inc(t_u.u);
    y1 := t_u.f;
    if y2 = Double(0.0) then
    begin
      Result := TanhDatabase(x, y0 + y1);
      Exit;
    end;
  end;
  Result := y0 + y1;
end;

function pcr_tanh(x: Double): Double;
var
  ax, v0, t, t0h, t1h, th, tl, t0l, t1l: Double;
  dx, dx2, p, rh, rl, e, lb, ub: Double;
  rqh, rql, ph, pl, qh, ql, qd, res: Double;
  x2, x3, x4, x8, p0, p1: Double;
  one, df, fsgn: Double;
  ix, jt, v_u, sp, lu: Tb64u64;
  aix: UInt64;
  i1, i0: Int32;
  ie: Int64;
  sh_sum, d_sum, sl_sum: Double;
begin
  ax := Abs(x);
  ix.f := ax;
  aix := ix.u;

  if aix >= cTanhAixLge then
  begin
    if aix > UInt64($7FF0000000000000) then begin Result := x + x; Exit; end;   // NaN
    if x >= Double(0.0) then fsgn := Double(1.0) else fsgn := -Double(1.0);
    if aix = UInt64($7FF0000000000000) then begin Result := fsgn; Exit; end;    // ±Inf
    if x >= Double(0.0) then df := cTanhP55.f else df := -cTanhP55.f;
    Result := fsgn - df;
    Exit;
  end;

  v0 := pcr_fma(ax, cTanhSBig.f, cTanhMagic.f);
  jt.f := v0;
  v_u.u := jt.u and cTanhMask27;
  t := v_u.f - cTanhP25.f;

  i1 := Int32((jt.u shr 27) and $3F);
  i0 := Int32((jt.u shr 33) and $3F);
  ie := SarInt64(Int64(jt.u shl 13), 52);
  sp.u := UInt64(Int64(1023) + ie) shl 52;

  t0h := cExpT0[i0, 1].f;
  t1h := cExpT1[i1, 1].f;
  th  := t0h * t1h;

  if aix < cTanhAixMed then
  begin
    // |x| < 3.683
    if aix < cTanhAixSml then
    begin
      // |x| < 0.25
      if aix < cTanhAixT30 then
      begin
        // |x| < 2^-30
        if aix < cTanhAixT32 then
        begin
          // |x| < 2^-32
          if aix = 0 then begin Result := x; Exit; end;
          Result := pcr_fma(x, -cTanhP55.f, x);
          Exit;
        end;
        x3 := x * x * x;
        Result := x - x3 / Double(3.0);
        Exit;
      end;

      x2 := x * x; x3 := x2 * x; x4 := x2 * x2; x8 := x4 * x4;
      p1 := (cTanhMedC[4].f + x2 * cTanhMedC[5].f)
            + x4 * (cTanhMedC[6].f + x2 * cTanhMedC[7].f);
      p0 := (cTanhMedC[0].f + x2 * cTanhMedC[1].f)
            + x4 * (cTanhMedC[2].f + x2 * cTanhMedC[3].f);
      p0 := p0 + x8 * p1;
      p0 := p0 * x3;

      pcr_fasttwosum(rh, rl, x, p0);
      e  := x3 * cTanh1a52.f;
      lb := rh + (rl - e);
      ub := rh + (rl + e);
      if lb = ub then begin Result := lb; Exit; end;
      Result := TanhZero(x);
      Exit;
    end;

    // 0.25 <= |x| < 3.683 — fast path with Ziv test
    t0l := cExpT0[i0, 0].f;
    t1l := cExpT1[i1, 0].f;
    tl  := t0h * t1l + t1h * t0l + pcr_fma(t0h, t1h, -th);
    th  := th * sp.f;
    tl  := tl * sp.f;

    dx  := (cTanhL2hM.f * t - ax) - cTanhL2lM.f * t;
    dx2 := dx * dx;
    p   := dx * ((cTanhChOuter[0].f + dx * cTanhChOuter[1].f)
                 + dx2 * (cTanhChOuter[2].f + dx * cTanhChOuter[3].f));
    rh  := th;
    rl  := tl + rh * p;
    pcr_fasttwosum(rh, rl, rh, rl);

    ph := rh; pl := rl;
    qh := rh; ql := rl;
    pcr_fasttwosum(qh, qd, Double(1.0), qh);
    ql := ql + qd;

    rqh := Double(1.0) / qh;
    rql := (ql * rqh + pcr_fma(rqh, qh, -Double(1.0))) * (-rqh);
    ph  := pcr_muldd(ph, pl, rqh, rql, pl);

    e  := rh * cTanh1p62.f;
    // fasttwosub(0.5, ph, &rl): s = 0.5 - ph; z = 0.5 - s; rl = z - ph;
    rh := Double(0.5) - ph;
    rl := (Double(0.5) - rh) - ph;
    rl := rl - pl;
    if x >= Double(0.0) then
    begin rh := rh * Double(2.0); rl := rl * Double(2.0); end
    else
    begin rh := rh * -Double(2.0); rl := rl * -Double(2.0); end;
    lb := rh + (rl - e);
    ub := rh + (rl + e);
    if lb = ub then begin Result := lb; Exit; end;
  end
  else
  begin
    // |x| >= 3.683 — fast fallback when tanh ≈ ±1
    dx  := pcr_fma(cTanhL2L.f, t, -ax);
    dx2 := dx * dx;
    p   := dx * ((cTanhChOuter[0].f + dx * cTanhChOuter[1].f)
                 + dx2 * (cTanhChOuter[2].f + dx * cTanhChOuter[3].f));
    rh  := th * sp.f;
    rh  := rh + (p + ((Double(2.0) * cTanh1p3_55.f) * ax)) * rh;
    e   := rh * cTanh11p49.f;
    rh  := (Double(2.0) * rh) / (Double(1.0) + rh);
    if x >= Double(0.0) then one := Double(1.0) else one := -Double(1.0);
    if x >= Double(0.0) then rh := rh else rh := -rh;
    lb  := one - (rh + e);
    ub  := one - (rh - e);
    if lb = ub then begin Result := lb; Exit; end;

    t0l := cExpT0[i0, 0].f;
    t1l := cExpT1[i1, 0].f;
    tl  := t0h * t1l + t1h * t0l + pcr_fma(t0h, t1h, -th);
    th  := th * sp.f;
    tl  := tl * sp.f;
  end;

  // Accurate (slow) path: shared by both medium and large branches
  rh := TanhExpAccurate(-Double(2.0) * ax, t, th, tl, rl);

  pcr_fasttwosum(qh, qd, Double(1.0), rh);
  ql := rl + qd;
  pcr_fasttwosum(qh, ql, qh, ql);

  rqh := Double(1.0) / qh;
  rql := (ql * rqh + pcr_fma(rqh, qh, -Double(1.0))) * (-rqh);
  ph  := pcr_muldd(rh, rl, rqh, rql, pl);

  // fasttwosub(0.5, ph, &rl)
  rh := Double(0.5) - ph;
  rl := (Double(0.5) - rh) - ph;
  rl := rl - pl;
  pcr_fasttwosum(rh, rl, rh, rl);

  if x >= Double(0.0) then
    res := Double(2.0) * rh + Double(2.0) * rl
  else
    res := -Double(2.0) * rh - Double(2.0) * rl;

  lu.f := rl;
  if (((lu.u + 32) and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) < 65) then
  begin
    Result := TanhDatabase(x, res);
    Exit;
  end;
  Result := res;
end;

// ---------------------------------------------------------------------------
// pcr_cospi — correctly-rounded binary64 cos(pi*x).
// Ported from core-math/src/binary64/cospi/cospi.c by Alexei Sibidanov.
// All helper routines (fasttwosum/muldd_acc/mulddd/polydd) inline; the
// accurate refinement is pure double-double (no TDInt64 path).
// Tables extracted programmatically from cospi.c (see tasklist64.md 1.07).
// ---------------------------------------------------------------------------

const
  // Main path constants sn[3], cn[2]
  cCospiSn: array[0..2] of Tb64u64 = (
    (u:$3B5921FB54442D18),(u:$B204ABBCE625BE51),(u:$289466BC6044BA16));
  cCospiCn: array[0..1] of Tb64u64 = (
    (u:$B6B3BD3CC9BE45DB),(u:$2D503C1F00186416));

  // as_cospi_zero tables
  cCospiZeroCh: array[0..1, 0..1] of Tb64u64 = (
    ((u:$C013BD3CC9BE45DE),(u:$BCB692B71366CC04)),
    ((u:$40103C1F081B5AC4),(u:$BCB32B33FDA9113C)));
  cCospiZeroCl: array[0..1] of Tb64u64 = (
    (u:$BFF55D3C7E3CBFF9),(u:$3FCE1F50604FA0FF));

  // Fast-path polynomial c[0..3] for |x| <= 2^-12 branch
  cCospiFastC: array[0..3] of Tb64u64 = (
    (u:$C013BD3CC9BE45DC),(u:$40103C1F081B0833),
    (u:$BFF55D3C6FC9AF15),(u:$3FCE1D3FF2AE3F9A));

  // as_sinpi_refine tables sh[3][2], ch[2][2]
  cSinpiRefSh: array[0..2, 0..1] of Tb64u64 = (
    ((u:$400921FB54442D18),(u:$3CA1A62633145C06)),
    ((u:$BE94ABBCE625BE53),(u:$3B305511CBC65743)),
    ((u:$3D0466BC6775AAE1),(u:$B8D9C3C168D990A0)));
  cSinpiRefCh: array[0..1, 0..1] of Tb64u64 = (
    ((u:$BE93BD3CC9BE45DE),(u:$BB3692B71366CC04)),
    ((u:$3D103C1F081B5AC4),(u:$B9B32B33FDA9113C)));

  // as_sinpi_refine scalar constants
  cSinpiRefSllK: Tb64u64 = (u:$BB632D2CC920DCB4); // -0x1.32d2cc920dcb4p-73
  cSinpiRefCll0: Tb64u64 = (u:$BB755D3C7E3CBFF9); // -0x1.55d3c7e3cbff9p-72
  cSinpiRefCll1: Tb64u64 = (u:$39CE1F50604FA0FF); //  0x1.e1f50604fa0ffp-99
  cSinpiRefScale: Tb64u64 = (u:$3C00000000000000); //  0x1p-63 (x = z*0x1p-63)
  cSinpiRefXlow: Tb64u64 = (u:$3F30000000000000); //  0x1p-12 (mulddd factor)
  cSinpiRefEr:   Tb64u64 = (u:$3840000000000000); //  0x1p-123

  // Fast-path thresholds
  cCospiAxT12:   UInt64  = $3F30000000000000;
  cCospiAxTiny:  UInt64  = $3E2CCF6429BE6621;
  cCospiEpsFast: Tb64u64 = (u:$3CFA000000000000); //  0x1.ap-48

  // 1ulp constants
  cCospiP55: Tb64u64 = (u:$3C80000000000000); //  0x1p-55

  // Database for as_sinpi_refine exceptions
  cSinpiDbIq: array[0..7] of Int32 = (903, 1029, 1078, 1217, 1025, 1026, 1033, 1235);
  cSinpiDbX:  array[0..7] of Tb64u64 = (
    (u:$BFDBDD02D1AD6000),(u:$BFCA4AD070549D00),(u:$3FCFBDB79CA3DA00),(u:$3FEC0CCEE4ADA200),
    (u:$3FDB536647B1FE98),(u:$BFDDC93EAAD12A18),(u:$3FEE78F0E592B360),(u:$3FDF13412A48D800));
  cSinpiDbR:  array[0..7] of Tb64u64 = (
    (u:$3FEF72C906962631),(u:$3FEFFFC4D2C6CA51),(u:$3FEFE3C8219054C3),(u:$3FEE99FD53791BCF),
    (u:$3FEFFFFC5DDD0738),(u:$3FEFFFF84B21C731),(u:$3FEFFF2270422604),(u:$3FEE55A7FA9A24C4));
  cSinpiDbD:  array[0..7] of Tb64u64 = (
    (u:$3C80000000000000),(u:$3C80000000000000),(u:$3C80000000000000),(u:$BC80000000000000),
    (u:$3C80000000000000),(u:$3C80000000000000),(u:$BC80000000000000),(u:$BC80000000000000));

  // sincosn tables (fast-path)
  cSincosN1_Sn: array[0..32, 0..1] of Tb64u64 = (
    ((u:$0000000000000000),(u:$0000000000000000)),
    ((u:$3FA91F6600000000),(u:$BE1DE44FD832257A)),
    ((u:$3FB917A6C0000000),(u:$BE0EB25EA0F138C7)),
    ((u:$3FC2C81060000000),(u:$3E3D1CC27444C003)),
    ((u:$3FC8F8B840000000),(u:$BE1CB2CFAA4DA337)),
    ((u:$3FCF19F980000000),(u:$BE237A839542DEEF)),
    ((u:$3FD2940630000000),(u:$BE12A60FA574A369)),
    ((u:$3FD58F9A70000000),(u:$3E36AC7F73F84090)),
    ((u:$3FD87DE2A0000000),(u:$3E3ABAA58B469891)),
    ((u:$3FDB5D1010000000),(u:$BE387A8CFF5264EA)),
    ((u:$3FDE2B5D40000000),(u:$BE3FE4271387C9DC)),
    ((u:$3FE0738798000000),(u:$3E222FFED9697FAF)),
    ((u:$3FE1C73B38000000),(u:$3E2AE68C86C9774A)),
    ((u:$3FE30FF800000000),(u:$BE38F47E58F7E631)),
    ((u:$3FE44CF328000000),(u:$BE37B7114F3FC4AF)),
    ((u:$3FE57D6938000000),(u:$BE3B989B02EAE413)),
    ((u:$3FE6A09E68000000),(u:$BE280C4336F74D05)),
    ((u:$3FE7B5DF20000000),(u:$3E33557D76F0AC85)),
    ((u:$3FE8BC8068000000),(u:$3E38A8BA05A743DA)),
    ((u:$3FE9B3E048000000),(u:$BDD8F17E98771434)),
    ((u:$3FEA9B6628000000),(u:$3E20EA1A3033EC62)),
    ((u:$3FEB728348000000),(u:$BE37348E1378D3E6)),
    ((u:$3FEC38B2F0000000),(u:$3E280BDB0D23E9D1)),
    ((u:$3FECED7AF8000000),(u:$BE3E19C46879EDAF)),
    ((u:$3FED906BD0000000),(u:$BE19AE573AEA067C)),
    ((u:$3FEE212108000000),(u:$BE384BC8DA0298EE)),
    ((u:$3FEE9F4158000000),(u:$BE239D225A27D387)),
    ((u:$3FEF0A7EF8000000),(u:$3E3C9186B952C7AE)),
    ((u:$3FEF6297D0000000),(u:$BDD1469FAA77A357)),
    ((u:$3FEFA75580000000),(u:$BE1EEB5D2BD05465)),
    ((u:$3FEFD88DA0000000),(u:$3E3E89292CF04139)),
    ((u:$3FEFF621E0000000),(u:$3E3BCB6BEF1D421F)),
    ((u:$3FF0000000000000),(u:$0000000000000000)));
  cSincosN1_Sm: array[0..32, 0..1] of Tb64u64 = (
    ((u:$0000000000000000),(u:$0000000000000000)),
    ((u:$3F59220000000000),(u:$BE354466E349EE53)),
    ((u:$3F6921F800000000),(u:$3E17D99497495D20)),
    ((u:$3F72D97800000000),(u:$3E017CCB5E27CB43)),
    ((u:$3F7921F000000000),(u:$3E2FCCE00E23572D)),
    ((u:$3F7F6A6400000000),(u:$3E3F9A2A3C5885CF)),
    ((u:$3F82D96C00000000),(u:$BE3E35ED1FA23D22)),
    ((u:$3F85FDA000000000),(u:$3E1BD602F04014C8)),
    ((u:$3F8921D200000000),(u:$BDD909C3DCCF0E28)),
    ((u:$3F8C460000000000),(u:$BE0E1B7526EBF9F2)),
    ((u:$3F8F6A2A00000000),(u:$BE32A8CD06AF94A2)),
    ((u:$3F91472700000000),(u:$3E0B5AE26618769E)),
    ((u:$3F92D93700000000),(u:$BE31073C40B25037)),
    ((u:$3F946B4400000000),(u:$BE3F80C5F9200607)),
    ((u:$3F95FD4D00000000),(u:$3E20FD5912EF3F57)),
    ((u:$3F978F5300000000),(u:$3E377727C10D8CA5)),
    ((u:$3F99215600000000),(u:$BE00B933040D8EB2)),
    ((u:$3F9AB35500000000),(u:$BE33ABEC0D92AE48)),
    ((u:$3F9C454F00000000),(u:$3E33394EC7229C28)),
    ((u:$3F9DD74600000000),(u:$BE3CE6D5319D653E)),
    ((u:$3F9F693700000000),(u:$3E28E8E7807F600B)),
    ((u:$3FA07D9200000000),(u:$BDC9EEA00D71246D)),
    ((u:$3FA1468600000000),(u:$BE325E9F40A5121A)),
    ((u:$3FA20F7700000000),(u:$3E19D6238D09F231)),
    ((u:$3FA2D86580000000),(u:$BE14D75465D2F213)),
    ((u:$3FA3A15100000000),(u:$BE137BCA0781D73A)),
    ((u:$3FA46A3980000000),(u:$BE20079E86EEC954)),
    ((u:$3FA5331F00000000),(u:$BE3E222F8A75DE87)),
    ((u:$3FA5FC0100000000),(u:$BE36B7995E4BB32D)),
    ((u:$3FA6C4DF80000000),(u:$BDF41523E139C56B)),
    ((u:$3FA78DBA80000000),(u:$3E32C3A342EB5F11)),
    ((u:$3FA8569200000000),(u:$3E35DA89E0B235C0)),
    ((u:$3FF0000000000000),(u:$0000000000000000)));
  cSincosN1_Cm: array[0..32, 0..1] of Tb64u64 = (
    ((u:$3FF0000000000000),(u:$0000000000000000)),
    ((u:$3FEFFFFD88000000),(u:$3E061BB991AF64F1)),
    ((u:$3FEFFFF620000000),(u:$3E2621D01D2A6063)),
    ((u:$3FEFFFE9C8000000),(u:$3E38F17465DE1773)),
    ((u:$3FEFFFD888000000),(u:$BE338BAB6D94C71D)),
    ((u:$3FEFFFC250000000),(u:$3E16BB5DD7625BCD)),
    ((u:$3FEFFFA730000000),(u:$BE3B439D8A459600)),
    ((u:$3FEFFF8718000000),(u:$3E237CE2EEC7251B)),
    ((u:$3FEFFF6218000000),(u:$BE2646D24A88970E)),
    ((u:$3FEFFF3828000000),(u:$BE39BB848BDB041E)),
    ((u:$3FEFFF0940000000),(u:$3E3E29DE85718CC2)),
    ((u:$3FEFFED570000000),(u:$3E3CC695B5E89E49)),
    ((u:$3FEFFE9CB8000000),(u:$BE3DA572F6A4BCCA)),
    ((u:$3FEFFE5F08000000),(u:$BE30D43929B71F74)),
    ((u:$3FEFFE1C68000000),(u:$3E0C32DDD89AA147)),
    ((u:$3FEFFDD4D8000000),(u:$3E3FBC7A9242CCF3)),
    ((u:$3FEFFD8860000000),(u:$3E1099A19765595D)),
    ((u:$3FEFFD36F8000000),(u:$BE2DBAFF3C93CDC4)),
    ((u:$3FEFFCE0A0000000),(u:$BE38EACC3AF7AC55)),
    ((u:$3FEFFC8558000000),(u:$BE3996F62D2DF41D)),
    ((u:$3FEFFC2520000000),(u:$BE3071603E8582DF)),
    ((u:$3FEFFBBFF8000000),(u:$3E07E5454B8225E9)),
    ((u:$3FEFFB55E8000000),(u:$BE3ED0128D24B027)),
    ((u:$3FEFFAE6E0000000),(u:$3E25569984BD1A7A)),
    ((u:$3FEFFA72F0000000),(u:$BDA08A362D33736D)),
    ((u:$3FEFF9FA10000000),(u:$3DFA441BA9901FD5)),
    ((u:$3FEFF97C40000000),(u:$3E304600A0A95596)),
    ((u:$3FEFF8F988000000),(u:$BE3387D3A589FA3D)),
    ((u:$3FEFF871D8000000),(u:$3E36DC0EF98B1C67)),
    ((u:$3FEFF7E540000000),(u:$3E301907C4C59658)),
    ((u:$3FEFF753B8000000),(u:$3E38DC8B1E83CCFF)),
    ((u:$3FEFF6BD48000000),(u:$BE2C4BBB2C348040)),
    ((u:$0000000000000000),(u:$0000000000000000)));

  // sincosn2 tables (accurate-path)
  cSincosN2_Sn: array[0..32, 0..1] of Tb64u64 = (
    ((u:$0000000000000000),(u:$0000000000000000)),
    ((u:$3FA91F65F10DD814),(u:$BC2912BD0D569A90)),
    ((u:$3FB917A6BC29B42C),(u:$BC3E2718D26ED688)),
    ((u:$3FC2C8106E8E613A),(u:$3C513000A89A11E0)),
    ((u:$3FC8F8B83C69A60B),(u:$BC626D19B9FF8D82)),
    ((u:$3FCF19F97B215F1B),(u:$BC642DEEF11DA2C4)),
    ((u:$3FD294062ED59F06),(u:$BC75D28DA2C4612D)),
    ((u:$3FD58F9A75AB1FDD),(u:$BC1EFDC0D58CF620)),
    ((u:$3FD87DE2A6AEA963),(u:$BC672CEDD3D5A610)),
    ((u:$3FDB5D1009E15CC0),(u:$3C65B362CB974183)),
    ((u:$3FDE2B5D3806F63B),(u:$3C5E0D891D3C6841)),
    ((u:$3FE073879922FFEE),(u:$BC8A5A014347406C)),
    ((u:$3FE1C73B39AE68C8),(u:$3C8B25DD267F6600)),
    ((u:$3FE30FF7FCE17035),(u:$BC6EFCC626F74A6F)),
    ((u:$3FE44CF325091DD6),(u:$3C68076A2CFDC6B3)),
    ((u:$3FE57D69348CECA0),(u:$BC875720992BFBB2)),
    ((u:$3FE6A09E667F3BCD),(u:$BC8BDD3413B26456)),
    ((u:$3FE7B5DF226AAFAF),(u:$BC70F537ACDF0AD7)),
    ((u:$3FE8BC806B151741),(u:$BC82C5E12ED1336D)),
    ((u:$3FE9B3E047F38741),(u:$BC830EE286712474)),
    ((u:$3FEA9B66290EA1A3),(u:$3C39F630E8B6DAC8)),
    ((u:$3FEB728345196E3E),(u:$BC8BC69F324E6D61)),
    ((u:$3FEC38B2F180BDB1),(u:$BC76E0B1757C8D07)),
    ((u:$3FECED7AF43CC773),(u:$BC5E7B6BB5AB58AE)),
    ((u:$3FED906BCF328D46),(u:$3C7457E610231AC2)),
    ((u:$3FEE212104F686E5),(u:$BC8014C76C126527)),
    ((u:$3FEE9F4156C62DDA),(u:$3C8760B1E2E3F81E)),
    ((u:$3FEF0A7EFB9230D7),(u:$3C752C7ADC6B4989)),
    ((u:$3FEF6297CFF75CB0),(u:$3C7562172A361FD3)),
    ((u:$3FEFA7557F08A517),(u:$BC87A0A8CA13571F)),
    ((u:$3FEFD88DA3D12526),(u:$BC887DF6378811C7)),
    ((u:$3FEFF621E3796D7E),(u:$BC6C57BC2E24AA15)),
    ((u:$3FF0000000000000),(u:$0000000000000000)));
  cSincosN2_Sm: array[0..31, 0..1] of Tb64u64 = (
    ((u:$0000000000000000),(u:$0000000000000000)),
    ((u:$3F5921FAAEE6472E),(u:$BBFEE52E284A9DF8)),
    ((u:$3F6921F8BECCA4BA),(u:$3C02BA407BCAB5B2)),
    ((u:$3F72D97822F996BC),(u:$3C13E5A15ED6AA3E)),
    ((u:$3F7921F0FE670071),(u:$3BFAB967FE6B7A9B)),
    ((u:$3F7F6A65F9A2A3C6),(u:$BC1DE8C48783F3AE)),
    ((u:$3F82D96B0E509703),(u:$BC01E9131FF52DC9)),
    ((u:$3F85FDA037AC05E1),(u:$BC2FF59BF4B574EE)),
    ((u:$3F8921D1FCDEC784),(u:$3C29878EBE836D9D)),
    ((u:$3F8C45FFE1E48AD9),(u:$3C04060E4BD32E79)),
    ((u:$3F8F6A296AB997CB),(u:$BC2F2943D8FE7033)),
    ((u:$3F9147270DAD7133),(u:$3C08769E00E01800)),
    ((u:$3F92D936BBE30EFD),(u:$3C2B5F91EE371D64)),
    ((u:$3F946B4381FCE81B),(u:$3C3FF9F89FB65BE3)),
    ((u:$3F95FD4D21FAB226),(u:$BC20C0A91C37851C)),
    ((u:$3F978F535DDC9F04),(u:$3C2B194AD9B1AA97)),
    ((u:$3F992155F7A3667E),(u:$BBFB1D63091A0130)),
    ((u:$3F9AB354B1504FCA),(u:$BC32AE47937CBDD3)),
    ((u:$3F9C454F4CE53B1D),(u:$BC3D63D7FEF0E36C)),
    ((u:$3F9DD7458C64AB3A),(u:$BC3D653DF3FCC281)),
    ((u:$3F9F693731D1CF01),(u:$BBD3FE9BC66286C7)),
    ((u:$3FA07D91FF984580),(u:$BC3AE248DA7A9007)),
    ((u:$3FA14685DB42C17F),(u:$BC42890D277CB974)),
    ((u:$3FA20F770CEB11C7),(u:$BC4EC1B9D46693B1)),
    ((u:$3FA2D865759455CD),(u:$3C2686F65BA93AC0)),
    ((u:$3FA3A150F6421AFC),(u:$3C3F8A318BA775FD)),
    ((u:$3FA46A396FF86179),(u:$3C2136AC00FA2DA9)),
    ((u:$3FA5331EC3BBA0EB),(u:$3C2442F2A9DAC128)),
    ((u:$3FA5FC00D290CD43),(u:$3C4A2669A693A8E1)),
    ((u:$3FA6C4DF7D7D5B84),(u:$BC339C56A9BD0A9B)),
    ((u:$3FA78DBAA5874686),(u:$BC34A0EF4035C29C)),
    ((u:$3FA856922BB513C1),(u:$3C491ADFD607CB2B)));
  cSincosN2_Cm: array[0..31, 0..1] of Tb64u64 = (
    ((u:$0000000000000000),(u:$0000000000000000)),
    ((u:$3EB3BD3C88CDCA13),(u:$3B5874628D2B6835)),
    ((u:$3ED3BD3BC5FC5AB4),(u:$BB48A1BEBF665CEF)),
    ((u:$3EE634E1D173443D),(u:$3B6198FF5804D8DC)),
    ((u:$3EF3BD38BAB6D94C),(u:$3B9C73BE2184804E)),
    ((u:$3EFED7A51288A277),(u:$BB9BCD47BF19555B)),
    ((u:$3F0634DA1CEC522D),(u:$BBA3FFF35C1BDB61)),
    ((u:$3F0E39B20C7444E3),(u:$3BAAE49F0BE31E8C)),
    ((u:$3F13BD2C8DA49511),(u:$3BA70DF810BCC0E2)),
    ((u:$3F18FB66EE122F6C),(u:$3B9079FF34BFFA7E)),
    ((u:$3F1ED7875885EA3A),(u:$BBA983278BAA11C2)),
    ((u:$3F22A8C672D4942F),(u:$BBBE4938F661AA14)),
    ((u:$3F2634BB4AE5ED49),(u:$3BCE64C904CC8156)),
    ((u:$3F2A0FA1A872536E),(u:$3BBF7469ECC1331F)),
    ((u:$3F2E3978F34889D9),(u:$3BC5EB897CDC1B0D)),
    ((u:$3F31592043856DBD),(u:$3BC9865D9CDC3744)),
    ((u:$3F33BCFBD9979A27),(u:$BBD595D548D9A586)),
    ((u:$3F36484EDD7F9E4A),(u:$BBB91DF31AAA6F7F)),
    ((u:$3F38FB18EACC3AF8),(u:$BBD4EAA508EEC2B7)),
    ((u:$3F3BD55996F62D2E),(u:$BBA7C51F7D8F9F71)),
    ((u:$3F3ED71071603E86),(u:$BBDF4827CCB50B62)),
    ((u:$3F41001E81ABAB48),(u:$BBD12F4A2E3616D5)),
    ((u:$3F42A86F68094692),(u:$3BE604E1F8F76C27)),
    ((u:$3F44647AAA599ED1),(u:$BBE1A79A86FB8FF5)),
    ((u:$3F46344004228D8B),(u:$3BE33736C96557C9)),
    ((u:$3F4817BF2DDF22B3),(u:$3BEFC0555A81D729)),
    ((u:$3F4A0EF7DCFFAFAB),(u:$3BE54D46B817BCA4)),
    ((u:$3F4C19E9C3E9D2C5),(u:$BB970A980659E790)),
    ((u:$3F4E389491F8833A),(u:$3BEC7313BEEAB883)),
    ((u:$3F50357BF9BE0ECF),(u:$BBF965827F33D907)),
    ((u:$3F515889C8DD385F),(u:$3BC98094FAB77B42)),
    ((u:$3F52857389776587),(u:$BBFBFDFCE09AEA7B)));

// ---------------------------------------------------------------------------
// CospiSincosN — fast-path sin/cos table lookup (port of sincosn in cospi.c)
// ---------------------------------------------------------------------------
procedure CospiSincosN(s: Int32; out sh, sl, ch, cl: Double);
var
  j, is_, ic, jm: Int32;
  ss, sc: Int32;
  sbh, sbl, cbh, cbl: Double;
  slh, sll, clh, cll: Double;
  sb, cb: Double;
  Ch_, Cl_, Sh_, Sl_: Double;
  tch, tcl, tsh, tsl: Double;
begin
  j := s and $3FF;
  if ((s shr 10) and 1) <> 0 then j := 1024 - j;
  is_ := j shr 5;
  ic  := $20 - is_;
  jm  := j and $1F;
  ss  := (s shr 11) and 1;
  sc  := (UInt32(s + 1024) shr 11) and 1;

  sbh := cSincosN1_Sn[is_, 0].f; sbl := cSincosN1_Sn[is_, 1].f;
  cbh := cSincosN1_Sn[ic,  0].f; cbl := cSincosN1_Sn[ic,  1].f;
  slh := cSincosN1_Sm[jm,  0].f; sll := cSincosN1_Sm[jm,  1].f;
  clh := cSincosN1_Cm[jm,  0].f; cll := cSincosN1_Cm[jm,  1].f;

  sb := sbh + sbl; cb := cbh + cbl;
  Ch_ := cbh*clh - sbh*slh;
  Cl_ := clh*cbl - slh*sbl + cb*cll - sb*sll;
  Sh_ := sbh*clh + cbh*slh;
  Sl_ := slh*cbl + clh*sbl + cb*sll + sb*cll;

  tch := Ch_ + Cl_; tcl := (Ch_ - tch) + Cl_;
  tsh := Sh_ + Sl_; tsl := (Sh_ - tsh) + Sl_;

  if sc = 0 then begin ch := tch; cl := tcl; end
  else             begin ch := -tch; cl := -tcl; end;
  if ss = 0 then begin sh := tsh; sl := tsl; end
  else             begin sh := -tsh; sl := -tsl; end;
end;

// ---------------------------------------------------------------------------
// CospiSincosN2 — accurate-path sin/cos (port of sincosn2 in cospi.c)
// ---------------------------------------------------------------------------
procedure CospiSincosN2(s: Int32; out sh, sl, ch, cl: Double);
var
  j, is_, ic, jm: Int32;
  ss, sc: Int32;
  sbh, sbl, cbh, cbl: Double;
  slh, sll, clh, cll: Double;
  cch, ccl, ssh, ssl, csh, csl, sch, scl: Double;
  tch, tcl, tsh, tsl, tcl2, tsl2: Double;
begin
  j := s and $3FF;
  if ((s shr 10) and 1) <> 0 then j := 1024 - j;
  is_ := j shr 5;
  ic  := $20 - is_;
  jm  := j and $1F;
  ss  := (s shr 11) and 1;
  sc  := (UInt32(s + 1024) shr 11) and 1;

  sbh := cSincosN2_Sn[is_, 0].f; sbl := cSincosN2_Sn[is_, 1].f;
  cbh := cSincosN2_Sn[ic,  0].f; cbl := cSincosN2_Sn[ic,  1].f;
  slh := cSincosN2_Sm[jm,  0].f; sll := cSincosN2_Sm[jm,  1].f;
  clh := cSincosN2_Cm[jm,  0].f; cll := cSincosN2_Cm[jm,  1].f;

  cch := pcr_muldd(clh, cll, cbh, cbl, ccl);
  ssh := pcr_muldd(slh, sll, sbh, sbl, ssl);
  csh := pcr_muldd(clh, cll, sbh, sbl, csl);
  sch := pcr_muldd(slh, sll, cbh, cbl, scl);

  // tch = fasttwosum(ssh, cch, &tcl); tcl += ccl + ssl;
  pcr_fasttwosum(tch, tcl, ssh, cch);
  tcl := tcl + ccl + ssl;
  // tsh = fasttwosum(-sch, csh, &tsl); tsl += csl - scl;
  pcr_fasttwosum(tsh, tsl, -sch, csh);
  tsl := tsl + csl - scl;

  // tch = fasttwosum(cbh, -tch, &tcl2); tcl = cbl - tcl + tcl2;
  pcr_fasttwosum(tch, tcl2, cbh, -tch);
  tcl := cbl - tcl + tcl2;
  // tsh = fasttwosum(sbh, -tsh, &tsl2); tsl = sbl - tsl + tsl2;
  pcr_fasttwosum(tsh, tsl2, sbh, -tsh);
  tsl := sbl - tsl + tsl2;

  if sc = 0 then begin ch := tch; cl := tcl; end
  else             begin ch := -tch; cl := -tcl; end;
  if ss = 0 then begin sh := tsh; sl := tsl; end
  else             begin sh := -tsh; sl := -tsl; end;
end;

// ---------------------------------------------------------------------------
// CospiAsZero — as_cospi_zero refinement (|x| near 0 slow path)
// ---------------------------------------------------------------------------
function CospiAsZero(x: Double): Double;
var
  x2, dx2, fl, fh: Double;
  chp, clp, thp, tlp: Double;
  y0, y1, y2: Double;
  s_tmp, z_tmp: Double;
  t_u, w_u: Tb64u64;
begin
  x2  := x * x;
  dx2 := pcr_fma(x, x, -x2);

  fl := x2 * (cCospiZeroCl[0].f + x2 * cCospiZeroCl[1].f);

  // polydd(x2, dx2, 2, ch, &fl) seeded
  chp := cCospiZeroCh[1,0].f + fl;
  clp := ((cCospiZeroCh[1,0].f - chp) + fl) + cCospiZeroCh[1,1].f;
  // i = 0
  chp := pcr_muldd(x2, dx2, chp, clp, clp);
  thp := chp + cCospiZeroCh[0,0].f;
  tlp := (cCospiZeroCh[0,0].f - thp) + chp;
  chp := thp;
  clp := clp + tlp + cCospiZeroCh[0,1].f;

  fh := chp; fl := clp;
  fh := pcr_muldd(x2, dx2, fh, fl, fl);

  // y0 = fasttwosum(1, fh, &y1);
  pcr_fasttwosum(y0, y1, Double(1.0), fh);
  // y1 = fasttwosum(y1, fl, &y2); — inlined, may violate |y1|>=|fl|? fasttwosum form
  s_tmp := y1 + fl;
  z_tmp := s_tmp - y1;
  y2 := fl - z_tmp;
  y1 := s_tmp;

  t_u.f := y1;
  if (t_u.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) = 0 then
  begin
    w_u.f := y2;
    if ((w_u.u xor t_u.u) shr 63) <> 0 then
      Dec(t_u.u)
    else
      Inc(t_u.u);
    y1 := t_u.f;
  end;
  Result := y0 + y1;
end;

// ---------------------------------------------------------------------------
// CospiSinpiRefine — as_sinpi_refine (accurate path for large args)
// ---------------------------------------------------------------------------
function CospiSinpiRefine(iq: Int32; z: Double): Double;
var
  x, x2, dx2: Double;
  sll, slh, cll, clh: Double;
  chp, clp, thp, tlp: Double;
  sbh, sbl, cbh, cbl: Double;
  csh, csl, sch, scl: Double;
  tsh, tsl, tsl2: Double;
  iqm: Int32;
  sgn, dbx: Double;
  i: Int32;
  t_u: Tb64u64;
begin
  x   := z * cSinpiRefScale.f;              // x = z * 2^-63
  x2  := x * x;
  dx2 := pcr_fma(x, x, -x2);

  // sll = -0x1.32d2cc920dcb4p-73 * x2
  sll := cSinpiRefSllK.f * x2;

  // polydd(x2, dx2, 3, sh, &sll) seeded
  chp := cSinpiRefSh[2,0].f + sll;
  clp := ((cSinpiRefSh[2,0].f - chp) + sll) + cSinpiRefSh[2,1].f;
  // i = 1
  chp := pcr_muldd(x2, dx2, chp, clp, clp);
  thp := chp + cSinpiRefSh[1,0].f;
  tlp := (cSinpiRefSh[1,0].f - thp) + chp;
  chp := thp;
  clp := clp + tlp + cSinpiRefSh[1,1].f;
  // i = 0
  chp := pcr_muldd(x2, dx2, chp, clp, clp);
  thp := chp + cSinpiRefSh[0,0].f;
  tlp := (cSinpiRefSh[0,0].f - thp) + chp;
  chp := thp;
  clp := clp + tlp + cSinpiRefSh[0,1].f;
  slh := chp; sll := clp;

  // slh = mulddd(slh, sll, x*0x1p-12, &sll)
  slh := pcr_mulddd_pd(slh, sll, x * cSinpiRefXlow.f, sll);

  // cll = x2*(-0x1.55d3c7e3cbff9p-72 + 0x1.e1f50604fa0ffp-99 * x2)
  cll := x2 * (cSinpiRefCll0.f + cSinpiRefCll1.f * x2);

  // polydd(x2, dx2, 2, ch, &cll) seeded
  chp := cSinpiRefCh[1,0].f + cll;
  clp := ((cSinpiRefCh[1,0].f - chp) + cll) + cSinpiRefCh[1,1].f;
  // i = 0
  chp := pcr_muldd(x2, dx2, chp, clp, clp);
  thp := chp + cSinpiRefCh[0,0].f;
  tlp := (cSinpiRefCh[0,0].f - thp) + chp;
  chp := thp;
  clp := clp + tlp + cSinpiRefCh[0,1].f;
  clh := chp; cll := clp;

  // clh = muldd_acc(clh, cll, x2, dx2, &cll)
  clh := pcr_muldd(clh, cll, x2, dx2, cll);

  // sincosn2(iq, &sbh, &sbl, &cbh, &cbl)
  CospiSincosN2(iq, sbh, sbl, cbh, cbl);

  // csh = muldd_acc(clh, cll, sbh, sbl, &csl)
  csh := pcr_muldd(clh, cll, sbh, sbl, csl);
  // sch = muldd_acc(slh, sll, cbh, cbl, &scl)
  sch := pcr_muldd(slh, sll, cbh, cbl, scl);
  // tsh = fasttwosum(sch, csh, &tsl); tsl += csl + scl
  pcr_fasttwosum(tsh, tsl, sch, csh);
  tsl := tsl + csl + scl;
  // tsh = fasttwosum(sbh, tsh, &tsl2); tsl = sbl + tsl + tsl2
  pcr_fasttwosum(tsh, tsl2, sbh, tsh);
  tsl := sbl + tsl + tsl2;

  // Exception database probe
  t_u.f := tsl;
  if ((t_u.u or (UInt64($FFF) shl 52)) = UInt64($FFFFFFFFFFFFFFFF))
     or ((t_u.u shl 12) = 0) then
  begin
    if iq > 2048 then sgn := -Double(1.0) else sgn := Double(1.0);
    iqm := iq and $7FF;
    for i := 0 to 7 do
    begin
      dbx := cSinpiDbX[i].f;
      if ((x = dbx) and (iqm = cSinpiDbIq[i]))
         or ((x = -dbx) and (iqm = 2048 - cSinpiDbIq[i])) then
      begin
        Result := sgn * cSinpiDbR[i].f + sgn * cSinpiDbD[i].f;
        Exit;
      end;
    end;
  end;
  Result := tsh + tsl;
end;

// ---------------------------------------------------------------------------
// pcr_cospi — entry point
// ---------------------------------------------------------------------------
function pcr_cospi(x: Double): Double;
var
  ix: Tb64u64;
  ax: UInt64;
  e, s_shift, si: Int32;
  m: Int64;
  iq_u: UInt64;
  iq: Int32;
  x2, x4, eps, p, lb, ub: Double;
  z, z2, fs, fc, er, r: Double;
  sh, sl, ch_, cl_: Double;
  k: Int64;
begin
  ix.f := x;
  ax   := ix.u and (UInt64($FFFFFFFFFFFFFFFF) shr 1);

  if ax = 0 then begin Result := Double(1.0); Exit; end;

  e := Int32(ax shr 52);
  m := Int64((ix.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) or (UInt64(1) shl 52));
  s_shift := 1063 - e;

  if s_shift < 0 then
  begin
    // |x| >= 2^41
    if e = $7FF then
    begin
      if (ix.u shl 12) = 0 then
      begin
        Result := cNaNDoublePos.f;
        Exit;
      end;
      Result := x + x; // NaN
      Exit;
    end;
    s_shift := -s_shift - 1; // 2^(41+s) <= |x| < 2^(42+s)
    if s_shift > 11 then begin Result := Double(1.0); Exit; end;
    iq_u := (UInt64(m) shl s_shift) + 1024;
    if (iq_u and 2047) = 0 then begin Result := Double(0.0); Exit; end;
    CospiSincosN(Int32(iq_u), sh, sl, ch_, cl_);
    Result := sh + sl;
    Exit;
  end;

  if ax <= cCospiAxT12 then
  begin
    // |x| <= 2^-12
    if ax <= cCospiAxTiny then
    begin
      Result := Double(1.0) - cCospiP55.f;
      Exit;
    end;
    x2 := x * x; x4 := x2 * x2;
    eps := x2 * cCospiEpsFast.f;
    p := x2 * ((cCospiFastC[0].f + x2 * cCospiFastC[1].f)
               + x4 * (cCospiFastC[2].f + x2 * cCospiFastC[3].f));
    lb := (p - eps) + Double(1.0);
    ub := (p + eps) + Double(1.0);
    if lb = ub then begin Result := lb; Exit; end;
    Result := CospiAsZero(x);
    Exit;
  end;

  // exact-zero detection: si = e - 1011; if (m << si) == 2^63 return 0.0
  si := e - 1011;
  if (si >= 0) and (si < 64)
     and ((UInt64(m) shl si) = UInt64($8000000000000000)) then
  begin
    Result := Double(0.0);
    Exit;
  end;

  iq_u := (UInt64(Int64(m) shr s_shift) + 2048) and 8191;
  iq_u := (iq_u + 1) shr 1;
  iq := Int32(iq_u);

  k := Int64(UInt64(m) shl (e - 1000));
  z := k; z2 := z * z;
  fs := cCospiSn[0].f + z2 * (cCospiSn[1].f + z2 * cCospiSn[2].f);
  fc := cCospiCn[0].f + z2 * cCospiCn[1].f;
  CospiSincosN(iq, sh, sl, ch_, cl_);
  er := z * cSinpiRefEr.f;  // z * 2^-123
  r  := sl + sh * (z2 * fc) + ch_ * (z * fs);
  lb := (r - er) + sh;
  ub := (r + er) + sh;
  if lb = ub then begin Result := lb; Exit; end;

  Result := CospiSinpiRefine(iq, z);
end;

// ---------------------------------------------------------------------------
// pcr_asin — correctly-rounded binary64 arc sine.
// Ported from core-math/src/binary64/asin/asin.c by Alexei Sibidanov.
// Reuses cAcosCC / cAcosSHi / cAcosSLo / cAcosCHi / cAcosCLo / cAcosCt:
// asin and acos share the same 33×8 polynomial table and the same
// 33-entry sin(pi/64*j) double-double table.
// ---------------------------------------------------------------------------

const
  cAsinSmallTh:  Tb64u64 = (u: $7CAE26E892247DEC); // ax<this → fma(2^-55,x,x)
  cAsinSmallC:   Tb64u64 = (u: $3C80000000000000); // 0x1p-55
  cAsinEps1:     Tb64u64 = (u: $3CB9620000000000); // 0x1.962p-52
  cAsinEps2:     Tb64u64 = (u: $3960000000000000); // 0x1p-100
  cAsinSignsMask: UInt32 = $1F73FFCB;

  // 29-entry rare-input database (xdb, ydb), per as_asin_database in asin.c.
  cAsinDbX: array[0..28] of UInt64 = (
    $3E57137449123EF6, $3E5D12ED0AF1A27E, $3E851C4B960778F5, $3E93CFC2A006A414,
    $3E9CBAA95DADB559, $3EBACD69F89AD8F1, $3EF2BFFFFFFC233B, $3EFFF0F3022B2E9D,
    $3F13217783D70D1D, $3F1C373FF4AAD79B, $3F6B3F28593CAD2F, $3F8E17B3F6BB5E6E,
    $3F941D60A76A82ED, $3F9921C0A0486537, $3F9A3A2919D6B19B, $3F9D6315F7EE7E01,
    $3F9EA6FDC56FC61A, $3FA69768DC89BB00, $3FAA4816B2066707, $3FAD77B117F230D6,
    $3FAFC7A07B2549AA, $3FB2DF0542154F1B, $3FB51CF5DB1B1956, $3FC9697CB602C582,
    $3FCD0EF799001BA9, $3FD4A8E1A96E38E3, $3FDDA4E0E6C717A5, $3FDEA8E8FDF47549,
    $3FE3B9994ABB81D4);
  cAsinDbY: array[0..28] of Tb64u64 = (
    (u:$3E57137449123EF7),(u:$3E5D12ED0AF1A27F),(u:$3E851C4B9607790D),(u:$3E93CFC2A006A465),
    (u:$3E9CBAA95DADB650),(u:$3EBACD69F89AE57A),(u:$3EF2C00000006DDD),(u:$3EFFF0F3024065E7),
    (u:$3F132177841FFBEF),(u:$3F1C373FF594D65B),(u:$3F6B3F2BA40DBC66),(u:$3F8E17FAEFAC7797),
    (u:$3F941DB571D96126),(u:$3F99226605233224),(u:$3F9A3AE514C9BEFA),(u:$3F9D641E6D5E769A),
    (u:$3F9EA829E3E988E5),(u:$3FA69949B3D51FB1),(u:$3FAA4B0BFB0454D5),(u:$3FAD7BDCD778049F),
    (u:$3FAFCCDC252CAD1F),(u:$3FB2E36813A98740),(u:$3FB5231B416BA885),(u:$3FC994FFB5DAF0F9),
    (u:$3FCD5064E6FE82C5),(u:$3FD50954B7BBF87B),(u:$3FDED25C5EB8C916),(u:$3FDFF92A8CA216CD),
    (u:$3FE540E24E5F33F3));

function AsinDatabase(x, f: Double): Double;
var
  t: Tb64u64;
  ax: UInt64;
  a, b, m: Int32;
  yt: Double;
begin
  Result := f;
  t.f := x;
  ax := t.u and (UInt64($7FFFFFFFFFFFFFFF));
  a := 0; b := High(cAsinDbX); m := (a + b) div 2;
  while a <= b do
  begin
    if cAsinDbX[m] < ax then a := m + 1
    else if cAsinDbX[m] = ax then
    begin
      // t.f = ydb[m]; t.u -= 54<<52; t.u |= ((signs>>m)&1) << 63
      t.u := cAsinDbY[m].u - (UInt64(54) shl 52);
      t.u := t.u or ((UInt64((cAsinSignsMask shr m) and 1)) shl 63);
      yt := cAsinDbY[m].f;        // ydb[m] is always positive
      if x >= Double(0.0) then
        Result :=  yt + t.f
      else
        Result := -yt - t.f;
      Exit;
    end
    else b := m - 1;
    m := (a + b) div 2;
  end;
end;

function AsinRefine(x, phi: Double): Double;
var
  s2, dx2, c2h, c2l, c2f, ch, cl: Double;
  jf: Int64;
  Ch_v, Cl_v, Sh_v, Sl_v, ax_r: Double;
  dsh, dsl, dch, dcl: Double;
  Sc, dSc, Cs, dCs, v, dv: Double;
  sgn: Double;
  jt: Int64;
  jtd: Double;
  v2, dv2: Double;
  fh, fl: Double;
  chp, clp, th, tl: Double;
  ph, pl, ps: Double;
  pl_new, ps_inner: Double;
  th_u, tl_u, tn_u: Tb64u64;
  dn, de: Int64;
  hard: Boolean;
  res: Double;
begin
  s2  := x * x;
  dx2 := pcr_fma(x, x, -s2);
  pcr_fasttwosum(c2h, c2l, Double(1.0), -s2);   // == fasttwosub(1, s2)
  c2l := c2l - dx2;
  pcr_fasttwosum(c2h, c2l, c2h, c2l);

  c2f := pcr_fma(x, -x, Double(1.0));
  ch  := Sqrt(c2f);
  cl  := (c2l - pcr_fma(ch, ch, -c2f)) * ((Double(0.5) / c2f) * ch);

  jf := Trunc(pcr_roundeven(Abs(phi) * cAcosRefScale.f));

  Ch_v := cAcosSHi[32 - jf].f;  Cl_v := cAcosSLo[32 - jf].f;
  Sh_v := cAcosSHi[jf].f;       Sl_v := cAcosSLo[jf].f;

  ax_r := Abs(x);
  dsh := ax_r - Sh_v;  dsl := -Sl_v;
  dch := ch   - Ch_v;  dcl := cl - Cl_v;

  Sc  := pcr_fma(Sh_v, dch, cAcosC2fK.f) - cAcosC2fK.f;
  dSc := pcr_fma(Sh_v, dch, -Sc);

  Cs  := pcr_fma(Ch_v, dsh, cAcosC2fK.f) - cAcosC2fK.f;
  dCs := pcr_fma(Ch_v, dsh, -Cs);

  v  := Cs - Sc;
  dv := (Ch_v * dsl + Cl_v * dsh) - (Sh_v * dcl + Sl_v * dch) - (dSc - dCs);
  pcr_fasttwosum(v, dv, v, dv);

  if x >= Double(0.0) then sgn := Double(1.0) else sgn := -Double(1.0);
  if x >= Double(0.0) then jt := jf else jt := -jf;     // jf*sgn, range [-32,32]
  jtd := jt;

  v2 := pcr_muldd(v, dv, v, dv, dv2);
  v  := v  * sgn;
  dv := dv * sgn;

  fl := v2 * (cAcosCt[0].f + v2 * (cAcosCt[1].f + v2 * cAcosCt[2].f));
  // fh = polydd(v2, dv2, 5, c, &fl) with incoming *l = fl (seeded variant)
  pcr_fasttwosum(chp, fl, cAcosCHi[4].f, fl);
  clp := cAcosCLo[4].f + fl;
  chp := pcr_muldd(v2, dv2, chp, clp, clp);
  pcr_fasttwosum(th, tl, cAcosCHi[3].f, chp);
  chp := th;  clp := (cAcosCLo[3].f + clp) + tl;
  chp := pcr_muldd(v2, dv2, chp, clp, clp);
  pcr_fasttwosum(th, tl, cAcosCHi[2].f, chp);
  chp := th;  clp := (cAcosCLo[2].f + clp) + tl;
  chp := pcr_muldd(v2, dv2, chp, clp, clp);
  pcr_fasttwosum(th, tl, cAcosCHi[1].f, chp);
  chp := th;  clp := (cAcosCLo[1].f + clp) + tl;
  chp := pcr_muldd(v2, dv2, chp, clp, clp);
  pcr_fasttwosum(th, tl, cAcosCHi[0].f, chp);
  chp := th;  clp := (cAcosCLo[0].f + clp) + tl;
  fh := chp;  fl := clp;

  fh := pcr_muldd(v, dv, fh, fl, fl);

  ph := jtd * cAcosPi64H.f;
  pl := cAcosPi64M.f * jtd;
  ps := cAcosPi64L.f * jtd;
  // pl = fastsum(fh, fl, pl, ps, &ps)  (fasttwosum form, not full twosum)
  pcr_fasttwosum(pl_new, ps_inner, fh, pl);
  ps := (fl + ps) + ps_inner;
  pl := pl_new;
  pcr_fasttwosum(ph, pl, ph, pl);
  pcr_fasttwosum(pl, ps, pl, ps);
  pcr_fasttwosum(ph, pl, ph, pl);
  pcr_fasttwosum(pl, ps, pl, ps);

  th_u.f := ph;  tl_u.f := pl;
  tn_u.u := (th_u.u and (UInt64($7FF) shl 52)) - (UInt64(53) shl 52);
  tl_u.u := tl_u.u and (UInt64($7FFFFFFFFFFFFFFF));   // |pl|
  dn := Int64(tl_u.u - tn_u.u);
  de := Int64((tn_u.u - tl_u.u) shr 52);
  hard := ((dn >= -2) and (dn <= 0)) or (de > 46);
  res := ph + pl;
  if hard then res := AsinDatabase(x, res);
  Result := res;
end;

function pcr_asin(x: Double): Double;
var
  ix: Tb64u64;
  ax: UInt64;
  k: Int64;
  j: Int64;
  f0h, f0l, t, z, zl, jd: Double;
  t2, d_poly, fh, fl, eps, lb, ub, sum_sh, fastsum_sl: Double;
begin
  ix.f := x;
  ax := ix.u shl 1;

  if ax > UInt64($7FC0000000000000) then
  begin
    // |x| > 0.5 branch
    k := ix.u shr 63;
    if k = 0 then
    begin
      f0h :=  cAcosPiHalfH.f;  f0l :=  cAcosPiHalfL.f;
    end
    else
    begin
      f0h := -cAcosPiHalfH.f;  f0l := -cAcosPiHalfL.f;
    end;
    if ax >= UInt64($7FE0000000000000) then
    begin
      if ax = UInt64($7FE0000000000000) then
      begin Result := f0h + f0l; Exit; end;        // |x| = 1
      if ax > UInt64($FFE0000000000000) then
      begin Result := x + x; Exit; end;            // NaN
      Result := Double(0.0) / Double(0.0); Exit;   // |x|>1: NaN (negative, like C)
    end;
    t  := Double(2.0) - Double(2.0) * Abs(x);
    jd := pcr_roundeven(t * cAcosP5.f);            // t * 32
    if x >= Double(0.0) then z := -Sqrt(t) else z := Sqrt(t);  // copysign(sqrt(t), -x)
    zl := pcr_fma(z, z, -t) * ((-Double(0.5) / t) * z);
    t  := Double(0.25) * t - jd * cAcosN7.f;       // 0.25*t - jd/128
  end
  else
  begin
    // |x| <= 0.5
    if ax < cAsinSmallTh.u then
    begin
      // tiny x: asin(x) ~ x. fma(2^-55, x, x) preserves sign of ±0 with
      // hardware FMA but the emulated pcr_fma_pascal can flip -0 → +0;
      // short-circuit for x=±0 to keep the sign bit intact.
      if ax = 0 then
        Result := x
      else
        Result := pcr_fma(cAsinSmallC.f, x, x);
      Exit;
    end;
    f0h := Double(0.0);
    f0l := Double(0.0);
    t   := x * x;
    jd  := pcr_roundeven(t * cAcosP7.f);           // t * 128
    t   := pcr_fma(x, x, -cAcosN7.f * jd);
    z   := x;
    zl  := Double(0.0);
  end;

  j := Trunc(jd);
  t2 := t * t;
  d_poly := t * ((cAcosCC[j,2].f + t * cAcosCC[j,3].f) +
                 t2 * ((cAcosCC[j,4].f + t * cAcosCC[j,5].f) +
                       t2 * (cAcosCC[j,6].f + t * cAcosCC[j,7].f)));
  fh := cAcosCC[j,0].f;
  fl := cAcosCC[j,1].f + d_poly;
  fh := pcr_muldd(z, zl, fh, fl, fl);
  // fastsum(f0h, f0l, fh, fl, &fl_out): sh = fasttwosum(f0h, fh, &sl); fl_out = (f0l+fl)+sl
  fastsum_sl := Double(0.0);
  pcr_fasttwosum(sum_sh, fastsum_sl, f0h, fh);
  fl := (f0l + fl) + fastsum_sl;
  fh := sum_sh;

  eps := Abs(z * t) * cAsinEps1.f + cAsinEps2.f;
  lb  := fh + (fl - eps);
  ub  := fh + (fl + eps);
  if lb <> ub then
    Result := AsinRefine(x, lb)
  else
    Result := lb;
end;

// ---------------------------------------------------------------------------
// pcr_exp — correctly-rounded binary64 exponential.
// Ported from core-math/src/binary64/exp/exp.c by Alexei Sibidanov.
// Reuses cExpT0 and cExpT1 (also used by tanh and the rest of the exp family).
// ---------------------------------------------------------------------------

const
  cExpS:        Tb64u64 = (u:$40B71547652B82FE); //  2^12/log(2)
  cExpL2H:      Tb64u64 = (u:$3F262E42FF000000); //  log(2) hi (29-bit exact)
  cExpL2L:      Tb64u64 = (u:$3D0718432A1B0E26); //  log(2) lo
  cExpL2LL:     Tb64u64 = (u:$3999FF0342542FC3); //  log(2) low-low (refine)
  cExpFastCh2:  Tb64u64 = (u:$3FC55555557E54FF); //  fast-path c2
  cExpFastCh3:  Tb64u64 = (u:$3FA55555553A12F4); //  fast-path c3
  cExpEpsFast:  Double  = 1.64e-19;              //  ub/lb bracket epsilon

  // Thresholds (compared against ix.u or aix = ix.u & 0x7fff..)
  cExpTinyAix:    UInt64 = $3C90000000000000; // |x| <= 0x1p-54 → 1+x
  cExpHugeAix:    UInt64 = $40862E42FEFA39F0; // |x| >= ln(DBL_MAX rounded)
  cExpUnderAix:   UInt64 = $40874910D52D3052; // x <= -1.74910...p+9 → underflow
  cExpSubnormalU: UInt64 = $C086232BDD7ABCD2; // ix.u > this (signed-bit set) → subnormal branch
  cExpAccTinyExp: UInt64 = $3C9;              // (ix>>52)&0x7ff < this → 1+x in refine
  cExpHugeMul:    Tb64u64 = (u:$7FE0000000000000); // 0x1p+1023
  cExpTinyMul:    Tb64u64 = (u:$0010000000000000); // 0x1p-1022

  // Accurate-path polynomial ch[7][2] (hi, lo)
  cExpAccCh: array[0..6, 0..1] of Tb64u64 = (
    ((u:$3FF0000000000000),(u:$0000000000000000)),
    ((u:$3FE0000000000000),(u:$39C712F72ECEC2CF)),
    ((u:$3FC5555555555555),(u:$3C65555555554D07)),
    ((u:$3FA5555555555555),(u:$3C455194D28275DA)),
    ((u:$3F81111111111111),(u:$3C012FAA0E1C0F7B)),
    ((u:$3F56C16C16DA6973),(u:$BBF4BA45AB25D2A3)),
    ((u:$3F2A01A019EB7F31),(u:$BBC9091D845ECD36)));

  // 51-entry exception database (sorted by ix.u of x).
  cExpDb: array[0..50] of UInt64 = (
    UInt64($3CAFFFFFFFFFFFFF), UInt64($3F1BA07D73250DE7), UInt64($3F76A4D1AF9CC989),
    UInt64($3F95A75293A5DCDA), UInt64($3FA42EA46949B3C7), UInt64($3FA7C8BB0CF5D160),
    UInt64($3FC0948D39A41695), UInt64($3FCA065FEFAE814F), UInt64($3FCF6E4C3CED7C72),
    UInt64($3FD1A0408712E00A), UInt64($3FDBCAB27D05ABDE), UInt64($3FE005AE04256BAB),
    UInt64($401273C188AA7B14), UInt64($40183D4BCDEBB3F4), UInt64($40308F51434652C3),
    UInt64($4031D5C2DAEBE367), UInt64($403C44CE0D716A1A), UInt64($404E07E71BFCF06F),
    UInt64($404F7216C4B435C9), UInt64($40654CD1FEA7663A), UInt64($407D6479EBA7C971),
    UInt64($BF1664716B68A409), UInt64($BF2A2FEFEFD580DF), UInt64($BF3CE3F638D0C742),
    UInt64($BF3CEFF32831E2C2), UInt64($BF433ACCAE78B371), UInt64($BF4D792B60084F92),
    UInt64($BF77FB235D76CCE7), UInt64($BF81FF9B8E8B38BE), UInt64($BF854511E930898C),
    UInt64($BF95C5ED0EC83666), UInt64($BF98C56FF5326197), UInt64($BF9A4187F2CA71F9),
    UInt64($BFBA8F783D749A8F), UInt64($BFBBD44FDAED819F), UInt64($BFBDAF693D64FADA),
    UInt64($BFC290EA09E36479), UInt64($BFC8AEB636F3CE35), UInt64($BFCD3F3799439415),
    UInt64($BFCEA16274B0109B), UInt64($BFE22E24FA3D5CF9), UInt64($BFE85068C07FBBF6),
    UInt64($BFEBDC7955D1482C), UInt64($BFF2A9CAD9998262), UInt64($BFFCC37EF7DE7501),
    UInt64($C0002393D5976769), UInt64($C0065061DAF79A78), UInt64($C02E8BDBFCD9144E),
    UInt64($C038F80E06F3A04C), UInt64($C0559F038076039C), UInt64($C06981587AD4542F));

// As-ldexp by raw exponent bits.
function ExpAsLdexp(x: Double; i: Int64): Double; inline;
var u: Tb64u64;
begin
  u.f := x;
  u.u := u.u + (UInt64(i) shl 52);
  Result := u.f;
end;

// Strip top 12 bits (sign+exponent) — used to build subnormal result.
function ExpAsToDenormal(x: Double): Double; inline;
var u: Tb64u64;
begin
  u.f := x;
  u.u := u.u and (UInt64($000FFFFFFFFFFFFF));
  Result := u.f;
end;

// as_exp_database (sorted binary search; on hit, ulp-corrects f via signs hash).
function ExpDatabase(x, f: Double): Double;
var
  ix: Tb64u64;
  a, b, m: Int32;
  jf, dr, r: Tb64u64;
  s: UInt64;
  s2: array[0..1] of UInt64;
  t: UInt64;
  k: Int64;
begin
  Result := f;
  ix.f := x;
  s2[0] := UInt64($57F5FE2E5BDE4075);
  s2[1] := UInt64($0000003C1F16B8ED);
  s := UInt64(333811522313371);
  a := 0; b := High(cExpDb); m := (a + b) div 2;
  while a <= b do
  begin
    if cExpDb[m] < ix.u then a := m + 1
    else if cExpDb[m] = ix.u then
    begin
      jf.f := f;
      dr.u := ((s shr m) shl 63) or UInt64($3C90000000000000);
      t := (s2[m shr 5] shr ((m shl 1) and 63)) and 3;
      for k := -1 to 1 do
      begin
        r.u := jf.u + UInt64(k);
        if (r.u and 3) = t then
        begin
          Result := r.f + dr.f;
          Exit;
        end;
      end;
      Exit;
    end
    else b := m - 1;
    m := (a + b) div 2;
  end;
end;

// as_exp_accurate: refinement when fast bracket fails.
function ExpRefine(x: Double): Double;
var
  ix, ixs: Tb64u64;
  t, dx, dxl, dxll, dxh, fh, fl, e: Double;
  th, tl, t0h, t0l, t1h, t1l: Double;
  jt, i0, i1, ie: Int64;
  ch_v, cl_v: Double;
  i: Int32;
  v: Tb64u64;
  d_kind: Int64;
begin
  ix.f := x;
  if ((ix.u shr 52) and $7FF) < cExpAccTinyExp then
  begin
    Result := Double(1.0) + x;
    Exit;
  end;
  t  := pcr_roundeven(x * cExpS.f);
  jt := Trunc(t);
  i0 := (jt shr 6) and $3F;
  i1 := jt and $3F;
  ie := SarInt64(jt, 12);
  t0h := cExpT0[i0, 1].f;  t0l := cExpT0[i0, 0].f;
  t1h := cExpT1[i1, 1].f;  t1l := cExpT1[i1, 0].f;
  th := pcr_muldd(t0h, t0l, t1h, t1l, tl);

  dx   := x - cExpL2H.f * t;
  dxl  := cExpL2L.f * t;
  dxll := cExpL2LL.f * t + pcr_fma(cExpL2L.f, t, -dxl);
  dxh  := dx + dxl;
  dxl  := (dx - dxh) + dxl + dxll;

  // opolydd unrolled (n=7, in/out l in cl_v)
  ch_v := cExpAccCh[6, 0].f;  cl_v := cExpAccCh[6, 1].f;
  for i := 5 downto 0 do
  begin
    ch_v := pcr_muldd(dxh, dxl, ch_v, cl_v, cl_v);
    // th_p = ch + c[i][0]; tl_p = (c[i][0] - th_p) + ch
    fh := ch_v + cExpAccCh[i, 0].f;
    fl := (cExpAccCh[i, 0].f - fh) + ch_v;
    ch_v := fh;
    cl_v := cl_v + fl + cExpAccCh[i, 1].f;
  end;
  fh := ch_v;  fl := cl_v;

  fh := pcr_muldd(dxh, dxl, fh, fl, fl);

  if ix.u > cExpSubnormalU then
  begin
    // x < -0x1.6232bdd7abcd2p+9 — subnormal branch
    ixs.u := UInt64(1 - ie) shl 52;
    fh := pcr_muldd(fh, fl, th, tl, fl);
    // fastsum(th, tl, fh, fl, &fl): sh = fasttwosum(th, fh, &sl); fl = (tl+fl)+sl
    pcr_fasttwosum(ch_v, cl_v, th, fh);
    fl := (tl + fl) + cl_v;
    fh := ch_v;
    pcr_fasttwosum(fh, e, ixs.f, fh);
    fl := fl + e;
    Result := ExpAsToDenormal(fh + fl);
  end
  else
  begin
    if th = Double(1.0) then
    begin
      pcr_fasttwosum(fh, e, th, fh);
      pcr_fasttwosum(fl, e, e, fl);
      v.f := fl;
      if (v.u and UInt64($000FFFFFFFFFFFFF)) = 0 then
      begin
        ixs.f := e;
        // C: d = ((sign(ix) XOR sign(v)) ? -1 : 1); ix.u += d (mod 2^64)
        if (v.u shr 63) <> (ixs.u shr 63) then v.u := v.u - 1 else v.u := v.u + 1;
        fl := v.f;
        d_kind := 0; // suppress hint
      end;
    end
    else
    begin
      fh := pcr_muldd(fh, fl, th, tl, fl);
      pcr_fasttwosum(ch_v, cl_v, th, fh);
      fl := (tl + fl) + cl_v;
      fh := ch_v;
    end;
    pcr_fasttwosum(fh, fl, fh, fl);
    v.f := fl;
    // d = (ix.u + 2) & (~0>>12)  (clear top 12 bits)
    if ((v.u + 2) and UInt64($000FFFFFFFFFFFFF)) <= 2 then
      fh := ExpDatabase(x, fh);
    Result := ExpAsLdexp(fh, ie);
  end;
end;

function pcr_exp(x: Double): Double;
var
  ix: Tb64u64;
  aix: UInt64;
  t, dx, dx2, p, fh, fl, tx, ub, lb, eps, e: Double;
  th, tl, t0h, t0l, t1h, t1l: Double;
  jt, i0, i1, ie: Int64;
  ixs: Tb64u64;
begin
  ix.f := x;
  aix := ix.u and UInt64($7FFFFFFFFFFFFFFF);
  if aix <= cExpTinyAix then
  begin
    Result := Double(1.0) + x;
    Exit;
  end;
  if aix >= cExpHugeAix then
  begin
    if aix > UInt64($7FF0000000000000) then begin Result := x + x; Exit; end;
    if aix = UInt64($7FF0000000000000) then
    begin
      if (ix.u shr 63) <> 0 then Result := Double(0.0)
      else Result := x;
      Exit;
    end;
    if (ix.u shr 63) = 0 then
    begin
      Result := cExpHugeMul.f * cExpHugeMul.f;
      Exit;
    end;
    if aix >= cExpUnderAix then
    begin
      Result := cExpTinyMul.f * cExpTinyMul.f;
      Exit;
    end;
  end;
  t  := pcr_roundeven(x * cExpS.f);
  jt := Trunc(t);
  i0 := (jt shr 6) and $3F;
  i1 := jt and $3F;
  ie := SarInt64(jt, 12);
  t0h := cExpT0[i0, 1].f;  t0l := cExpT0[i0, 0].f;
  t1h := cExpT1[i1, 1].f;  t1l := cExpT1[i1, 0].f;
  th := pcr_muldd(t0h, t0l, t1h, t1l, tl);

  dx  := (x - cExpL2H.f * t) + cExpL2L.f * t;
  dx2 := dx * dx;
  p   := (Double(1.0) + dx * Double(0.5)) + dx2 * (cExpFastCh2.f + dx * cExpFastCh3.f);
  fh  := th;
  tx  := th * dx;
  fl  := tl + tx * p;
  eps := cExpEpsFast;

  if ix.u > cExpSubnormalU then
  begin
    ixs.u := UInt64(1 - ie) shl 52;
    pcr_fasttwosum(fh, e, ixs.f, fh);
    fl := fl + e;
    ub := fh + (fl + eps);
    lb := fh + (fl - eps);
    if ub <> lb then begin Result := ExpRefine(x); Exit; end;
    Result := ExpAsToDenormal(lb);
  end
  else
  begin
    ub := fh + (fl + eps);
    lb := fh + (fl - eps);
    if ub <> lb then begin Result := ExpRefine(x); Exit; end;
    Result := ExpAsLdexp(lb, ie);
  end;
end;

// ---------------------------------------------------------------------------
// pcr_exp2 — correctly-rounded binary64 2^x.
// Ported from core-math/src/binary64/exp2/exp2.c. Reuses cExpT0/cExpT1.
// ---------------------------------------------------------------------------

const
  cExp2Cd: array[0..5, 0..1] of Tb64u64 = (
    ((u:$3F262E42FEFA39EF),(u:$3BBABC9E3B39873E)),
    ((u:$3E4EBFBDFF82C58F),(u:$BAE5E43A53E44950)),
    ((u:$3D6C6B08D704A0C0),(u:$BA0D3A15710D3D83)),
    ((u:$3C83B2AB6FBA4E77),(u:$3914DD5D2A5E025A)),
    ((u:$3B95D87FE7A66459),(u:$B83DC47E47BEB9DD)),
    ((u:$3AA430912F9FB79D),(u:$B744FCD51FCB7640)));

  cExp2FastC: array[0..3] of Tb64u64 = (
    (u:$3F262E42FEFA39EF),(u:$3E4EBFBDFF82C58F),
    (u:$3D6C6B08D73B3E01),(u:$3C83B2AB6FDDA001));

  cExp2Db: array[0..42] of UInt64 = (
    UInt64($3F5E4596526BF94D), UInt64($3F5E76049073067F), UInt64($3F6755AA6FA428CD),
    UInt64($3F679015CE2843D7), UInt64($3F7F99AFEFA30D65), UInt64($3F98D040898B73F5),
    UInt64($3FB673A7779D5293), UInt64($3FB8859F5E252908), UInt64($3FBFA18DFAD6E466),
    UInt64($3FC6C4175EA0C6E1), UInt64($3FC926961243BABA), UInt64($3FE3E34FA6AB969E),
    UInt64($3FEB32A6C92D1185), UInt64($3FF9F1A7D355CB4F), UInt64($BF243C1CEA9BD4D9),
    UInt64($BF277970470A37ED), UInt64($BF27D44C7C8229A6), UInt64($BF395A914543EAB7),
    UInt64($BF399BE01D01064A), UInt64($BF413F898B1E4F28), UInt64($BF468E7A49000B1C),
    UInt64($BF486D2A6E5E8368), UInt64($BF5120D3BDFB6ED8), UInt64($BF53EC814D260D02),
    UInt64($BF647B667916C4B2), UInt64($BF6899E0474BA2D5), UInt64($BF6BA84C6EBFB038),
    UInt64($BF7111BC29CCDBB1), UInt64($BF71FB57E1996E26), UInt64($BF772E40977492C3),
    UInt64($BF7EBF8CF367FCB8), UInt64($BF807F812303F10A), UInt64($BF9234ADA2403885),
    UInt64($BF935DD739305031), UInt64($BFA526CE079B05A5), UInt64($BFB3EA95A5C16E4A),
    UInt64($BFC33564DB4BB9EC), UInt64($BFCD4854D9F87FCA), UInt64($BFCFE89353E31CBF),
    UInt64($BFD83960B2A8D2C4), UInt64($BFDE242801B45D0D), UInt64($BFECEF4C143B5ADF),
    UInt64($BFF60E582CAA34B1));

  cExp2TinyAx:   UInt64 = $792E2A8ECA5705FC; // ax <= this → 1 + copysign(2^-54, x)
  cExp2HugeAx:   UInt64 = $8120000000000000; // ax >= this → |x| >= 1024
  cExp2UnderIxU: UInt64 = $C090CC0000000000; // ix.u >= this → underflow
  cExp2SubIxU:   UInt64 = $C08FF00000000000; // ix.u <= this → normal (> subnormal)
  cExp2P54:      Tb64u64 = (u:$3C90000000000000); // 0x1p-54
  cExp2TinyLo:   Tb64u64 = (u:$BC971547652B82FE); // -0x1.71547652b82fep-54
  cExp2TinyHi:   Tb64u64 = (u:$3CA71547652B82FD); //  0x1.71547652b82fdp-53

function Exp2Database(x, f: Double): Double;
var
  ix, jf, dy, y: Tb64u64;
  s2: array[0..1] of UInt64;
  k: UInt64;
  p: UInt64;
  a, b, m: Int32;
  i: Int64;
begin
  Result := f;
  ix.f := x;
  s2[0] := UInt64($3B216FBD5FD7665F);
  s2[1] := UInt64($000000000034C797);
  k := UInt64(8677191773140);
  a := 0; b := High(cExp2Db); m := (a + b) div 2;
  while a <= b do
  begin
    if cExp2Db[m] < ix.u then a := m + 1
    else if cExp2Db[m] = ix.u then
    begin
      p := (s2[m shr 5] shr ((m * 2) and 63)) and 3;
      jf.f := f;
      dy.u := (UInt64($3C90) or ((k shr m) shl 15)) shl 48;
      for i := -1 to 1 do
      begin
        y.u := jf.u + UInt64(i);
        if (y.u and 3) = p then
        begin
          Result := y.f + dy.f;
          Exit;
        end;
      end;
      Exit;
    end
    else b := m - 1;
    m := (a + b) div 2;
  end;
end;

function Exp2Refine(x: Double): Double;
var
  ix, ixs, v: Tb64u64;
  sx, fx, z, t, th, tl, t0h, t0l, t1h, t1l, fh, fl, e, ch_v, cl_v: Double;
  k, i0, i1, ie: Int64;
  i: Int32;
begin
  ix.f := x;
  sx := Double(4096.0) * x;
  fx := pcr_roundeven(sx);
  z  := sx - fx;
  k  := Trunc(fx);
  i0 := (k shr 6) and $3F;
  i1 := k and $3F;
  ie := k shr 12;
  t0h := cExpT0[i0, 1].f;  t0l := cExpT0[i0, 0].f;
  t1h := cExpT1[i1, 1].f;  t1l := cExpT1[i1, 0].f;
  th := pcr_muldd(t0h, t0l, t1h, t1l, tl);

  // polydd(z, 6, cd, &fl): scalar*dd polynomial
  ch_v := cExp2Cd[5, 0].f;  cl_v := cExp2Cd[5, 1].f;
  for i := 4 downto 0 do
  begin
    ch_v := pcr_mulddd_pd(ch_v, cl_v, z, cl_v);
    fh := ch_v + cExp2Cd[i, 0].f;
    fl := (cExp2Cd[i, 0].f - fh) + ch_v;
    ch_v := fh;
    cl_v := cl_v + fl + cExp2Cd[i, 1].f;
  end;
  fh := ch_v; fl := cl_v;
  fh := pcr_mulddd_pd(fh, fl, z, fl);

  if ix.u <= cExp2SubIxU then
  begin
    // tiny-x guard: if -0x1.71547652b82fep-54 <= x <= 0x1.71547652b82fdp-53, exp2(x)=fma(x,0.5,1)
    if (x >= cExp2TinyLo.f) and (x <= cExp2TinyHi.f) then
    begin
      Result := pcr_fma(x, Double(0.5), Double(1.0));
      Exit;
    end;
    if (k and $FFF) = 0 then
    begin
      // 4096*x rounds to 4096*integer → z=0 means 2^x near exact
      pcr_fasttwosum(fh, e, th, fh);
      pcr_fasttwosum(fl, e, e, fl);
      v.f := fl;
      if (v.u and UInt64($000FFFFFFFFFFFFF)) = 0 then
      begin
        if ((v.u shr 52) and $7FF) <> 0 then
        begin
          ixs.f := e;
          if (v.u shr 63) <> (ixs.u shr 63) then v.u := v.u - 1 else v.u := v.u + 1;
          fl := v.f;
        end;
      end;
    end
    else
    begin
      fh := pcr_muldd(fh, fl, th, tl, fl);
      pcr_fasttwosum(ch_v, cl_v, th, fh);
      fl := (tl + fl) + cl_v;
      fh := ch_v;
    end;
    pcr_fasttwosum(fh, fl, fh, fl);
    v.f := fl;
    if ((v.u + 2) and UInt64($000FFFFFFFFFFFFF)) <= 2 then
      fh := Exp2Database(x, fh);
    Result := ExpAsLdexp(fh, ie);
  end
  else
  begin
    ixs.u := UInt64(1 - ie) shl 52;
    fh := pcr_muldd(fh, fl, th, tl, fl);
    pcr_fasttwosum(ch_v, cl_v, th, fh);
    fl := (tl + fl) + cl_v;
    fh := ch_v;
    pcr_fasttwosum(fh, e, ixs.f, fh);
    fl := fl + e;
    Result := ExpAsToDenormal(fh + fl);
  end;
end;

function pcr_exp2(x: Double): Double;
var
  ix, ixs: Tb64u64;
  ax, m_bits, ex, frac: UInt64;
  sx, fx, z, z2, th, tl, t0h, t0l, t1h, t1l, fh, fl, tz, eps, ub, lb, e, signed_p54: Double;
  k, i0, i1, ie: Int64;
begin
  ix.f := x;
  ax := ix.u shl 1;
  if ax = 0 then begin Result := Double(1.0); Exit; end;
  if ax >= cExp2HugeAx then
  begin
    if ax > UInt64($FFE0000000000000) then begin Result := x + x; Exit; end;
    if ax = UInt64($FFE0000000000000) then
    begin
      if (ix.u shr 63) <> 0 then Result := Double(0.0) else Result := x;
      Exit;
    end;
    if (ix.u shr 63) <> 0 then
    begin
      if ix.u >= cExp2UnderIxU then
      begin
        Result := cExpTinyMul.f * cExpTinyMul.f;
        Exit;
      end;
    end
    else
    begin
      Result := cExpHugeMul.f * x;
      Exit;
    end;
  end;
  if ax <= cExp2TinyAx then
  begin
    if (ix.u shr 63) <> 0 then signed_p54 := -cExp2P54.f else signed_p54 := cExp2P54.f;
    Result := Double(1.0) + signed_p54;
    Exit;
  end;

  m_bits := ix.u shl 12;
  ex     := (ax shr 53) - UInt64($3FF);
  frac   := (ex shr 63) or (m_bits shl (ex and 63));

  sx := Double(4096.0) * x;
  fx := pcr_roundeven(sx);
  z  := sx - fx;
  z2 := z * z;
  k  := Trunc(fx);
  i0 := (k shr 6) and $3F;
  i1 := k and $3F;
  ie := k shr 12;
  t0h := cExpT0[i0, 1].f;  t0l := cExpT0[i0, 0].f;
  t1h := cExpT1[i1, 1].f;  t1l := cExpT1[i1, 0].f;
  th := pcr_muldd(t0h, t0l, t1h, t1l, tl);

  tz := th * z;
  fh := th;
  fl := tz * ((cExp2FastC[0].f + z * cExp2FastC[1].f)
              + z2 * (cExp2FastC[2].f + z * cExp2FastC[3].f)) + tl;
  eps := cExpEpsFast;

  if ix.u <= cExp2SubIxU then
  begin
    if frac <> 0 then
    begin
      ub := fh + (fl + eps);
      fh := fh + (fl - eps);
      if ub <> fh then begin Result := Exp2Refine(x); Exit; end;
    end;
    Result := ExpAsLdexp(fh, ie);
  end
  else
  begin
    ixs.u := UInt64(1 - ie) shl 52;
    pcr_fasttwosum(fh, e, ixs.f, fh);
    fl := fl + e;
    if frac <> 0 then
    begin
      ub := fh + (fl + eps);
      fh := fh + (fl - eps);
      if ub <> fh then begin Result := Exp2Refine(x); Exit; end;
    end;
    Result := ExpAsToDenormal(fh);
  end;
  lb := 0; // suppress hint
end;

// ---------------------------------------------------------------------------
// pcr_exp10 — correctly-rounded binary64 10^x.
// Ported from core-math/src/binary64/exp10/exp10.c. Reuses cExpT0/cExpT1.
// ---------------------------------------------------------------------------

const
  cExp10AccC: array[0..5, 0..1] of Tb64u64 = (
    ((u:$40026BB1BBB55516),(u:$BCAF48AD494EA102)),
    ((u:$40053524C73CEA69),(u:$BCAE2BFAB318D399)),
    ((u:$4000470591DE2CA4),(u:$3CA81F50779E162B)),
    ((u:$3FF2BD7609FD98C4),(u:$3C931A5CC5D3D313)),
    ((u:$3FE1429FFD336AA3),(u:$3C8910DE8C68A0C2)),
    ((u:$3FCA7ED7086882B4),(u:$BC605E703D496537)));

  cExp10FastCh: array[0..3] of Tb64u64 = (
    (u:$40026BB1BBB55516),(u:$40053524C73CEA69),
    (u:$4000470591FD74E1),(u:$3FF2BD760A1F32A5));

  cExp10Scale: Tb64u64 = (u:$40CA934F0979A371); // 0x1.a934f0979a371p+13
  cExp10L0:    Tb64u64 = (u:$3F13441350800000); // 0x1.34413508p-14
  cExp10L1:    Tb64u64 = (u:$3D1F79FEF311F12B); // 0x1.f79fef311f12bp-46
  cExp10L2:    Tb64u64 = (u:$39AAC0B7C917826B); // 0x1.ac0b7c917826bp-101
  cExp10EpsFast: Double = 1.63e-19;

  cExp10HugeAix:  UInt64 = $40734413509F79FE; // aix > this → |x| > 0x1.34413509f79fep+8
  cExp10UnderAix: UInt64 = $407439B746E36B52; // aix > this → underflow
  cExp10TinyAix:  UInt64 = $3C7BCB7B1526E50E; // aix <= this → return 1+x
  cExp10SubIxU:   UInt64 = $C0733A7146F72A42; // ix.u < this → normal branch
  cExp10UnderHi:  Tb64u64 = (u:$0018000000000000); // 0x1.8p-1022
  cExp10UnderLo:  Tb64u64 = (u:$3C80000000000000); // 0x1p-55

  cExp10Db: array[0..48] of UInt64 = (
    UInt64($3F4821E0F2AFB970), UInt64($3F57C3DDD23AC8CA), UInt64($3F5A2D7C1699E82D),
    UInt64($3F7EC65645EDC394), UInt64($3F890D7373B3A546), UInt64($3F97E3C84F2CB9B5),
    UInt64($3FA25765968ECD68), UInt64($3FA9AA6FD4D21A47), UInt64($3FAE7B525705EDEF),
    UInt64($3FD12E02AA997AF2), UInt64($3FDC414AA8BD83B1), UInt64($3FDD7D271AB4EEB4),
    UInt64($3FE1FE5F30572361), UInt64($3FE522C9F19CC202), UInt64($3FF1DAF94CF0BD01),
    UInt64($3FF75F49C6AD3BAD), UInt64($3FFA3C782D4F54FC), UInt64($3FFCC30B915EC8C4),
    UInt64($400EE9674267E65F), UInt64($4012D5494EB1DD13), UInt64($40389063309F3004),
    UInt64($4052A59B82B6FC5E), UInt64($406CDE37694F4D10), UInt64($BF045DDB10382E3F),
    UInt64($BF0485426A688467), UInt64($BF06506061AAE6F7), UInt64($BF0898A8C3990624),
    UInt64($BF117362E953393B), UInt64($BF1E40231E216CAD), UInt64($BF27A7F33CC3FD0B),
    UInt64($BF363DF14C04AB23), UInt64($BF3A1B18D3A28957), UInt64($BF3E12494018E44C),
    UInt64($BF44C7A2BE09B10E), UInt64($BF4DE686910F4F52), UInt64($BF5EBB11D32C9493),
    UInt64($BF7F6F96F005FD47), UInt64($BF8B44E17164CE91), UInt64($BF93B95082297EA7),
    UInt64($BF95B25114A07A72), UInt64($BFBA9CF11E5ADBC5), UInt64($BFCC360CDDE773F7),
    UInt64($BFD56FF305822F26), UInt64($BFDC03419F51B93E), UInt64($BFE1416C72A588A6),
    UInt64($BFED18176754AAC7), UInt64($C01AA5575135E2D3), UInt64($C034CD4AF2FCA2B4),
    UInt64($C05DA5B10D8689FD));

function Exp10Database(x, f: Double): Double;
var
  ix, jf, d, probe: Tb64u64;
  s: UInt64;
  s2: array[0..1] of UInt64;
  p: UInt64;
  a, b, m: Int32;
  i: Int32;
begin
  Result := f;
  ix.f := x;
  s2[0] := UInt64($7EB37EF5AC3FE7C6);
  s2[1] := UInt64($00000003781B19E1);
  s := UInt64(371470981966157);
  a := 0; b := High(cExp10Db); m := (a + b) div 2;
  while a <= b do
  begin
    if cExp10Db[m] < ix.u then a := m + 1
    else if cExp10Db[m] = ix.u then
    begin
      d.u := (((s shr m) and 1) shl 63) or UInt64($3C90000000000000);
      p := s2[m shr 5] shr (2 * (m and 31));
      jf.f := f;
      for i := 0 to 2 do
      begin
        case i of
          0: probe.u := jf.u;
          1: probe.u := jf.u - 1;
          2: probe.u := jf.u + 1;
        end;
        if ((probe.u xor p) and 3) = 0 then
        begin
          Result := probe.f + d.f;
          Exit;
        end;
      end;
      Exit;
    end
    else b := m - 1;
    m := (a + b) div 2;
  end;
end;

function Exp10Refine(x: Double): Double;
var
  ix, v, l: Tb64u64;
  t, dx, dxl, dxll, dxh, th, tl, t0h, t0l, t1h, t1l, fh, fl, ch_v, cl_v: Double;
  jt, i0, i1, ie: Int64;
  i: Int32;
  sfh_differ: Boolean;
  delta: UInt64;
begin
  ix.f := x;
  t  := pcr_roundeven(cExp10Scale.f * x);
  jt := Trunc(t);
  i0 := (jt shr 6) and $3F;
  i1 := jt and $3F;
  ie := SarInt64(jt, 12);
  t0h := cExpT0[i0, 1].f;  t0l := cExpT0[i0, 0].f;
  t1h := cExpT1[i1, 1].f;  t1l := cExpT1[i1, 0].f;
  th := pcr_muldd(t0h, t0l, t1h, t1l, tl);

  // Accurate path uses l1 = -0x1.f79fef311f12bp-46, l2 = -0x1.ac0b7c917826bp-101.
  dx   := x - cExp10L0.f * t;
  dxl  := -cExp10L1.f * t;
  dxll := -cExp10L2.f * t + pcr_fma(-cExp10L1.f, t, -dxl);
  dxh  := dx + dxl;
  dxl  := ((dx - dxh) + dxl) + dxll;

  // opolydd(dxh, dxl, 6, c, &fl)
  ch_v := cExp10AccC[5, 0].f;  cl_v := cExp10AccC[5, 1].f;
  for i := 4 downto 0 do
  begin
    ch_v := pcr_muldd(dxh, dxl, ch_v, cl_v, cl_v);
    fh := ch_v + cExp10AccC[i, 0].f;
    fl := (cExp10AccC[i, 0].f - fh) + ch_v;
    ch_v := fh;
    cl_v := cl_v + fl + cExp10AccC[i, 1].f;
  end;
  fh := ch_v; fl := cl_v;
  fh := pcr_muldd(dxh, dxl, fh, fl, fl);

  if ix.u < cExp10SubIxU then
  begin
    if (jt and $FFF) = 0 then
    begin
      pcr_fasttwosum(fh, fl, fh, fl);
      pcr_fasttwosum(th, fh, th, fh);
      pcr_fasttwosum(fh, fl, fh, fl);
      v.f := fh;
      if (v.u shl 12) = 0 then
      begin
        l.f := fl;
        sfh_differ := (v.u shr 63) <> (l.u shr 63);
        if sfh_differ then delta := UInt64($FFF7FFFFFFFFFFFF)
        else delta := UInt64($0008000000000000);
        v.u := v.u + delta;
      end;
      fh := th + v.f;
    end
    else
    begin
      fh := pcr_muldd(fh, fl, th, tl, fl);
      pcr_fasttwosum(ch_v, cl_v, th, fh);
      fl := (tl + fl) + cl_v;
      fh := ch_v;
      pcr_fasttwosum(fh, fl, fh, fl);
      v.f := fl;
      if (((v.u + 4) and UInt64($000FFFFFFFFFFFFF)) <= 4) or (((v.u shr 52) and $7FF) < 918) then
        fh := Exp10Database(x, fh);
    end;
    Result := ExpAsLdexp(fh, ie);
  end
  else
  begin
    v.u := UInt64(1 - ie) shl 52;
    fh := pcr_muldd(fh, fl, th, tl, fl);
    pcr_fasttwosum(ch_v, cl_v, th, fh);
    fl := (tl + fl) + cl_v;
    fh := ch_v;
    pcr_fasttwosum(fh, tl, v.f, fh);
    fl := fl + tl;
    Result := ExpAsToDenormal(fh + fl);
  end;
end;

function pcr_exp10(x: Double): Double;
var
  ix, v: Tb64u64;
  aix: UInt64;
  t, th, tl, t0h, t0l, t1h, t1l, dx, dx2, p, fh, fx_poly, fl, eps, ub, lb: Double;
  jt, i0, i1, ie: Int64;
begin
  ix.f := x;
  aix := ix.u and UInt64($7FFFFFFFFFFFFFFF);
  if aix > cExp10HugeAix then
  begin
    if aix > UInt64($7FF0000000000000) then begin Result := x + x; Exit; end;
    if aix = UInt64($7FF0000000000000) then
    begin
      if (ix.u shr 63) <> 0 then Result := Double(0.0) else Result := x;
      Exit;
    end;
    if (ix.u shr 63) = 0 then
    begin
      Result := cExpHugeMul.f * Double(2.0); Exit;
    end;
    if aix > cExp10UnderAix then
    begin
      Result := cExp10UnderHi.f * cExp10UnderLo.f; Exit;
    end;
  end;
  if aix <= cExp10TinyAix then
  begin
    Result := Double(1.0) + x;
    Exit;
  end;

  t  := pcr_roundeven(cExp10Scale.f * x);
  jt := Trunc(t);
  i0 := (jt shr 6) and $3F;
  i1 := jt and $3F;
  ie := SarInt64(jt, 12);
  t0h := cExpT0[i0, 1].f;  t0l := cExpT0[i0, 0].f;
  t1h := cExpT1[i1, 1].f;  t1l := cExpT1[i1, 0].f;
  th := pcr_muldd(t0h, t0l, t1h, t1l, tl);

  dx  := (x - cExp10L0.f * t) - cExp10L1.f * t;
  dx2 := dx * dx;
  p   := (cExp10FastCh[0].f + dx * cExp10FastCh[1].f)
       + dx2 * (cExp10FastCh[2].f + dx * cExp10FastCh[3].f);
  fh := th;
  fx_poly := th * dx;
  fl := tl + fx_poly * p;
  eps := cExp10EpsFast;

  if ix.u < cExp10SubIxU then
  begin
    ub := fh + (fl + eps);
    lb := fh + (fl - eps);
    if ub <> lb then begin Result := Exp10Refine(x); Exit; end;
    Result := ExpAsLdexp(fh + fl, ie);
  end
  else
  begin
    v.u := UInt64(1 - ie) shl 52;
    pcr_fasttwosum(fh, tl, v.f, fh);
    fl := fl + tl;
    ub := fh + (fl + eps);
    lb := fh + (fl - eps);
    if ub <> lb then begin Result := Exp10Refine(x); Exit; end;
    Result := ExpAsToDenormal(fh + fl);
  end;
end;

// ---------------------------------------------------------------------------
// pcr_expm1 — correctly-rounded binary64 e^x - 1.
// Ported from core-math/src/binary64/expm1/expm1.c. For |x|>=0.25 reuses
// cExpT0/cExpT1 and the cExpAccCh / cExpFastCh* polynomials. The small-x
// branch (|x|<0.25) carries its own tz table and polynomials.
// ---------------------------------------------------------------------------

const
  // tz[i+32] for i in [-32, 32]. Each entry is (lo, hi) for 2^(i/128) - 1
  // approximation used by the small-x fast path.
  cExpm1Tz: array[0..64, 0..1] of Tb64u64 = (
    ((u:$BC6797D4686C5393),(u:$BFCC5041854DF7D4)),
    ((u:$BC8EA1CB9D163339),(u:$BFCB881A23AEBB48)),
    ((u:$3C8F483A3E8CD60F),(u:$BFCABE60E1F21838)),
    ((u:$3C7DFFD920F493DB),(u:$BFC9F3129931FAB0)),
    ((u:$BC851BFDBB129094),(u:$BFC9262C1C3430A0)),
    ((u:$3C8CD3E5225E2206),(u:$BFC857AA375DB4E4)),
    ((u:$3C5E3A6BDAECE8F9),(u:$BFC78789B0A5E0C0)),
    ((u:$BC8DAF2AE0C2D3D4),(u:$BFC6B5C7478983D8)),
    ((u:$BC7FD36226FADD44),(u:$BFC5E25FB4FDE210)),
    ((u:$3C7D887CD0341AB0),(u:$BFC50D4FAB639758)),
    ((u:$BC8676A52A1A618B),(u:$BFC43693D679612C)),
    ((u:$3C79776B420AD283),(u:$BFC35E28DB4ECD9C)),
    ((u:$3C73D5FD7D70A5ED),(u:$BFC2840B5836CF68)),
    ((u:$3C5A94AD2C8FA0BF),(u:$BFC1A837E4BA3760)),
    ((u:$3C26AD4C353465B0),(u:$BFC0CAAB118A1278)),
    ((u:$BC78BBA170E59B65),(u:$BFBFD6C2D0E3D910)),
    ((u:$BC8E1E0A76CB0685),(u:$BFBE14AED893EEF0)),
    ((u:$3C8FE131F55E75F8),(u:$BFBC4F1331D22D40)),
    ((u:$BC8B5BEEE8BCEE31),(u:$BFBA85E8C62D9C10)),
    ((u:$BC77FE9B02C25E9B),(u:$BFB8B92870FA2B58)),
    ((u:$BC832AE7BDAF1116),(u:$BFB6E8CAFF341FE8)),
    ((u:$3C7A6CFE58CBD73B),(u:$BFB514C92F634788)),
    ((u:$3C68798DE3138A56),(u:$BFB33D1BB17DF2E8)),
    ((u:$BC3589321A7EF10B),(u:$BFB161BB26CBB590)),
    ((u:$BC78D0E700FCFB65),(u:$BFAF0540438FD5C0)),
    ((u:$3C8473EF07D5DD3B),(u:$BFAB3F864C080000)),
    ((u:$BC838E62149C16E2),(u:$BFA7723950130400)),
    ((u:$BC508BB6309BD394),(u:$BFA39D4A1A77E050)),
    ((u:$BC8BAD3FD501A227),(u:$BF9F8152AEE94500)),
    ((u:$3C63D27AC39ED253),(u:$BF97B88F290230E0)),
    ((u:$BC8B60BBD08AAC55),(u:$BF8FC055004416C0)),
    ((u:$BC4A00D03B3359DE),(u:$BF7FE0154AAEED80)),
    ((u:$0000000000000000),(u:$0000000000000000)),
    ((u:$3C8861931C15E39B),(u:$3F80100AB00222C0)),
    ((u:$3C77AB864B3E9045),(u:$3F90202AD5778E40)),
    ((u:$3C74E5659D75E95B),(u:$3F984890D9043740)),
    ((u:$3C78E0BD083ABA81),(u:$3FA040AC0224FD90)),
    ((u:$3C345CC1CF959B1B),(u:$3FA465509D383EB0)),
    ((u:$BC8EB6980CE14DA7),(u:$3FA89246D053D180)),
    ((u:$3C77324137D6C342),(u:$3FACC79F4F5613A0)),
    ((u:$BC45272FF30EED1B),(u:$3FB082B577D34ED8)),
    ((u:$BC81280F19DACE1C),(u:$3FB2A5DD543CCC50)),
    ((u:$BC8D550AF31C8EC3),(u:$3FB4CD4FC989CD68)),
    ((u:$3C87923B72AA582D),(u:$3FB6F91575870690)),
    ((u:$BC776C2E732457F1),(u:$3FB92937074E0CD8)),
    ((u:$3C881F5C92A5200F),(u:$3FBB5DBD3F681220)),
    ((u:$3C8E8AC7A4D3206C),(u:$3FBD96B0EFF0E790)),
    ((u:$BC712DB6F4BBE33B),(u:$3FBFD41AFCBA45E8)),
    ((u:$BC58C4A5DF1EC7E5),(u:$3FC10B022DB7AE68)),
    ((u:$BC6BD4B1C37EA8A2),(u:$3FC22E3B09DC54D8)),
    ((u:$3C85AEB9860044D0),(u:$3FC353BC9FB00B20)),
    ((u:$BC64C26602C63FDA),(u:$3FC47B8B853AAFEC)),
    ((u:$BC87F644C1F9D314),(u:$3FC5A5AC59B963CC)),
    ((u:$3C8F5AA8EC61FC2D),(u:$3FC6D223C5B10638)),
    ((u:$3C27AB912C69FFEB),(u:$3FC800F67B00D7B8)),
    ((u:$BC5B3564BC0EC9CD),(u:$3FC9322934F54148)),
    ((u:$3C86A7062465BE33),(u:$3FCA65C0B85AC1A8)),
    ((u:$BC885718D2FF1BF4),(u:$3FCB9BC1D3910094)),
    ((u:$BC8045CB0C685E08),(u:$3FCCD4315E9E0834)),
    ((u:$BC16E7FB859D5055),(u:$3FCE0F143B41A554)),
    ((u:$3C851BBDEE020603),(u:$3FCF4C6F5508EE5C)),
    ((u:$3C6E17611AFC42C5),(u:$3FD04623D0B0F8C8)),
    ((u:$BC71C5B2E8735A43),(u:$3FD0E7510FD7C564)),
    ((u:$BC825FE139C4CFFD),(u:$3FD189C1ECAEB084)),
    ((u:$BC789843C4964554),(u:$3FD22D78F0FA061A)));

  // Fast path small-x c[6]
  cExpm1FastC: array[0..5] of Tb64u64 = (
    (u:$3F80000000000000),(u:$3F00000000000000),
    (u:$3E755555555551AD),(u:$3DE555555555599C),
    (u:$3D511111AD1AD69D),(u:$3CB6C16C168B1FB5));

  // Accurate small-x cl[6]
  cExpm1AccCl: array[0..5] of Tb64u64 = (
    (u:$3DA93974A8CA5354),(u:$3D6AE7F3E71E4908),
    (u:$3D2AE7F357341648),(u:$3CE952C7F96664CB),
    (u:$3CA686F8CE633AAE),(u:$3C62F49B2FBFB5B6));

  // Accurate small-x ch[11][2]
  cExpm1AccCh: array[0..10, 0..1] of Tb64u64 = (
    ((u:$3FC5555555555555),(u:$3C65555555555554)),
    ((u:$3FA5555555555555),(u:$3C45555555555123)),
    ((u:$3F81111111111111),(u:$3C01111111118167)),
    ((u:$3F56C16C16C16C17),(u:$BBEF49F49E220CEA)),
    ((u:$3F2A01A01A01A01A),(u:$3B6A019EFF6F919C)),
    ((u:$3EFA01A01A01A01A),(u:$3B39FCFF48A75B41)),
    ((u:$3EC71DE3A556C734),(u:$BB6C14F73758CD7F)),
    ((u:$3E927E4FB7789F5C),(u:$3B3DFCE97931018F)),
    ((u:$3E5AE64567F544E3),(u:$3AFC513DA9E4C9C5)),
    ((u:$3E21EED8EFF8D831),(u:$3ACCA00AF84F2B60)),
    ((u:$3DE6124613A86E8F),(u:$3A8F27AC6000898F)));

  cExpm1Db: array[0..37] of UInt64 = (
    UInt64($3FBE923C188EA79B), UInt64($3FD1A0408712E00A), UInt64($3FD1C38132777B26),
    UInt64($3FD27F4980D511FF), UInt64($3FD8172A0E02F90E), UInt64($3FD8BBE2FB45C151),
    UInt64($3FDBCAB27D05ABDE), UInt64($3FE005AE04256BAB), UInt64($3FEACCFBE46B4EF0),
    UInt64($3FED086543694C5A), UInt64($401273C188AA7B14), UInt64($40183D4BCDEBB3F4),
    UInt64($40308F51434652C3), UInt64($4031D5C2DAEBE367), UInt64($403C44CE0D716A1A),
    UInt64($4042EE70220FB1C5), UInt64($40489D56A0C38E6F), UInt64($4057A60EE15E3E9D),
    UInt64($4061F0DA93354198), UInt64($40654CD1FEA7663A), UInt64($406556C678D5E976),
    UInt64($4072DA9E5E6AF0B0), UInt64($4079E7B643238A14), UInt64($407D6479EBA7C971),
    UInt64($4080BC04AF1B09F5), UInt64($BFBAB86CB1743B75), UInt64($BFD119AAE6072D39),
    UInt64($BFD175693A03B590), UInt64($BFD474D4DE7C14BB), UInt64($BFD789D025948EFA),
    UInt64($BFD82B5DFAF59B4C), UInt64($BFD9D871E078EBCE), UInt64($BFE1397ADD4538AC),
    UInt64($BFE22E24FA3D5CF9), UInt64($BFEDC2B5DF1F7D3D), UInt64($BFF0A54D87783D6F),
    UInt64($BFF2A9CAD9998262), UInt64($BFFE42A2ABB1BF0F));

  cExpm1Quarter:  UInt64 = $3FD0000000000000; // |x| < this → small-x branch
  cExpm1TinyAix:  UInt64 = $3CA0000000000000; // |x| < 2^-53 → fma path
  cExpm1HugeAix:  UInt64 = $40862E42FEFA39F0; // |x| >= ln(huge) → over/underflow handling
  cExpm1UnderIxA: UInt64 = $C0425E4F7B2737FA; // ix.u >= → soft underflow approx
  cExpm1UnderIxB: UInt64 = $C042B708872320E2; // ix.u >= → return -1 + 2^-55
  cExpm1Eps0:     Tb64u64 = (u:$3BEA000000000000); //  0x1.ap-65
  cExpm1EpsAdd:   Tb64u64 = (u:$3970000000000000); //  0x1p-104

  // Soft-underflow approx constants
  cExpm1UA: Tb64u64 = (u:$40425E4F7B2737FA); // 0x1.25e4f7b2737fap+5
  cExpm1UB: Tb64u64 = (u:$3CC8486612173C69); // 0x1.8486612173c69p-51
  cExpm1UC: Tb64u64 = (u:$3C971547652B82FE); // 0x1.71547652b82fep-54
  cExpm1UD: Tb64u64 = (u:$BFEFFFFFFFFFFFFF); // -0x1.fffffffffffffp-1

function Expm1Database(x, f: Double): Double;
var
  ix, jf, dr, r: Tb64u64;
  s: UInt64;
  s2: array[0..1] of UInt64;
  t: UInt64;
  a, b, m: Int32;
  k: Int64;
begin
  Result := f;
  ix.f := x;
  s2[0] := UInt64($76F58B0D65BD5553);
  s2[1] := UInt64($0000000000000C06);
  s := UInt64($300E81651C);
  a := 0; b := High(cExpm1Db); m := (a + b) div 2;
  while a <= b do
  begin
    if cExpm1Db[m] < ix.u then a := m + 1
    else if cExpm1Db[m] = ix.u then
    begin
      jf.f := f;
      dr.u := ((s shr m) shl 63) or
              ((((jf.u shr 52) and $7FF) - 54) shl 52);
      t := (s2[m shr 5] shr ((m shl 1) and 63)) and 3;
      for k := -1 to 1 do
      begin
        r.u := jf.u + UInt64(k);
        if (r.u and 3) = t then
        begin
          Result := r.f + dr.f;
          Exit;
        end;
      end;
      Exit;
    end
    else b := m - 1;
    m := (a + b) div 2;
  end;
end;

function Expm1RefineSmall(x: Double): Double;
var
  fl, fh, hx, x2h, x2l, v0, v1, v2, e: Double;
  ch_v, cl_v, fh_t, fl_t: Double;
  i: Int32;
  v_u, v2_u: Tb64u64;
  delta: UInt64;
begin
  // fl seed = polynomial in cl[6]
  fl := x*(cExpm1AccCl[0].f + x*(cExpm1AccCl[1].f + x*(cExpm1AccCl[2].f
       + x*(cExpm1AccCl[3].f + x*(cExpm1AccCl[4].f + x*cExpm1AccCl[5].f)))));

  // opolyddd(x, 11, ch, &fl)
  pcr_fasttwosum(ch_v, fl, cExpm1AccCh[10, 0].f, fl);
  cl_v := cExpm1AccCh[10, 1].f + fl;
  for i := 9 downto 0 do
  begin
    ch_v := pcr_mulddd_pd(ch_v, cl_v, x, cl_v);
    // fastsum(c[i][0], c[i][1], ch, cl, &cl)
    pcr_fasttwosum(fh_t, fl_t, cExpm1AccCh[i, 0].f, ch_v);
    cl_v := (cExpm1AccCh[i, 1].f + cl_v) + fl_t;
    ch_v := fh_t;
  end;
  fl := cl_v;
  fh := pcr_mulddd_pd(ch_v, fl, x, fl);
  fh := pcr_mulddd_pd(fh, fl, x, fl);
  fh := pcr_mulddd_pd(fh, fl, x, fl);

  hx := Double(0.5) * x;
  x2h := x * hx;
  x2l := pcr_fma(x, hx, -x2h);

  // fastsum(x2h, x2l, fh, fl, &fl)
  pcr_fasttwosum(fh_t, fl_t, x2h, fh);
  fl := (x2l + fl) + fl_t;
  fh := fh_t;

  pcr_fasttwosum(v0, v2, x, fh);
  pcr_fasttwosum(v1, v2, v2, fl);
  pcr_fasttwosum(v0, v1, v0, v1);
  pcr_fasttwosum(v1, v2, v1, v2);

  v_u.f := v1;
  if (v_u.u and UInt64($000FFFFFFFFFFFFF)) = 0 then
  begin
    if (v_u.u shl 1) = 0 then
    begin
      Result := Expm1Database(x, v0);
      Exit;
    end;
    v2_u.f := v2;
    // d = ((sign(v1) ^ sign(v2)) ? -1 : +1)
    if (v_u.u shr 63) <> (v2_u.u shr 63) then delta := UInt64($FFFFFFFFFFFFFFFF)
    else delta := UInt64(1);
    v_u.u := v_u.u + delta;
    v1 := v_u.f;
  end;
  Result := v0 + v1;
  e := 0; // suppress hint
end;

function Expm1RefineLarge(x: Double): Double;
var
  ix, off, v: Tb64u64;
  t, dx, dxl, dxll, dxh, fh, fl, e, th, tl, t0h, t0l, t1h, t1l: Double;
  ch_v, cl_v, fh_t, fl_t: Double;
  jt, i0, i1, ie: Int64;
  i: Int32;
begin
  ix.f := x;
  t  := pcr_roundeven(x * cExpS.f);
  jt := Trunc(t);
  i0 := (jt shr 6) and $3F;
  i1 := jt and $3F;
  ie := SarInt64(jt, 12);
  t0h := cExpT0[i0, 1].f;  t0l := cExpT0[i0, 0].f;
  t1h := cExpT1[i1, 1].f;  t1l := cExpT1[i1, 0].f;
  th := pcr_muldd(t0h, t0l, t1h, t1l, tl);

  dx   := x - cExpL2H.f * t;
  dxl  := cExpL2L.f * t;
  dxll := cExpL2LL.f * t + pcr_fma(cExpL2L.f, t, -dxl);
  dxh  := dx + dxl;
  dxl  := (dx - dxh) + dxl + dxll;

  ch_v := cExpAccCh[6, 0].f;  cl_v := cExpAccCh[6, 1].f;
  for i := 5 downto 0 do
  begin
    ch_v := pcr_muldd(dxh, dxl, ch_v, cl_v, cl_v);
    fh_t := ch_v + cExpAccCh[i, 0].f;
    fl_t := (cExpAccCh[i, 0].f - fh_t) + ch_v;
    ch_v := fh_t;
    cl_v := cl_v + fl_t + cExpAccCh[i, 1].f;
  end;
  fh := ch_v; fl := cl_v;
  fh := pcr_muldd(dxh, dxl, fh, fl, fl);
  fh := pcr_muldd(fh, fl, th, tl, fl);

  pcr_fasttwosum(fh_t, fl_t, th, fh);
  fl := (tl + fl) + fl_t;
  fh := fh_t;

  off.u := UInt64(2048 + 1023 - ie) shl 52;
  if ie < 53 then
    pcr_fasttwosum(fh, e, off.f, fh)
  else if ie < 104 then
    pcr_fasttwosum(fh, e, fh, off.f)
  else
    e := Double(0.0);
  fl := fl + e;
  pcr_fasttwosum(fh, fl, fh, fl);
  v.f := fl;
  if ((v.u + 8) and UInt64($000FFFFFFFFFFFFF)) <= 8 then
    fh := Expm1Database(x, fh);
  Result := ExpAsLdexp(fh, ie);
end;

function Expm1Refine(x: Double): Double;
begin
  if Abs(x) < Double(0.25) then
    Result := Expm1RefineSmall(x)
  else
    Result := Expm1RefineLarge(x);
end;

function pcr_expm1(x: Double): Double;
var
  ix, off: Tb64u64;
  aix: UInt64;
  sx, fx, z, z2, th, tl, fh, fl, eps, ub, lb, e0, rh, rl: Double;
  t, t0h, t0l, t1h, t1l, dx, dx2, p, tx: Double;
  jt, i0, i1, ie: Int64;
  i_idx: Int64;
  e: Double;
  ch_v: Double;
begin
  ix.f := x;
  aix := ix.u and UInt64($7FFFFFFFFFFFFFFF);
  if aix < cExpm1Quarter then
  begin
    if aix < cExpm1TinyAix then
    begin
      // For |x| < 2^-53, expm1(x) rounds to x: the increment 2^-54*|x|
      // is at most 0.5 ULP(x) (specifically m/4 ULP for normal x with
      // m in [1,2)) so round-to-nearest never changes the bits. The C
      // code calls fma to set the inexact flag; we omit that since the
      // numerical result is identical and the emulated pcr_fma rounds
      // some subnormals incorrectly.
      Result := x;
      Exit;
    end;
    sx := Double(128.0) * x;
    fx := pcr_roundeven(sx);
    z  := sx - fx;
    z2 := z * z;
    i_idx := Trunc(fx);
    th := cExpm1Tz[i_idx + 32, 1].f;
    tl := cExpm1Tz[i_idx + 32, 0].f;
    fh := z * cExpm1FastC[0].f;
    fl := z2 * ((cExpm1FastC[1].f + z * cExpm1FastC[2].f)
              + z2 * (cExpm1FastC[3].f + z * (cExpm1FastC[4].f + z * cExpm1FastC[5].f)));
    e0 := cExpm1Eps0.f;
    eps := z2 * e0 + cExpm1EpsAdd.f;
    pcr_fasttwosum(rh, rl, th, fh);
    rl := rl + tl + fl;
    fh := pcr_muldd(th, tl, fh, fl, fl);
    pcr_fasttwosum(ch_v, ub, rh, fh);
    fl := (rl + fl) + ub;
    fh := ch_v;
    ub := fh + (fl + eps);
    lb := fh + (fl - eps);
    if ub <> lb then begin Result := Expm1Refine(x); Exit; end;
    Result := lb;
    Exit;
  end;

  if aix >= cExpm1HugeAix then
  begin
    if aix > UInt64($7FF0000000000000) then begin Result := x + x; Exit; end;
    if aix = UInt64($7FF0000000000000) then
    begin
      if (ix.u shr 63) <> 0 then Result := -Double(1.0) else Result := x;
      Exit;
    end;
    if (ix.u shr 63) = 0 then
    begin
      Result := cExpHugeMul.f * cExpHugeMul.f;
      Exit;
    end;
  end;
  if ix.u >= cExpm1UnderIxA then
  begin
    if ix.u >= cExpm1UnderIxB then
    begin
      Result := -Double(1.0) + cExp10UnderLo.f;
      Exit;
    end;
    Result := (cExpm1UA.f + x + cExpm1UB.f) * cExpm1UC.f + cExpm1UD.f;
    Exit;
  end;

  t  := pcr_roundeven(x * cExpS.f);
  jt := Trunc(t);
  i0 := (jt shr 6) and $3F;
  i1 := jt and $3F;
  ie := SarInt64(jt, 12);
  t0h := cExpT0[i0, 1].f;  t0l := cExpT0[i0, 0].f;
  t1h := cExpT1[i1, 1].f;  t1l := cExpT1[i1, 0].f;
  th := pcr_muldd(t0h, t0l, t1h, t1l, tl);

  dx  := (x - cExpL2H.f * t) + cExpL2L.f * t;
  dx2 := dx * dx;
  p   := (Double(1.0) + dx * Double(0.5)) + dx2 * (cExpFastCh2.f + dx * cExpFastCh3.f);
  fh := th;
  tx := th * dx;
  fl := tl + tx * p;
  eps := cExpEpsFast * th;

  off.u := UInt64(2048 + 1023 - ie) shl 52;
  if ie < 53 then
    pcr_fasttwosum(fh, e, off.f, fh)
  else if ie < 75 then
    pcr_fasttwosum(fh, e, fh, off.f)
  else
    e := Double(0.0);
  fl := fl + e;
  ub := fh + (fl + eps);
  lb := fh + (fl - eps);
  if ub <> lb then begin Result := Expm1Refine(x); Exit; end;
  Result := ExpAsLdexp(lb, ie);
end;

{$I cos_port.inc}
{$I sin_port.inc}
{$I tan_port.inc}
{$I sincos_port.inc}
{$I sinpi_port.inc}
{$I atanpi_port.inc}
{$I cosh_port.inc}
{$I tanpi_port.inc}
{$I sinh_port.inc}
{$I acosh_port.inc}
{$I atanh_port.inc}
{$I asinh_port.inc}
{$I log1p_port.inc}
{$I log_port.inc}
{$I log10_port.inc}
{$I erf_port.inc}
{$I exp2m1_port.inc}
{$I exp10m1_port.inc}
{$I acospi_port.inc}
{$I tgamma_port.inc}
{$I erfc_port.inc}
{$I hypot_port.inc}
{$I lgamma_port.inc}
{$I log10p1_port.inc}
{$I log2p1_port.inc}
{$I atan2_port.inc}
{$I atan2pi_port.inc}
{$I asinpi_port.inc}
{$I pow_port.inc}

// ---------------------------------------------------------------------------
// Stubs — delegate to C reference until each function is ported.
// Replace each stub body with the real Pascal port as phases 1-5 progress.
// ---------------------------------------------------------------------------
// pcr_acos — ported above
// pcr_acosh — ported above
// pcr_acospi — ported above
// pcr_asin — ported above
// pcr_asinh — ported above
// pcr_asinpi — ported above
// pcr_atan — ported above
// pcr_atanh — ported above
// pcr_atanpi — ported above
// pcr_cbrt — ported above
// pcr_cos — ported above
// pcr_cosh — ported above
// pcr_cospi — ported above
// pcr_erf — ported above
// pcr_erfc — ported above
// pcr_exp — ported above
// pcr_exp10 — ported above
// pcr_exp10m1 — ported above
// pcr_exp2 — ported above
// pcr_exp2m1 — ported above
// pcr_expm1 — ported above
function  pcr_lgamma(x: Double): Double;  begin Result := pcr_lgamma_pas(x); end;
// pcr_log — ported above
// pcr_log10 — ported above
// pcr_log10p1 — ported above
// pcr_log1p — ported above
// pcr_log2 — ported above
// pcr_log2p1 — ported above
// pcr_sin — ported above
// pcr_sinh — ported above
// pcr_sinpi — ported above
// pcr_tan — ported above
// pcr_tanpi — ported above
function  pcr_tgamma(x: Double): Double;  begin Result := pcr_tgamma_pas(x); end;
function  pcr_atan2(y, x: Double): Double;  begin Result := pcr_atan2_pas(y, x); end;
function  pcr_atan2pi(y, x: Double): Double; begin Result := pcr_atan2pi_pas(y, x); end;
function  pcr_hypot(x, y: Double): Double;  begin Result := pcr_hypot_pas(x, y); end;
function  pcr_pow(x, y: Double): Double;    begin Result := PcrPowPas(x, y);  end;
// pcr_sincos — ported above

end.
