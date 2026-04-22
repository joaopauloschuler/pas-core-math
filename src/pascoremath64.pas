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

// ── Stub functions (delegate to C reference until ported) ────────────────────
// pcr_acos declared in ported section above
function  pcr_acosh(x: Double): Double; inline;
function  pcr_acospi(x: Double): Double; inline;
function  pcr_asin(x: Double): Double; inline;
function  pcr_asinh(x: Double): Double; inline;
function  pcr_asinpi(x: Double): Double; inline;
// pcr_atan declared in ported section above
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
// pcr_log2 declared in ported section above
function  pcr_log2p1(x: Double): Double; inline;
function  pcr_sin(x: Double): Double; inline;
function  pcr_sinh(x: Double): Double; inline;
function  pcr_sinpi(x: Double): Double; inline;
function  pcr_tan(x: Double): Double; inline;
// pcr_tanh declared in ported section above
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
  cTanhT0: array[0..63, 0..1] of Tb64u64 = (
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
  cTanhT1: array[0..63, 0..1] of Tb64u64 = (
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

  t0h := cTanhT0[i0, 1].f;
  t1h := cTanhT1[i1, 1].f;
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
    t0l := cTanhT0[i0, 0].f;
    t1l := cTanhT1[i1, 0].f;
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

    t0l := cTanhT0[i0, 0].f;
    t1l := cTanhT1[i1, 0].f;
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
// Stubs — delegate to C reference until each function is ported.
// Replace each stub body with the real Pascal port as phases 1-5 progress.
// ---------------------------------------------------------------------------
// pcr_acos — ported above
function  pcr_acosh(x: Double): Double;   begin Result := cr_acosh(x);   end;
function  pcr_acospi(x: Double): Double;  begin Result := cr_acospi(x);  end;
function  pcr_asin(x: Double): Double;    begin Result := cr_asin(x);    end;
function  pcr_asinh(x: Double): Double;   begin Result := cr_asinh(x);   end;
function  pcr_asinpi(x: Double): Double;  begin Result := cr_asinpi(x);  end;
// pcr_atan — ported above
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
// pcr_log2 — ported above
function  pcr_log2p1(x: Double): Double;  begin Result := cr_log2p1(x);  end;
function  pcr_sin(x: Double): Double;     begin Result := cr_sin(x);     end;
function  pcr_sinh(x: Double): Double;    begin Result := cr_sinh(x);    end;
function  pcr_sinpi(x: Double): Double;   begin Result := cr_sinpi(x);   end;
function  pcr_tan(x: Double): Double;     begin Result := cr_tan(x);     end;
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
