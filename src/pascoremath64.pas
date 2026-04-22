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

// ── Stub functions (delegate to C reference until ported) ────────────────────
function  pcr_acos(x: Double): Double; inline;
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
// Stubs — delegate to C reference until each function is ported.
// Replace each stub body with the real Pascal port as phases 1-5 progress.
// ---------------------------------------------------------------------------
function  pcr_acos(x: Double): Double;    begin Result := cr_acos(x);    end;
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
