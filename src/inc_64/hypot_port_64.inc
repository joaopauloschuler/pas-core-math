// ---------------------------------------------------------------------------
// pcr_hypot — correctly-rounded binary64 sqrt(x^2 + y^2).
// Ported from core-math/src/binary64/hypot/hypot.c (Alexei Sibidanov, 2022-2025).
// ---------------------------------------------------------------------------
//
// Notes on infrastructure:
//   - SSE MXCSR get/set is provided by FPC's Math unit (GetMXCSR/SetMXCSR).
//   - The single u128 use site (lm*lm in as_hypot_hard) is open-coded as a
//     64x64->128 multiply, since both halves are needed only for one shifted
//     read and a sticky-bit test.  No general TUInt128 type required.

const
  cHypEMsk:    UInt64 = $7FF0000000000000;  // exponent mask (0x7ff << 52)
  cHypEx1023:  UInt64 = $3FF0000000000000;  // 0x3ff << 52
  cHypGap27:   UInt64 = $01B0000000000000;  // 27 << 52
  cHypAidrLo:  UInt64 = $39B0000000000000;  // lower edge for accurate dispatch
  cHypAidrHi:  UInt64 = $3C9FFFFFFFFFFF80;  // upper edge for accurate dispatch
  cHypAidrMid: UInt64 = $3C90000000000000;  // mid base for ((aidr-base+16)>>5)
  cHypInfMant: UInt64 = $7FF0000000000000;  // 0x7ff << 52
  cHypTwo53:   UInt64 = UInt64(1) shl 53;
  cHyp1pM27:   Tb64u64 = (u: $3E40000000000000);  // 0x1p-27
  cHyp1pM54:   Tb64u64 = (u: $3C90000000000000);  // 0x1p-54
  cHypMaxNorm: Tb64u64 = (u: $7FEFFFFFFFFFFFFF);  // 0x1.fffffffffffffp+1023

// ---------------------------------------------------------------------------
// Local helpers
// ---------------------------------------------------------------------------

// 64x64 -> 128-bit unsigned multiply.  Implements one limb of __int128 math.
procedure HypMul64x64(a, b: UInt64; out hi, lo: UInt64); inline;
var
  ah, al, bh, bl, m1, m2, c: UInt64;
begin
  ah := a shr 32;  al := a and $FFFFFFFF;
  bh := b shr 32;  bl := b and $FFFFFFFF;
  lo := al * bl;
  m1 := ah * bl;
  m2 := al * bh;
  hi := ah * bh;
  // Add (m1 + m2) << 32 to lo|hi, propagating carries.
  c  := lo + (m1 shl 32);
  if c < lo then Inc(hi);
  lo := c;
  c  := lo + (m2 shl 32);
  if c < lo then Inc(hi);
  lo := c;
  hi := hi + (m1 shr 32) + (m2 shr 32);
end;

procedure HypFastTwoSum(out s, e: Double; x, y: Double); inline;
var ss: Double;
begin
  ss := x + y;
  e  := y - (ss - x);
  s  := ss;
end;

// ---------------------------------------------------------------------------
// as_hypot_denorm:  both x and y are subnormal, encoded as a (>= b) integer
// representations.  Returns the correctly-rounded Double sqrt(x^2+y^2).
// ---------------------------------------------------------------------------

function HypotDenorm(a, b: UInt64): Double;
var
  af, bf, op, om, dummy: Double;
  rm: UInt64;
  tm, D: Int64;
  rb, rb2, sb: Boolean;
  underflow: Boolean;
  res: Tb64u64;
  trig_uf: Double;
begin
  af := Double(Int64(a));
  bf := Double(Int64(b));
  underflow := False;
  // a, b are x, y multiplied by 2^1074
  a := a shl 1;
  b := b shl 1;
  rm := UInt64(Trunc(Sqrt(af*af + bf*bf)));
  tm := Int64(rm shl 1);
  // D = a*a + b*b - tm*tm using modular u64 arithmetic; signed difference fits i64
  D := Int64(a*a + b*b - UInt64(tm)*UInt64(tm));
  while D > 2 * tm do begin
    D := D - (2 * tm + 1);
    Inc(tm);
  end;
  while D < 0 do begin
    D := D + (2 * tm - 1);
    Dec(tm);
  end;
  rb  := (tm and 1) <> 0;
  rb2 := D >= tm;
  sb  := D <> 0;
  rm  := UInt64(tm) shr 1;
  underflow := rm < $0010000000000000;  // < 2^52
  if rb or sb then begin
    op := 1.0 + cHyp1pM54.f;
    om := 1.0 - cHyp1pM54.f;
    if op = om then begin // round-to-nearest
      if sb then begin
        if rb then Inc(rm);
        if ((rm shr 52) <> 0) and rb2 then underflow := False;
      end else begin
        // sticky=0, round=1: round-half-to-even
        if (rm and 1) <> 0 then Inc(rm);
      end;
    end else if op > 1.0 then begin // round upwards
      Inc(rm);
      if ((rm shr 52) <> 0) and ((tm and 1) <> 0) then underflow := False;
    end;
    if underflow then begin
      // trigger the underflow exception explicitly.
      trig_uf := Double(Int64($0010000000000000));  // 2^-1022
      dummy := trig_uf * trig_uf;
      if dummy = 0 then ;  // suppress "unused" warning
    end;
  end;
  res.u := rm;
  Result := res.f;
end;

// ---------------------------------------------------------------------------
// as_hypot_hard:  Newton-refined sqrt with 128-bit comparison.
// flag is the entry MXCSR; restored only when D=0 (exact result).
// ---------------------------------------------------------------------------

function HypotHard(x, y: Double; flag: LongWord): Double;
const
  bs = 2;
var
  op, om: Double;
  xi_, yi_, ri, xi_out: Tb64u64;
  bm, lm, m2, lm_lo, lm_hi, rm, tm, rm2, mask: UInt64;
  be, le, re, de, ls, k, i: Int32;
  D, e: Int64;
begin
  op := 1.0 + cHyp1pM54.f;
  om := 1.0 - cHyp1pM54.f;
  xi_.f := x;  yi_.f := y;
  bm := (xi_.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) or (UInt64(1) shl 52);
  lm := (yi_.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) or (UInt64(1) shl 52);
  be := Int32(xi_.u shr 52);
  le := Int32(yi_.u shr 52);
  ri.f := Sqrt(x*x + y*y);
  rm := ri.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12);
  re := Int32(ri.u shr 52) - $3FF;
  rm := rm or (UInt64(1) shl 52);
  for i := 0 to 2 do begin
    if rm = (UInt64(1) shl 52) then begin
      rm := UInt64($FFFFFFFFFFFFFFFF) shr 11;
      Dec(re);
    end else
      Dec(rm);
  end;
  bm := bm shl bs;
  m2 := bm * bm;
  de := be - le;
  ls := bs - de;
  if ls >= 0 then begin
    lm := lm shl ls;
    m2 := m2 + lm * lm;
  end else begin
    HypMul64x64(lm, lm, lm_hi, lm_lo);
    ls := ls * 2;            // ls < 0; shift by -ls
    // m2 += lm2 >> -ls   (low 64 bits)
    if -ls < 64 then
      m2 := m2 + ((lm_hi shl (64 + ls)) or (lm_lo shr (-ls)))
    else if -ls = 64 then
      m2 := m2 + lm_hi
    else
      m2 := m2 + (lm_hi shr (-ls - 64));
    // sticky bit: m2 |= !!(lm2 << (128+ls))  ==  any nonzero in low (-ls) bits of lm2
    if -ls <= 0 then begin
      // unreachable in this branch
    end else if -ls < 64 then begin
      mask := (UInt64(1) shl (-ls)) - 1;
      if (lm_lo and mask) <> 0 then m2 := m2 or 1;
    end else begin
      // -ls in [64..127] won't occur (de bounded), but be safe
      if lm_lo <> 0 then m2 := m2 or 1
      else if -ls > 64 then begin
        mask := (UInt64(1) shl (-ls - 64)) - 1;
        if (lm_hi and mask) <> 0 then m2 := m2 or 1;
      end;
    end;
  end;
  k := bs + re;
  // Find the smallest rm with rm^2 << (2k) >= m2
  repeat
    if rm >= cHypTwo53 then Inc(rm, 2) else Inc(rm);
    tm := rm shl k;
    rm2 := tm * tm;
    D := Int64(m2 - rm2);
  until D <= 0;
  if D = 0 then begin
    SetMXCSR(flag);
  end else begin
    if op = om then begin
      // tm = (rm << k) - (1 << (k - (rm <= 2^53 ? 1 : 0)))
      if rm <= cHypTwo53 then
        tm := (rm shl k) - (UInt64(1) shl (k - 1))
      else
        tm := (rm shl k) - (UInt64(1) shl k);
      D := Int64(m2 - tm*tm);
      if D <> 0 then begin
        if D < 0 then Inc(rm);  // D>>63 == 1 means negative, so rm += 1; matches C semantics for D<0
      end else begin
        rm := rm and (not UInt64(1));  // rm -= rm & 1
      end;
    end else begin
      if op = 1.0 then begin
        if rm > cHypTwo53 then Dec(rm, 2) else Dec(rm);
      end;
    end;
  end;
  if rm >= (UInt64(1) shl 53) then begin
    rm := rm shr 1;
    Inc(re);
  end;
  e := Int64(be) - 1 + Int64(re);
  xi_out.u := UInt64(e shl 52) + rm;
  Result := xi_out.f;
end;

// ---------------------------------------------------------------------------
// as_hypot_overflow:  produce Inf with correct flags for hypot >= 2^1024.
// ---------------------------------------------------------------------------

function HypotOverflow: Double;
var z, f: Double;
begin
  z := cHypMaxNorm.f;
  f := z + z;
  Result := f;
end;

// ---------------------------------------------------------------------------
// pcr_hypot — main entry.
// ---------------------------------------------------------------------------

function pcr_hypot(x, y: Double): Double;
var
  flag: LongWord;
  xi, yi, xd, yd, thd: Tb64u64;
  emsk, ex, ey, wx, wy, wm, de, off_u: UInt64;
  off: Int64;
  ninf, nqnn: Boolean;
  u, v, x2, dx2, y2, dy2, r2, ir2, dr2, th, rsqrt, dz, tl: Double;
  xa, ya, r: Double;
  nz: Int32;
  aidr, mid: UInt64;
begin
  flag := GetMXCSR;
  xi.f := x;  yi.f := y;
  emsk := cHypEMsk;
  ex := xi.u and emsk;
  ey := yi.u and emsk;
  xa := Abs(x);
  ya := Abs(y);
  if (ex = emsk) or (ey = emsk) then begin
    wx := xi.u shl 1;
    wy := yi.u shl 1;
    wm := emsk shl 1;
    ninf := (wx = wm) xor (wy = wm);
    nqnn := ((wx shr 52) = $FFF) xor ((wy shr 52) = $FFF);
    if ninf and nqnn then begin
      if wx = wm then Result := xa * xa
      else            Result := ya * ya;
      Exit;
    end;
    Result := xa + ya;
    Exit;
  end;
  if xa >= ya then begin u := xa; v := ya; end
  else             begin u := ya; v := xa; end;
  xd.f := u;  yd.f := v;
  ey := yd.u;
  if (ey shr 52) = 0 then begin
    if yd.u = 0 then begin Result := xd.f; Exit; end;
    ex := xd.u;
    if (ex shr 52) = 0 then begin
      if ex = 0 then begin Result := 0.0; Exit; end;
      Result := HypotDenorm(ex, ey);
      Exit;
    end;
    nz := pcr_clzll(ey);
    ey := ey shl (nz - 11);
    ey := ey and (UInt64($FFFFFFFFFFFFFFFF) shr 12);
    ey := ey - (UInt64(nz - 12) shl 52);
    yd.u := ey;
  end;
  de := xd.u - yd.u;
  if de > cHypGap27 then begin
    r := pcr_fma(cHyp1pM27.f, v, u);
    Result := r;
    Exit;
  end;
  off := Int64(cHypEx1023) - Int64(xd.u and emsk);
  off_u := UInt64(off);
  xd.u := xd.u + off_u;
  yd.u := yd.u + off_u;
  x := xd.f;
  y := yd.f;
  x2  := x*x;   dx2 := pcr_fma(x, x, -x2);
  y2  := y*y;   dy2 := pcr_fma(y, y, -y2);
  r2  := x2 + y2;
  ir2 := 0.5 / r2;
  dr2 := ((x2 - r2) + y2) + (dx2 + dy2);
  th  := Sqrt(r2);
  rsqrt := th * ir2;
  dz  := dr2 - pcr_fma(th, th, -r2);
  tl  := rsqrt * dz;
  HypFastTwoSum(th, tl, th, tl);
  thd.f := th;
  ey    := UInt64(0);
  if tl >= 0 then begin xi.f := tl; ey := xi.u; end
  else            begin xi.f := -tl; ey := xi.u; end;
  ex := thd.u and emsk;
  // aidr = ey + (0x3fe<<52) - ex
  aidr := ey + UInt64($3FE0000000000000) - ex;
  mid := (aidr - cHypAidrMid + 16) shr 5;
  if (mid = 0) or (aidr < cHypAidrLo) or (aidr > cHypAidrHi) then
    thd.f := HypotHard(x, y, flag);
  thd.u := thd.u - off_u;
  if thd.u >= cHypInfMant then begin
    Result := HypotOverflow;
    Exit;
  end;
  Result := thd.f;
end;
