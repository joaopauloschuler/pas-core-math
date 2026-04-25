// TestQInt64.pas — sanity check for the 256-bit TQInt64 arithmetic in
// pascoremathtypes.pas (port of core-math/src/binary64/pow/qint.h).
{$I pascoremath.inc}
program TestQInt64;
uses
  SysUtils, Math, pascoremathtypes;

var
  fail: Integer = 0;

procedure Check(cond: Boolean; const msg: string);
begin
  if not cond then begin
    WriteLn('FAIL: ', msg);
    Inc(fail);
  end;
end;

// Convert a small integer to TQInt64 (n in 1..2^53). qint convention:
// mantissa in [1, 2), so r0's MSB = 1 always represents 1.0.
// value = (r0 / 2^63) * 2^ex.
procedure QFromUInt(out a: TQInt64; n: UInt64; sgn: Byte = 0);
var k: Integer;
begin
  if n = 0 then begin a := QINT_ZERO; Exit; end;
  k := BsrQWord(n);            // ex of the value: 2^k <= n < 2^(k+1)
  a.r0 := n shl (63 - k);      // shift n so MSB lands at bit 63
  a.r1 := 0; a.r2 := 0; a.r3 := 0;
  a.ex := k;
  a.sgn := sgn;
end;

// Extract the small-integer value (works when 0 <= ex <= 62).
function QToUInt(const a: TQInt64): Int64;
var v: UInt64;
  shift: Integer;
begin
  if (a.r0 = 0) and (a.r1 = 0) then begin Result := 0; Exit; end;
  shift := 63 - Integer(a.ex);
  if shift < 0 then v := a.r0 shl (-shift)
  else if shift >= 64 then v := 0
  else v := a.r0 shr shift;
  if a.sgn <> 0 then Result := -Int64(v) else Result := Int64(v);
end;

procedure TestCmp;
var a, b: TQInt64;
begin
  WriteLn('-- TestCmp --');
  QFromUInt(a, 7);
  QFromUInt(b, 5);
  Check(CmpQIntAbs(a, b) = 1, 'cmp 7,5');
  Check(CmpQIntAbs(b, a) = -1, 'cmp 5,7');
  Check(CmpQIntAbs(a, a) = 0, 'cmp 7,7');
  Check(CmpQIntAbs22(a, b) = 1, 'cmp22 7,5');
end;

procedure TestAdd;
var a, b, r: TQInt64;
begin
  WriteLn('-- TestAdd --');
  QFromUInt(a, 3); QFromUInt(b, 5);
  AddQInt(r, a, b);
  Check(QToUInt(r) = 8, '3+5=8');

  QFromUInt(a, 100); QFromUInt(b, 100, 1);
  AddQInt(r, a, b);
  Check(QIntZeroP(r), '100 + (-100) = 0');

  QFromUInt(a, 7); QFromUInt(b, 3, 1);
  AddQInt(r, a, b);
  Check(QToUInt(r) = 4, '7 + (-3) = 4');

  QFromUInt(a, 3, 1); QFromUInt(b, 7);
  AddQInt(r, a, b);
  Check(QToUInt(r) = 4, '(-3) + 7 = 4');

  // |a| = |b| with same sign: result = 2*a
  QFromUInt(a, 5); QFromUInt(b, 5);
  AddQInt(r, a, b);
  Check(QToUInt(r) = 10, '5+5 = 10 (Sterbenz exact)');

  // alias: r aliases a
  QFromUInt(a, 11); QFromUInt(b, 4);
  AddQInt(a, a, b);
  Check(QToUInt(a) = 15, '11+4 alias');

  // 22 variant: only upper 128 bits, behaves like Add for small ints
  QFromUInt(a, 9); QFromUInt(b, 16);
  AddQInt22(r, a, b);
  Check(QToUInt(r) = 25, '9+16 (22 variant)');

  // ZERO
  AddQInt(r, QINT_ZERO, QINT_ONE);
  Check((r.r0 = QINT_ONE.r0) and (r.ex = QINT_ONE.ex), '0+1 = 1');
  AddQInt(r, QINT_ONE, QINT_ZERO);
  Check((r.r0 = QINT_ONE.r0) and (r.ex = QINT_ONE.ex), '1+0 = 1');
end;

procedure TestMul;
var a, b, r: TQInt64;
begin
  WriteLn('-- TestMul --');
  QFromUInt(a, 3); QFromUInt(b, 7);
  MulQInt(r, a, b);
  Check(QToUInt(r) = 21, '3*7=21');

  QFromUInt(a, 12345); QFromUInt(b, 67890);
  MulQInt(r, a, b);
  Check(QToUInt(r) = 12345 * 67890, '12345*67890');

  QFromUInt(a, 5, 1); QFromUInt(b, 8);
  MulQInt(r, a, b);
  Check(QToUInt(r) = -40, '(-5)*8 = -40');

  // ONE * x = x
  QFromUInt(a, 12345);
  MulQInt(r, QINT_ONE, a);
  Check(QToUInt(r) = 12345, '1 * 12345 = 12345');

  // M_ONE * x = -x
  MulQInt(r, QINT_M_ONE, a);
  Check(QToUInt(r) = -12345, '(-1) * 12345 = -12345');

  // mul_qint_11 (only top limbs)
  QFromUInt(a, 17); QFromUInt(b, 19);
  MulQInt11(r, a, b);
  Check(QToUInt(r) = 17 * 19, '17*19 (mul11)');

  // mul_qint_21
  MulQInt21(r, a, b);
  Check(QToUInt(r) = 17 * 19, '17*19 (mul21)');

  // mul_qint_22, 31, 33, 41 — for inputs with no low limbs they all equal mul_qint
  MulQInt22(r, a, b);
  Check(QToUInt(r) = 17 * 19, '17*19 (mul22)');
  MulQInt31(r, a, b);
  Check(QToUInt(r) = 17 * 19, '17*19 (mul31)');
  MulQInt33(r, a, b);
  Check(QToUInt(r) = 17 * 19, '17*19 (mul33)');
  MulQInt41(r, a, b);
  Check(QToUInt(r) = 17 * 19, '17*19 (mul41)');

  // mul_qint_2 (integer scaling)
  QFromUInt(a, 1234);
  MulQIntInt(r, 7, a);
  Check(QToUInt(r) = 1234 * 7, 'mul_qint_2: 1234*7');
  MulQIntInt(r, -3, a);
  Check(QToUInt(r) = -1234 * 3, 'mul_qint_2: 1234*(-3)');
  MulQIntInt(r, 0, a);
  Check(QIntZeroP(r), 'mul_qint_2: 0');
  MulQIntInt(r, 1, a);
  Check(QToUInt(r) = 1234, 'mul_qint_2: *1');
end;

// Verify identity: ONE_Q + M_ONE_Q = 0
procedure TestONEPlusMONE;
var r: TQInt64;
begin
  WriteLn('-- TestONE+M_ONE --');
  AddQInt(r, QINT_ONE, QINT_M_ONE);
  Check(QIntZeroP(r), '1 + (-1) = 0');
end;

// LOG2_Q * LOG2_INV_Q ≈ 2^12 ⇒ result.ex = 12 and r0 ≈ 0x8000... within tight ulp
procedure TestLog2Identity;
var r: TQInt64;
begin
  WriteLn('-- TestLog2 * LOG2_INV --');
  MulQInt(r, QINT_LOG2, QINT_LOG2_INV);
  // log(2) * 2^12/log(2) = 2^12 → ex = 13 with r0 = 0x8000... after normalisation
  // (because the qint normalised representation puts MSB at bit 63, so 2^12 = (0.5)*2^13)
  Check(r.ex = 12, Format('LOG2 * LOG2_INV: ex=%d (want 12)', [r.ex]));
  // Top limb should be very close to 0x8000000000000000
  Check((r.r0 >= UInt64($8000000000000000) - 16) and
        (r.r0 <= UInt64($8000000000000000) + 16),
        Format('LOG2*LOG2_INV r0=%x', [r.r0]));
end;

// Spot-check a 256-bit product where the low limbs matter.
// Compute (1 + 2^-64) * (1 + 2^-64) = 1 + 2^-63 + 2^-128
// In normalised qint form: r0 = 0x8000... , r1 = 0x8000... (= 2^-128 in low limb position),
// after the *2 from the not(t6>>127)=0 branch... actually with both factors close to 1
// we have ex=1 and the leading bit fits.
procedure TestLowLimbMul;
var
  a, b, r: TQInt64;
begin
  WriteLn('-- TestLowLimbMul --');
  // a = 1 + 2^-64. qint mantissa in [1,2), so r0 = 0x8000... (the leading 1)
  // and r1 = 0x8000... (the 2^-64 bit at MSB of second limb), ex = 0.
  a.r0 := UInt64($8000000000000000); a.r1 := UInt64($8000000000000000);
  a.r2 := 0; a.r3 := 0; a.ex := 0; a.sgn := 0;
  b := a;
  MulQInt(r, a, b);
  // (1 + 2^-64)^2 = 1 + 2^-63 + 2^-128. ex = 0.
  // r0 = 0x8000000000000001 (1.0 at bit 63, 2^-63 at bit 0)
  // r2 = 0x8000000000000000 (2^-128 at bit 63 of r2)
  Check(r.ex = 0, Format('(1+2^-64)^2 ex=%d (want 0)', [r.ex]));
  Check(r.r0 = UInt64($8000000000000001),
        Format('r0=%x want 0x8000000000000001', [r.r0]));
  Check(r.r1 = 0, Format('r1=%x want 0', [r.r1]));
  Check(r.r2 = UInt64($8000000000000000),
        Format('r2=%x want 0x8000000000000000', [r.r2]));
end;

begin
  TestCmp;
  TestAdd;
  TestMul;
  TestONEPlusMONE;
  TestLog2Identity;
  TestLowLimbMul;
  WriteLn;
  if fail = 0 then WriteLn('ALL OK')
  else WriteLn('FAILURES: ', fail);
  Halt(fail);
end.
