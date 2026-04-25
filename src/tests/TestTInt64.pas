// TestTInt64.pas — sanity check for the 192-bit TInt64 arithmetic in
// pascoremathtypes.pas (port of core-math/src/binary64/atan2/tint.h).
{$I pascoremath.inc}
program TestTInt64;
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

procedure CheckClose(got, want: Double; const msg: string; relTol: Double = 1e-15);
var diff: Double;
begin
  if want = 0 then diff := Abs(got)
  else diff := Abs((got - want) / want);
  if not (diff <= relTol) then begin
    WriteLn(Format('FAIL: %s  got=%.17e want=%.17e reldiff=%.3e', [msg, got, want, diff]));
    Inc(fail);
  end;
end;

procedure TestFromTo;
var
  a: TInt64;
  vals: array[0..9] of Double = (
    1.0, -1.0, 0.5, 2.0, 3.141592653589793, -1.5e-100, 7.0, 1e300, 1e-300, 12345.6789);
  i: Integer;
  back: Double;
begin
  WriteLn('-- TestFromTo --');
  for i := 0 to High(vals) do begin
    TIntFromD(a, vals[i]);
    back := TIntToD(a, 0, 0.0, 0.0);
    CheckClose(back, vals[i], 'fromd/tod for ' + FloatToStr(vals[i]), 0.0);
  end;
end;

procedure TestAdd;
var a, b, r: TInt64; got: Double;
begin
  WriteLn('-- TestAdd --');
  TIntFromD(a, 1.0);
  TIntFromD(b, 2.0);
  AddTInt(r, a, b);
  got := TIntToD(r, 0, 0.0, 0.0);
  CheckClose(got, 3.0, '1+2', 0.0);

  TIntFromD(a, 1.0);
  TIntFromD(b, -1.0);
  AddTInt(r, a, b);
  got := TIntToD(r, 0, 0.0, 0.0);
  CheckClose(got, 0.0, '1-1', 0.0);

  TIntFromD(a, 1e16);
  TIntFromD(b, 1.0);
  AddTInt(r, a, b);
  got := TIntToD(r, 0, 0.0, 0.0);
  CheckClose(got, 1e16 + 1.0, '1e16+1', 1e-16);

  // Subtraction with cancellation: 1 - (1 - 2^-53) = exactly 2^-53
  TIntFromD(a, 1.0);
  TIntFromD(b, 1.0 - 1.1102230246251565e-16);  // 1 - 2^-53
  b.sgn := 1; // negate
  AddTInt(r, a, b);
  got := TIntToD(r, 0, 0.0, 0.0);
  CheckClose(got, 1.1102230246251565e-16, '1-(1-2^-53)', 0.0);
end;

procedure TestMul;
var a, b, r: TInt64; got: Double;
begin
  WriteLn('-- TestMul --');
  TIntFromD(a, 3.0);
  TIntFromD(b, 7.0);
  MulTInt(r, a, b);
  got := TIntToD(r, 0, 0.0, 0.0);
  CheckClose(got, 21.0, '3*7', 0.0);

  TIntFromD(a, 1.5);
  TIntFromD(b, 2.5);
  MulTInt(r, a, b);
  got := TIntToD(r, 0, 0.0, 0.0);
  CheckClose(got, 3.75, '1.5*2.5', 0.0);

  TIntFromD(a, -2.0);
  TIntFromD(b, 5.0);
  MulTInt(r, a, b);
  got := TIntToD(r, 0, 0.0, 0.0);
  CheckClose(got, -10.0, '-2*5', 0.0);
end;

procedure TestPi;
var pi_d: Double;
begin
  WriteLn('-- TestPi --');
  // TINT_PI's leading 53 bits should round to math.Pi.
  pi_d := TIntToD(TINT_PI, 1, 0.0, 0.0);
  CheckClose(pi_d, Pi, 'TINT_PI -> double', 1e-16);
  pi_d := TIntToD(TINT_PI2, 1, 0.0, 0.0);
  CheckClose(pi_d, Pi/2, 'TINT_PI2 -> double', 1e-16);
end;

procedure TestDiv;
var r: TInt64; got, b, c: Double;
begin
  WriteLn('-- TestDiv --');
  b := 1.0; c := 3.0;
  DivTIntD(r, b, c);
  got := TIntToD(r, 0, 0.0, 0.0);
  CheckClose(got, b/c, '1/3', 1e-16);

  b := 22.0; c := 7.0;
  DivTIntD(r, b, c);
  got := TIntToD(r, 0, 0.0, 0.0);
  CheckClose(got, b/c, '22/7', 1e-16);

  b := 2.0; c := 1.0;
  DivTIntD(r, b, c);
  got := TIntToD(r, 0, 0.0, 0.0);
  CheckClose(got, b/c, '2/1', 1e-16);

  b := -1.0; c := 2.0;
  DivTIntD(r, b, c);
  got := TIntToD(r, 0, 0.0, 0.0);
  CheckClose(got, b/c, '-1/2', 0.0);
end;

procedure TestShifts;
var a, r: TInt64;
begin
  WriteLn('-- TestShifts --');
  TIntFromD(a, 1.0);
  // After fromd: h=$8000...0, m=0, l=0, ex=1
  Check(a.h = UInt64($8000000000000000), 'fromd(1.0): h MSB');
  Check(a.m = 0, 'fromd(1.0): m=0');
  Check(a.l = 0, 'fromd(1.0): l=0');
  Check(a.ex = 1, 'fromd(1.0): ex=1');

  RShiftTInt(r, a, 1);
  Check(r.h = UInt64($4000000000000000), 'rshift 1: h');
  RShiftTInt(r, a, 64);
  Check((r.h = 0) and (r.m = UInt64($8000000000000000)) and (r.l = 0), 'rshift 64');
  RShiftTInt(r, a, 128);
  Check((r.h = 0) and (r.m = 0) and (r.l = UInt64($8000000000000000)), 'rshift 128');

  LShiftTInt(r, a, 1);
  Check((r.h = 0) and (r.m = 0) and (r.l = 0), 'lshift 1 from 1.0: top bit shifted out');

  // shift 1.0 (h=0x80...) right by 65, bit lands in m
  RShiftTInt(r, a, 65);
  Check(r.h = 0, 'rshift 65: h=0');
  Check(r.m = UInt64($4000000000000000), 'rshift 65: m');
end;

begin
  TestShifts;
  TestFromTo;
  TestAdd;
  TestMul;
  TestPi;
  TestDiv;
  WriteLn;
  if fail = 0 then WriteLn('ALL OK')
  else WriteLn('FAILURES: ', fail);
  Halt(fail);
end.
