// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//                                                                                                                                                                                                      
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
{$I pascoremath.inc}
unit hexfloat;

// Utility to convert C99 hex float literals to Pascal Double/Single values.
// Handles strings like '0x1.8p+0', '-0x1.fffffep-2', '0x0p+0', '0x1p127', etc.

interface

function HexToDouble(const s: string): Double;
function HexToSingle(const s: string): Single;

implementation

uses Math, SysUtils;

function HexToDouble(const s: string): Double;
var
  p: Int32;
  neg: Boolean;
  mantissa: Double;
  exp: Int32;
  negexp: Boolean;
  c: Char;
  digit: Int32;
  intpart: Double;
  fracpart: Double;
  fracscale: Double;
begin
  p := 1;
  Result := 0.0;

  // Skip leading whitespace
  while (p <= Length(s)) and (s[p] = ' ') do
    Inc(p);

  // Parse optional sign
  neg := False;
  if (p <= Length(s)) and (s[p] = '-') then
  begin
    neg := True;
    Inc(p);
  end
  else if (p <= Length(s)) and (s[p] = '+') then
    Inc(p);

  // Expect '0x' or '0X'
  if (p + 1 > Length(s)) or (s[p] <> '0') or not (s[p+1] in ['x', 'X']) then
    raise EConvertError.CreateFmt('HexToDouble: expected 0x prefix in "%s"', [s]);
  Inc(p, 2);

  // Parse integer part of hex mantissa
  intpart := 0.0;
  while p <= Length(s) do
  begin
    c := s[p];
    if (c >= '0') and (c <= '9') then digit := Ord(c) - Ord('0')
    else if (c >= 'a') and (c <= 'f') then digit := Ord(c) - Ord('a') + 10
    else if (c >= 'A') and (c <= 'F') then digit := Ord(c) - Ord('A') + 10
    else break;
    intpart := intpart * 16.0 + digit;
    Inc(p);
  end;

  // Parse optional fractional part
  fracpart := 0.0;
  if (p <= Length(s)) and (s[p] = '.') then
  begin
    Inc(p);
    fracscale := 1.0 / 16.0;
    while p <= Length(s) do
    begin
      c := s[p];
      if (c >= '0') and (c <= '9') then digit := Ord(c) - Ord('0')
      else if (c >= 'a') and (c <= 'f') then digit := Ord(c) - Ord('a') + 10
      else if (c >= 'A') and (c <= 'F') then digit := Ord(c) - Ord('A') + 10
      else break;
      fracpart := fracpart + digit * fracscale;
      fracscale := fracscale / 16.0;
      Inc(p);
    end;
  end;

  mantissa := intpart + fracpart;

  // Expect 'p' or 'P' (binary exponent)
  exp := 0;
  if (p <= Length(s)) and (s[p] in ['p', 'P']) then
  begin
    Inc(p);
    negexp := False;
    if (p <= Length(s)) and (s[p] = '-') then
    begin
      negexp := True;
      Inc(p);
    end
    else if (p <= Length(s)) and (s[p] = '+') then
      Inc(p);

    while p <= Length(s) do
    begin
      c := s[p];
      if (c >= '0') and (c <= '9') then
      begin
        exp := exp * 10 + (Ord(c) - Ord('0'));
        Inc(p);
      end
      else
        break;
    end;

    if negexp then exp := -exp;
  end;

  // Result = sign * mantissa * 2^exp
  if mantissa = 0.0 then
    Result := 0.0
  else
    Result := mantissa * Power(2.0, exp);

  if neg then Result := -Result;
end;

function HexToSingle(const s: string): Single;
begin
  Result := Single(HexToDouble(s));
end;

end.
