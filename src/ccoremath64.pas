// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
{$I pascoremath.inc}
unit ccoremath64;

interface

{$linklib coremath64}
{$linklib m}

// Univariate functions
function cr_acos(x: Double): Double; cdecl; external 'coremath';
function cr_acosh(x: Double): Double; cdecl; external 'coremath';
function cr_acospi(x: Double): Double; cdecl; external 'coremath';
function cr_asin(x: Double): Double; cdecl; external 'coremath';
function cr_asinh(x: Double): Double; cdecl; external 'coremath';
function cr_asinpi(x: Double): Double; cdecl; external 'coremath';
function cr_atan(x: Double): Double; cdecl; external 'coremath';
function cr_atanh(x: Double): Double; cdecl; external 'coremath';
function cr_atanpi(x: Double): Double; cdecl; external 'coremath';
function cr_cbrt(x: Double): Double; cdecl; external 'coremath';
function cr_cos(x: Double): Double; cdecl; external 'coremath';
function cr_cosh(x: Double): Double; cdecl; external 'coremath';
function cr_cospi(x: Double): Double; cdecl; external 'coremath';
function cr_erf(x: Double): Double; cdecl; external 'coremath';
function cr_erfc(x: Double): Double; cdecl; external 'coremath';
function cr_exp(x: Double): Double; cdecl; external 'coremath';
function cr_exp10(x: Double): Double; cdecl; external 'coremath';
function cr_exp10m1(x: Double): Double; cdecl; external 'coremath';
function cr_exp2(x: Double): Double; cdecl; external 'coremath';
function cr_exp2m1(x: Double): Double; cdecl; external 'coremath';
function cr_expm1(x: Double): Double; cdecl; external 'coremath';
function cr_lgamma(x: Double): Double; cdecl; external 'coremath';
function cr_log(x: Double): Double; cdecl; external 'coremath';
function cr_log10(x: Double): Double; cdecl; external 'coremath';
function cr_log10p1(x: Double): Double; cdecl; external 'coremath';
function cr_log1p(x: Double): Double; cdecl; external 'coremath';
function cr_log2(x: Double): Double; cdecl; external 'coremath';
function cr_log2p1(x: Double): Double; cdecl; external 'coremath';
function cr_rsqrt(x: Double): Double; cdecl; external 'coremath';
function cr_sin(x: Double): Double; cdecl; external 'coremath';
function cr_sinh(x: Double): Double; cdecl; external 'coremath';
function cr_sinpi(x: Double): Double; cdecl; external 'coremath';
function cr_tan(x: Double): Double; cdecl; external 'coremath';
function cr_tanh(x: Double): Double; cdecl; external 'coremath';
function cr_tanpi(x: Double): Double; cdecl; external 'coremath';
function cr_tgamma(x: Double): Double; cdecl; external 'coremath';

// Bivariate functions
function cr_atan2(y: Double; x: Double): Double; cdecl; external 'coremath';
function cr_atan2pi(y: Double; x: Double): Double; cdecl; external 'coremath';
function cr_hypot(x: Double; y: Double): Double; cdecl; external 'coremath';
function cr_pow(x: Double; y: Double): Double; cdecl; external 'coremath';

// sincos: void cr_sincos(double x, double *sout, double *cout)
procedure cr_sincos(x: Double; sout: PDouble; cout: PDouble); cdecl; external 'coremath';

implementation

end.
