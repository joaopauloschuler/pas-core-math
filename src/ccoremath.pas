// pas-core-math - Pascal port of CORE-MATH
// https://github.com/joaopauloschuler/pas-core-math
//                                                                                                                                                                                                      
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
{$I pascoremath.inc}
unit ccoremath;

interface

{$linklib coremath}
{$linklib m}

// Univariate functions
function cr_acosf(x: Single): Single; cdecl; external 'coremath';
function cr_acoshf(x: Single): Single; cdecl; external 'coremath';
function cr_acospif(x: Single): Single; cdecl; external 'coremath';
function cr_asinf(x: Single): Single; cdecl; external 'coremath';
function cr_asinhf(x: Single): Single; cdecl; external 'coremath';
function cr_asinpif(x: Single): Single; cdecl; external 'coremath';
function cr_atanf(x: Single): Single; cdecl; external 'coremath';
function cr_atanhf(x: Single): Single; cdecl; external 'coremath';
function cr_atanpif(x: Single): Single; cdecl; external 'coremath';
function cr_cbrtf(x: Single): Single; cdecl; external 'coremath';
function cr_cosf(x: Single): Single; cdecl; external 'coremath';
function cr_coshf(x: Single): Single; cdecl; external 'coremath';
function cr_cospif(x: Single): Single; cdecl; external 'coremath';
function cr_erff(x: Single): Single; cdecl; external 'coremath';
function cr_erfcf(x: Single): Single; cdecl; external 'coremath';
function cr_expf(x: Single): Single; cdecl; external 'coremath';
function cr_exp10f(x: Single): Single; cdecl; external 'coremath';
function cr_exp10m1f(x: Single): Single; cdecl; external 'coremath';
function cr_exp2f(x: Single): Single; cdecl; external 'coremath';
function cr_exp2m1f(x: Single): Single; cdecl; external 'coremath';
function cr_expm1f(x: Single): Single; cdecl; external 'coremath';
function cr_lgammaf(x: Single): Single; cdecl; external 'coremath';
function cr_logf(x: Single): Single; cdecl; external 'coremath';
function cr_log10f(x: Single): Single; cdecl; external 'coremath';
function cr_log10p1f(x: Single): Single; cdecl; external 'coremath';
function cr_log1pf(x: Single): Single; cdecl; external 'coremath';
function cr_log2f(x: Single): Single; cdecl; external 'coremath';
function cr_log2p1f(x: Single): Single; cdecl; external 'coremath';
function cr_rsqrtf(x: Single): Single; cdecl; external 'coremath';
function cr_sinf(x: Single): Single; cdecl; external 'coremath';
function cr_sinhf(x: Single): Single; cdecl; external 'coremath';
function cr_sinpif(x: Single): Single; cdecl; external 'coremath';
function cr_tanf(x: Single): Single; cdecl; external 'coremath';
function cr_tanhf(x: Single): Single; cdecl; external 'coremath';
function cr_tanpif(x: Single): Single; cdecl; external 'coremath';
function cr_tgammaf(x: Single): Single; cdecl; external 'coremath';

// Bivariate functions
function cr_atan2f(y: Single; x: Single): Single; cdecl; external 'coremath';
function cr_atan2pif(y: Single; x: Single): Single; cdecl; external 'coremath';
function cr_hypotf(x: Single; y: Single): Single; cdecl; external 'coremath';
function cr_powf(x: Single; y: Single): Single; cdecl; external 'coremath';
function cr_compoundf(x: Single; y: Single): Single; cdecl; external 'coremath';

// sincosf: void cr_sincosf(float x, float *sout, float *cout)
procedure cr_sincosf(x: Single; sout: PSingle; cout: PSingle); cdecl; external 'coremath';

implementation

end.
