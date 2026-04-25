// ---------------------------------------------------------------------------
// pcr_pow — correctly-rounded binary64 power.
// Ported from core-math/src/binary64/pow/{pow.c, pow.h, dint.h, qint.h} by
// Paul Zimmermann, Tom Hubrecht and Claude Etienne Rohrbacher.
//
// Three-tier Ziv strategy:
//   * fast path: pure double-double, ~2^-67 absolute error
//     (uses pow.h tables: cPowInverse, cPowLogInv*, cPowT1*, cPowT2*, cPowP1, cPowQ1)
//   * slow path: TDInt64, ~2^-122 absolute error
//     (uses dint tables: cPowInverse2*, cPowLogInv2*, cPowT12, cPowT22, cPowP2, cPowQ2)
//   * very-slow path: TQInt64, ~2^-242 absolute error
//     (uses qint tables: cPowInverse3*, cPowLogInv3*, cPowT13, cPowT23, cPowP3, cPowQ3)
//
// Conventions:
//   * TDInt64 .ex shifted by +1 vs pow/dint.h (mantissa in [0.5,1) here vs [1,2) there)
//   * TQInt64 .ex matches pow/qint.h exactly (no offset)
//   * Tables are pre-shifted at generation time; see tools/pow_gen.py.
// ---------------------------------------------------------------------------

{$I pow_const.inc}

// --- 5.02e: q_1 / p_1 — pure-double polynomial helpers --------------------

// Given |z| < 2^-12.905, compute (qh, ql) ≈ exp(z) with relative error < 2^-64.9.
// Mirrors pow.c:120 (q_1).
procedure PowQ1(out qh, ql: Double; z: Double);
var
  q, h0, h1, l1: Double;
begin
  q  := pcr_fma(cPowQ1_4.f, z, cPowQ1_3.f);
  q  := pcr_fma(q,           z, cPowQ1_2.f);
  h0 := pcr_fma(q,           z, cPowQ1_1.f);
  pcr_a_mul(h1, l1, z, h0);                 // exact: h1 + l1 = z * h0
  pcr_fasttwosum(qh, ql, cPowQ1_0.f, h1);   // |Q_1[0]| = 1 dominates
  ql := ql + l1;
end;

// Given |z| < 32*2^-13, compute (ph, pl) ≈ log(1+z) - z + z^2/2 such that
//   |ph + pl - (log(1+z) - z)| < 2^-75.5
// (ph holds -z^2/2's leading part; the polynomial tail is folded into pl).
// Mirrors pow.c:323 (p_1).
procedure PowP1(out ph, pl: Double; z: Double);
var
  wh, wl, t, u, v: Double;
begin
  pcr_a_mul(wh, wl, z, z);                  // exact: wh + wl = z^2
  t := pcr_fma(cPowP1_5.f, z, cPowP1_4.f);
  u := pcr_fma(cPowP1_3.f, z, cPowP1_2.f);
  v := pcr_fma(cPowP1_1.f, z, cPowP1_0.f);
  u := pcr_fma(t, wh, u);
  v := pcr_fma(u, wh, v);
  u := v * wh;
  ph := Double(-0.5) * wh;
  pl := pcr_fma(u, z, Double(-0.5) * wl);
end;
