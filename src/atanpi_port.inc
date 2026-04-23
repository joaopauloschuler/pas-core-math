// ---------------------------------------------------------------------------
// pcr_atanpi — correctly-rounded atan(x)/pi (binary64).
// Ported from CORE-MATH cr_atanpi (Sibidanov & Zimmermann, 2023-2025).
// Reuses cAtanAHi/Lo, cAtanC0/1/2, cAtanIdHi/Lo/Lo2/Lo3, cAtanPhiScale,
// cAtanRefCh*/Cl*, cAtanCh*/Ch2*, cAtanPiHalfH/L from the atan port above.
// ---------------------------------------------------------------------------

const
  cAtanpiOneOverPiH: Tb64u64 = (u:$3FD45F306DC9C883); //  0x1.45f306dc9c883p-2
  cAtanpiOneOverPiL: Tb64u64 = (u:$BC76B01EC5417056); // -0x1.6b01ec5417056p-56
  cAtanpiOneOver3Pi: Tb64u64 = (u:$3FBB2995E7B7B604); //  0x1.b2995e7b7b604p-4
  cAtanpiEFast:      Tb64u64 = (u:$3CB4100000000000); //  0x1.41p-52
  cAtanpiETiny:      Tb64u64 = (u:$3BB6F00000000000); //  0x1.6fp-68
  cAtanpiThrBig:     Tb64u64 = (u:$413BE00000000000); //  0x1.bep20
  cAtanpiThrTiny:    Tb64u64 = (u:$3E4C700000000000); //  0x1.c7p-27
  cAtanpiP600:       Tb64u64 = (u:$1A70000000000000); //  0x1p-600
  cAtanpiPp106:      Tb64u64 = (u:$4690000000000000); //  0x1p+106
  cAtanpiPm106:      Tb64u64 = (u:$3950000000000000); //  0x1p-106
  cAtanpiSmallIn:    Tb64u64 = (u:$0015CBA89AF1F855); //  0x1.5cba89af1f855p-1022
  cAtanpiSmallOut:   Tb64u64 = (u:$0006F00F7CD3A40B); //  0x1.bc03df34e902cp-1024
  cAtanpiAs1:        Tb64u64 = (u:$3C80000000000001); //  0x1.0000000000001p-55
  cAtanpiAs2:        Tb64u64 = (u:$3C7FFFFFFFFFFFFF); //  0x1.fffffffffffffp-56
  cAtanpiHalfP55:    Tb64u64 = (u:$3C80000000000000); //  0x1p-55
  cAtanpiPm56:       Tb64u64 = (u:$3C70000000000000); //  0x1p-56

  cAtanpiTinyExc: array[0..19, 0..2] of Tb64u64 = (
    ((u:$3DA9B46CFB0F10CA),(u:$3D905D3AA746E163),(u:$B6A94A704A0E7980)),
    ((u:$3E0EAC33F52580C1),(u:$3DF386DE1F01225A),(u:$B742BE95E5723A66)),
    ((u:$3CA56A4AA740A5A7),(u:$3C8B44453E2404E7),(u:$B5BC5376C33C7A50)),
    ((u:$3D99D47B2E2A4DEE),(u:$3D8071A2E1E2E7A0),(u:$B6DE948955ED137A)),
    ((u:$3D3932DDAEECDABE),(u:$3D200ABFB23406C5),(u:$B9CFFFFFFFFFFFFF)),
    ((u:$3D1E0B9044E5A1CF),(u:$3D032099F7C94E82),(u:$B6722752CF9B0BC5)),
    ((u:$3CD13B412EA2182D),(u:$3CB5F095755008F0),(u:$B5FA8AD4502B5F10)),
    ((u:$3D4FF878736A1BCD),(u:$3D345A65501D94D3),(u:$B677CB6951D7806D)),
    ((u:$3D94F03BD4CEDD49),(u:$3D7AA8DC8DA24739),(u:$B6DFE11160EC4EE9)),
    ((u:$3DB0F2BF03E38A9E),(u:$3D9594435F8C0326),(u:$B6E8EEF034A3526E)),
    ((u:$3CB508ECB38F52F9),(u:$3C9AC84C88F979A2),(u:$B5B63773CBD64D5C)),
    ((u:$3D58DDCD3207FE60),(u:$3D3FA930B4EB436C),(u:$39E0000000000000)),
    ((u:$3CF1331B44547444),(u:$3CD5E63596A62718),(u:$B634A9F744D7CEA7)),
    ((u:$3E0F7681FA730E73),(u:$3DF407A8AC1C768A),(u:$B753C6A80C812C54)),
    ((u:$3DB09DF162E4049C),(u:$3D952849D0B1E523),(u:$BA3FFFFFFFFFFFFF)),
    ((u:$3C959AF9A1194EFE),(u:$3C7B824198B94A89),(u:$391FFFFFFFFFFFFF)),
    ((u:$3D165F203D6BB18E),(u:$3CFC7C00D31634A6),(u:$B62D4AF6699E375F)),
    ((u:$3DD6F89F6DF8258D),(u:$3DBD3F710390DB87),(u:$B6E936061594FAD7)),
    ((u:$3CC3E4D2D87B5AEF),(u:$3CA954626979D7D3),(u:$B5E6F8DBD4F38CDB)),
    ((u:$3DBF46E9E25A51D1),(u:$3DA3E95C097E08E0),(u:$B6F8E5F60920B571)));

  cAtanpiRefExc: array[0..55, 0..2] of Tb64u64 = (
    ((u:$40116265F83E8121),(u:$3FDB648E5B0804DA),(u:$3929DB921C686893)),
    ((u:$4058A5C063607FBA),(u:$3FDFCB1A9F45A781),(u:$BC7FFFFFFFFFFFFF)),
    ((u:$3F91A5E744FF1806),(u:$3F7677CDC7DEA093),(u:$38C417D05C2A4506)),
    ((u:$3F82EEE26511C903),(u:$3F681B11AC013890),(u:$38B508C83C3DB5E3)),
    ((u:$3F7BB0BA10084DE6),(u:$3F61A0C1411A3DC1),(u:$38CA854B560FE4DB)),
    ((u:$3FD55166A4D7D836),(u:$3FBA33741ABEA942),(u:$B9154A4C868A708E)),
    ((u:$3FC5FA4B5E7CE940),(u:$3FABB6715A0B1F81),(u:$39000037A0BFDFAE)),
    ((u:$4032A04FA64D3600),(u:$3FDEE846D3680301),(u:$3903015B69C79EF4)),
    ((u:$4030FB7E1298320C),(u:$3FDECD42A2EF8667),(u:$B9421D6585C13D10)),
    ((u:$402B6DAA66352E75),(u:$3FDE84651D4C9934),(u:$39334E0B5A33A1DB)),
    ((u:$4005717DF9ABCCD3),(u:$3FD8B9BE647E02A9),(u:$B926C5194938E976)),
    ((u:$3FBE0AC5527E7F0D),(u:$3FA309CE21202CCA),(u:$B90141BE39E1919D)),
    ((u:$402932F168987539),(u:$3FDE62F20ABA165C),(u:$B92839C5638718F2)),
    ((u:$4022E342CA2440EE),(u:$3FDDD9D1ECFFC9A0),(u:$B9253DEDDDE1DDE1)),
    ((u:$3FAB37EF5831A537),(u:$3F914FC0E0C2DACE),(u:$B8EF6B67EF93B398)),
    ((u:$3F7D8EA551D1D795),(u:$3F62D0FE21DF78D0),(u:$B8D0F4D438D8DD24)),
    ((u:$3FCBADB08CD8B1C5),(u:$3FB15A7DB07DAA97),(u:$BC5FFFFFFFFFFFFF)),
    ((u:$40618071548E23A6),(u:$3FDFDAC0C22A86C5),(u:$3C7FFFFFFFFFFFFF)),
    ((u:$405E8B9C1FE81CA3),(u:$3FDFD55124D78891),(u:$BC7FFFFFFFFFFFFE)),
    ((u:$3F9E2FAAB41052E7),(u:$3F833629443D32F5),(u:$3C2FFFFFFFFFFFFF)),
    ((u:$3F803F709107440D),(u:$3F64AFDA1E054259),(u:$3C0FFFFFFFFFFFFE)),
    ((u:$3F0928A1CA94D1DC),(u:$3EF0043BD1FD5BAE),(u:$B8512E31518BAD19)),
    ((u:$40ADEDA0DDB0795C),(u:$3FDFFEA37CAE8F37),(u:$BC7FFFFFFFFFFFFE)),
    ((u:$3F472F4CB85B8404),(u:$3F2D850E9850AB33),(u:$B88014404207EB42)),
    ((u:$3EA4B3C41E10A1D4),(u:$3E8A5BDF2B61615B),(u:$B7C180C42F2C36BB)),
    ((u:$3F64AD2F60D5240B),(u:$3F4A537A6DE092A1),(u:$B8AF2E112F36ECA4)),
    ((u:$4131151841058E68),(u:$3FDFFFFECEB437E6),(u:$B932DD2F50460971)),
    ((u:$3EB8B012D530C5C8),(u:$3E9F6EF7B3EE090C),(u:$B7F9D79D9FB153C4)),
    ((u:$4080FA189169E75E),(u:$3FDFF666764F1751),(u:$B93DA2823E32E295)),
    ((u:$40646C64E8717CCF),(u:$3FDFE014CAC09E73),(u:$3C7FFFFFFFFFFFFF)),
    ((u:$40B5A9EE4B581471),(u:$3FDFFF0F445BEE08),(u:$B934256EAE14A448)),
    ((u:$40B6AB7E3B9E35CD),(u:$3FDFFF19F36D9901),(u:$B94246ECD835D1C3)),
    ((u:$3F2842E540B07AC1),(u:$3F0EE3F52B4AF0DC),(u:$B868C285BFA8C0EC)),
    ((u:$4115D2C01B65F659),(u:$3FDFFFFC441A900D),(u:$B948137CC73DD777)),
    ((u:$4109D18F49824D6A),(u:$3FDFFFF9B00D35D1),(u:$391FE1E36287467B)),
    ((u:$413400B2B748DA28),(u:$3FDFFFFEFB46ADE3),(u:$B9010CA50F3BFEF4)),
    ((u:$41327B786B2F1BAE),(u:$3FDFFFFEE5D3F843),(u:$B938AFEE3BDFB07B)),
    ((u:$3EBBBC4A4251B4B9),(u:$3EA1A82F5F7F86BD),(u:$3B4FFFFFFFFFFFFF)),
    ((u:$3F472A01AEF844E5),(u:$3F2D7E514F3FFE99),(u:$3887C2D5BE2B88D8)),
    ((u:$40C95E2B447C5ED5),(u:$3FDFFF99356BF911),(u:$3C7FFFFFFFFFFFFF)),
    ((u:$41288D8109598821),(u:$3FDFFFFE572F576A),(u:$3935430C0B185E51)),
    ((u:$3E83EF59F63222CC),(u:$3E6961C9F0BC5519),(u:$3B0FFFFFFFFFFFFF)),
    ((u:$41172F2198034E17),(u:$3FDFFFFC7C36896E),(u:$B92ADFC30E218AAF)),
    ((u:$3E4DF1C8E36BBCC0),(u:$3E331030ABF0398F),(u:$3ADFFFFFFFFFFFFE)),
    ((u:$408AD5781AB29BBA),(u:$3FDFF9ED31438683),(u:$3C7FFFFFFFFFFFFF)),
    ((u:$40AEDE23D6624A3A),(u:$3FDFFEAE182F312E),(u:$B9426EEE3B0A2CCB)),
    ((u:$4124A7E234CDCE9B),(u:$3FDFFFFE0709E4E1),(u:$BC7FFFFFFFFFFFFD)),
    ((u:$3E81C404899CA420),(u:$3E669EB74CE9E7ED),(u:$37CD126A89BF70E8)),
    ((u:$4077422A83423830),(u:$3FDFF1FC5DEB1C5C),(u:$B930CF8B48B700A5)),
    ((u:$3F0F37BE8C1F5DF0),(u:$3EF3DFB3D39744CF),(u:$B84CAE8F0C659368)),
    ((u:$4138C9F8999242B7),(u:$3FDFFFFF2D9DCA7B),(u:$3C7FFFFFFFFFFFFF)),
    ((u:$412427F51DCF1D68),(u:$3FDFFFFDFA850495),(u:$3C7FFFFFFFFFFFFD)),
    ((u:$413ABB7B457CDF84),(u:$3FDFFFFF3CE932E2),(u:$38F6F17A18F91149)),
    ((u:$3F3BE630D1090FE6),(u:$3F21C2DC0A9340EF),(u:$B868147B0A8785F5)),
    ((u:$40C27F67F799CCBD),(u:$3FDFFF730800F741),(u:$3C7FFFFFFFFFFFFF)),
    ((u:$4083F91B5B198CDD),(u:$3FDFF7D71D1D656B),(u:$B91C90DDE6F07A40)));

// Hardware FMA (AVX2) — for AtanpiSmall's subnormal-result fmas where
// Dekker emulation isn't correctly rounded.
function AtanpiHwFma(a, b, c: Double): Double; assembler;
asm
  vfmadd213sd xmm0, xmm1, xmm2
end;

// --- atanpi asymptotic: |x| > 0x1.bep20 ---
function AtanpiAsympt(x: Double): Double;
var h, m, l, yh, yl: Double;
    m_u: Tb64u64;
begin
  if x >= Double(0.0) then h := Double(0.5) else h := -Double(0.5);
  yh := Double(1.0) / x;
  yl := yh * pcr_fma(yh, -x, Double(1.0));
  m := pcr_muldd(yh, yl, cAtanpiOneOverPiH.f, cAtanpiOneOverPiL.f, l);
  m := -m;
  l := pcr_fma(cAtanpiOneOver3Pi.f * yh, yh * yh, -l);
  pcr_fasttwosum(h, m, h, m);
  pcr_fasttwosum(m, l, m, l);
  m_u.f := m;
  if (m_u.u and UInt64($7FFFFFFFFFFFFFFF)) = cAtanpiHalfP55.u then
  begin
    if (m * l) > Double(0.0) then
    begin
      if m >= Double(0.0) then m :=  cAtanpiAs1.f else m := -cAtanpiAs1.f;
    end
    else
    begin
      if m >= Double(0.0) then m :=  cAtanpiAs2.f else m := -cAtanpiAs2.f;
    end;
  end;
  Result := h + m;
end;

// --- atanpi tiny: |x| < 0x1.c7p-27 ---
function AtanpiTiny(x: Double): Double;
var h, l: Double; ax: Double; i: Integer;
begin
  h := x * cAtanpiOneOverPiH.f;
  l := pcr_fma(x, cAtanpiOneOverPiH.f, -h);
  l := pcr_fma(x, cAtanpiOneOverPiL.f, l);
  l := pcr_fma(-cAtanpiOneOver3Pi.f * x, x * x, l);
  ax := Abs(x);
  for i := 0 to 19 do
    if ax = cAtanpiTinyExc[i, 0].f then
    begin
      if x >= Double(0.0) then
      begin
        h :=  cAtanpiTinyExc[i, 1].f;
        l :=  cAtanpiTinyExc[i, 2].f;
      end
      else
      begin
        h := -cAtanpiTinyExc[i, 1].f;
        l := -cAtanpiTinyExc[i, 2].f;
      end;
      Break;
    end;
  Result := h + l;
end;

// --- atanpi small: |x| < 2^-54 (including 0) ---
function AtanpiSmall(x: Double): Double;
var v: Tb64u64; e: UInt64; h, corr: Double;
begin
  if x = Double(0.0) then begin Result := x; Exit; end;
  if Abs(x) = cAtanpiSmallIn.f then
  begin
    if x > Double(0.0) then
      Result := AtanpiHwFma(-cAtanpiP600.f, cAtanpiP600.f,  cAtanpiSmallOut.f)
    else
      Result := AtanpiHwFma( cAtanpiP600.f, cAtanpiP600.f, -cAtanpiSmallOut.f);
    Exit;
  end;
  v.f := x;
  if (v.u and UInt64($FFFFFFFFFFFFF)) = UInt64($59AF9A1194EFE) then
  begin
    e := v.u shr 52;
    if (e and UInt64($7FF)) > UInt64(2) then
    begin
      v.u := ((e - UInt64(2)) shl 52) or UInt64($B824198B94A89);
      if x > Double(0.0) then
        Result := AtanpiHwFma( cAtanpiP600.f, cAtanpiP600.f, v.f)
      else
        Result := AtanpiHwFma(-cAtanpiP600.f, cAtanpiP600.f, v.f);
      Exit;
    end;
  end;
  h := x * cAtanpiOneOverPiH.f;
  corr := AtanpiHwFma(x * cAtanpiPp106.f, cAtanpiOneOverPiH.f, -h * cAtanpiPp106.f);
  corr := AtanpiHwFma(x * cAtanpiPp106.f, cAtanpiOneOverPiL.f, corr);
  Result := AtanpiHwFma(corr, cAtanpiPm106.f, h);
end;

// --- atanpi refine (slow path) ---
function AtanpiRefine(x, a: Double): Double;
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
  j: Integer;
begin
  ax := Abs(x);
  if ax > cAtanpiThrBig.f then begin Result := AtanpiAsympt(x); Exit; end;
  if ax < cAtanpiThrTiny.f then begin Result := AtanpiTiny(x);   Exit; end;

  for j := 0 to 55 do
    if ax = cAtanpiRefExc[j, 0].f then
    begin
      if x >= Double(0.0) then
        Result :=  cAtanpiRefExc[j, 1].f +  cAtanpiRefExc[j, 2].f
      else
        Result := -cAtanpiRefExc[j, 1].f + -cAtanpiRefExc[j, 2].f;
      Exit;
    end;

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

  // Multiply v0+v1 (which approximates atan(x)) by 1/pi.
  v0 := pcr_muldd(v0, v1, cAtanpiOneOverPiH.f, cAtanpiOneOverPiL.f, v1);
  Result := v1 + v0;
end;

function pcr_atanpi(x: Double): Double;
var
  t: Tb64u64;
  at_u: UInt64;
  i_row: Int64;
  ii: Int64;
  u, ut, ut2: UInt64;
  ta_u, id_u: Tb64u64;
  h, ah, al_v, f_v: Double;
  h2, h4: Double;
  e, e0, ub, lb, ub0: Double;
  x2, x3, x4, fs: Double;
begin
  t.f  := x;
  at_u := t.u and UInt64($7FFFFFFFFFFFFFFF);
  i_row := Int64(at_u shr 51) - Int64(2030);

  if at_u < UInt64($3F7B21C475E6362A) then
  begin
    if at_u < UInt64($3C90000000000000) then
    begin
      Result := AtanpiSmall(x);
      Exit;
    end;
    if x = Double(0.0) then begin Result := x; Exit; end;

    x2 := x * x;
    x3 := x * x2;
    x4 := x2 * x2;
    fs := x3 * ((cAtanCh2_0.f + x2 * cAtanCh2_1.f) +
                 x4 * (cAtanCh2_2.f + x2 * cAtanCh2_3.f));
    // (x+fs)/pi; fast-path error bound 0x1.6fp-68 * x
    h := pcr_muldd(x, fs, cAtanpiOneOverPiH.f, cAtanpiOneOverPiL.f, al_v);
    ub := h + pcr_fma( cAtanpiETiny.f, x, al_v);
    lb := h + pcr_fma(-cAtanpiETiny.f, x, al_v);
    if ub = lb then begin Result := ub; Exit; end;
    // for refine: ub0 = (fs + fs*0x4.8p-52) + x (atan(x) fast approximation)
    ub := (fs + fs * cAtanFmaUb.f) + x;
    Result := AtanpiRefine(x, ub);
    Exit;
  end;

  if at_u > UInt64($4062DED8E34A9035) then
  begin
    if at_u >= UInt64($43445F306DC9C883) then
    begin
      if at_u >= (UInt64($7FF) shl 52) then
      begin
        if at_u = (UInt64($7FF) shl 52) then
        begin
          if x > Double(0.0) then Result :=  Double(0.5)
          else                    Result := -Double(0.5);
          Exit;
        end;
        Result := x + x; // NaN
        Exit;
      end;
      if x > Double(0.0) then
        Result :=  Double(0.5) - cAtanpiPm56.f
      else
        Result := -Double(0.5) + cAtanpiPm56.f;
      Exit;
    end;
    h    := -Double(1.0) / x;
    ah   := cAtanPiHalfH.f;
    al_v := cAtanPiHalfL.f;
    if x < Double(0.0) then begin ah := -ah; al_v := -al_v; end;
  end
  else
  begin
    if at_u = UInt64($3FF0000000000000) then
    begin
      if x > Double(0.0) then Result :=  Double(0.25) else Result := -Double(0.25);
      Exit;
    end;
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
  // atan-space estimate for refine
  e0   := h * cAtanEFactor.f;
  ub0  := (al_v + e0) + ah;
  // ah + al_v ~ atan(x); normalize and multiply by 1/pi
  pcr_fasttwosum(ah, al_v, ah, al_v);
  ah   := pcr_muldd(ah, al_v, cAtanpiOneOverPiH.f, cAtanpiOneOverPiL.f, al_v);
  e    := h * cAtanpiEFast.f;
  ub   := (al_v + e) + ah;
  lb   := (al_v - e) + ah;
  if ub = lb then begin Result := ub; Exit; end;
  Result := AtanpiRefine(x, ub0);
end;
