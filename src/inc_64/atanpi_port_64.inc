// ---------------------------------------------------------------------------
// pcr_atanpi — correctly-rounded atan(x)/pi (binary64).
// Ported from CORE-MATH cr_atanpi (Sibidanov & Zimmermann, 2023-2025).
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
// Reuses cAtanAHi/Lo, cAtanC0/1/2, cAtanIdHi/Lo/Lo2/Lo3, cAtanPhiScale,
// cAtanRefCh*/Cl*, cAtanCh*/Ch2*, cAtanPiHalfH/L from the atan port above.
// ---------------------------------------------------------------------------

const
  cAtanpiOneOverPiH: Tb64u64 = (u:$3FD45F306DC9C883); //  0x1.45f306dc9c883p-2
  cAtanpiOneOverPiL: Tb64u64 = (u:QWord($BC76B01EC5417056)); // -0x1.6b01ec5417056p-56
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

  // Pillar B (Phase 6.3): cAtanpiTinyExc lifted to named scalars; the
  // for i := 0 to 19 exception scan in AtanpiTiny is unrolled.
  cAtanpiTinyExc_0_0: Tb64u64 = (u:$3DA9B46CFB0F10CA);
  cAtanpiTinyExc_0_1: Tb64u64 = (u:$3D905D3AA746E163);
  cAtanpiTinyExc_0_2: Tb64u64 = (u:QWord($B6A94A704A0E7980));
  cAtanpiTinyExc_1_0: Tb64u64 = (u:$3E0EAC33F52580C1);
  cAtanpiTinyExc_1_1: Tb64u64 = (u:$3DF386DE1F01225A);
  cAtanpiTinyExc_1_2: Tb64u64 = (u:QWord($B742BE95E5723A66));
  cAtanpiTinyExc_2_0: Tb64u64 = (u:$3CA56A4AA740A5A7);
  cAtanpiTinyExc_2_1: Tb64u64 = (u:$3C8B44453E2404E7);
  cAtanpiTinyExc_2_2: Tb64u64 = (u:QWord($B5BC5376C33C7A50));
  cAtanpiTinyExc_3_0: Tb64u64 = (u:$3D99D47B2E2A4DEE);
  cAtanpiTinyExc_3_1: Tb64u64 = (u:$3D8071A2E1E2E7A0);
  cAtanpiTinyExc_3_2: Tb64u64 = (u:QWord($B6DE948955ED137A));
  cAtanpiTinyExc_4_0: Tb64u64 = (u:$3D3932DDAEECDABE);
  cAtanpiTinyExc_4_1: Tb64u64 = (u:$3D200ABFB23406C5);
  cAtanpiTinyExc_4_2: Tb64u64 = (u:QWord($B9CFFFFFFFFFFFFF));
  cAtanpiTinyExc_5_0: Tb64u64 = (u:$3D1E0B9044E5A1CF);
  cAtanpiTinyExc_5_1: Tb64u64 = (u:$3D032099F7C94E82);
  cAtanpiTinyExc_5_2: Tb64u64 = (u:QWord($B6722752CF9B0BC5));
  cAtanpiTinyExc_6_0: Tb64u64 = (u:$3CD13B412EA2182D);
  cAtanpiTinyExc_6_1: Tb64u64 = (u:$3CB5F095755008F0);
  cAtanpiTinyExc_6_2: Tb64u64 = (u:QWord($B5FA8AD4502B5F10));
  cAtanpiTinyExc_7_0: Tb64u64 = (u:$3D4FF878736A1BCD);
  cAtanpiTinyExc_7_1: Tb64u64 = (u:$3D345A65501D94D3);
  cAtanpiTinyExc_7_2: Tb64u64 = (u:QWord($B677CB6951D7806D));
  cAtanpiTinyExc_8_0: Tb64u64 = (u:$3D94F03BD4CEDD49);
  cAtanpiTinyExc_8_1: Tb64u64 = (u:$3D7AA8DC8DA24739);
  cAtanpiTinyExc_8_2: Tb64u64 = (u:QWord($B6DFE11160EC4EE9));
  cAtanpiTinyExc_9_0: Tb64u64 = (u:$3DB0F2BF03E38A9E);
  cAtanpiTinyExc_9_1: Tb64u64 = (u:$3D9594435F8C0326);
  cAtanpiTinyExc_9_2: Tb64u64 = (u:QWord($B6E8EEF034A3526E));
  cAtanpiTinyExc_10_0: Tb64u64 = (u:$3CB508ECB38F52F9);
  cAtanpiTinyExc_10_1: Tb64u64 = (u:$3C9AC84C88F979A2);
  cAtanpiTinyExc_10_2: Tb64u64 = (u:QWord($B5B63773CBD64D5C));
  cAtanpiTinyExc_11_0: Tb64u64 = (u:$3D58DDCD3207FE60);
  cAtanpiTinyExc_11_1: Tb64u64 = (u:$3D3FA930B4EB436C);
  cAtanpiTinyExc_11_2: Tb64u64 = (u:$39E0000000000000);
  cAtanpiTinyExc_12_0: Tb64u64 = (u:$3CF1331B44547444);
  cAtanpiTinyExc_12_1: Tb64u64 = (u:$3CD5E63596A62718);
  cAtanpiTinyExc_12_2: Tb64u64 = (u:QWord($B634A9F744D7CEA7));
  cAtanpiTinyExc_13_0: Tb64u64 = (u:$3E0F7681FA730E73);
  cAtanpiTinyExc_13_1: Tb64u64 = (u:$3DF407A8AC1C768A);
  cAtanpiTinyExc_13_2: Tb64u64 = (u:QWord($B753C6A80C812C54));
  cAtanpiTinyExc_14_0: Tb64u64 = (u:$3DB09DF162E4049C);
  cAtanpiTinyExc_14_1: Tb64u64 = (u:$3D952849D0B1E523);
  cAtanpiTinyExc_14_2: Tb64u64 = (u:QWord($BA3FFFFFFFFFFFFF));
  cAtanpiTinyExc_15_0: Tb64u64 = (u:$3C959AF9A1194EFE);
  cAtanpiTinyExc_15_1: Tb64u64 = (u:$3C7B824198B94A89);
  cAtanpiTinyExc_15_2: Tb64u64 = (u:$391FFFFFFFFFFFFF);
  cAtanpiTinyExc_16_0: Tb64u64 = (u:$3D165F203D6BB18E);
  cAtanpiTinyExc_16_1: Tb64u64 = (u:$3CFC7C00D31634A6);
  cAtanpiTinyExc_16_2: Tb64u64 = (u:QWord($B62D4AF6699E375F));
  cAtanpiTinyExc_17_0: Tb64u64 = (u:$3DD6F89F6DF8258D);
  cAtanpiTinyExc_17_1: Tb64u64 = (u:$3DBD3F710390DB87);
  cAtanpiTinyExc_17_2: Tb64u64 = (u:QWord($B6E936061594FAD7));
  cAtanpiTinyExc_18_0: Tb64u64 = (u:$3CC3E4D2D87B5AEF);
  cAtanpiTinyExc_18_1: Tb64u64 = (u:$3CA954626979D7D3);
  cAtanpiTinyExc_18_2: Tb64u64 = (u:QWord($B5E6F8DBD4F38CDB));
  cAtanpiTinyExc_19_0: Tb64u64 = (u:$3DBF46E9E25A51D1);
  cAtanpiTinyExc_19_1: Tb64u64 = (u:$3DA3E95C097E08E0);
  cAtanpiTinyExc_19_2: Tb64u64 = (u:QWord($B6F8E5F60920B571));

  // Pillar B (Phase 6.3): cAtanpiRefExc lifted to named scalars; the
  // for j := 0 to 55 exception scan in AtanpiRefine is unrolled.
  cAtanpiRefExc_0_0: Tb64u64 = (u:$40116265F83E8121);
  cAtanpiRefExc_0_1: Tb64u64 = (u:$3FDB648E5B0804DA);
  cAtanpiRefExc_0_2: Tb64u64 = (u:$3929DB921C686893);
  cAtanpiRefExc_1_0: Tb64u64 = (u:$4058A5C063607FBA);
  cAtanpiRefExc_1_1: Tb64u64 = (u:$3FDFCB1A9F45A781);
  cAtanpiRefExc_1_2: Tb64u64 = (u:QWord($BC7FFFFFFFFFFFFF));
  cAtanpiRefExc_2_0: Tb64u64 = (u:$3F91A5E744FF1806);
  cAtanpiRefExc_2_1: Tb64u64 = (u:$3F7677CDC7DEA093);
  cAtanpiRefExc_2_2: Tb64u64 = (u:$38C417D05C2A4506);
  cAtanpiRefExc_3_0: Tb64u64 = (u:$3F82EEE26511C903);
  cAtanpiRefExc_3_1: Tb64u64 = (u:$3F681B11AC013890);
  cAtanpiRefExc_3_2: Tb64u64 = (u:$38B508C83C3DB5E3);
  cAtanpiRefExc_4_0: Tb64u64 = (u:$3F7BB0BA10084DE6);
  cAtanpiRefExc_4_1: Tb64u64 = (u:$3F61A0C1411A3DC1);
  cAtanpiRefExc_4_2: Tb64u64 = (u:$38CA854B560FE4DB);
  cAtanpiRefExc_5_0: Tb64u64 = (u:$3FD55166A4D7D836);
  cAtanpiRefExc_5_1: Tb64u64 = (u:$3FBA33741ABEA942);
  cAtanpiRefExc_5_2: Tb64u64 = (u:QWord($B9154A4C868A708E));
  cAtanpiRefExc_6_0: Tb64u64 = (u:$3FC5FA4B5E7CE940);
  cAtanpiRefExc_6_1: Tb64u64 = (u:$3FABB6715A0B1F81);
  cAtanpiRefExc_6_2: Tb64u64 = (u:$39000037A0BFDFAE);
  cAtanpiRefExc_7_0: Tb64u64 = (u:$4032A04FA64D3600);
  cAtanpiRefExc_7_1: Tb64u64 = (u:$3FDEE846D3680301);
  cAtanpiRefExc_7_2: Tb64u64 = (u:$3903015B69C79EF4);
  cAtanpiRefExc_8_0: Tb64u64 = (u:$4030FB7E1298320C);
  cAtanpiRefExc_8_1: Tb64u64 = (u:$3FDECD42A2EF8667);
  cAtanpiRefExc_8_2: Tb64u64 = (u:QWord($B9421D6585C13D10));
  cAtanpiRefExc_9_0: Tb64u64 = (u:$402B6DAA66352E75);
  cAtanpiRefExc_9_1: Tb64u64 = (u:$3FDE84651D4C9934);
  cAtanpiRefExc_9_2: Tb64u64 = (u:$39334E0B5A33A1DB);
  cAtanpiRefExc_10_0: Tb64u64 = (u:$4005717DF9ABCCD3);
  cAtanpiRefExc_10_1: Tb64u64 = (u:$3FD8B9BE647E02A9);
  cAtanpiRefExc_10_2: Tb64u64 = (u:QWord($B926C5194938E976));
  cAtanpiRefExc_11_0: Tb64u64 = (u:$3FBE0AC5527E7F0D);
  cAtanpiRefExc_11_1: Tb64u64 = (u:$3FA309CE21202CCA);
  cAtanpiRefExc_11_2: Tb64u64 = (u:QWord($B90141BE39E1919D));
  cAtanpiRefExc_12_0: Tb64u64 = (u:$402932F168987539);
  cAtanpiRefExc_12_1: Tb64u64 = (u:$3FDE62F20ABA165C);
  cAtanpiRefExc_12_2: Tb64u64 = (u:QWord($B92839C5638718F2));
  cAtanpiRefExc_13_0: Tb64u64 = (u:$4022E342CA2440EE);
  cAtanpiRefExc_13_1: Tb64u64 = (u:$3FDDD9D1ECFFC9A0);
  cAtanpiRefExc_13_2: Tb64u64 = (u:QWord($B9253DEDDDE1DDE1));
  cAtanpiRefExc_14_0: Tb64u64 = (u:$3FAB37EF5831A537);
  cAtanpiRefExc_14_1: Tb64u64 = (u:$3F914FC0E0C2DACE);
  cAtanpiRefExc_14_2: Tb64u64 = (u:QWord($B8EF6B67EF93B398));
  cAtanpiRefExc_15_0: Tb64u64 = (u:$3F7D8EA551D1D795);
  cAtanpiRefExc_15_1: Tb64u64 = (u:$3F62D0FE21DF78D0);
  cAtanpiRefExc_15_2: Tb64u64 = (u:QWord($B8D0F4D438D8DD24));
  cAtanpiRefExc_16_0: Tb64u64 = (u:$3FCBADB08CD8B1C5);
  cAtanpiRefExc_16_1: Tb64u64 = (u:$3FB15A7DB07DAA97);
  cAtanpiRefExc_16_2: Tb64u64 = (u:QWord($BC5FFFFFFFFFFFFF));
  cAtanpiRefExc_17_0: Tb64u64 = (u:$40618071548E23A6);
  cAtanpiRefExc_17_1: Tb64u64 = (u:$3FDFDAC0C22A86C5);
  cAtanpiRefExc_17_2: Tb64u64 = (u:$3C7FFFFFFFFFFFFF);
  cAtanpiRefExc_18_0: Tb64u64 = (u:$405E8B9C1FE81CA3);
  cAtanpiRefExc_18_1: Tb64u64 = (u:$3FDFD55124D78891);
  cAtanpiRefExc_18_2: Tb64u64 = (u:QWord($BC7FFFFFFFFFFFFE));
  cAtanpiRefExc_19_0: Tb64u64 = (u:$3F9E2FAAB41052E7);
  cAtanpiRefExc_19_1: Tb64u64 = (u:$3F833629443D32F5);
  cAtanpiRefExc_19_2: Tb64u64 = (u:$3C2FFFFFFFFFFFFF);
  cAtanpiRefExc_20_0: Tb64u64 = (u:$3F803F709107440D);
  cAtanpiRefExc_20_1: Tb64u64 = (u:$3F64AFDA1E054259);
  cAtanpiRefExc_20_2: Tb64u64 = (u:$3C0FFFFFFFFFFFFE);
  cAtanpiRefExc_21_0: Tb64u64 = (u:$3F0928A1CA94D1DC);
  cAtanpiRefExc_21_1: Tb64u64 = (u:$3EF0043BD1FD5BAE);
  cAtanpiRefExc_21_2: Tb64u64 = (u:QWord($B8512E31518BAD19));
  cAtanpiRefExc_22_0: Tb64u64 = (u:$40ADEDA0DDB0795C);
  cAtanpiRefExc_22_1: Tb64u64 = (u:$3FDFFEA37CAE8F37);
  cAtanpiRefExc_22_2: Tb64u64 = (u:QWord($BC7FFFFFFFFFFFFE));
  cAtanpiRefExc_23_0: Tb64u64 = (u:$3F472F4CB85B8404);
  cAtanpiRefExc_23_1: Tb64u64 = (u:$3F2D850E9850AB33);
  cAtanpiRefExc_23_2: Tb64u64 = (u:QWord($B88014404207EB42));
  cAtanpiRefExc_24_0: Tb64u64 = (u:$3EA4B3C41E10A1D4);
  cAtanpiRefExc_24_1: Tb64u64 = (u:$3E8A5BDF2B61615B);
  cAtanpiRefExc_24_2: Tb64u64 = (u:QWord($B7C180C42F2C36BB));
  cAtanpiRefExc_25_0: Tb64u64 = (u:$3F64AD2F60D5240B);
  cAtanpiRefExc_25_1: Tb64u64 = (u:$3F4A537A6DE092A1);
  cAtanpiRefExc_25_2: Tb64u64 = (u:QWord($B8AF2E112F36ECA4));
  cAtanpiRefExc_26_0: Tb64u64 = (u:$4131151841058E68);
  cAtanpiRefExc_26_1: Tb64u64 = (u:$3FDFFFFECEB437E6);
  cAtanpiRefExc_26_2: Tb64u64 = (u:QWord($B932DD2F50460971));
  cAtanpiRefExc_27_0: Tb64u64 = (u:$3EB8B012D530C5C8);
  cAtanpiRefExc_27_1: Tb64u64 = (u:$3E9F6EF7B3EE090C);
  cAtanpiRefExc_27_2: Tb64u64 = (u:QWord($B7F9D79D9FB153C4));
  cAtanpiRefExc_28_0: Tb64u64 = (u:$4080FA189169E75E);
  cAtanpiRefExc_28_1: Tb64u64 = (u:$3FDFF666764F1751);
  cAtanpiRefExc_28_2: Tb64u64 = (u:QWord($B93DA2823E32E295));
  cAtanpiRefExc_29_0: Tb64u64 = (u:$40646C64E8717CCF);
  cAtanpiRefExc_29_1: Tb64u64 = (u:$3FDFE014CAC09E73);
  cAtanpiRefExc_29_2: Tb64u64 = (u:$3C7FFFFFFFFFFFFF);
  cAtanpiRefExc_30_0: Tb64u64 = (u:$40B5A9EE4B581471);
  cAtanpiRefExc_30_1: Tb64u64 = (u:$3FDFFF0F445BEE08);
  cAtanpiRefExc_30_2: Tb64u64 = (u:QWord($B934256EAE14A448));
  cAtanpiRefExc_31_0: Tb64u64 = (u:$40B6AB7E3B9E35CD);
  cAtanpiRefExc_31_1: Tb64u64 = (u:$3FDFFF19F36D9901);
  cAtanpiRefExc_31_2: Tb64u64 = (u:QWord($B94246ECD835D1C3));
  cAtanpiRefExc_32_0: Tb64u64 = (u:$3F2842E540B07AC1);
  cAtanpiRefExc_32_1: Tb64u64 = (u:$3F0EE3F52B4AF0DC);
  cAtanpiRefExc_32_2: Tb64u64 = (u:QWord($B868C285BFA8C0EC));
  cAtanpiRefExc_33_0: Tb64u64 = (u:$4115D2C01B65F659);
  cAtanpiRefExc_33_1: Tb64u64 = (u:$3FDFFFFC441A900D);
  cAtanpiRefExc_33_2: Tb64u64 = (u:QWord($B948137CC73DD777));
  cAtanpiRefExc_34_0: Tb64u64 = (u:$4109D18F49824D6A);
  cAtanpiRefExc_34_1: Tb64u64 = (u:$3FDFFFF9B00D35D1);
  cAtanpiRefExc_34_2: Tb64u64 = (u:$391FE1E36287467B);
  cAtanpiRefExc_35_0: Tb64u64 = (u:$413400B2B748DA28);
  cAtanpiRefExc_35_1: Tb64u64 = (u:$3FDFFFFEFB46ADE3);
  cAtanpiRefExc_35_2: Tb64u64 = (u:QWord($B9010CA50F3BFEF4));
  cAtanpiRefExc_36_0: Tb64u64 = (u:$41327B786B2F1BAE);
  cAtanpiRefExc_36_1: Tb64u64 = (u:$3FDFFFFEE5D3F843);
  cAtanpiRefExc_36_2: Tb64u64 = (u:QWord($B938AFEE3BDFB07B));
  cAtanpiRefExc_37_0: Tb64u64 = (u:$3EBBBC4A4251B4B9);
  cAtanpiRefExc_37_1: Tb64u64 = (u:$3EA1A82F5F7F86BD);
  cAtanpiRefExc_37_2: Tb64u64 = (u:$3B4FFFFFFFFFFFFF);
  cAtanpiRefExc_38_0: Tb64u64 = (u:$3F472A01AEF844E5);
  cAtanpiRefExc_38_1: Tb64u64 = (u:$3F2D7E514F3FFE99);
  cAtanpiRefExc_38_2: Tb64u64 = (u:$3887C2D5BE2B88D8);
  cAtanpiRefExc_39_0: Tb64u64 = (u:$40C95E2B447C5ED5);
  cAtanpiRefExc_39_1: Tb64u64 = (u:$3FDFFF99356BF911);
  cAtanpiRefExc_39_2: Tb64u64 = (u:$3C7FFFFFFFFFFFFF);
  cAtanpiRefExc_40_0: Tb64u64 = (u:$41288D8109598821);
  cAtanpiRefExc_40_1: Tb64u64 = (u:$3FDFFFFE572F576A);
  cAtanpiRefExc_40_2: Tb64u64 = (u:$3935430C0B185E51);
  cAtanpiRefExc_41_0: Tb64u64 = (u:$3E83EF59F63222CC);
  cAtanpiRefExc_41_1: Tb64u64 = (u:$3E6961C9F0BC5519);
  cAtanpiRefExc_41_2: Tb64u64 = (u:$3B0FFFFFFFFFFFFF);
  cAtanpiRefExc_42_0: Tb64u64 = (u:$41172F2198034E17);
  cAtanpiRefExc_42_1: Tb64u64 = (u:$3FDFFFFC7C36896E);
  cAtanpiRefExc_42_2: Tb64u64 = (u:QWord($B92ADFC30E218AAF));
  cAtanpiRefExc_43_0: Tb64u64 = (u:$3E4DF1C8E36BBCC0);
  cAtanpiRefExc_43_1: Tb64u64 = (u:$3E331030ABF0398F);
  cAtanpiRefExc_43_2: Tb64u64 = (u:$3ADFFFFFFFFFFFFE);
  cAtanpiRefExc_44_0: Tb64u64 = (u:$408AD5781AB29BBA);
  cAtanpiRefExc_44_1: Tb64u64 = (u:$3FDFF9ED31438683);
  cAtanpiRefExc_44_2: Tb64u64 = (u:$3C7FFFFFFFFFFFFF);
  cAtanpiRefExc_45_0: Tb64u64 = (u:$40AEDE23D6624A3A);
  cAtanpiRefExc_45_1: Tb64u64 = (u:$3FDFFEAE182F312E);
  cAtanpiRefExc_45_2: Tb64u64 = (u:QWord($B9426EEE3B0A2CCB));
  cAtanpiRefExc_46_0: Tb64u64 = (u:$4124A7E234CDCE9B);
  cAtanpiRefExc_46_1: Tb64u64 = (u:$3FDFFFFE0709E4E1);
  cAtanpiRefExc_46_2: Tb64u64 = (u:QWord($BC7FFFFFFFFFFFFD));
  cAtanpiRefExc_47_0: Tb64u64 = (u:$3E81C404899CA420);
  cAtanpiRefExc_47_1: Tb64u64 = (u:$3E669EB74CE9E7ED);
  cAtanpiRefExc_47_2: Tb64u64 = (u:$37CD126A89BF70E8);
  cAtanpiRefExc_48_0: Tb64u64 = (u:$4077422A83423830);
  cAtanpiRefExc_48_1: Tb64u64 = (u:$3FDFF1FC5DEB1C5C);
  cAtanpiRefExc_48_2: Tb64u64 = (u:QWord($B930CF8B48B700A5));
  cAtanpiRefExc_49_0: Tb64u64 = (u:$3F0F37BE8C1F5DF0);
  cAtanpiRefExc_49_1: Tb64u64 = (u:$3EF3DFB3D39744CF);
  cAtanpiRefExc_49_2: Tb64u64 = (u:QWord($B84CAE8F0C659368));
  cAtanpiRefExc_50_0: Tb64u64 = (u:$4138C9F8999242B7);
  cAtanpiRefExc_50_1: Tb64u64 = (u:$3FDFFFFF2D9DCA7B);
  cAtanpiRefExc_50_2: Tb64u64 = (u:$3C7FFFFFFFFFFFFF);
  cAtanpiRefExc_51_0: Tb64u64 = (u:$412427F51DCF1D68);
  cAtanpiRefExc_51_1: Tb64u64 = (u:$3FDFFFFDFA850495);
  cAtanpiRefExc_51_2: Tb64u64 = (u:$3C7FFFFFFFFFFFFD);
  cAtanpiRefExc_52_0: Tb64u64 = (u:$413ABB7B457CDF84);
  cAtanpiRefExc_52_1: Tb64u64 = (u:$3FDFFFFF3CE932E2);
  cAtanpiRefExc_52_2: Tb64u64 = (u:$38F6F17A18F91149);
  cAtanpiRefExc_53_0: Tb64u64 = (u:$3F3BE630D1090FE6);
  cAtanpiRefExc_53_1: Tb64u64 = (u:$3F21C2DC0A9340EF);
  cAtanpiRefExc_53_2: Tb64u64 = (u:QWord($B868147B0A8785F5));
  cAtanpiRefExc_54_0: Tb64u64 = (u:$40C27F67F799CCBD);
  cAtanpiRefExc_54_1: Tb64u64 = (u:$3FDFFF730800F741);
  cAtanpiRefExc_54_2: Tb64u64 = (u:$3C7FFFFFFFFFFFFF);
  cAtanpiRefExc_55_0: Tb64u64 = (u:$4083F91B5B198CDD);
  cAtanpiRefExc_55_1: Tb64u64 = (u:$3FDFF7D71D1D656B);
  cAtanpiRefExc_55_2: Tb64u64 = (u:QWord($B91C90DDE6F07A40));

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
var h, l: Double; ax: Double;
begin
  h := x * cAtanpiOneOverPiH.f;
  l := pcr_fma(x, cAtanpiOneOverPiH.f, -h);
  l := pcr_fma(x, cAtanpiOneOverPiL.f, l);
  l := pcr_fma(-cAtanpiOneOver3Pi.f * x, x * x, l);
  ax := Abs(x);
  // Pillar A (Phase 6.3): exception-db scan unrolled (20 entries).
  if x >= Double(0.0) then
  begin
    if ax = cAtanpiTinyExc_0_0.f then begin h := cAtanpiTinyExc_0_1.f; l := cAtanpiTinyExc_0_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_1_0.f then begin h := cAtanpiTinyExc_1_1.f; l := cAtanpiTinyExc_1_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_2_0.f then begin h := cAtanpiTinyExc_2_1.f; l := cAtanpiTinyExc_2_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_3_0.f then begin h := cAtanpiTinyExc_3_1.f; l := cAtanpiTinyExc_3_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_4_0.f then begin h := cAtanpiTinyExc_4_1.f; l := cAtanpiTinyExc_4_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_5_0.f then begin h := cAtanpiTinyExc_5_1.f; l := cAtanpiTinyExc_5_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_6_0.f then begin h := cAtanpiTinyExc_6_1.f; l := cAtanpiTinyExc_6_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_7_0.f then begin h := cAtanpiTinyExc_7_1.f; l := cAtanpiTinyExc_7_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_8_0.f then begin h := cAtanpiTinyExc_8_1.f; l := cAtanpiTinyExc_8_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_9_0.f then begin h := cAtanpiTinyExc_9_1.f; l := cAtanpiTinyExc_9_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_10_0.f then begin h := cAtanpiTinyExc_10_1.f; l := cAtanpiTinyExc_10_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_11_0.f then begin h := cAtanpiTinyExc_11_1.f; l := cAtanpiTinyExc_11_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_12_0.f then begin h := cAtanpiTinyExc_12_1.f; l := cAtanpiTinyExc_12_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_13_0.f then begin h := cAtanpiTinyExc_13_1.f; l := cAtanpiTinyExc_13_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_14_0.f then begin h := cAtanpiTinyExc_14_1.f; l := cAtanpiTinyExc_14_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_15_0.f then begin h := cAtanpiTinyExc_15_1.f; l := cAtanpiTinyExc_15_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_16_0.f then begin h := cAtanpiTinyExc_16_1.f; l := cAtanpiTinyExc_16_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_17_0.f then begin h := cAtanpiTinyExc_17_1.f; l := cAtanpiTinyExc_17_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_18_0.f then begin h := cAtanpiTinyExc_18_1.f; l := cAtanpiTinyExc_18_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_19_0.f then begin h := cAtanpiTinyExc_19_1.f; l := cAtanpiTinyExc_19_2.f; Result := h + l; Exit; end;
  end
  else
  begin
    if ax = cAtanpiTinyExc_0_0.f then begin h := -cAtanpiTinyExc_0_1.f; l := -cAtanpiTinyExc_0_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_1_0.f then begin h := -cAtanpiTinyExc_1_1.f; l := -cAtanpiTinyExc_1_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_2_0.f then begin h := -cAtanpiTinyExc_2_1.f; l := -cAtanpiTinyExc_2_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_3_0.f then begin h := -cAtanpiTinyExc_3_1.f; l := -cAtanpiTinyExc_3_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_4_0.f then begin h := -cAtanpiTinyExc_4_1.f; l := -cAtanpiTinyExc_4_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_5_0.f then begin h := -cAtanpiTinyExc_5_1.f; l := -cAtanpiTinyExc_5_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_6_0.f then begin h := -cAtanpiTinyExc_6_1.f; l := -cAtanpiTinyExc_6_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_7_0.f then begin h := -cAtanpiTinyExc_7_1.f; l := -cAtanpiTinyExc_7_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_8_0.f then begin h := -cAtanpiTinyExc_8_1.f; l := -cAtanpiTinyExc_8_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_9_0.f then begin h := -cAtanpiTinyExc_9_1.f; l := -cAtanpiTinyExc_9_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_10_0.f then begin h := -cAtanpiTinyExc_10_1.f; l := -cAtanpiTinyExc_10_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_11_0.f then begin h := -cAtanpiTinyExc_11_1.f; l := -cAtanpiTinyExc_11_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_12_0.f then begin h := -cAtanpiTinyExc_12_1.f; l := -cAtanpiTinyExc_12_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_13_0.f then begin h := -cAtanpiTinyExc_13_1.f; l := -cAtanpiTinyExc_13_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_14_0.f then begin h := -cAtanpiTinyExc_14_1.f; l := -cAtanpiTinyExc_14_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_15_0.f then begin h := -cAtanpiTinyExc_15_1.f; l := -cAtanpiTinyExc_15_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_16_0.f then begin h := -cAtanpiTinyExc_16_1.f; l := -cAtanpiTinyExc_16_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_17_0.f then begin h := -cAtanpiTinyExc_17_1.f; l := -cAtanpiTinyExc_17_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_18_0.f then begin h := -cAtanpiTinyExc_18_1.f; l := -cAtanpiTinyExc_18_2.f; Result := h + l; Exit; end;
    if ax = cAtanpiTinyExc_19_0.f then begin h := -cAtanpiTinyExc_19_1.f; l := -cAtanpiTinyExc_19_2.f; Result := h + l; Exit; end;
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
      Result := pcr_fma(-cAtanpiP600.f, cAtanpiP600.f,  cAtanpiSmallOut.f)
    else
      Result := pcr_fma( cAtanpiP600.f, cAtanpiP600.f, -cAtanpiSmallOut.f);
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
        Result := pcr_fma( cAtanpiP600.f, cAtanpiP600.f, v.f)
      else
        Result := pcr_fma(-cAtanpiP600.f, cAtanpiP600.f, v.f);
      Exit;
    end;
  end;
  h := x * cAtanpiOneOverPiH.f;
  corr := pcr_fma(x * cAtanpiPp106.f, cAtanpiOneOverPiH.f, -h * cAtanpiPp106.f);
  corr := pcr_fma(x * cAtanpiPp106.f, cAtanpiOneOverPiL.f, corr);
  Result := pcr_fma(corr, cAtanpiPm106.f, h);
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
begin
  ax := Abs(x);
  if ax > cAtanpiThrBig.f then begin Result := AtanpiAsympt(x); Exit; end;
  if ax < cAtanpiThrTiny.f then begin Result := AtanpiTiny(x);   Exit; end;

  // Pillar A (Phase 6.3): exception-db scan unrolled (56 entries).
  if x >= Double(0.0) then
  begin
    if ax = cAtanpiRefExc_0_0.f then begin Result := cAtanpiRefExc_0_1.f + cAtanpiRefExc_0_2.f; Exit; end;
    if ax = cAtanpiRefExc_1_0.f then begin Result := cAtanpiRefExc_1_1.f + cAtanpiRefExc_1_2.f; Exit; end;
    if ax = cAtanpiRefExc_2_0.f then begin Result := cAtanpiRefExc_2_1.f + cAtanpiRefExc_2_2.f; Exit; end;
    if ax = cAtanpiRefExc_3_0.f then begin Result := cAtanpiRefExc_3_1.f + cAtanpiRefExc_3_2.f; Exit; end;
    if ax = cAtanpiRefExc_4_0.f then begin Result := cAtanpiRefExc_4_1.f + cAtanpiRefExc_4_2.f; Exit; end;
    if ax = cAtanpiRefExc_5_0.f then begin Result := cAtanpiRefExc_5_1.f + cAtanpiRefExc_5_2.f; Exit; end;
    if ax = cAtanpiRefExc_6_0.f then begin Result := cAtanpiRefExc_6_1.f + cAtanpiRefExc_6_2.f; Exit; end;
    if ax = cAtanpiRefExc_7_0.f then begin Result := cAtanpiRefExc_7_1.f + cAtanpiRefExc_7_2.f; Exit; end;
    if ax = cAtanpiRefExc_8_0.f then begin Result := cAtanpiRefExc_8_1.f + cAtanpiRefExc_8_2.f; Exit; end;
    if ax = cAtanpiRefExc_9_0.f then begin Result := cAtanpiRefExc_9_1.f + cAtanpiRefExc_9_2.f; Exit; end;
    if ax = cAtanpiRefExc_10_0.f then begin Result := cAtanpiRefExc_10_1.f + cAtanpiRefExc_10_2.f; Exit; end;
    if ax = cAtanpiRefExc_11_0.f then begin Result := cAtanpiRefExc_11_1.f + cAtanpiRefExc_11_2.f; Exit; end;
    if ax = cAtanpiRefExc_12_0.f then begin Result := cAtanpiRefExc_12_1.f + cAtanpiRefExc_12_2.f; Exit; end;
    if ax = cAtanpiRefExc_13_0.f then begin Result := cAtanpiRefExc_13_1.f + cAtanpiRefExc_13_2.f; Exit; end;
    if ax = cAtanpiRefExc_14_0.f then begin Result := cAtanpiRefExc_14_1.f + cAtanpiRefExc_14_2.f; Exit; end;
    if ax = cAtanpiRefExc_15_0.f then begin Result := cAtanpiRefExc_15_1.f + cAtanpiRefExc_15_2.f; Exit; end;
    if ax = cAtanpiRefExc_16_0.f then begin Result := cAtanpiRefExc_16_1.f + cAtanpiRefExc_16_2.f; Exit; end;
    if ax = cAtanpiRefExc_17_0.f then begin Result := cAtanpiRefExc_17_1.f + cAtanpiRefExc_17_2.f; Exit; end;
    if ax = cAtanpiRefExc_18_0.f then begin Result := cAtanpiRefExc_18_1.f + cAtanpiRefExc_18_2.f; Exit; end;
    if ax = cAtanpiRefExc_19_0.f then begin Result := cAtanpiRefExc_19_1.f + cAtanpiRefExc_19_2.f; Exit; end;
    if ax = cAtanpiRefExc_20_0.f then begin Result := cAtanpiRefExc_20_1.f + cAtanpiRefExc_20_2.f; Exit; end;
    if ax = cAtanpiRefExc_21_0.f then begin Result := cAtanpiRefExc_21_1.f + cAtanpiRefExc_21_2.f; Exit; end;
    if ax = cAtanpiRefExc_22_0.f then begin Result := cAtanpiRefExc_22_1.f + cAtanpiRefExc_22_2.f; Exit; end;
    if ax = cAtanpiRefExc_23_0.f then begin Result := cAtanpiRefExc_23_1.f + cAtanpiRefExc_23_2.f; Exit; end;
    if ax = cAtanpiRefExc_24_0.f then begin Result := cAtanpiRefExc_24_1.f + cAtanpiRefExc_24_2.f; Exit; end;
    if ax = cAtanpiRefExc_25_0.f then begin Result := cAtanpiRefExc_25_1.f + cAtanpiRefExc_25_2.f; Exit; end;
    if ax = cAtanpiRefExc_26_0.f then begin Result := cAtanpiRefExc_26_1.f + cAtanpiRefExc_26_2.f; Exit; end;
    if ax = cAtanpiRefExc_27_0.f then begin Result := cAtanpiRefExc_27_1.f + cAtanpiRefExc_27_2.f; Exit; end;
    if ax = cAtanpiRefExc_28_0.f then begin Result := cAtanpiRefExc_28_1.f + cAtanpiRefExc_28_2.f; Exit; end;
    if ax = cAtanpiRefExc_29_0.f then begin Result := cAtanpiRefExc_29_1.f + cAtanpiRefExc_29_2.f; Exit; end;
    if ax = cAtanpiRefExc_30_0.f then begin Result := cAtanpiRefExc_30_1.f + cAtanpiRefExc_30_2.f; Exit; end;
    if ax = cAtanpiRefExc_31_0.f then begin Result := cAtanpiRefExc_31_1.f + cAtanpiRefExc_31_2.f; Exit; end;
    if ax = cAtanpiRefExc_32_0.f then begin Result := cAtanpiRefExc_32_1.f + cAtanpiRefExc_32_2.f; Exit; end;
    if ax = cAtanpiRefExc_33_0.f then begin Result := cAtanpiRefExc_33_1.f + cAtanpiRefExc_33_2.f; Exit; end;
    if ax = cAtanpiRefExc_34_0.f then begin Result := cAtanpiRefExc_34_1.f + cAtanpiRefExc_34_2.f; Exit; end;
    if ax = cAtanpiRefExc_35_0.f then begin Result := cAtanpiRefExc_35_1.f + cAtanpiRefExc_35_2.f; Exit; end;
    if ax = cAtanpiRefExc_36_0.f then begin Result := cAtanpiRefExc_36_1.f + cAtanpiRefExc_36_2.f; Exit; end;
    if ax = cAtanpiRefExc_37_0.f then begin Result := cAtanpiRefExc_37_1.f + cAtanpiRefExc_37_2.f; Exit; end;
    if ax = cAtanpiRefExc_38_0.f then begin Result := cAtanpiRefExc_38_1.f + cAtanpiRefExc_38_2.f; Exit; end;
    if ax = cAtanpiRefExc_39_0.f then begin Result := cAtanpiRefExc_39_1.f + cAtanpiRefExc_39_2.f; Exit; end;
    if ax = cAtanpiRefExc_40_0.f then begin Result := cAtanpiRefExc_40_1.f + cAtanpiRefExc_40_2.f; Exit; end;
    if ax = cAtanpiRefExc_41_0.f then begin Result := cAtanpiRefExc_41_1.f + cAtanpiRefExc_41_2.f; Exit; end;
    if ax = cAtanpiRefExc_42_0.f then begin Result := cAtanpiRefExc_42_1.f + cAtanpiRefExc_42_2.f; Exit; end;
    if ax = cAtanpiRefExc_43_0.f then begin Result := cAtanpiRefExc_43_1.f + cAtanpiRefExc_43_2.f; Exit; end;
    if ax = cAtanpiRefExc_44_0.f then begin Result := cAtanpiRefExc_44_1.f + cAtanpiRefExc_44_2.f; Exit; end;
    if ax = cAtanpiRefExc_45_0.f then begin Result := cAtanpiRefExc_45_1.f + cAtanpiRefExc_45_2.f; Exit; end;
    if ax = cAtanpiRefExc_46_0.f then begin Result := cAtanpiRefExc_46_1.f + cAtanpiRefExc_46_2.f; Exit; end;
    if ax = cAtanpiRefExc_47_0.f then begin Result := cAtanpiRefExc_47_1.f + cAtanpiRefExc_47_2.f; Exit; end;
    if ax = cAtanpiRefExc_48_0.f then begin Result := cAtanpiRefExc_48_1.f + cAtanpiRefExc_48_2.f; Exit; end;
    if ax = cAtanpiRefExc_49_0.f then begin Result := cAtanpiRefExc_49_1.f + cAtanpiRefExc_49_2.f; Exit; end;
    if ax = cAtanpiRefExc_50_0.f then begin Result := cAtanpiRefExc_50_1.f + cAtanpiRefExc_50_2.f; Exit; end;
    if ax = cAtanpiRefExc_51_0.f then begin Result := cAtanpiRefExc_51_1.f + cAtanpiRefExc_51_2.f; Exit; end;
    if ax = cAtanpiRefExc_52_0.f then begin Result := cAtanpiRefExc_52_1.f + cAtanpiRefExc_52_2.f; Exit; end;
    if ax = cAtanpiRefExc_53_0.f then begin Result := cAtanpiRefExc_53_1.f + cAtanpiRefExc_53_2.f; Exit; end;
    if ax = cAtanpiRefExc_54_0.f then begin Result := cAtanpiRefExc_54_1.f + cAtanpiRefExc_54_2.f; Exit; end;
    if ax = cAtanpiRefExc_55_0.f then begin Result := cAtanpiRefExc_55_1.f + cAtanpiRefExc_55_2.f; Exit; end;
  end
  else
  begin
    if ax = cAtanpiRefExc_0_0.f then begin Result := -cAtanpiRefExc_0_1.f + -cAtanpiRefExc_0_2.f; Exit; end;
    if ax = cAtanpiRefExc_1_0.f then begin Result := -cAtanpiRefExc_1_1.f + -cAtanpiRefExc_1_2.f; Exit; end;
    if ax = cAtanpiRefExc_2_0.f then begin Result := -cAtanpiRefExc_2_1.f + -cAtanpiRefExc_2_2.f; Exit; end;
    if ax = cAtanpiRefExc_3_0.f then begin Result := -cAtanpiRefExc_3_1.f + -cAtanpiRefExc_3_2.f; Exit; end;
    if ax = cAtanpiRefExc_4_0.f then begin Result := -cAtanpiRefExc_4_1.f + -cAtanpiRefExc_4_2.f; Exit; end;
    if ax = cAtanpiRefExc_5_0.f then begin Result := -cAtanpiRefExc_5_1.f + -cAtanpiRefExc_5_2.f; Exit; end;
    if ax = cAtanpiRefExc_6_0.f then begin Result := -cAtanpiRefExc_6_1.f + -cAtanpiRefExc_6_2.f; Exit; end;
    if ax = cAtanpiRefExc_7_0.f then begin Result := -cAtanpiRefExc_7_1.f + -cAtanpiRefExc_7_2.f; Exit; end;
    if ax = cAtanpiRefExc_8_0.f then begin Result := -cAtanpiRefExc_8_1.f + -cAtanpiRefExc_8_2.f; Exit; end;
    if ax = cAtanpiRefExc_9_0.f then begin Result := -cAtanpiRefExc_9_1.f + -cAtanpiRefExc_9_2.f; Exit; end;
    if ax = cAtanpiRefExc_10_0.f then begin Result := -cAtanpiRefExc_10_1.f + -cAtanpiRefExc_10_2.f; Exit; end;
    if ax = cAtanpiRefExc_11_0.f then begin Result := -cAtanpiRefExc_11_1.f + -cAtanpiRefExc_11_2.f; Exit; end;
    if ax = cAtanpiRefExc_12_0.f then begin Result := -cAtanpiRefExc_12_1.f + -cAtanpiRefExc_12_2.f; Exit; end;
    if ax = cAtanpiRefExc_13_0.f then begin Result := -cAtanpiRefExc_13_1.f + -cAtanpiRefExc_13_2.f; Exit; end;
    if ax = cAtanpiRefExc_14_0.f then begin Result := -cAtanpiRefExc_14_1.f + -cAtanpiRefExc_14_2.f; Exit; end;
    if ax = cAtanpiRefExc_15_0.f then begin Result := -cAtanpiRefExc_15_1.f + -cAtanpiRefExc_15_2.f; Exit; end;
    if ax = cAtanpiRefExc_16_0.f then begin Result := -cAtanpiRefExc_16_1.f + -cAtanpiRefExc_16_2.f; Exit; end;
    if ax = cAtanpiRefExc_17_0.f then begin Result := -cAtanpiRefExc_17_1.f + -cAtanpiRefExc_17_2.f; Exit; end;
    if ax = cAtanpiRefExc_18_0.f then begin Result := -cAtanpiRefExc_18_1.f + -cAtanpiRefExc_18_2.f; Exit; end;
    if ax = cAtanpiRefExc_19_0.f then begin Result := -cAtanpiRefExc_19_1.f + -cAtanpiRefExc_19_2.f; Exit; end;
    if ax = cAtanpiRefExc_20_0.f then begin Result := -cAtanpiRefExc_20_1.f + -cAtanpiRefExc_20_2.f; Exit; end;
    if ax = cAtanpiRefExc_21_0.f then begin Result := -cAtanpiRefExc_21_1.f + -cAtanpiRefExc_21_2.f; Exit; end;
    if ax = cAtanpiRefExc_22_0.f then begin Result := -cAtanpiRefExc_22_1.f + -cAtanpiRefExc_22_2.f; Exit; end;
    if ax = cAtanpiRefExc_23_0.f then begin Result := -cAtanpiRefExc_23_1.f + -cAtanpiRefExc_23_2.f; Exit; end;
    if ax = cAtanpiRefExc_24_0.f then begin Result := -cAtanpiRefExc_24_1.f + -cAtanpiRefExc_24_2.f; Exit; end;
    if ax = cAtanpiRefExc_25_0.f then begin Result := -cAtanpiRefExc_25_1.f + -cAtanpiRefExc_25_2.f; Exit; end;
    if ax = cAtanpiRefExc_26_0.f then begin Result := -cAtanpiRefExc_26_1.f + -cAtanpiRefExc_26_2.f; Exit; end;
    if ax = cAtanpiRefExc_27_0.f then begin Result := -cAtanpiRefExc_27_1.f + -cAtanpiRefExc_27_2.f; Exit; end;
    if ax = cAtanpiRefExc_28_0.f then begin Result := -cAtanpiRefExc_28_1.f + -cAtanpiRefExc_28_2.f; Exit; end;
    if ax = cAtanpiRefExc_29_0.f then begin Result := -cAtanpiRefExc_29_1.f + -cAtanpiRefExc_29_2.f; Exit; end;
    if ax = cAtanpiRefExc_30_0.f then begin Result := -cAtanpiRefExc_30_1.f + -cAtanpiRefExc_30_2.f; Exit; end;
    if ax = cAtanpiRefExc_31_0.f then begin Result := -cAtanpiRefExc_31_1.f + -cAtanpiRefExc_31_2.f; Exit; end;
    if ax = cAtanpiRefExc_32_0.f then begin Result := -cAtanpiRefExc_32_1.f + -cAtanpiRefExc_32_2.f; Exit; end;
    if ax = cAtanpiRefExc_33_0.f then begin Result := -cAtanpiRefExc_33_1.f + -cAtanpiRefExc_33_2.f; Exit; end;
    if ax = cAtanpiRefExc_34_0.f then begin Result := -cAtanpiRefExc_34_1.f + -cAtanpiRefExc_34_2.f; Exit; end;
    if ax = cAtanpiRefExc_35_0.f then begin Result := -cAtanpiRefExc_35_1.f + -cAtanpiRefExc_35_2.f; Exit; end;
    if ax = cAtanpiRefExc_36_0.f then begin Result := -cAtanpiRefExc_36_1.f + -cAtanpiRefExc_36_2.f; Exit; end;
    if ax = cAtanpiRefExc_37_0.f then begin Result := -cAtanpiRefExc_37_1.f + -cAtanpiRefExc_37_2.f; Exit; end;
    if ax = cAtanpiRefExc_38_0.f then begin Result := -cAtanpiRefExc_38_1.f + -cAtanpiRefExc_38_2.f; Exit; end;
    if ax = cAtanpiRefExc_39_0.f then begin Result := -cAtanpiRefExc_39_1.f + -cAtanpiRefExc_39_2.f; Exit; end;
    if ax = cAtanpiRefExc_40_0.f then begin Result := -cAtanpiRefExc_40_1.f + -cAtanpiRefExc_40_2.f; Exit; end;
    if ax = cAtanpiRefExc_41_0.f then begin Result := -cAtanpiRefExc_41_1.f + -cAtanpiRefExc_41_2.f; Exit; end;
    if ax = cAtanpiRefExc_42_0.f then begin Result := -cAtanpiRefExc_42_1.f + -cAtanpiRefExc_42_2.f; Exit; end;
    if ax = cAtanpiRefExc_43_0.f then begin Result := -cAtanpiRefExc_43_1.f + -cAtanpiRefExc_43_2.f; Exit; end;
    if ax = cAtanpiRefExc_44_0.f then begin Result := -cAtanpiRefExc_44_1.f + -cAtanpiRefExc_44_2.f; Exit; end;
    if ax = cAtanpiRefExc_45_0.f then begin Result := -cAtanpiRefExc_45_1.f + -cAtanpiRefExc_45_2.f; Exit; end;
    if ax = cAtanpiRefExc_46_0.f then begin Result := -cAtanpiRefExc_46_1.f + -cAtanpiRefExc_46_2.f; Exit; end;
    if ax = cAtanpiRefExc_47_0.f then begin Result := -cAtanpiRefExc_47_1.f + -cAtanpiRefExc_47_2.f; Exit; end;
    if ax = cAtanpiRefExc_48_0.f then begin Result := -cAtanpiRefExc_48_1.f + -cAtanpiRefExc_48_2.f; Exit; end;
    if ax = cAtanpiRefExc_49_0.f then begin Result := -cAtanpiRefExc_49_1.f + -cAtanpiRefExc_49_2.f; Exit; end;
    if ax = cAtanpiRefExc_50_0.f then begin Result := -cAtanpiRefExc_50_1.f + -cAtanpiRefExc_50_2.f; Exit; end;
    if ax = cAtanpiRefExc_51_0.f then begin Result := -cAtanpiRefExc_51_1.f + -cAtanpiRefExc_51_2.f; Exit; end;
    if ax = cAtanpiRefExc_52_0.f then begin Result := -cAtanpiRefExc_52_1.f + -cAtanpiRefExc_52_2.f; Exit; end;
    if ax = cAtanpiRefExc_53_0.f then begin Result := -cAtanpiRefExc_53_1.f + -cAtanpiRefExc_53_2.f; Exit; end;
    if ax = cAtanpiRefExc_54_0.f then begin Result := -cAtanpiRefExc_54_1.f + -cAtanpiRefExc_54_2.f; Exit; end;
    if ax = cAtanpiRefExc_55_0.f then begin Result := -cAtanpiRefExc_55_1.f + -cAtanpiRefExc_55_2.f; Exit; end;
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
