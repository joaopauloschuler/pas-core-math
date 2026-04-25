// ===========================================================================
// pcr_acosh — correctly-rounded binary64 acosh(x).
// Ported from core-math/src/binary64/acosh/acosh.c by Alexei Sibidanov.
// https://github.com/joaopauloschuler/pas-core-math
//
// Copyright (c) 2024-2026 Joao Paulo Schwarz Schuler and contributors.
// Refer to the git commit history for individual authorship.
// SPDX-License-Identifier: MIT
// ===========================================================================

type
  TAcoshBEntry = record c0: UInt16; c1: Int16; end;

const
  cAcoshB: array[0..31] of TAcoshBEntry = (
    (c0:   301; c1:  27565), (c0:  7189; c1:  24786), (c0: 13383; c1:  22167),
    (c0: 18923; c1:  19696), (c0: 23845; c1:  17361), (c0: 28184; c1:  15150),
    (c0: 31969; c1:  13054), (c0: 35231; c1:  11064), (c0: 37996; c1:   9173),
    (c0: 40288; c1:   7372), (c0: 42129; c1:   5657), (c0: 43542; c1:   4020),
    (c0: 44546; c1:   2457), (c0: 45160; c1:    962), (c0: 45399; c1:   -468),
    (c0: 45281; c1:  -1838), (c0: 44821; c1:  -3151), (c0: 44032; c1:  -4412),
    (c0: 42929; c1:  -5622), (c0: 41522; c1:  -6786), (c0: 39825; c1:  -7905),
    (c0: 37848; c1:  -8982), (c0: 35602; c1: -10020), (c0: 33097; c1: -11020),
    (c0: 30341; c1: -11985), (c0: 27345; c1: -12916), (c0: 24115; c1: -13816),
    (c0: 20661; c1: -14685), (c0: 16989; c1: -15526), (c0: 13107; c1: -16339),
    (c0:  9022; c1: -17126), (c0:  4740; c1: -17889));

  cAcoshR1: array[0..32] of Tb64u64 = (
    (u:$3FF0000000000000), (u:$3FEF507600000000), (u:$3FEEA4B000000000),
    (u:$3FEDFC9800000000), (u:$3FED581800000000), (u:$3FECB72000000000),
    (u:$3FEC199C00000000), (u:$3FEB7F7600000000), (u:$3FEAE8A000000000),
    (u:$3FEA550400000000), (u:$3FE9C49200000000), (u:$3FE9373800000000),
    (u:$3FE8ACE600000000), (u:$3FE8258A00000000), (u:$3FE7A11400000000),
    (u:$3FE71F7600000000), (u:$3FE6A09E00000000), (u:$3FE6247E00000000),
    (u:$3FE5AB0800000000), (u:$3FE5342C00000000), (u:$3FE4BFDA00000000),
    (u:$3FE44E0800000000), (u:$3FE3DEA600000000), (u:$3FE371A800000000),
    (u:$3FE306FE00000000), (u:$3FE29E9E00000000), (u:$3FE2387A00000000),
    (u:$3FE1D48800000000), (u:$3FE172B800000000), (u:$3FE1130200000000),
    (u:$3FE0B55800000000), (u:$3FE059B000000000), (u:$3FE0000000000000));

  cAcoshR2: array[0..32] of Tb64u64 = (
    (u:$3FF0000000000000), (u:$3FEFFA7400000000), (u:$3FEFF4EA00000000),
    (u:$3FEFEF6200000000), (u:$3FEFE9DA00000000), (u:$3FEFE45200000000),
    (u:$3FEFDECC00000000), (u:$3FEFD94600000000), (u:$3FEFD3C200000000),
    (u:$3FEFCE3E00000000), (u:$3FEFC8BC00000000), (u:$3FEFC33A00000000),
    (u:$3FEFBDBA00000000), (u:$3FEFB83A00000000), (u:$3FEFB2BC00000000),
    (u:$3FEFAD3E00000000), (u:$3FEFA7C200000000), (u:$3FEFA24600000000),
    (u:$3FEF9CCA00000000), (u:$3FEF975000000000), (u:$3FEF91D800000000),
    (u:$3FEF8C6000000000), (u:$3FEF86E800000000), (u:$3FEF817200000000),
    (u:$3FEF7BFE00000000), (u:$3FEF768A00000000), (u:$3FEF711600000000),
    (u:$3FEF6BA400000000), (u:$3FEF663200000000), (u:$3FEF60C200000000),
    (u:$3FEF5B5200000000), (u:$3FEF55E400000000), (u:$3FEF507600000000));

  // l1[i][0] = low, l1[i][1] = high (C indexing: [0]=lo, [1]=hi)
  cAcoshL1: array[0..32, 0..1] of Tb64u64 = (
    ((u:$0000000000000000),(u:$0000000000000000)),
    ((u:QWord($BD1269E2038315B3)),(u:$3F962E4EACD40000)),
    ((u:QWord($BD23F2558BDDFC47)),(u:$3FA62E3CE7218000)),
    ((u:$3D207EA13C34EFB5),(u:$3FB0A2AB6D3EC000)),
    ((u:$3D38F3E77084D3BA),(u:$3FB62E4A86D8C000)),
    ((u:QWord($BD18D92A005F1A7E)),(u:$3FBBB9DB7062C000)),
    ((u:$3D358239E799BFE5),(u:$3FC0A2B1A22CC000)),
    ((u:QWord($BD3A93FCF5F593B7)),(u:$3FC3687F0A298000)),
    ((u:QWord($BD1DB4CAC32FD2B5)),(u:$3FC62E4116B64000)),
    ((u:QWord($BD10E65A92EE0F3B)),(u:$3FC8F409E4DF6000)),
    ((u:QWord($BD38261383D475F1)),(u:$3FCBB9D15001C000)),
    ((u:QWord($BD3359886207513B)),(u:$3FCE7F9A8C940000)),
    ((u:$3D3811F87496CEB7),(u:$3FD0A2B052DDB000)),
    ((u:$3D34991EC6CB435C),(u:$3FD205955EF73000)),
    ((u:QWord($BD34581ABFEB8927)),(u:$3FD3687BD9121000)),
    ((u:$3D3CAB48F6942703),(u:$3FD4CB5E8F2B5000)),
    ((u:QWord($BD0DF2C452FDE132)),(u:$3FD62E4420E20000)),
    ((u:$3D26109F4FDB74BD),(u:$3FD791292C46A000)),
    ((u:QWord($BD36B95FBDAC7696)),(u:$3FD8F40AF84E7000)),
    ((u:$3D17394FA880CBDA),(u:$3FDA56ED8F865000)),
    ((u:QWord($BD150B06A94ECCAB)),(u:$3FDBB9D6505B4000)),
    ((u:QWord($BD3BE2ABF0B38989)),(u:$3FDD1CB91E728000)),
    ((u:QWord($BD37D6BF1E34DA04)),(u:$3FDE7F9D139E2000)),
    ((u:QWord($BD3423C1E14DE6ED)),(u:$3FDFE27DB9B0E000)),
    ((u:$3D3C46F1A0EFBBC2),(u:$3FE0A2B25060A800)),
    ((u:$3D2834FE4E3E6018),(u:$3FE154244482A000)),
    ((u:$3D16A03D0F02B650),(u:$3FE2059731298800)),
    ((u:$3D3D437056526F30),(u:$3FE2B707145DE000)),
    ((u:QWord($BD2A0233728405C5)),(u:$3FE3687B0E0B2800)),
    ((u:QWord($BD24DBDDA10D2BF1)),(u:$3FE419EC5D3F6800)),
    ((u:$3D3F7D0A25D154F2),(u:$3FE4CB5F9FC02000)),
    ((u:$3D315EDE4D803B18),(u:$3FE57CD28421A800)),
    ((u:$3D2EF35793C76730),(u:$3FE62E42FEFA3800)));

  cAcoshL2: array[0..32, 0..1] of Tb64u64 = (
    ((u:$0000000000000000),(u:$0000000000000000)),
    ((u:$3D35ABDAC3638E99),(u:$3F4631EC81E00000)),
    ((u:QWord($BD216B8BE9BBE239)),(u:$3F562FD812700000)),
    ((u:QWord($BD3364C6315542EB)),(u:$3F60A25205080000)),
    ((u:$3D2734ABE459C900),(u:$3F662DADC1D00000)),
    ((u:$3D30CF8A761431BF),(u:$3F6BB9FF94D00000)),
    ((u:$3D2DA2718EB78708),(u:$3F70A2A2DEF80000)),
    ((u:$3D334ADA62C59B93),(u:$3F7368C0FAE40000)),
    ((u:$3D3D09AB376682D4),(u:$3F762E58E4F80000)),
    ((u:QWord($BD23CB7B94329211)),(u:$3F78F46BD28C0000)),
    ((u:QWord($BD2EEC5C297C41D0)),(u:$3F7BB9F831200000)),
    ((u:QWord($BD36411B9395D150)),(u:$3F7E7FFF8F300000)),
    ((u:QWord($BD31C0E59A43053C)),(u:$3F80A2C0006E0000)),
    ((u:$3D16506596E077B6),(u:$3F8205BDB6F00000)),
    ((u:$3D3E256BCE6FAA27),(u:$3F836877C86E0000)),
    ((u:$3CCBD42467B0C8D1),(u:$3F84CB6F55780000)),
    ((u:QWord($BD3C4F92132FF0F0)),(u:$3F862E230E8C0000)),
    ((u:QWord($BD380BE08BFAB390)),(u:$3F87911440F60000)),
    ((u:QWord($BD3F0B1319CEB1F7)),(u:$3F88F443020A0000)),
    ((u:$3D2A65FCFB8DE99B),(u:$3F8A572DBEF40000)),
    ((u:$3D14233885D3779C),(u:$3F8BB9D449A60000)),
    ((u:$3D3F46A59E646EDB),(u:$3F8D1CB8491C0000)),
    ((u:QWord($BD3C3D2F11C11446)),(u:$3F8E7FD9D2AA0000)),
    ((u:$3D27763F78A1E0CC),(u:$3F8FE2B6F9780000)),
    ((u:$3D3B4C37FC60C043),(u:$3F90A2A7C7A50000)),
    ((u:QWord($BD15B8A822859BE3)),(u:$3F915412CA860000)),
    ((u:QWord($BD3F2D8C9FC06400)),(u:$3F92059C90050000)),
    ((u:QWord($BD3E80E79C20378D)),(u:$3F92B703F49B0000)),
    ((u:$3D368256E4329BDB),(u:$3F93688A1A8D0000)),
    ((u:$3D37E9741DA248C3),(u:$3F9419EDC7BA0000)),
    ((u:$3D2E330DCCCE602B),(u:$3F94CB7034FA0000)),
    ((u:$3CE2F32B5D18EEFB),(u:$3F957CD011870000)),
    ((u:QWord($BD1269E2038315B3)),(u:$3F962E4EACD40000)));

  cAcoshC: array[0..4] of Tb64u64 = (
    (u:QWord($BFE0000000000000)), (u:$3FD5555555555530), (u:QWord($BFCFFFFFFFFFFFA0)),
    (u:$3FC99999E33A6366), (u:QWord($BFC555559EF9525F)));

  // log(2) split (matches cLog values elsewhere)
  cAcoshL2h: Tb64u64 = (u:$3FE62E42FEFA3800); // 0x1.62e42fefa38p-1
  cAcoshL2l: Tb64u64 = (u:$3D2EF35793C76730); // 0x1.ef35793c7673p-45

  // Fast-path (|x-1| small) inline cl[9]
  cAcoshOneFastCl: array[0..8] of Tb64u64 = (
    (u:QWord($BFB5555555555555)), (u:$3F93333333332F95), (u:QWord($BF76DB6DB6D5534C)),
    (u:$3F5F1C71C1E04356), (u:QWord($BF46E8B8E3E40D58)), (u:$3F31C4BA825AC4FE),
    (u:QWord($BF1C9045534E6D9E)), (u:$3F071FEDAE26A76B), (u:QWord($BEEF1F4F8CC65342)));

  // Medium-range polynomials for g
  cAcoshMed4: array[0..3] of Tb64u64 = (
    (u:$3BD5C4B6148816E2), (u:QWord($BFD000000000005C)), (u:QWord($BFB7FFFFFEBF3E6C)),
    (u:QWord($BFAAAB6691F2BAE7)));
  cAcoshMed3: array[0..2] of Tb64u64 = (
    (u:QWord($BBC7F77C8429C6C6)), (u:QWord($BFCFFFFFFFFFF214)), (u:QWord($BFB8000268641BFE)));
  cAcoshMed2: array[0..1] of Tb64u64 = (
    (u:$3BC7A0ED2EFFDD10), (u:QWord($BFD000000017D048)));

  // as_acosh_one ch[10][2], cl[6]
  cAcoshOneCh: array[0..9, 0..1] of Tb64u64 = (
    ((u:QWord($BFB5555555555555)),(u:QWord($BC55555555554AF1))),
    ((u:$3F93333333333333),(u:$3C29999998933F0E)),
    ((u:QWord($BF76DB6DB6DB6DB7)),(u:$3C024929B16EC6B7)),
    ((u:$3F5F1C71C71C71C7),(u:$3BDC56D45E265E2C)),
    ((u:QWord($BF46E8BA2E8BA2E9)),(u:$3BE6D50CE7188D3D)),
    ((u:$3F31C4EC4EC4EC43),(u:$3BDC6791D1CF399A)),
    ((u:QWord($BF1C99999999914F)),(u:$3BBEE0D9408A2E2A)),
    ((u:$3F07A878787648E2),(u:QWord($BBA1CEA281E08012))),
    ((u:QWord($BEF3FDE50D0CB4B9)),(u:$3B70335101403D9D)),
    ((u:$3EE12EF3BF8A0A74),(u:$3AFF9C6B51787043)));
  cAcoshOneCl: array[0..5] of Tb64u64 = (
    (u:QWord($BECDF3B9D1296EA9)), (u:$3EBA681D7D2298EB), (u:QWord($BEA77EAD7B1CA449)),
    (u:$3E94EDD2DDB3721F), (u:QWord($BE81BF173531EE23)), (u:$3E6613229230E255));

  // Refine path: t1 (17), t2/t3/t4 (16 each)
  cAcoshT1: array[0..16] of Tb64u64 = (
    (u:$3FF0000000000000), (u:$3FEEA4AFA0000000), (u:$3FED5818E0000000),
    (u:$3FEC199BE0000000), (u:$3FEAE89F98000000), (u:$3FE9C49180000000),
    (u:$3FE8ACE540000000), (u:$3FE7A11470000000), (u:$3FE6A09E68000000),
    (u:$3FE5AB07E0000000), (u:$3FE4BFDAD8000000), (u:$3FE3DEA650000000),
    (u:$3FE306FE08000000), (u:$3FE2387A70000000), (u:$3FE172B840000000),
    (u:$3FE0B55870000000), (u:$3FE0000000000000));
  cAcoshT2: array[0..15] of Tb64u64 = (
    (u:$3FF0000000000000), (u:$3FEFE9D968000000), (u:$3FEFD3C228000000),
    (u:$3FEFBDBA38000000), (u:$3FEFA7C180000000), (u:$3FEF91D800000000),
    (u:$3FEF7BFDB0000000), (u:$3FEF663278000000), (u:$3FEF507658000000),
    (u:$3FEF3AC948000000), (u:$3FEF252B38000000), (u:$3FEF0F9C20000000),
    (u:$3FEEFA1BF0000000), (u:$3FEEE4AAA0000000), (u:$3FEECF4830000000),
    (u:$3FEEB9F488000000));
  cAcoshT3: array[0..15] of Tb64u64 = (
    (u:$3FF0000000000000), (u:$3FEFFE9D20000000), (u:$3FEFFD3A58000000),
    (u:$3FEFFBD798000000), (u:$3FEFFA74E8000000), (u:$3FEFF91248000000),
    (u:$3FEFF7AFB8000000), (u:$3FEFF64D38000000), (u:$3FEFF4EAC8000000),
    (u:$3FEFF38868000000), (u:$3FEFF22618000000), (u:$3FEFF0C3D0000000),
    (u:$3FEFEF61A0000000), (u:$3FEFEDFF78000000), (u:$3FEFEC9D68000000),
    (u:$3FEFEB3B60000000));
  cAcoshT4: array[0..15] of Tb64u64 = (
    (u:$3FF0000000000000), (u:$3FEFFFE9D0000000), (u:$3FEFFFD3A0000000),
    (u:$3FEFFFBD78000000), (u:$3FEFFFA748000000), (u:$3FEFFF9118000000),
    (u:$3FEFFF7AE8000000), (u:$3FEFFF64C0000000), (u:$3FEFFF4E90000000),
    (u:$3FEFFF3860000000), (u:$3FEFFF2238000000), (u:$3FEFFF0C08000000),
    (u:$3FEFFEF5D8000000), (u:$3FEFFEDFA8000000), (u:$3FEFFEC980000000),
    (u:$3FEFFEB350000000));

  // LL[4][17][3] triples
  cAcoshLL: array[0..3, 0..16, 0..2] of Tb64u64 = (
    ( ((u:$0000000000000000),(u:$0000000000000000),(u:$0000000000000000)),
      ((u:$3F962E432B240000),(u:QWord($BD5745AF34BB54B8)),(u:QWord($B9E17E3EC05CDE70))),
      ((u:$3FA62E42E4A80000),(u:$3D3111A4EADF3120),(u:$3A2CFF3027ABB119)),
      ((u:$3FB0A2B233F10000),(u:QWord($BD588AC4EC78AF80)),(u:$3A24FA087CA75DFD)),
      ((u:$3FB62E43056C0000),(u:$3D16BD65E8B0B700),(u:QWord($BA0B18E160362C24))),
      ((u:$3FBBB9D3CBD60000),(u:$3D5DE14AA55EC2B0),(u:QWord($BA1C6AC3F1862A6B))),
      ((u:$3FC0A2B244DA0000),(u:$3D594DEF487FEA70),(u:QWord($BA1DEAD1A4581ACF))),
      ((u:$3FC3687AA9B78000),(u:$3D49CEC9A50DB220),(u:$3A234A70684F8E0E)),
      ((u:$3FC62E42FABA0000),(u:QWord($BD3D69047A3AEB00)),(u:QWord($BA04E061F79144E2))),
      ((u:$3FC8F40B56D28000),(u:$3D5DE7D755FD2E20),(u:$3A1BDC7ECF001489)),
      ((u:$3FCBB9D3B61F0000),(u:$3D1C14F1445B1200),(u:$3A2A1D78CBDC5B58)),
      ((u:$3FCE7F9C11F08000),(u:QWord($BD46E3E0000DAE70)),(u:$3A16A4559FADDE98)),
      ((u:$3FD0A2B242EC4000),(u:$3D5BB7CF852A5FE8),(u:$3A2A6AEF11EE43BD)),
      ((u:$3FD205966C764000),(u:$3D2AD3A5F2142940),(u:$3A25CC344FA10652)),
      ((u:$3FD3687A98AAC000),(u:$3D21623671842F00),(u:QWord($BA10B428FE1F9E43))),
      ((u:$3FD4CB5EC93F4000),(u:$3D53D50980EA5130),(u:$3A267F0EA083B1C4)),
      ((u:$3FD62E42FEFA4000),(u:QWord($BD38432A1B0E2640)),(u:$3A2803F2F6AF40F3)) ),
    ( ((u:$0000000000000000),(u:$0000000000000000),(u:$0000000000000000)),
      ((u:$3F562E462B400000),(u:$3D5061D003B97318),(u:$3A2D7FAEE66A2E1E)),
      ((u:$3F662E44C9200000),(u:$3D595A7BFF5E2390),(u:QWord($BA0F7E788A871350))),
      ((u:$3F70A2B1E3300000),(u:$3D42A3A1A65AA3A0),(u:QWord($BA254599C9605442))),
      ((u:$3F762E4367C00000),(u:QWord($BD24A995B6D9DDC0)),(u:QWord($B9B56BB79B254F33))),
      ((u:$3F7BB9D449A00000),(u:$3D58A119C42E9BC0),(u:QWord($BA28ECF7D8D661F1))),
      ((u:$3F80A2B1F1900000),(u:$3D58863771BD10A8),(u:$3A1E9731DE7F0155)),
      ((u:$3F83687AD1100000),(u:$3D5E026A347CA1C8),(u:$39EFADC62522444D)),
      ((u:$3F862E436F280000),(u:$3D525B84F71B70B8),(u:QWord($B9FFCB3F98612D27))),
      ((u:$3F88F40B7B380000),(u:QWord($BD462A0A4FD47580)),(u:$3A23CB3C35D9F6A1)),
      ((u:$3F8BB9D3ABB00000),(u:QWord($BD50EC48F94D7860)),(u:QWord($BA26B47D410E4CC7))),
      ((u:$3F8E7F9BB2300000),(u:$3D4E4415CBC97A00),(u:QWord($BA23729FDB677231))),
      ((u:$3F90A2B224780000),(u:QWord($BD5CB73F4505B030)),(u:QWord($BA21B3B3A3BC370A))),
      ((u:$3F92059691E80000),(u:QWord($BD4ABCC3412F2640)),(u:QWord($BA0FE6E998E48673))),
      ((u:$3F93687A76800000),(u:QWord($BD543901E5C97A90)),(u:$39FB54CDD52A5D88)),
      ((u:$3F94CB5EB5D80000),(u:QWord($BD58F106F00F13B8)),(u:QWord($BA28F793F5FCE148))),
      ((u:$3F962E432B240000),(u:QWord($BD5745AF34BB54B8)),(u:QWord($B9E17E3EC05CDE70))) ),
    ( ((u:$0000000000000000),(u:$0000000000000000),(u:$0000000000000000)),
      ((u:$3F162E7B00000000),(u:QWord($BD3868625640A680)),(u:QWord($BA234BF0DB910F65))),
      ((u:$3F262E35F6000000),(u:QWord($BD42EE3D96B696A0)),(u:$3A1A2948CD558655)),
      ((u:$3F30A2B4B2000000),(u:$3D053EDBCF116500),(u:QWord($B9ECFC26CCF6D0E4))),
      ((u:$3F362E4BE1000000),(u:$3CB783E334614000),(u:QWord($BA204B96DA30E63A))),
      ((u:$3F3BB9E085000000),(u:QWord($BD460785F20ACB20)),(u:QWord($B9FF33369BF7DFF1))),
      ((u:$3F40A2B94D000000),(u:$3D5FD4B3A2733530),(u:QWord($B9F685A35575EFF1))),
      ((u:$3F4368810F800000),(u:$3D07DED26DC81300),(u:QWord($B9F4C4D1ABCA79BF))),
      ((u:$3F462E4787800000),(u:$3D57D2BEE9A1F630),(u:$3A2860233B7AD130)),
      ((u:$3F48F40CB4800000),(u:QWord($BD5AF034EAF471C0)),(u:$3A1AE748822D57B7)),
      ((u:$3F4BB9D094000000),(u:QWord($BD57A223013A20F0)),(u:QWord($BA21E499087075B6))),
      ((u:$3F4E7FA32C800000),(u:QWord($BD4B2E67B1B59BD0)),(u:QWord($BA254A41EDA30FA6))),
      ((u:$3F50A2B237000000),(u:QWord($BD37AD97FF4AC7A0)),(u:$3A2F932DA91371DD)),
      ((u:$3F52059A33800000),(u:QWord($BD396422D90DF400)),(u:QWord($BA190800FBBF2ED3))),
      ((u:$3F53687982400000),(u:$3D30F90540018120),(u:$3A29567E01E48F9A)),
      ((u:$3F54CB602C000000),(u:QWord($BD40D709A5EC0B50)),(u:$3A1253DFD44635D2)),
      ((u:$3F562E462B400000),(u:$3D5061D003B97318),(u:$3A2D7FAEE66A2E1E)) ),
    ( ((u:$0000000000000000),(u:$0000000000000000),(u:$0000000000000000)),
      ((u:$3ED63007C0000000),(u:QWord($BD4DB0E38E5AAAA0)),(u:$3A2259A7B94815B9)),
      ((u:$3EE6300F60000000),(u:$3D32B1C755804380),(u:$3A278CABBA01E3E4)),
      ((u:$3EF0A21150000000),(u:QWord($BD55FF2237307590)),(u:$3A08074FEACFE49D)),
      ((u:$3EF62E1EC0000000),(u:QWord($BD285D6F6487CE40)),(u:$3A205485074B9276)),
      ((u:$3EFBBA3010000000),(u:QWord($BD4AF5D58A7C9210)),(u:QWord($BA230A8C0FD2FF5F))),
      ((u:$3F00A32298000000),(u:$3D4590FAA0883BD0),(u:$3A295E9BDA999947)),
      ((u:$3F03682F10000000),(u:$3D5F0224376EFAF8),(u:QWord($BA25843C0DB50D10))),
      ((u:$3F062E3D80000000),(u:QWord($BD4142C13DAED4A0)),(u:$3A2C68A61183CE87)),
      ((u:$3F08F44DD8000000),(u:QWord($BD4AA489F3999310)),(u:$3A111C5C376854EA)),
      ((u:$3F0BB96010000000),(u:$3D59904D8B6A3638),(u:$3A28C89554493C8F)),
      ((u:$3F0E7F7440000000),(u:$3D55785DDBE7CBA8),(u:$3A1E7FF3CDE7D70C)),
      ((u:$3F10A2C530000000),(u:QWord($BD46D9E8780D0D50)),(u:$3A1AD9C178106693)),
      ((u:$3F1205D134000000),(u:QWord($BD4214A2E893FCC0)),(u:$3A2548A9500C9822)),
      ((u:$3F13685E28000000),(u:$3D4E235886461030),(u:$3A12A97B26DA2D88)),
      ((u:$3F14CB6C18000000),(u:$3D52B7CFCEA9E0D8),(u:QWord($BA25095048A6B824))),
      ((u:$3F162E7B00000000),(u:QWord($BD3868625640A680)),(u:QWord($BA234BF0DB910F65))) ));

  cAcoshRefineCh: array[0..2, 0..1] of Tb64u64 = (
    ((u:$3FE0000000000000),(u:$39024B67EE516E3B)),
    ((u:QWord($BFD0000000000000)),(u:QWord($B91932CE43199A8D))),
    ((u:$3FC5555555555555),(u:$3C655540C15CF91F)));
  cAcoshRefineCl: array[0..2] of Tb64u64 = (
    (u:QWord($BFC0000000000000)), (u:$3FB9999999A0754F), (u:QWord($BFB55555555C3157)));

  // log(2) triple split for refine path
  cAcoshL20: Tb64u64 = (u:$3FD62E42FEFA3800); // 0x1.62e42fefa38p-2
  cAcoshL21: Tb64u64 = (u:$3D1EF35793C76800); // 0x1.ef35793c768p-46
  cAcoshL22: Tb64u64 = (u:QWord($BA49FF0342542FC3)); // -0x1.9ff0342542fc3p-91

  cAcoshDb: array[0..9, 0..2] of Tb64u64 = (
    ((u:$3FF5BFF041B260FE),(u:$3FEA6031CD5F93BA),(u:$3C80000000000000)),
    ((u:$3FF9EFDCA62B700A),(u:$3FF104B648F113A1),(u:$3C90000000000000)),
    ((u:$3FF9EFDCA62B700A),(u:$3FF104B648F113A1),(u:$3C90000000000000)),
    ((u:$3FFA5BF3ACFDE4B2),(u:$3FF1585720F35CD9),(u:QWord($BC90000000000000))),
    ((u:$400D888DD2101D93),(u:$3FFFAF8B7A12CF9F),(u:QWord($BC90000000000000))),
    ((u:$4040151DEF34C2B8),(u:$4010A7B6E3FED720),(u:$3CB0000000000000)),
    ((u:$40645EA160DDC71F),(u:$401725811DCF6782),(u:$3CB0000000000000)),
    ((u:$40813570067ACC9F),(u:$401C04672343DCCF),(u:QWord($BCB0000000000000))),
    ((u:$4092A686E4B567CE),(u:$401F1C928E7F1E65),(u:$3CB0000000000000)),
    ((u:$40ECB62EEC26BD78),(u:$402759A2AD4C4D56),(u:$3CC0000000000000)));

  cAcoshBr0: UInt64 = $3FF1E83E425AEE63; // 0x1.1e83e425aee63p+0
  cAcoshBr1: UInt64 = $405BF00000000000; // 0x1.bfp+6 upper limit for main path
  cAcoshBr2: UInt64 = $4087100000000000; // 0x1.71p+9
  cAcoshBr3: UInt64 = $40E0100000000000; // 0x1.01p+15
  cAcoshBr4: UInt64 = $41EA000000000000; // 0x1.ap+31
  cAcoshRefineBr1: UInt64 = $4190000000000000; // 0x1p+26
  cAcoshRefineBr2: UInt64 = $4330000000000000; // 0x1p+52
  cAcoshInvLog2:   Tb64u64 = (u:$3FF71547652B82FE); // 0x1.71547652b82fep+0
  cAcoshVOffset:   Tb64u64 = (u:$3FF0000800000000); // 0x1.00008p+0
  cAcoshFastEps1:  Tb64u64 = (u:$3CCFC00000000000); // 0x1.fcp-51
  cAcoshFastEps2:  Tb64u64 = (u:$3970000000000000); // 0x1p-104
  cAcoshMainEps:   Double  = 2.8e-19;

// Local helper: adddd — double-double + double-double addition.
function AcoshAddDD(xh, xl, ch, cl: Double; out l: Double): Double; inline;
var s, d: Double;
begin
  s := xh + ch;
  d := s - xh;
  l := ((ch - d) + (xh + (d - s))) + (xl + cl);
  Result := s;
end;

function AcoshDatabase(x, f: Double): Double;
var a, b, m: Int32;
begin
  a := 0; b := High(cAcoshDb);
  m := (a + b) div 2;
  while a <= b do
  begin
    if cAcoshDb[m, 0].f < x then a := m + 1
    else if cAcoshDb[m, 0].f = x then
    begin
      f := cAcoshDb[m, 1].f + cAcoshDb[m, 2].f;
      Break;
    end
    else b := m - 1;
    m := (a + b) div 2;
  end;
  Result := f;
end;

// as_acosh_one: near-one fast-refine path
function AcoshOne(x, sh, sl: Double): Double;
var
  y0, y1, y2: Double;
  ch, cl, th, tl, xl0: Double;
  i: Int32;
begin
  y2 := x * (cAcoshOneCl[0].f + x * (cAcoshOneCl[1].f + x * (cAcoshOneCl[2].f
      + x * (cAcoshOneCl[3].f + x * (cAcoshOneCl[4].f + x * cAcoshOneCl[5].f)))));
  // Inlined seeded polydd: n=10, xh=x, xl=0, seed *l = y2
  xl0 := Double(0.0);
  ch := cAcoshOneCh[9, 0].f + y2;
  cl := ((cAcoshOneCh[9, 0].f - ch) + y2) + cAcoshOneCh[9, 1].f;
  for i := 8 downto 0 do
  begin
    ch := pcr_muldd(x, xl0, ch, cl, cl);
    th := ch + cAcoshOneCh[i, 0].f;
    tl := (cAcoshOneCh[i, 0].f - th) + ch;
    ch := th;
    cl := cl + tl + cAcoshOneCh[i, 1].f;
  end;
  y2 := cl;
  y1 := ch;
  y1 := pcr_mulddd_pd(y1, y2, x, y2);
  pcr_fasttwosum(y0, y1, Double(1.0), y1);
  y1 := y1 + y2;
  y0 := pcr_muldd(y0, y1, sh, sl, y1);
  Result := y0 + y1;
end;

function AcoshRefine(x, a: Double): Double;
var
  ix, t, v, t0, w: Tb64u64;
  zh, zl, x2h, x2l, wh, wl, sh, ish, sl: Double;
  ex, e, i1, i2, i3, i4: Int32;
  iv: UInt64;
  ed, el0, el1, el2: Double;
  L0, L1, L2: Double;
  t12, t34, th, tl, dh, dl, shm, slm: Double;
  xh, xl: Double;
 xh2, sh_: Double;
  sh2, sl2: Double;
  v0, v1, v2, vtmp: Double;
  er, de: UInt64;
  res: Double;
begin
  ix.f := x;
  if ix.u < cAcoshRefineBr1 then
  begin
    x2h := x * x;
    x2l := pcr_fma(x, x, -x2h);
    wh  := x2h - Double(1.0);
    pcr_fasttwosum(wh, wl, wh, x2l);
    sh  := pcr_sqrt(wh);
    ish := Double(0.5) / wh;
    sl  := (ish * sh) * (wl - pcr_fma(sh, sh, -wh));
    pcr_fasttwosum(zh, zl, x, sh);
    zl  := zl + sl;
    pcr_fasttwosum(zh, zl, zh, zl);
  end
  else if ix.u < cAcoshRefineBr2 then
  begin
    zh := Double(2.0) * x;
    zl := Double(-0.5) / x;
  end
  else
  begin
    zh := x;
    zl := Double(0.0);
  end;

  t.f := zh;
  ex := (t.u shr 52) and $7FF;
  if zl = Double(0.0) then
    e := ex - $3FF + 1
  else
    e := ex - $3FF;
  t.u := t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12);
  t.u := t.u or (UInt64($3FF) shl 52);
  ed := e;

  v.f := a - ed + cAcoshVOffset.f;
  iv := (v.u - (UInt64($3FF) shl 52)) shr (52 - 16);
  i1 := Int32((iv shr 12) and $1F);
  i2 := Int32((iv shr  8) and $F);
  i3 := Int32((iv shr  4) and $F);
  i4 := Int32(iv and $F);

  el2 := cAcoshL22.f * ed;
  el1 := cAcoshL21.f * ed;
  el0 := cAcoshL20.f * ed;

  L0 := cAcoshLL[0, i1, 0].f + cAcoshLL[1, i2, 0].f
      + (cAcoshLL[2, i3, 0].f + cAcoshLL[3, i4, 0].f);
  L1 := cAcoshLL[0, i1, 1].f + cAcoshLL[1, i2, 1].f
      + (cAcoshLL[2, i3, 1].f + cAcoshLL[3, i4, 1].f);
  L2 := cAcoshLL[0, i1, 2].f + cAcoshLL[1, i2, 2].f
      + (cAcoshLL[2, i3, 2].f + cAcoshLL[3, i4, 2].f);
  L0 := L0 + el0;

  t12 := cAcoshT1[i1].f * cAcoshT2[i2].f;
  t34 := cAcoshT3[i3].f * cAcoshT4[i4].f;
  th  := t12 * t34;
  tl  := pcr_fma(t12, t34, -th);
  dh  := th * t.f;
  dl  := pcr_fma(th, t.f, -dh);
  shm := tl * t.f;
  slm := pcr_fma(tl, t.f, -shm);

  pcr_fasttwosum(xh, xl, dh - Double(1.0), dl);
  if zl <> Double(0.0) then
  begin
    t.f := zl;
    t.u := UInt64(Int64(t.u) - (Int64(e) shl 52));
    xl := xl + th * t.f;
  end;
  xh := AcoshAddDD(xh, xl, shm, slm, xl);

  // sl = xh*(cl[0] + xh*(cl[1] + xh*cl[2])); sh = polydd(xh, xl, 3, ch, &sl);
  sl2 := xh * (cAcoshRefineCl[0].f + xh * (cAcoshRefineCl[1].f + xh * cAcoshRefineCl[2].f));
  // Inlined seeded polydd: n=3, seed *l = sl2
  xh2 := cAcoshRefineCh[2, 0].f + sl2;
  sh2 := ((cAcoshRefineCh[2, 0].f - xh2) + sl2) + cAcoshRefineCh[2, 1].f;
  for i1 := 1 downto 0 do   // reuse i1 var; safe
  begin
    xh2 := pcr_muldd(xh, xl, xh2, sh2, sh2);
    th  := xh2 + cAcoshRefineCh[i1, 0].f;
    tl  := (cAcoshRefineCh[i1, 0].f - th) + xh2;
    xh2 := th;
    sh2 := sh2 + tl + cAcoshRefineCh[i1, 1].f;
  end;
  sh_ := pcr_muldd(xh, xl, xh2, sh2, sl2);
  sh_ := AcoshAddDD(sh_, sl2, el1, el2, sl2);
  sh_ := AcoshAddDD(sh_, sl2, L1, L2, sl2);
  pcr_fasttwosum(v0, vtmp, L0, sh_);
  pcr_fasttwosum(v1, v2, vtmp, sl2);
  v0 := v0 * Double(2.0);
  v1 := v1 * Double(2.0);
  v2 := v2 * Double(2.0);

  t.f := v1;
  if (t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12)) = 0 then
  begin
    w.f := v2;
    if ((w.u xor t.u) shr 63) <> 0 then t.u := t.u - 1
    else t.u := t.u + 1;
    v1 := t.f;
  end;

  t0.f := v0;
  er := ((t.u + 7) and (UInt64($FFFFFFFFFFFFFFFF) shr 12));
  de := ((t0.u shr 52) and $7FF) - ((t.u shr 52) and $7FF);
  res := v0 + v1;
  if (de > 102) or (er < 15) then
    Result := AcoshDatabase(x, res)
  else
    Result := res;
end;

function pcr_acosh(x: Double): Double;
var
 ix, t: Tb64u64;
  aix, iuint, jshift: UInt64;
  off: Int32;
  z, iz, zt, sh, sl: Double;
  z2, z4, ds, eps, lb, ub: Double;
  g, x2h, wh, wl, ish: Double;
  tl, th: Double;
  ex, e, i1, i2: Int32;
  dsigned: Int64;
  ed, r, dx, dx2, f: Double;
  lh, ll: Double;
begin
  ix.f := x;
  // x < 0 or NaN/Inf
  if ix.u >= UInt64($7FF0000000000000) then
  begin
    aix := ix.u shl 1;
    if (ix.u = UInt64($7FF0000000000000)) or (aix > (UInt64($7FF) shl 53)) then
    begin
      Result := x + x;
      Exit;
    end;
    Result := Double(0.0) / Double(0.0);
    Exit;
  end;
  // x <= 1
  if Int64(ix.u) <= Int64($3FF0000000000000) then
  begin
    if ix.u = UInt64($3FF0000000000000) then begin Result := Double(0.0); Exit; end;
    Result := Double(0.0) / Double(0.0);
    Exit;
  end;

  g := Double(0.0);
  off := $3FE;
  t.u := ix.u;

  if ix.u < cAcoshBr0 then
  begin
    // 1 < x < 0x1.1e83e425aee63p+0
    z  := x - Double(1.0);
    iz := Double(-0.25) / z;
    zt := Double(2.0) * z;
    sh := pcr_sqrt(zt);
    sl := pcr_fma(sh, sh, -zt) * (sh * iz);
    z2 := z * z;
    z4 := z2 * z2;
    ds := pcr_fma(sh * z,
           (cAcoshOneFastCl[0].f + z * (((cAcoshOneFastCl[1].f + z * cAcoshOneFastCl[2].f)
             + z2 * (cAcoshOneFastCl[3].f + z * cAcoshOneFastCl[4].f))
             + z4 * ((cAcoshOneFastCl[5].f + z * cAcoshOneFastCl[6].f)
                    + z2 * (cAcoshOneFastCl[7].f + z * cAcoshOneFastCl[8].f)))),
           sl);
    eps := ds * cAcoshFastEps1.f - cAcoshFastEps2.f * sh;
    lb := sh + (ds - eps);
    ub := sh + (ds + eps);
    if lb = ub then begin Result := lb; Exit; end;
    Result := AcoshOne(z, sh, sl);
    Exit;
  end
  else if ix.u < cAcoshBr1 then
  begin
    off := $3FF;
    x2h := x * x;
    wh  := x2h - Double(1.0);
    wl  := pcr_fma(x, x, -x2h);
    sh  := pcr_sqrt(wh);
    ish := Double(0.5) / wh;
    sl  := (wl - pcr_fma(sh, sh, -wh)) * (sh * ish);
    pcr_fasttwosum(th, tl, x, sh);
    tl := tl + sl;
    t.f := th;
    g := tl / th;
  end
  else if ix.u < cAcoshBr2 then
  begin
    z := Double(1.0) / (x * x);
    g := cAcoshMed4[0].f + z * (cAcoshMed4[1].f + z * (cAcoshMed4[2].f + z * cAcoshMed4[3].f));
  end
  else if ix.u < cAcoshBr3 then
  begin
    z := Double(1.0) / (x * x);
    g := cAcoshMed3[0].f + z * (cAcoshMed3[1].f + z * cAcoshMed3[2].f);
  end
  else if ix.u < cAcoshBr4 then
  begin
    z := Double(1.0) / (x * x);
    g := cAcoshMed2[0].f + z * cAcoshMed2[1].f;
  end
  else
  begin
    g := Double(0.0);
  end;

  ex := (t.u shr 52) and $7FF;
  e  := ex - off;
  t.u := t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 12);
  ed := e;
  iuint := t.u shr (52 - 5);
  dsigned := Int64(t.u and (UInt64($FFFFFFFFFFFFFFFF) shr 17));
  jshift := UInt64(
              Int64(t.u)
              + (Int64(cAcoshB[iuint].c0) shl 33)
              + (Int64(cAcoshB[iuint].c1) * (dsigned shr 16))
            ) shr (52 - 10);
  t.u := t.u or (UInt64($3FF) shl 52);
  i1 := Int32(jshift shr 5);
  i2 := Int32(jshift and $1F);
  r := cAcoshR1[i1].f * cAcoshR2[i2].f;
  dx := pcr_fma(r, t.f, -Double(1.0));
  dx2 := dx * dx;
  f := dx2 * ((cAcoshC[0].f + dx * cAcoshC[1].f)
            + dx2 * ((cAcoshC[2].f + dx * cAcoshC[3].f) + dx2 * cAcoshC[4].f));
  lh := (cAcoshL1[i1, 1].f + cAcoshL2[i2, 1].f) + cAcoshL2h.f * ed;
  ll := dx + cAcoshL2l.f * ed;
  ll := ll + g;
  ll := ll + cAcoshL1[i1, 0].f + cAcoshL2[i2, 0].f;
  ll := ll + f;
  eps := cAcoshMainEps;
  lb := lh + (ll - eps);
  ub := lh + (ll + eps);
  if lb = ub then begin Result := lb; Exit; end;
  Result := AcoshRefine(x, cAcoshInvLog2.f * lb);
end;
