#!/usr/bin/env python3
import struct

def enc(s):
    # hex-float literal like "0x1.2p+3" or "-0x1.2p-3" or "0x0p+0"
    v = float.fromhex(s)
    return struct.unpack('<Q', struct.pack('<d', v))[0]

def tb64(label, s):
    u = enc(s)
    print(f"    ((u:${u:016X}),(u:${enc('0')})),  // {label}")

def arr3(db):
    for i, row in enumerate(db):
        a,b,c = row
        print(f"    ((u:${enc(a):016X}),(u:${enc(b):016X}),(u:${enc(c):016X})){','  if i < len(db)-1 else ''}  // {i}")

# Fast-path c[5]
fast = ['0x1.5555555555555p-3', '0x1.111111111151ep-7', '0x1.a01a019d0c767p-13', '0x1.71de444a96e11p-19', '0x1.ae8465375242p-26']
for i, s in enumerate(fast):
    print(f"  cSinhFastC{i}: Tb64u64 = (u:${enc(s):016X}); // {s}")

print()
# as_sinh_zero ch[5][2]
zch = [
    ('0x1.5555555555555p-3', '0x1.555555555552fp-57'),
    ('0x1.1111111111111p-7', '0x1.11111115cf00dp-63'),
    ('0x1.a01a01a01a01ap-13', '0x1.a0011c925b85cp-73'),
    ('0x1.71de3a556c734p-19', '-0x1.b4e2835532bcdp-73'),
    ('0x1.ae64567f54482p-26', '-0x1.defcf17a6ab79p-81'),
]
print("  cSinhZeroCh: array[0..4, 0..1] of Tb64u64 = (")
for i,(a,b) in enumerate(zch):
    sep = ',' if i < len(zch)-1 else ''
    print(f"    ((u:${enc(a):016X}),(u:${enc(b):016X})){sep}")
print("  );")

print()
# as_sinh_zero cl[3]
cl = ['0x1.6124613aef206p-33', '0x1.ae7f36beea815p-41', '0x1.95785063cd974p-49']
for i, s in enumerate(cl):
    print(f"  cSinhZeroCl{i}: Tb64u64 = (u:${enc(s):016X}); // {s}")

print()
# eps constants
print(f"  cSinhFastEps: Tb64u64 = (u:${enc('0x1.cp-53'):016X}); // 0x1.cp-53 (fast-path)")
print(f"  cSinhErrLrg:  Tb64u64 = (u:${enc('0x1.1b6p-63'):016X}); // 0x1.1b6p-63")
# cSinhErrMed = 0x1.202p-63 same as cCoshErrMed
# cSinhErrSml = 0x1.c0ap-62 same as cCoshErrSml

# thresholds
print(f"  cSinhAxTiny:   UInt64 = $3FD0000000000000; // 0.25")
print(f"  cSinhAxTinier: UInt64 = $3E57137449123EF7; // 0x1.7137449123ef7p-26")

print()
# Database 50 entries of (x, h, l)
db = [
 ('0x1.364303e1ad8f6p-2','0x1.3b07e0c779ddap-2','-0x1.bcp-106'),
 ('0x1.4169f234f23b9p-2','0x1.46b7b3b358f99p-2','-0x1p-56'),
 ('0x1.616cc75d49226p-2','0x1.687bd068c1c1ep-2','0x1.ap-111'),
 ('0x1.ae3773250e7d2p-2','0x1.bafc3479fc9ccp-2','-0x1p-105'),
 ('0x1.b7efa91915c95p-2','0x1.c59869f17b483p-2','-0x1p-104'),
 ('0x1.d68039861ab53p-2','0x1.e73b46abb01e1p-2','-0x1.2p-109'),
 ('0x1.e90f16eb88c09p-2','0x1.fbdd4a37760b7p-2','-0x1.f8p-108'),
 ('0x1.a3fc7e4dd47d1p-1','0x1.d4b21ebf542fp-1','0x1.ep-107'),
 ('0x1.aa3b649a96091p-1','0x1.dd32c5ed1e93p-1','0x1.8ap-106'),
 ('0x1.c13876341b62ep-1','0x1.fd1d7f1c8170cp-1','0x1.72p-105'),
 ('0x1.2f5d3b178914ap+0','0x1.7b8516ffd2406p+0','-0x1.28p-104'),
 ('0x1.3ffc12b81cbc2p+0','0x1.9a0ff413a1af3p+0','0x1.cp-107'),
 ('0x1.44f65dff00782p+0','0x1.a38a3c3227609p+0','-0x1p-103'),
 ('0x1.7346e3c591a14p+0','0x1.01e9cfa77b855p+1','0x1.p-102'),
 ('0x1.b6e2c73f41415p+0','0x1.57e377b3f0b4bp+1','-0x1p-102'),
 ('0x1.dc5059d4e507dp+0','0x1.9168c60ed5256p+1','0x1.c6p-104'),
 ('0x1.f737f1e8378c7p+0','0x1.bffd3f94f40fbp+1','0x1.a4p-104'),
 ('0x1.3359640329982p+1','0x1.5e40df3f985bep+2','0x1.97p-102'),
 ('0x1.58a4ff5adac35p+1','0x1.d671928665bddp+2','0x1p-102'),
 ('0x1.8c0a26d055288p+1','0x1.6056b06a21918p+3','-0x1.bep-102'),
 ('0x1.bc3c2d0c95f52p+1','0x1.00fef7383a978p+4','0x1.61p-100'),
 ('0x1.0a19aebb51e9p+3','0x1.fee8f69c4cd25p+10','0x1.48p-95'),
 ('0x1.3eb8f61734227p+3','0x1.4ab1cf45e4e26p+13','0x1p-90'),
 ('0x1.43a81752eabe7p+3','0x1.81d364845ecfap+13','-0x1p-90'),
 ('0x1.16369cd53bb69p+4','0x1.0fbc6c02b1c9p+24','-0x1.9p-81'),
 ('0x1.20e29ea8b51e2p+4','0x1.08b8abba28abcp+25','0x1.9bp-79'),
 ('0x1.92a5c27afbe82p+4','0x1.3c81f9a247253p+35','0x1p-67'),
 ('0x1.a1e4f11b513d7p+4','0x1.9a65b6c2e2185p+36','-0x1.bcp-70'),
 ('0x1.c089fcf166171p+4','0x1.5c452e0e37569p+39','0x1.4p-69'),
 ('0x1.e42a98b3a0be5p+4','0x1.938768ca4f8aap+42','0x1.6dp-62'),
 ('0x1.04db52248cbb8p+5','0x1.0794072349523p+46','0x1.0e8p-57'),
 ('0x1.21bc021eeb97ep+5','0x1.3065064a170fbp+51','0x1.088p-52'),
 ('0x1.39fc4d3bb711p+5','0x1.8a4e90733b95ep+55','0x1.6ep-50'),
 ('0x1.3c895d86e96c9p+5','0x1.0f33837882a6p+56','-0x1.28p-49'),
 ('0x1.e07e71bfcf06fp+5','0x1.91ec4412c344fp+85','0x1p-24'),
 ('0x1.f7216c4b435c9p+5','0x1.a97e7be23e65ap+89','-0x1p-15'),
 ('0x1.6474c604cc0d7p+6','0x1.7a8f65ad009bdp+127','-0x1.08p+20'),
 ('0x1.7a60ee15e3e9dp+6','0x1.62e4dc3bbf53fp+135','0x1.bp+29'),
 ('0x1.1f0da93354198p+7','0x1.0bd73b73fc74cp+206','0x1.59p+102'),
 ('0x1.54cd1fea7663ap+7','0x1.c90810d354618p+244','0x1.2p+135'),
 ('0x1.556c678d5e976p+7','0x1.37e7ac4e7f9b3p+245','0x1.02p+141'),
 ('0x1.7945e34b18a9ap+7','0x1.1b0e4936a8c9bp+271','-0x1.fap+166'),
 ('0x1.2da9e5e6af0bp+8','0x1.27d6fe867d6f6p+434','0x1.0ap+329'),
 ('0x1.54ceba01331d5p+8','0x1.9a86785b5ef3ep+490','-0x1.22p+386'),
 ('0x1.9e7b643238a14p+8','0x1.f5da7fe652978p+596','0x1p+493'),
 ('0x1.c7206c1b753e4p+8','0x1.8670de0b68cadp+655','-0x1.78p+548'),
 ('0x1.d6479eba7c971p+8','0x1.62a88613629b6p+677','-0x1.4p+568'),
 ('0x1.eb9914d4ac1c8p+8','0x1.2b67eff65dce8p+708','-0x1.02p+603'),
 ('0x1.0bc04af1b09f5p+9','0x1.7b1d97c902985p+771','0x1.56p+666'),
 ('0x1.26ee1a46d8c8bp+9','0x1.fbe20477df4a7p+849','-0x1.55p+745'),
 ('0x1.4a869881f72acp+9','0x1.9ea7540a3d1f9p+952','-0x1.2dp+848'),
]
assert len(db)==51
print(f"  cSinhDb: array[0..50, 0..2] of Tb64u64 = (")
for i,(a,b,c) in enumerate(db):
    sep = ',' if i<50 else ''
    print(f"    ((u:${enc(a):016X}),(u:${enc(b):016X}),(u:${enc(c):016X})){sep}")
print("  );")
