import struct
def enc(f):
    u = int.from_bytes(struct.pack('<d', f), 'little')
    return f'${u:016X}'

def H(h): return float.fromhex(h)

# core constants
core = [
 ('ONE_OVER_PIH  0x1.45f306dc9c883p-2',  H('0x1.45f306dc9c883p-2')),
 ('ONE_OVER_PIL -0x1.6b01ec5417056p-56', -H('0x1.6b01ec5417056p-56')),
 ('ONE_OVER_3PI  0x1.b2995e7b7b604p-4',  H('0x1.b2995e7b7b604p-4')),
 ('Asympt_Pos    0x1.0000000000001p-55', H('0x1.0000000000001p-55')),
 ('Asympt_Neg    0x1.fffffffffffffp-56', H('0x1.fffffffffffffp-56')),
 ('Halfp55       0x1p-55',               H('0x1p-55')),
 ('Small_X       0x1.5cba89af1f855p-1022', H('0x1.5cba89af1f855p-1022')),
 ('Small_Out     0x1.bc03df34e902cp-1024', H('0x1.bc03df34e902cp-1024')),
 ('P600   0x1p-600',  H('0x1p-600')),
 ('Pp106  0x1p106',   H('0x1p+106')),
 ('Pm106  0x1p-106',  H('0x1p-106')),
 ('Pm56   0x1p-56',   H('0x1p-56')),
 ('EFact  0x1.41p-52', H('0x1.41p-52')),
 ('EFactSm 0x1.6fp-68', H('0x1.6fp-68')),
 ('Thresh_Big 0x1.bep20', H('0x1.bep20')),
 ('Thresh_Tiny 0x1.c7p-27', H('0x1.c7p-27')),
 ('PiHalfH 0x1.921fb54442d18p+0', H('0x1.921fb54442d18p+0')),
 ('PiHalfL 0x1.1a62633145c07p-54', H('0x1.1a62633145c07p-54')),
 ('Half 0x1p-1', H('0x1p-1')),
 ('Quarter 0x1p-2', H('0x1p-2')),
 ('Underflow_lim 0x1.921fb54442d17p-1021', H('0x1.921fb54442d17p-1021')),
 ('Pm1022 0x1p-1022', H('0x1p-1022')),
]
for n,v in core:
    print(f'  {n:42s} {enc(v)}')

print('\n--- TINY EXCEPTIONS (20 rows, 3 cols) ---')
tiny = [
 ('0x1.9b46cfb0f10cap-37',  '0x1.05d3aa746e163p-38',  '-0x1.94a704a0e798p-149'),
 ('0x1.eac33f52580c1p-31',  '0x1.386de1f01225ap-32',  '-0x1.2be95e5723a66p-139'),
 ('0x1.56a4aa740a5a7p-53',  '0x1.b44453e2404e7p-55',  '-0x1.c5376c33c7a5p-164'),
 ('0x1.9d47b2e2a4deep-38',  '0x1.071a2e1e2e7ap-39',   '-0x1.e948955ed137ap-146'),
 ('0x1.932ddaeecdabep-44',  '0x1.00abfb23406c5p-45',  '-0x1.fffffffffffffp-99'),
 ('0x1.e0b9044e5a1cfp-46',  '0x1.32099f7c94e82p-47',  '-0x1.22752cf9b0bc5p-152'),
 ('0x1.13b412ea2182dp-50',  '0x1.5f095755008fp-52',   '-0x1.a8ad4502b5f1p-160'),
 ('0x1.ff878736a1bcdp-43',  '0x1.45a65501d94d3p-44',  '-0x1.7cb6951d7806dp-152'),
 ('0x1.4f03bd4cedd49p-38',  '0x1.aa8dc8da24739p-40',  '-0x1.fe11160ec4ee9p-146'),
 ('0x1.0f2bf03e38a9ep-36',  '0x1.594435f8c0326p-38',  '-0x1.8eef034a3526ep-145'),
 ('0x1.508ecb38f52f9p-52',  '0x1.ac84c88f979a2p-54',  '-0x1.63773cbd64d5cp-164'),
 ('0x1.8ddcd3207fe6p-42',   '0x1.fa930b4eb436cp-44',   '0x1p-97'),
 ('0x1.1331b44547444p-48',  '0x1.5e63596a62718p-50',  '-0x1.4a9f744d7cea7p-156'),
 ('0x1.f7681fa730e73p-31',  '0x1.407a8ac1c768ap-32',  '-0x1.3c6a80c812c54p-138'),
 ('0x1.09df162e4049cp-36',  '0x1.52849d0b1e523p-38',  '-0x1.fffffffffffffp-92'),
 ('0x1.59af9a1194efep-54',  '0x1.b824198b94a89p-56',   '0x1.fffffffffffffp-110'),
 ('0x1.65f203d6bb18ep-46',  '0x1.c7c00d31634a6p-48',  '-0x1.d4af6699e375fp-157'),
 ('0x1.6f89f6df8258dp-34',  '0x1.d3f710390db87p-36',  '-0x1.936061594fad7p-145'),
 ('0x1.3e4d2d87b5aefp-51',  '0x1.954626979d7d3p-53',  '-0x1.6f8dbd4f38cdbp-161'),
 ('0x1.f46e9e25a51d1p-36',  '0x1.3e95c097e08ep-37',   '-0x1.8e5f60920b571p-144'),
]
print('cAtanpiTinyExc:')
for row in tiny:
    cells = []
    for col in row:
        col = col.strip()
        sign = 1.0
        if col.startswith('-'):
            sign = -1.0
            col = col[1:]
        v = sign * H(col)
        cells.append(enc(v))
    print(f'    (({cells[0]}),({cells[1]}),({cells[2]})),')

print('\n--- REFINE EXCEPTIONS (56 rows, 3 cols) ---')
refine = [
 ('0x1.16265f83e8121p+2', '0x1.b648e5b0804dap-2', '0x1.9db921c686893p-109'),
 ('0x1.8a5c063607fbap+6', '0x1.fcb1a9f45a781p-2', '-0x1.fffffffffffffp-56'),
 ('0x1.1a5e744ff1806p-6', '0x1.677cdc7dea093p-8', '0x1.417d05c2a4506p-115'),
 ('0x1.2eee26511c903p-7', '0x1.81b11ac01389p-9', '0x1.508c83c3db5e3p-116'),
 ('0x1.bb0ba10084de6p-8', '0x1.1a0c1411a3dc1p-9', '0x1.a854b560fe4dbp-115'),
 ('0x1.55166a4d7d836p-2', '0x1.a33741abea942p-4', '-0x1.54a4c868a708ep-110'),
 ('0x1.5fa4b5e7ce94p-3', '0x1.bb6715a0b1f81p-5', '0x1.00037a0bfdfaep-111'),
 ('0x1.2a04fa64d36p+4', '0x1.ee846d3680301p-2', '0x1.3015b69c79ef4p-111'),
 ('0x1.0fb7e1298320cp+4', '0x1.ecd42a2ef8667p-2', '-0x1.21d6585c13d1p-107'),
 ('0x1.b6daa66352e75p+3', '0x1.e84651d4c9934p-2', '0x1.34e0b5a33a1dbp-108'),
 ('0x1.5717df9abccd3p+1', '0x1.8b9be647e02a9p-2', '-0x1.6c5194938e976p-109'),
 ('0x1.e0ac5527e7f0dp-4', '0x1.309ce21202ccap-5', '-0x1.141be39e1919dp-111'),
 ('0x1.932f168987539p+3', '0x1.e62f20aba165cp-2', '-0x1.839c5638718f2p-109'),
 ('0x1.2e342ca2440eep+3', '0x1.dd9d1ecffc9ap-2', '-0x1.53deddde1dde1p-109'),
 ('0x1.b37ef5831a537p-5', '0x1.14fc0e0c2dacep-6', '-0x1.f6b67ef93b398p-113'),
 ('0x1.d8ea551d1d795p-8', '0x1.2d0fe21df78dp-9', '-0x1.0f4d438d8dd24p-114'),
 ('0x1.badb08cd8b1c5p-3', '0x1.15a7db07daa97p-4', '-0x1.fffffffffffffp-58'),
 ('0x1.18071548e23a6p+7', '0x1.fdac0c22a86c5p-2', '0x1.fffffffffffffp-56'),
 ('0x1.e8b9c1fe81ca3p+6', '0x1.fd55124d78891p-2', '-0x1.ffffffffffffep-56'),
 ('0x1.e2faab41052e7p-6', '0x1.33629443d32f5p-7', '0x1.fffffffffffffp-61'),
 ('0x1.03f709107440dp-7', '0x1.4afda1e054259p-9', '0x1.ffffffffffffep-63'),
 ('0x1.928a1ca94d1dcp-15', '0x1.0043bd1fd5baep-16', '-0x1.12e31518bad19p-122'),
 ('0x1.deda0ddb0795cp+11', '0x1.ffea37cae8f37p-2', '-0x1.ffffffffffffep-56'),
 ('0x1.72f4cb85b8404p-11', '0x1.d850e9850ab33p-13', '-0x1.014404207eb42p-119'),
 ('0x1.4b3c41e10a1d4p-21', '0x1.a5bdf2b61615bp-23', '-0x1.180c42f2c36bbp-131'),
 ('0x1.4ad2f60d5240bp-9', '0x1.a537a6de092a1p-11', '-0x1.f2e112f36eca4p-117'),
 ('0x1.1151841058e68p+20', '0x1.ffffeceb437e6p-2', '-0x1.2dd2f50460971p-108'),
 ('0x1.8b012d530c5c8p-20', '0x1.f6ef7b3ee090cp-22', '-0x1.9d79d9fb153c4p-128'),
 ('0x1.0fa189169e75ep+9', '0x1.ff666764f1751p-2', '-0x1.da2823e32e295p-108'),
 ('0x1.46c64e8717ccfp+7', '0x1.fe014cac09e73p-2', '0x1.fffffffffffffp-56'),
 ('0x1.5a9ee4b581471p+12', '0x1.fff0f445bee08p-2', '-0x1.4256eae14a448p-108'),
 ('0x1.6ab7e3b9e35cdp+12', '0x1.fff19f36d9901p-2', '-0x1.246ecd835d1c3p-107'),
 ('0x1.842e540b07ac1p-13', '0x1.ee3f52b4af0dcp-15', '-0x1.8c285bfa8c0ecp-121'),
 ('0x1.5d2c01b65f659p+18', '0x1.ffffc441a900dp-2', '-0x1.8137cc73dd777p-107'),
 ('0x1.9d18f49824d6ap+17', '0x1.ffff9b00d35d1p-2', '0x1.fe1e36287467bp-110'),
 ('0x1.400b2b748da28p+20', '0x1.ffffefb46ade3p-2', '-0x1.10ca50f3bfef4p-111'),
 ('0x1.27b786b2f1baep+20', '0x1.ffffee5d3f843p-2', '-0x1.8afee3bdfb07bp-108'),
 ('0x1.bbc4a4251b4b9p-20', '0x1.1a82f5f7f86bdp-21', '0x1.fffffffffffffp-75'),
 ('0x1.72a01aef844e5p-11', '0x1.d7e514f3ffe99p-13', '0x1.7c2d5be2b88d8p-119'),
 ('0x1.95e2b447c5ed5p+13', '0x1.fff99356bf911p-2', '0x1.fffffffffffffp-56'),
 ('0x1.88d8109598821p+19', '0x1.ffffe572f576ap-2', '0x1.5430c0b185e51p-108'),
 ('0x1.3ef59f63222ccp-23', '0x1.961c9f0bc5519p-25', '0x1.fffffffffffffp-79'),
 ('0x1.72f2198034e17p+18', '0x1.ffffc7c36896ep-2', '-0x1.adfc30e218aafp-109'),
 ('0x1.df1c8e36bbccp-27', '0x1.31030abf0398fp-28', '0x1.ffffffffffffep-82'),
 ('0x1.ad5781ab29bbap+9', '0x1.ff9ed31438683p-2', '0x1.fffffffffffffp-56'),
 ('0x1.ede23d6624a3ap+11', '0x1.ffeae182f312ep-2', '-0x1.26eee3b0a2ccbp-107'),
 ('0x1.4a7e234cdce9bp+19', '0x1.ffffe0709e4e1p-2', '-0x1.ffffffffffffdp-56'),
 ('0x1.1c404899ca42p-23', '0x1.69eb74ce9e7edp-25', '0x1.d126a89bf70e8p-131'),
 ('0x1.7422a8342383p+8', '0x1.ff1fc5deb1c5cp-2', '-0x1.0cf8b48b700a5p-108'),
 ('0x1.f37be8c1f5dfp-15', '0x1.3dfb3d39744cfp-16', '-0x1.cae8f0c659368p-123'),
 ('0x1.8c9f8999242b7p+20', '0x1.fffff2d9dca7bp-2', '0x1.fffffffffffffp-56'),
 ('0x1.427f51dcf1d68p+19', '0x1.ffffdfa850495p-2', '0x1.ffffffffffffdp-56'),
 ('0x1.abb7b457cdf84p+20', '0x1.fffff3ce932e2p-2', '0x1.6f17a18f91149p-112'),
 ('0x1.be630d1090fe6p-12', '0x1.1c2dc0a9340efp-13', '-0x1.8147b0a8785f5p-121'),
 ('0x1.27f67f799ccbdp+13', '0x1.fff730800f741p-2', '0x1.fffffffffffffp-56'),
 ('0x1.3f91b5b198cddp+9', '0x1.ff7d71d1d656bp-2', '-0x1.c90dde6f07a4p-110'),
]
print('cAtanpiRefExc:')
for row in refine:
    cells = []
    for col in row:
        col = col.strip()
        sign = 1.0
        if col.startswith('-'):
            sign = -1.0
            col = col[1:]
        v = sign * H(col)
        cells.append(enc(v))
    print(f'    (({cells[0]}),({cells[1]}),({cells[2]})),')
