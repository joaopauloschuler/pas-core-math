import struct
def enc(f):
    u = int.from_bytes(struct.pack('<d', f), 'little')
    return f'${u:016X}'

vals = [
 ('Ph 0x1.921fb54442d18p+1', float.fromhex('0x1.921fb54442d18p+1')),
 ('Pl 0x1.1a62633145c07p-53', float.fromhex('0x1.1a62633145c07p-53')),
 ('FastC[0] -0x1.4abbce625be51p+2', -float.fromhex('0x1.4abbce625be51p+2')),
 ('FastC[1] 0x1.466bc67754b46p+1', float.fromhex('0x1.466bc67754b46p+1')),
 ('FastC[2] -0x1.32d2cc12a51f4p-1', -float.fromhex('0x1.32d2cc12a51f4p-1')),
 ('FastC[3] 0x1.5060540058476p-4', float.fromhex('0x1.5060540058476p-4')),
 ('FastEps1 0x1p-47', float.fromhex('0x1p-47')),
 ('FastEps2 0x1p-102', float.fromhex('0x1p-102')),
 ('FastEr 5.5e-19', 5.5e-19),
 ('2p106', float.fromhex('0x1p+106')),
 ('2pm106', float.fromhex('0x1p-106')),
 ('Pi0 0x1.92p+1', float.fromhex('0x1.92p+1')),
 ('Pi1 0x1.fb54442d1846ap-11', float.fromhex('0x1.fb54442d1846ap-11')),
 ('Pi2 -0x1.d9cceba3f91f2p-65', -float.fromhex('0x1.d9cceba3f91f2p-65')),
 ('ZeroCh[0][0] -0x1.4abbce625be53p+2', -float.fromhex('0x1.4abbce625be53p+2')),
 ('ZeroCh[0][1] 0x1.05511c68477bep-52', float.fromhex('0x1.05511c68477bep-52')),
 ('ZeroCh[1][0] 0x1.466bc6775aae2p+1', float.fromhex('0x1.466bc6775aae2p+1')),
 ('ZeroCh[1][1] -0x1.6dc0cbefae1dap-54', -float.fromhex('0x1.6dc0cbefae1dap-54')),
 ('ZeroCh[2][0] -0x1.32d2cce62bd86p-1', -float.fromhex('0x1.32d2cce62bd86p-1')),
 ('ZeroCh[2][1] 0x1.066bd54973829p-55', float.fromhex('0x1.066bd54973829p-55')),
 ('ZeroCh[3][0] 0x1.50783487ee781p-4', float.fromhex('0x1.50783487ee781p-4')),
 ('ZeroCh[3][1] 0x1.832989f39a743p-58', float.fromhex('0x1.832989f39a743p-58')),
 ('ZeroCl[0] -0x1.e3074fde861fp-8', -float.fromhex('0x1.e3074fde861fp-8')),
 ('ZeroCl[1] 0x1.e8f4344534da6p-12', float.fromhex('0x1.e8f4344534da6p-12')),
 ('ZeroCl[2] -0x1.6f9cd7b8cb9dbp-16', -float.fromhex('0x1.6f9cd7b8cb9dbp-16')),
 ('Db4_X[0] -0x1.276b3fef466p-2', -float.fromhex('0x1.276b3fef466p-2')),
 ('Db4_X[1] -0x1.33caea0f24cp-2', -float.fromhex('0x1.33caea0f24cp-2')),
 ('Db4_X[2] -0x1.8a1e8a3e82cp-1', -float.fromhex('0x1.8a1e8a3e82cp-1')),
 ('Db4_X[3] -0x1.bdd02d1ad60p-2', -float.fromhex('0x1.bdd02d1ad60p-2')),
 ('Db4_R[0] 0x1.db8a79a80c3a0p-4', float.fromhex('0x1.db8a79a80c3a0p-4')),
 ('Db4_R[1] 0x1.5146c0bc45bcep-3', float.fromhex('0x1.5146c0bc45bcep-3')),
 ('Db4_R[2] 0x1.5d6561936b699p-1', float.fromhex('0x1.5d6561936b699p-1')),
 ('Db4_R[3] 0x1.f72c906962631p-1', float.fromhex('0x1.f72c906962631p-1')),
 ('Db4_D[0] 0x1p-110', float.fromhex('0x1p-110')),
 ('Db4_D[1] 0x1p-109', float.fromhex('0x1p-109')),
 ('Db4_D[2] -0x1p-55', -float.fromhex('0x1p-55')),
 ('Db4_D[3] 0x1p-55', float.fromhex('0x1p-55')),
]
for n,v in vals:
    print(f'  {n:45s} {enc(v)}')
