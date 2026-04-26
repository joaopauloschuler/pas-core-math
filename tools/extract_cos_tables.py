#!/usr/bin/env python3
"""Extract cos.c tables as Pascal const literals.

Reads core-math/src/binary64/cos/cos.c and emits Pascal-syntax constants for:
  T[20]   : uint64 2/(2pi) reduction table
  S[256]  : dint64_t sin2pi(i/2^11)
  C[256]  : dint64_t cos2pi(i/2^11)
  PSfast  : Double polynomial (6 entries)
  PCfast  : Double polynomial (6 entries)
  PS[6]   : dint64_t polynomial for sin
  PC[6]   : dint64_t polynomial for cos
  SC[256][3] : Double triples
Doubles are emitted as Tb64u64 records with exact IEEE-754 bit patterns.
"""
import re, struct, sys

COS_C = "/home/bpsa/app/core-math/src/binary64/cos/cos.c"

def hex_to_double(s):
    # C99 hex float literal -> Python float
    return float.fromhex(s)

def double_to_u64hex(d):
    return struct.unpack('<Q', struct.pack('<d', d))[0]

def parse_double(tok):
    tok = tok.strip().rstrip(',')
    # Accept forms like -0x1.8p+0
    return hex_to_double(tok)

def as_u64(d):
    return f"(u:${double_to_u64hex(d):016X})"

def read():
    with open(COS_C) as f:
        return f.read()

def extract_T(src):
    m = re.search(r'static const uint64_t T\[20\] = \{([^}]*)\};', src, re.S)
    body = m.group(1)
    # Entries are hex numbers separated by commas, with // comments
    vals = []
    for line in body.splitlines():
        line = re.sub(r'//.*', '', line).strip()
        for tok in line.split(','):
            tok = tok.strip()
            if tok:
                vals.append(int(tok, 16))
    assert len(vals) == 20, len(vals)
    return vals

def extract_dint_array(src, name, count):
    # match: static const dint64_t NAME[...] = { ... };
    pat = re.compile(r'static const dint64_t ' + re.escape(name) + r'\[[^\]]*\] = \{(.*?)\n\};', re.S)
    m = pat.search(src)
    if not m:
        raise RuntimeError(f"{name} not found")
    body = m.group(1)
    # Each entry: {.hi = 0x..., .lo = 0x..., .ex = N, .sgn = N},
    entry_pat = re.compile(
        r'\{\s*\.hi\s*=\s*(0x[0-9a-fA-F]+|0)\s*,\s*'
        r'\.lo\s*=\s*(0x[0-9a-fA-F]+|0)\s*,\s*'
        r'\.ex\s*=\s*(-?\d+)\s*,\s*'
        r'\.sgn\s*=\s*(-?\d+|0x[0-9a-fA-F]+)\s*\}'
    )
    entries = []
    for em in entry_pat.finditer(body):
        hi = int(em.group(1), 0)
        lo = int(em.group(2), 0)
        ex = int(em.group(3))
        sgn = int(em.group(4), 0)
        entries.append((hi, lo, ex, sgn))
    assert len(entries) == count, f"{name}: got {len(entries)}, want {count}"
    return entries

def extract_doubles_array(src, name, count, ctype='double'):
    pat = re.compile(r'static const\s+' + re.escape(ctype) + r'\s+' + re.escape(name) + r'\[\] = \{(.*?)\n\};', re.S)
    m = pat.search(src)
    if not m:
        raise RuntimeError(f"{name} not found")
    body = m.group(1)
    # strip comments
    body = re.sub(r'/\*.*?\*/', '', body, flags=re.S)
    body = re.sub(r'//.*', '', body)
    toks = [t.strip() for t in body.split(',') if t.strip()]
    vals = [hex_to_double(t) for t in toks]
    assert len(vals) == count, f"{name}: {len(vals)} != {count}"
    return vals

def extract_SC(src):
    pat = re.compile(r'static const double SC\[256\]\[3\] = \{(.*?)\n\};', re.S)
    m = pat.search(src)
    body = m.group(1)
    # Entries: {x, y, z}, /* idx */
    entry_pat = re.compile(r'\{([^}]*)\}')
    rows = []
    for em in entry_pat.finditer(body):
        parts = [p.strip() for p in em.group(1).split(',')]
        assert len(parts) == 3
        rows.append(tuple(hex_to_double(p) for p in parts))
    assert len(rows) == 256, len(rows)
    return rows

def emit_dint_table(name, entries):
    lines = [f'  {name}: array[0..{len(entries)-1}] of TDInt64 = (']
    for i, (hi, lo, ex, sgn) in enumerate(entries):
        sep = ',' if i < len(entries)-1 else ');'
        lines.append(f'    (hi:${hi:016X}; lo:${lo:016X}; ex:{ex}; sgn:{sgn}){sep}')
    return '\n'.join(lines)

def emit_u64_table(name, vals):
    lines = [f'  {name}: array[0..{len(vals)-1}] of UInt64 = (']
    for i, v in enumerate(vals):
        sep = ',' if i < len(vals)-1 else ');'
        lines.append(f'    UInt64(${v:016X}){sep}')
    return '\n'.join(lines)

def emit_double_table(name, vals):
    lines = [f'  {name}: array[0..{len(vals)-1}] of Tb64u64 = (']
    for i, d in enumerate(vals):
        sep = ',' if i < len(vals)-1 else ');'
        lines.append(f'    {as_u64(d)}{sep}')
    return '\n'.join(lines)

def emit_SC_table(name, rows):
    lines = [f'  {name}: array[0..{len(rows)-1}, 0..2] of Tb64u64 = (']
    for i, (a,b,c) in enumerate(rows):
        sep = ',' if i < len(rows)-1 else ');'
        lines.append(f'    ({as_u64(a)}, {as_u64(b)}, {as_u64(c)}){sep}')
    return '\n'.join(lines)

def main():
    src = read()
    T    = extract_T(src)
    S    = extract_dint_array(src, 'S', 256)
    C    = extract_dint_array(src, 'C', 256)
    PS   = extract_dint_array(src, 'PS', 6)
    PC   = extract_dint_array(src, 'PC', 6)
    PSf  = extract_doubles_array(src, 'PSfast', 5)
    PCf  = extract_doubles_array(src, 'PCfast', 5)
    SC   = extract_SC(src)

    print('// Auto-generated from core-math/src/binary64/cos/cos.c — do not edit by hand.')
    print('const')
    print(emit_u64_table('cSincosT', T))
    print()
    print(emit_dint_table('cSincosS', S))
    print()
    print(emit_dint_table('cSincosC', C))
    print()
    print(emit_dint_table('cCosPS', PS))
    print()
    print(emit_dint_table('cCosPC', PC))
    print()
    print(emit_double_table('cCosPSfast', PSf))
    print()
    print(emit_double_table('cCosPCfast', PCf))
    print()
    print(emit_SC_table('cCosSC', SC))

if __name__ == '__main__':
    main()
