#!/usr/bin/env python3
"""
x87_audit.py — flag potential x87-FPU-promoting patterns in Pascal math source.

Phase 6 of tasklist64.md: ensure FPC's code generator stays on SSE2 by
(a) typing every numeric literal explicitly via Double(...) typecasts,
(b) replacing literal-indexed constant-array reads in hot polynomial
chains with named scalar consts, and (c) unrolling fixed-trip-count
polynomial loops over coefficient arrays so Pillar B becomes possible.

What this script flags
----------------------
Pillar A — Loop-indexed `cFoo[i].f` reads inside `for ... do ... end`
           bodies (loop variable -> Pillar A unroll candidate).
Pillar B — Literal-indexed `cFoo[3].f` reads (Pillar B lift candidate).
Pillar C — Untyped numeric literals (`1.0`, `0.5`, `-2.0e-3`, ...) that
           are NOT preceded by `Double(`, NOT inside a `const ... :
           Double = ...;` declaration, and NOT a comment / hex-float bit
           pattern.

Exit code is 0 when a file is clean, 1 otherwise — handy for CI / per-
file iteration. The audit is heuristic: it errs on the side of reporting
borderline cases. Manual triage is part of the per-file Phase 6 pass.

Usage
-----
    python3 tmp/x87_audit.py src/inc_64/atan_port_64.inc
    python3 tmp/x87_audit.py src/inc_64/*.inc
    python3 tmp/x87_audit.py --pillar c src/pascoremath64.pas
    python3 tmp/x87_audit.py --summary src/inc_64/*.inc
"""

import argparse
import re
import sys
from pathlib import Path

# Pillar C: a numeric literal containing a decimal point or scientific
# exponent. We match either explicitly to keep it readable.
LITERAL_RE = re.compile(
    r"""
    (?<![A-Za-z0-9_.])           # not in middle of identifier / fp suffix
    -?                           # optional sign
    (?:
        \d+\.\d+(?:[eE][+-]?\d+)?  # 1.0 / 1.5e-3
      | \d+[eE][+-]?\d+            # 1e10
    )
    """,
    re.VERBOSE,
)

# A literal is "typed" when wrapped by Double(...). Look back ~16 chars.
TYPED_PREFIX_RE = re.compile(r"Double\s*\(\s*-?\s*$")

# Pillar B / A: cFoo[index].f where index is either a literal or an ident.
ARRAY_READ_RE = re.compile(
    r"""
    \b(c[A-Za-z][A-Za-z0-9_]*)   # 1: array name (convention: starts with c)
    \[                           #
    \s*([^\]]+?)\s*              # 2: index expression
    \]
    \.f\b                        # .f field
    """,
    re.VERBOSE,
)

# Track whether we're inside a const block so we can suppress literal
# warnings on declarations like `cFoo: Double = 1.5;`.
CONST_RE = re.compile(r"^\s*const\b", re.IGNORECASE)
# Heuristic: a typed-decl line has `: Double = ` or `: Tb64u64 = `.
TYPED_DECL_RE = re.compile(r":\s*[A-Za-z_][A-Za-z0-9_]*\s*=", re.IGNORECASE)
FOR_RE = re.compile(r"^\s*for\s+([A-Za-z_]\w*)\s*:=", re.IGNORECASE)
END_RE = re.compile(r"^\s*end\s*[;.]?\s*$", re.IGNORECASE)


def strip_comments(line: str) -> str:
    """Remove `//`, `{...}`, `(*...*)` comments at the line level."""
    out = []
    i, n = 0, len(line)
    in_brace = in_paren = False
    while i < n:
        c = line[i]
        c2 = line[i:i+2]
        if not in_brace and not in_paren:
            if c2 == "//":
                break
            if c == "{":
                in_brace = True
                i += 1
                continue
            if c2 == "(*":
                in_paren = True
                i += 2
                continue
            out.append(c)
        else:
            if in_brace and c == "}":
                in_brace = False
            elif in_paren and c2 == "*)":
                in_paren = False
                i += 2
                continue
        i += 1
    return "".join(out)


def is_int_literal(text: str) -> bool:
    return bool(re.fullmatch(r"\d+", text))


def audit_file(path: Path, want: set) -> list:
    findings = []
    try:
        src = path.read_text(encoding="utf-8", errors="replace")
    except OSError as e:
        return [(0, "?", f"cannot read: {e}")]

    loop_vars: list = []  # stack of (var, indent) — heuristic; we only
                          # track for-loops at procedure-body depth.
    in_const = False

    for lineno, raw in enumerate(src.splitlines(), start=1):
        code = strip_comments(raw)
        stripped = code.strip()
        if not stripped:
            continue

        # Track const block — naive: const ... until next begin/var/...
        if CONST_RE.match(code):
            in_const = True
        elif re.match(r"^\s*(begin|var|type|function|procedure|implementation|interface|uses)\b", code, re.IGNORECASE):
            in_const = False

        # Track for-loop variable (heuristic — pop on bare `end`).
        m = FOR_RE.match(code)
        if m:
            loop_vars.append(m.group(1).lower())
        elif END_RE.match(code) and loop_vars:
            loop_vars.pop()

        # Pillar B/A — array reads.
        if "b" in want or "a" in want:
            for m in ARRAY_READ_RE.finditer(code):
                arr, idx = m.group(1), m.group(2).strip()
                if is_int_literal(idx):
                    if "b" in want:
                        findings.append((lineno, "B",
                            f"literal-indexed read {arr}[{idx}].f — lift to named scalar"))
                else:
                    idx_low = idx.lower()
                    if "a" in want and any(lv == idx_low for lv in loop_vars):
                        findings.append((lineno, "A",
                            f"loop-indexed read {arr}[{idx}].f inside for-loop — unroll candidate"))

        # Pillar C — untyped float literals.
        if "c" in want and not in_const:
            # Skip lines that look like typed const declarations even
            # outside a `const` block (e.g. record-style initialisers).
            if TYPED_DECL_RE.search(code):
                continue
            for m in LITERAL_RE.finditer(code):
                lit = m.group(0)
                start = m.start()
                # Look back: is there a `Double(` opener within ~12 chars?
                back = code[max(0, start - 16):start]
                if TYPED_PREFIX_RE.search(back):
                    continue
                # Negative literal in `-Double(...)` form is also fine —
                # but our regex already matched a `-` only when attached
                # to digits, so plain `-` operators don't trigger here.
                findings.append((lineno, "C",
                    f"untyped literal `{lit}` — wrap as Double({lit})"))

    return findings


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+", type=Path)
    ap.add_argument("--pillar", default="abc",
                    help="subset of 'abc' (default all three pillars)")
    ap.add_argument("--summary", action="store_true",
                    help="one-line per file: counts per pillar")
    args = ap.parse_args()

    want = set(args.pillar.lower())
    if not want <= {"a", "b", "c"}:
        ap.error("--pillar must be a subset of 'abc'")

    any_dirty = False
    for path in args.paths:
        findings = audit_file(path, want)
        if not findings:
            if not args.summary:
                print(f"{path}: clean")
            continue
        any_dirty = True
        if args.summary:
            counts = {"A": 0, "B": 0, "C": 0}
            for _, p, _ in findings:
                counts[p] = counts.get(p, 0) + 1
            print(f"{path}: A={counts['A']} B={counts['B']} C={counts['C']}")
        else:
            for lineno, pillar, msg in findings:
                print(f"{path}:{lineno}: [{pillar}] {msg}")

    return 1 if any_dirty else 0


if __name__ == "__main__":
    sys.exit(main())
