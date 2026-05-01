#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SRC_DIR="$SCRIPT_DIR/.."
ROOT_DIR="$SRC_DIR/.."
BIN_DIR="$ROOT_DIR/bin"
CMATH_SRC32="$ROOT_DIR/../core-math/src/binary32"
CMATH_SRC64="$ROOT_DIR/../core-math/src/binary64"

mkdir -p "$BIN_DIR"

# ---- Step 1a: Build libcoremath.so (binary32) if not present ----
SOFILE="$SRC_DIR/libcoremath.so"
if [ ! -f "$SOFILE" ]; then
  echo "Building libcoremath.so (binary32) ..."

  FUNCS32=(
    "acos/acosf.c"      "acosh/acoshf.c"    "acospi/acospif.c"
    "asin/asinf.c"      "asinh/asinhf.c"    "asinpi/asinpif.c"
    "atan/atanf.c"      "atan2/atan2f.c"    "atan2pi/atan2pif.c"
    "atanh/atanhf.c"    "atanpi/atanpif.c"  "cbrt/cbrtf.c"
    "compound/compoundf.c"
    "cos/cosf.c"        "cosh/coshf.c"      "cospi/cospif.c"
    "erf/erff.c"        "erfc/erfcf.c"
    "exp/expf.c"        "exp10/exp10f.c"    "exp10m1/exp10m1f.c"
    "exp2/exp2f.c"      "exp2m1/exp2m1f.c"  "expm1/expm1f.c"
    "hypot/hypotf.c"    "lgamma/lgammaf.c"
    "log/logf.c"        "log10/log10f.c"    "log10p1/log10p1f.c"
    "log1p/log1pf.c"    "log2/log2f.c"      "log2p1/log2p1f.c"
    "pow/powf.c"        "rsqrt/rsqrtf.c"
    "sin/sinf.c"        "sincos/sincosf.c"  "sinh/sinhf.c"
    "sinpi/sinpif.c"    "tan/tanf.c"        "tanh/tanhf.c"
    "tanpi/tanpif.c"    "tgamma/tgammaf.c"
  )

  OBJS=()
  for src in "${FUNCS32[@]}"; do
    name="$(basename "$src" .c)"
    out="/tmp/cr_${name}.o"
    echo -n "  Compiling $src ... "
    if gcc -O2 -fPIC -c "$CMATH_SRC32/$src" -o "$out" 2>/tmp/err_${name}.txt; then
      echo "OK"
      OBJS+=("$out")
    else
      cat "/tmp/err_${name}.txt"
      echo "  SKIPPED $src"
    fi
  done

  echo "  Linking ${#OBJS[@]} objects -> $SOFILE"
  gcc -shared -O2 -o "$SOFILE" "${OBJS[@]}" -lm
  echo "  libcoremath.so built OK"
else
  echo "libcoremath.so already present, skipping binary32 C build."
fi

# ---- Step 1b: Build libcoremath64.so (binary64) if not present ----
# Only attempt if the binary64 source directory exists.
SOFILE64="$SRC_DIR/libcoremath64.so"
if [ -d "$CMATH_SRC64" ] && [ ! -f "$SOFILE64" ]; then
  echo "Building libcoremath64.so (binary64) ..."

  FUNCS64=(
    "acos/acos.c"       "acosh/acosh.c"     "acospi/acospi.c"
    "asin/asin.c"       "asinh/asinh.c"     "asinpi/asinpi.c"
    "atan/atan.c"       "atan2/atan2.c"     "atan2pi/atan2pi.c"
    "atanh/atanh.c"     "atanpi/atanpi.c"   "cbrt/cbrt.c"
    "cos/cos.c"         "cosh/cosh.c"       "cospi/cospi.c"
    "erf/erf.c"         "erfc/erfc.c"
    "exp/exp.c"         "exp10/exp10.c"     "exp10m1/exp10m1.c"
    "exp2/exp2.c"       "exp2m1/exp2m1.c"   "expm1/expm1.c"
    "hypot/hypot.c"     "lgamma/lgamma.c"
    "log/log.c"         "log10/log10.c"     "log10p1/log10p1.c"
    "log1p/log1p.c"     "log2/log2.c"       "log2p1/log2p1.c"
    "pow/pow.c"         "rsqrt/rsqrt.c"
    "sin/sin.c"         "sincos/sincos.c"   "sinh/sinh.c"
    "sinpi/sinpi.c"     "tan/tan.c"         "tanh/tanh.c"
    "tanpi/tanpi.c"     "tgamma/tgamma.c"
  )

  OBJS64=()
  for src in "${FUNCS64[@]}"; do
    name="$(basename "$src" .c)"
    out="/tmp/cr64_${name}.o"
    echo -n "  Compiling $src ... "
    if gcc -O2 -fPIC -c "$CMATH_SRC64/$src" -o "$out" 2>/tmp/err64_${name}.txt; then
      echo "OK"
      OBJS64+=("$out")
    else
      cat "/tmp/err64_${name}.txt"
      echo "  SKIPPED $src"
    fi
  done

  if [ ${#OBJS64[@]} -gt 0 ]; then
    echo "  Linking ${#OBJS64[@]} objects -> $SOFILE64"
    # Note: FPC {$linklib coremath64} in ccoremath64.pas expects libcoremath64.so
    gcc -shared -O2 -o "$SOFILE64" "${OBJS64[@]}" -lm
    echo "  libcoremath64.so built OK"
  else
    echo "  No objects compiled; libcoremath64.so not built."
  fi
elif [ -f "$SOFILE64" ]; then
  echo "libcoremath64.so already present, skipping binary64 C build."
else
  echo "binary64 source not found at $CMATH_SRC64 — skipping libcoremath64.so."
fi

FPC_FLAGS="-O3 -Fi.. -Fu.. -FE$BIN_DIR -Fl$SRC_DIR $@"

# ---- Clean compiled Pascal artifacts ----
find "$SRC_DIR" -maxdepth 1 \( -name '*.ppu' -o -name '*.o' -o -name '*.compiled' -o -name '*.s' \) -delete
find "$BIN_DIR" -maxdepth 1 \( -name '*.ppu' -o -name '*.o' -o -name '*.compiled' -o -name '*.s' \) -delete

# ---- Step 2: Compile TestHarness32 ----
echo
echo "Compiling TestHarness32.pas ..."
fpc $FPC_FLAGS "$SCRIPT_DIR/TestHarness32.pas"
echo "TestHarness32 compiled -> $BIN_DIR/TestHarness32"

# ---- Step 3: Compile Benchmark32 ----
echo
echo "Compiling Benchmark32.pas ..."
fpc $FPC_FLAGS "$SCRIPT_DIR/Benchmark32.pas"
echo "Benchmark32 compiled -> $BIN_DIR/Benchmark32"

# ---- Step 4: Compile BenchmarkFPC32 ----
echo
echo "Compiling BenchmarkFPC32.pas ..."
fpc $FPC_FLAGS "$SCRIPT_DIR/BenchmarkFPC32.pas"
echo "BenchmarkFPC32 compiled -> $BIN_DIR/BenchmarkFPC32"

# ---- Step 5: Compile FixedTest32 ----
echo
echo "Compiling FixedTest32.pas ..."
fpc $FPC_FLAGS "$SCRIPT_DIR/FixedTest32.pas"
echo "FixedTest32 compiled -> $BIN_DIR/FixedTest32"

# ---- Step 6: Compile binary64 test programs (only if libcoremath64.so exists) ----
if [ -f "$SOFILE64" ]; then
  echo
  echo "Compiling TestHarness64.pas ..."
  fpc $FPC_FLAGS "$SCRIPT_DIR/TestHarness64.pas" && \
    echo "TestHarness64 compiled -> $BIN_DIR/TestHarness64" || \
    echo "TestHarness64 skipped (functions not yet ported)"

  echo
  echo "Compiling Benchmark64.pas ..."
  fpc $FPC_FLAGS "$SCRIPT_DIR/Benchmark64.pas" && \
    echo "Benchmark64 compiled -> $BIN_DIR/Benchmark64" || \
    echo "Benchmark64 skipped (functions not yet ported)"

  echo
  echo "Compiling BenchmarkFPC64.pas ..."
  fpc $FPC_FLAGS "$SCRIPT_DIR/BenchmarkFPC64.pas" && \
    echo "BenchmarkFPC64 compiled -> $BIN_DIR/BenchmarkFPC64" || \
    echo "BenchmarkFPC64 skipped (functions not yet ported)"

  echo
  echo "Compiling FixedTest64.pas ..."
  fpc $FPC_FLAGS "$SCRIPT_DIR/FixedTest64.pas" && \
    echo "FixedTest64 compiled -> $BIN_DIR/FixedTest64" || \
    echo "FixedTest64 skipped (functions not yet ported)"

  echo
  echo "Compiling TestQInt64.pas ..."
  fpc $FPC_FLAGS "$SCRIPT_DIR/TestQInt64.pas" && \
    echo "TestQInt64 compiled -> $BIN_DIR/TestQInt64" || \
    echo "TestQInt64 skipped"

  echo
  echo "Compiling FmaCompare64.pas ..."
  fpc $FPC_FLAGS "$SCRIPT_DIR/FmaCompare64.pas" && \
    echo "FmaCompare64 compiled -> $BIN_DIR/FmaCompare64" || \
    echo "FmaCompare64 skipped"
else
  echo
  echo "Skipping binary64 test programs (libcoremath64.so not present)."
  echo "Run 'bash build.sh' after porting functions to enable binary64 tests."
fi

echo
echo "Build complete."
echo
echo "Run binary32 tests with:"
echo "  LD_LIBRARY_PATH=$SRC_DIR $BIN_DIR/TestHarness32"
echo "  LD_LIBRARY_PATH=$SRC_DIR $BIN_DIR/Benchmark32"
echo "  LD_LIBRARY_PATH=$SRC_DIR $BIN_DIR/BenchmarkFPC32"
echo "  LD_LIBRARY_PATH=$SRC_DIR $BIN_DIR/FixedTest32"
echo
echo "Run binary64 tests with (after porting functions):"
echo "  LD_LIBRARY_PATH=$SRC_DIR $BIN_DIR/TestHarness64"
echo "  taskset -c 1 env LD_LIBRARY_PATH=$SRC_DIR $BIN_DIR/Benchmark64"
echo "  taskset -c 1 env LD_LIBRARY_PATH=$SRC_DIR $BIN_DIR/BenchmarkFPC64"
echo "  LD_LIBRARY_PATH=$SRC_DIR $BIN_DIR/FixedTest64"
