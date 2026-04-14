#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SRC_DIR="$SCRIPT_DIR/.."
ROOT_DIR="$SRC_DIR/.."
BIN_DIR="$ROOT_DIR/bin"
CMATH_SRC="$ROOT_DIR/../core-math/src/binary32"

mkdir -p "$BIN_DIR"

# ---- Step 1: Build libcoremath.so if not present ----
SOFILE="$SRC_DIR/libcoremath.so"
if [ ! -f "$SOFILE" ]; then
  echo "Building libcoremath.so ..."

  FUNCS=(
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
  for src in "${FUNCS[@]}"; do
    name="$(basename "$src" .c)"
    out="/tmp/cr_${name}.o"
    echo -n "  Compiling $src ... "
    if gcc -O2 -fPIC -c "$CMATH_SRC/$src" -o "$out" 2>/tmp/err_${name}.txt; then
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
  echo "libcoremath.so already present, skipping C build."
fi

FPC_FLAGS="-O3 -Fi.. -Fu.. -FE$BIN_DIR -Fl$SRC_DIR $@"

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

# ---- Clean compiled Pascal artifacts ----
find "$SRC_DIR" -maxdepth 1 \( -name '*.ppu' -o -name '*.o' -o -name '*.compiled' \) -delete
find "$BIN_DIR" -maxdepth 1 \( -name '*.ppu' -o -name '*.o' -o -name '*.compiled' \) -delete

echo
echo "Build complete."
echo
echo "Run tests with:"
echo "  LD_LIBRARY_PATH=$SRC_DIR $BIN_DIR/TestHarness32"
echo "  LD_LIBRARY_PATH=$SRC_DIR $BIN_DIR/Benchmark32"
echo "  LD_LIBRARY_PATH=$SRC_DIR $BIN_DIR/BenchmarkFPC32"
echo "  LD_LIBRARY_PATH=$SRC_DIR $BIN_DIR/FixedTest32"
