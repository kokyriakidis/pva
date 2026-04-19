#!/usr/bin/env bash
# install.sh
# Installs/builds external dependencies, builds pva, and prepares tests for the
# pangenome VNTR alignment pipeline.
#
# What this script installs/builds:
#   bin/vg            - variation graph toolkit (prebuilt binary)
#   bin/gbwt          - GBWT library (already cloned, built here)
#   bin/sdsl-lite     - SDSL library (cloned and built here, required by gbwt)
#   bin/gbwtgraph     - GBWTGraph library (cloned and built here)
#   bin/gbz-base      - GBZ database tools (already cloned, built here via Cargo)
#   bin/pva + bin/libpva.so - main CLI and shared library
#   tests/unit/test_unit    - bundled unit test binary

set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
BIN="$ROOT/bin"
mkdir -p "$BIN"

echo "=== Pipeline installation ==="
echo "Root: $ROOT"

# ── 0. system dependencies ────────────────────────────────────────────────────
echo ""
echo "--- Checking system dependencies ---"
echo "    (if missing, run: sudo apt install -y build-essential libomp-dev cmake git curl)"
for cmd in git cmake cargo g++; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: '$cmd' not found. Install it before running this script."
        exit 1
    fi
done
# Check for libomp (required by gbwtgraph's -fopenmp)
if ! dpkg -s libomp-dev &>/dev/null 2>&1; then
    echo "WARNING: libomp-dev not found. gbwtgraph may fail to link."
    echo "         Run: sudo apt install -y libomp-dev"
fi
echo "OK: system dependencies present"

# ── 1. vg (prebuilt binary) ───────────────────────────────────────────────────
echo ""
echo "--- vg ---"
if [ -x "$BIN/vg" ]; then
    echo "OK: vg already present at bin/vg"
else
    echo "Downloading vg binary..."
    VG_VERSION="v1.70.0"
    curl -fsSL -o "$BIN/vg" \
        "https://github.com/vgteam/vg/releases/download/${VG_VERSION}/vg"
    chmod +x "$BIN/vg"
    echo "OK: vg ${VG_VERSION} downloaded"
fi

# ── 2. sdsl-lite (vgteam fork) ────────────────────────────────────────────────
echo ""
echo "--- sdsl-lite ---"
if [ -f "$BIN/sdsl-lite/lib/libsdsl.a" ]; then
    echo "OK: sdsl-lite already built"
else
    if [ ! -d "$BIN/sdsl-lite" ]; then
        echo "Cloning sdsl-lite..."
        git clone https://github.com/vgteam/sdsl-lite "$BIN/sdsl-lite"
    else
        echo "sdsl-lite source already present, rebuilding..."
    fi
    echo "Building sdsl-lite..."
    cmake -S "$BIN/sdsl-lite" \
          -B "$BIN/sdsl-lite/build" \
          -DCMAKE_INSTALL_PREFIX="$BIN/sdsl-lite" \
          -DCMAKE_BUILD_TYPE=Release
    cmake --build "$BIN/sdsl-lite/build" --parallel "$(nproc)"
    cmake --install "$BIN/sdsl-lite/build"
    echo "OK: sdsl-lite built"
fi

# ── 3. libhandlegraph ────────────────────────────────────────────────────────
echo ""
echo "--- libhandlegraph ---"
if [ -f "$BIN/sdsl-lite/lib/libhandlegraph.a" ]; then
    echo "OK: libhandlegraph already built and installed"
else
    if [ ! -d "$BIN/handlegraph" ]; then
        echo "Cloning libhandlegraph..."
        git clone https://github.com/vgteam/libhandlegraph "$BIN/handlegraph"
    fi
    echo "Building libhandlegraph..."
    cmake -S "$BIN/handlegraph" \
          -B "$BIN/handlegraph/build" \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX="$BIN/sdsl-lite"
    cmake --build "$BIN/handlegraph/build" --parallel "$(nproc)"
    cmake --install "$BIN/handlegraph/build"
    echo "OK: libhandlegraph built and installed into sdsl-lite prefix"
fi

# ── 4. gbwt ──────────────────────────────────────────────────────────────────
echo ""
echo "--- gbwt ---"
if [ ! -d "$BIN/gbwt" ]; then
    echo "Cloning gbwt..."
    git clone https://github.com/jltsiren/gbwt "$BIN/gbwt"
fi
if [ -f "$BIN/sdsl-lite/lib/libgbwt.a" ]; then
    echo "OK: gbwt already built and installed"
else
    echo "Building gbwt..."
    # gbwt Makefile expects sdsl-lite at ../sdsl-lite relative to itself
    make -C "$BIN/gbwt" SDSL_DIR="$BIN/sdsl-lite" -j "$(nproc)"
    # Install gbwt headers and library into sdsl-lite prefix so gbwtgraph
    # finds everything under INC_DIR/LIB_DIR (set by sdsl-lite/Make.helper)
    echo "Installing gbwt into sdsl-lite prefix..."
    cp -r "$BIN/gbwt/include/gbwt" "$BIN/sdsl-lite/include/"
    cp "$BIN/gbwt/lib/libgbwt.a"   "$BIN/sdsl-lite/lib/"
    echo "OK: gbwt built and installed"
fi

# ── 4. gbwtgraph ──────────────────────────────────────────────────────────────
echo ""
echo "--- gbwtgraph ---"
if [ ! -d "$BIN/gbwtgraph" ]; then
    echo "Cloning gbwtgraph..."
    git clone https://github.com/jltsiren/gbwtgraph "$BIN/gbwtgraph"
fi
if [ -f "$BIN/gbwtgraph/lib/libgbwtgraph.a" ]; then
    echo "OK: gbwtgraph already built"
else
    echo "Building gbwtgraph..."
    # gbwtgraph picks up gbwt + sdsl via INC_DIR/LIB_DIR from Make.helper
    make -C "$BIN/gbwtgraph" SDSL_DIR="$BIN/sdsl-lite" -j "$(nproc)"
    echo "OK: gbwtgraph built"
fi

# ── 5. gbz-base (Rust) ───────────────────────────────────────────────────────
echo ""
echo "--- gbz-base ---"
if [ -x "$BIN/gbz-base/target/release/gbz2db" ] && \
   [ -x "$BIN/gbz-base/target/release/query" ]; then
    echo "OK: gbz-base already built"
else
    if [ ! -d "$BIN/gbz-base" ]; then
        echo "Cloning gbz-base..."
        git clone https://github.com/jltsiren/gbz-base "$BIN/gbz-base"
    fi
    echo "Building gbz-base..."
    cargo build --release --manifest-path "$BIN/gbz-base/Cargo.toml"
    echo "OK: gbz-base built"
fi

# ── 6. centrolign ─────────────────────────────────────────────────────────────
echo ""
echo "--- centrolign ---"
CENTROLIGN_BUILD="$BIN/centrolign-src/build"
CENTROLIGN_INC="$BIN/centrolign-src/include"
LIBCENTROLIGN="$CENTROLIGN_BUILD/libcentrolign.so"

if [ -x "$BIN/centrolign" ] && [ -f "$LIBCENTROLIGN" ]; then
    echo "OK: centrolign already present at bin/centrolign"
else
    echo "Building centrolign from source..."
    if [ ! -d "$BIN/centrolign-src" ]; then
        git clone --recurse-submodules \
            https://github.com/jeizenga/centrolign "$BIN/centrolign-src"
    fi
    cmake -S "$BIN/centrolign-src" \
          -B "$CENTROLIGN_BUILD" \
          -DCMAKE_BUILD_TYPE=Release \
          -DBUILD_SHARED_LIBS=ON
    cmake --build "$CENTROLIGN_BUILD" --parallel "$(nproc)"
    cp "$CENTROLIGN_BUILD/centrolign" "$BIN/centrolign"
    echo "OK: centrolign built (binary + libcentrolign.so)"
fi

# ── 7. Python dependencies ────────────────────────────────────────────────────
echo ""
echo "--- Python dependencies (scikit-bio) ---"
if python3 -c "import skbio" &>/dev/null 2>&1; then
    echo "OK: scikit-bio already installed"
else
    echo "Installing scikit-bio..."
    pip install --break-system-packages scikit-bio
    echo "OK: scikit-bio installed"
fi

# ── 8. pva (libpva.so + pva binary) ─────────────────────────────────────────
echo ""
echo "--- pva ---"

# Sources compiled into libpva.so
# pva_align.cpp requires centrolign headers + libcentrolign.so
LIBPVA_SOURCES="$ROOT/pva_query.cpp $ROOT/pva_snarls.cpp $ROOT/pva_boundary.cpp $ROOT/pva_prep_guide_tree.cpp $ROOT/pva_snarls_dist_filter.cpp $ROOT/pva_infer_guide_tree.cpp $ROOT/pva_align.cpp"

# pva_trace_haplotypes.cpp links against gbwtgraph static libs — compiled directly into binary
PVA_HEADER="$ROOT/pva.h"
LIBPVA="$BIN/libpva.so"
PVA_BIN="$BIN/pva"

# Rebuild if any source or header is newer than the target
needs_rebuild() {
    local target=$1; shift
    [ ! -f "$target" ] && return 0
    for src in "$@"; do
        [ "$src" -nt "$target" ] && return 0
    done
    return 1
}

if needs_rebuild "$LIBPVA" $LIBPVA_SOURCES "$PVA_HEADER" "$ROOT/pva_utils.h" "$LIBCENTROLIGN"; then
    echo "Building libpva.so..."
    g++ -O2 -std=c++17 -fPIC -shared \
        -I"$CENTROLIGN_INC" \
        -o "$LIBPVA" \
        $LIBPVA_SOURCES \
        -L"$CENTROLIGN_BUILD" -lcentrolign \
        -Wl,-rpath,"$CENTROLIGN_BUILD"
    echo "OK: libpva.so built"
else
    echo "OK: libpva.so up to date"
fi

if needs_rebuild "$PVA_BIN" \
       "$ROOT/pva_main.cpp" "$ROOT/pva_trace_haplotypes.cpp" "$PVA_HEADER" "$ROOT/pva_utils.h" "$LIBPVA"; then
    echo "Building pva binary..."
    g++ -O2 -std=c++17 -fopenmp \
        -I"$BIN/gbwtgraph/include" \
        -I"$BIN/gbwt/include" \
        -I"$BIN/sdsl-lite/include" \
        -o "$PVA_BIN" \
        "$ROOT/pva_main.cpp" \
        "$ROOT/pva_trace_haplotypes.cpp" \
        -L"$BIN" -lpva \
        -Wl,-rpath,"$BIN" \
        -Wl,-rpath,"$CENTROLIGN_BUILD" \
        "$BIN/gbwtgraph/lib/libgbwtgraph.a" \
        "$BIN/sdsl-lite/lib/libgbwt.a" \
        "$BIN/sdsl-lite/lib/libhandlegraph.a" \
        "$BIN/sdsl-lite/lib/libsdsl.a" \
        -lcrypto -lzstd -lz -ldivsufsort -ldivsufsort64 -fopenmp
    echo "OK: pva built"
else
    echo "OK: pva up to date"
fi

# ── 9. tests ─────────────────────────────────────────────────────────────────
echo ""
echo "--- tests ---"
TEST_UNIT="$ROOT/tests/unit/test_unit"
if needs_rebuild "$TEST_UNIT" "$ROOT/tests/unit/test_unit.cpp" "$ROOT/pva_utils.h"; then
    echo "Building unit tests..."
    g++ -O0 -g -std=c++17 \
        -I"$ROOT" \
        -o "$TEST_UNIT" \
        "$ROOT/tests/unit/test_unit.cpp"
    echo "OK: unit tests built"
else
    echo "OK: unit tests up to date"
fi

echo "Running unit tests..."
"$TEST_UNIT"

chmod +x "$ROOT/tests/test_pva.sh"
echo "OK: integration test script ready (run: tests/test_pva.sh)"

# ── done ──────────────────────────────────────────────────────────────────────
echo ""
echo "=== All done ==="
echo ""
echo "Binaries available in bin/:"
echo "  bin/vg"
echo "  bin/gbz-base/target/release/gbz2db"
echo "  bin/gbz-base/target/release/query"
echo "  bin/centrolign"
echo "  bin/pva  (+ bin/libpva.so)"
echo ""
echo "Tests:"
echo "  tests/unit/test_unit       (unit tests — run automatically above)"
echo "  tests/test_pva.sh          (integration tests — run manually)"
echo "  TEST_DATA_DIR=<path> tests/test_pva.sh  (with real graph data)"
