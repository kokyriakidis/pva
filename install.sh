#!/usr/bin/env bash
# install.sh
# Sets up all tools and libraries needed for the pangenome VNTR alignment pipeline.
#
# What this installs:
#   bin/vg            - variation graph toolkit (prebuilt binary)
#   bin/gbwt          - GBWT library (already cloned, built here)
#   bin/sdsl-lite     - SDSL library (cloned and built here, required by gbwt)
#   bin/gbwtgraph     - GBWTGraph library (cloned and built here)
#   bin/gbz-base      - GBZ database tools (already cloned, built here via Cargo)
#   bin/filter_paths  - path filtering tool (built from filter_paths.cpp)
#   bin/trace_haplotypes - haplotype tracing tool (built from trace_haplotypes.cpp)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
BIN="$ROOT/bin"

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

# ── 6. filter_paths ──────────────────────────────────────────────────────────
echo ""
echo "--- filter_paths ---"
if [ -x "$BIN/filter_paths" ] && \
   [ "$ROOT/filter_paths.cpp" -ot "$BIN/filter_paths" ]; then
    echo "OK: filter_paths up to date"
else
    echo "Building filter_paths..."
    g++ -O3 -std=c++17 -o "$BIN/filter_paths" "$ROOT/filter_paths.cpp"
    echo "OK: filter_paths built"
fi

# ── 7. trace_haplotypes ───────────────────────────────────────────────────────
echo ""
echo "--- trace_haplotypes ---"
if [ -x "$BIN/trace_haplotypes" ] && \
   [ "$ROOT/trace_haplotypes.cpp" -ot "$BIN/trace_haplotypes" ]; then
    echo "OK: trace_haplotypes up to date"
else
    echo "Building trace_haplotypes..."
    g++ -O2 -std=c++17 -fopenmp \
        -I"$BIN/gbwtgraph/include" \
        -I"$BIN/gbwt/include" \
        -I"$BIN/sdsl-lite/include" \
        -o "$BIN/trace_haplotypes" \
        "$ROOT/trace_haplotypes.cpp" \
        -L"$BIN/gbwtgraph/lib" \
        -L"$BIN/gbwt/lib" \
        -L"$BIN/sdsl-lite/lib" \
        "$BIN/gbwtgraph/lib/libgbwtgraph.a" \
        "$BIN/sdsl-lite/lib/libgbwt.a" \
        "$BIN/sdsl-lite/lib/libhandlegraph.a" \
        "$BIN/sdsl-lite/lib/libsdsl.a" \
        -lcrypto -lzstd -lz -ldivsufsort -ldivsufsort64 -fopenmp
    echo "OK: trace_haplotypes built"
fi

# ── 8. centrolign ─────────────────────────────────────────────────────────────
echo ""
echo "--- centrolign ---"
if [ -x "$BIN/centrolign" ]; then
    echo "OK: centrolign already present at bin/centrolign"
else
    echo "Building centrolign from source..."
    if [ ! -d "$BIN/centrolign-src" ]; then
        git clone --recurse-submodules \
            https://github.com/jeizenga/centrolign "$BIN/centrolign-src"
    fi
    cmake -S "$BIN/centrolign-src" \
          -B "$BIN/centrolign-src/build" \
          -DCMAKE_BUILD_TYPE=Release
    cmake --build "$BIN/centrolign-src/build" --parallel "$(nproc)"
    cp "$BIN/centrolign-src/build/centrolign" "$BIN/centrolign"
    echo "OK: centrolign built"
fi

# ── 9. Python dependencies ────────────────────────────────────────────────────
echo ""
echo "--- Python dependencies (scikit-bio) ---"
if python3 -c "import skbio" &>/dev/null 2>&1; then
    echo "OK: scikit-bio already installed"
else
    echo "Installing scikit-bio..."
    pip install --user scikit-bio
    echo "OK: scikit-bio installed"
fi

# ── done ──────────────────────────────────────────────────────────────────────
echo ""
echo "=== All done ==="
echo ""
echo "Binaries available in bin/:"
echo "  bin/vg"
echo "  bin/gbz-base/target/release/gbz2db"
echo "  bin/gbz-base/target/release/query"
echo "  bin/filter_paths"
echo "  bin/trace_haplotypes"
echo "  bin/centrolign"
