#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "$PROJECT_ROOT"

module load gcc/9.2.0 || true
module load eigen/eigen-3.4.0-gcc-9.2.0 || true

g++ -fopenmp -O3 -std=c++0x \
    -Iinclude \
    src/load_data.cpp src/align.cpp src/buildmap.cpp src/msa.cpp apps/BPmap.cpp \
    $(pkg-config --cflags RNAlib2) $(pkg-config --libs RNAlib2) \
    -o bin/BPmap
