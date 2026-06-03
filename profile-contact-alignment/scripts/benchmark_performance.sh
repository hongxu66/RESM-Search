#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
THREADS="${THREADS:-8}"

cd "$PROJECT_ROOT"

echo "=== Alignment Algorithm Performance Benchmark ==="
echo
echo "Build: Release"
echo "Threads: ${THREADS}"
echo

cmake . -DCMAKE_BUILD_TYPE=Release
cmake --build . -j "${THREADS}"

echo
echo "Run a benchmark with your data, for example:"
echo "  /usr/bin/time -v ./bin/RESM_Search_align -a query.dbn -b target.fa -threads ${THREADS}"
echo
echo "For a baseline comparison:"
echo "  cmake . -DCMAKE_BUILD_TYPE=Debug"
echo "  cmake --build . -j ${THREADS}"
echo "  /usr/bin/time -v ./bin/RESM_Search_align -a query.dbn -b target.fa -threads ${THREADS}"
