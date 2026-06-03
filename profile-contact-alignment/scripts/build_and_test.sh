#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
THREADS="${THREADS:-8}"

cd "$PROJECT_ROOT"

echo "=== RESM-Search Release Build ==="
echo "Threads: ${THREADS}"

rm -rf CMakeCache.txt CMakeFiles/
cmake . -DCMAKE_BUILD_TYPE=Release
cmake --build . -j "${THREADS}"

echo
echo "Build successful. Executables:"
ls -lh bin/
echo
echo "BPmap usage:"
echo "  ./bin/BPmap -a query.map -b target.map"
