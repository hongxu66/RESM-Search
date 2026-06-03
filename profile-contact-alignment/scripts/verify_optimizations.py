#!/usr/bin/env python
"""Lightweight checks for the RESM-Search alignment optimizations."""

from __future__ import print_function
import os
import re
import sys


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def check_file(path, checks):
    with open(path, "rb") as handle:
        text = handle.read().decode("utf-8", "ignore")
    ok = True
    print("\nChecking {}".format(os.path.basename(path)))
    for name, pattern in checks.items():
        found = re.search(pattern, text, re.MULTILINE) is not None
        print("  {} {}".format("OK " if found else "ERR", name))
        ok = ok and found
    return ok


def main():
    checks = [
        check_file(
            os.path.join(ROOT, "src", "align.cpp"),
            {
                "parallel parameter sweep": r"#pragma\s+omp\s+parallel\s+for\s+collapse\(3\)",
                "runtime thread count": r"num_threads\(threads\)",
                "nested OpenMP guard": r"omp_in_parallel\(\)",
                "thread-local DP buffers": r"static\s+thread_local\s+vector<",
                "flat alignment DP": r"size_t\s+cells\s*=\s*\(size_t\)\(rows\s*\+\s*1\)\s*\*\s*stride",
                "one-row Falign DP": r"double\s+diag\s*=\s*0",
                "single-task outer loop guard": r"#pragma\s+omp\s+parallel\s+for\s+if\(",
                "PAM profile index fix": r"PAMatrix\(k1,k2\)",
            },
        ),
        check_file(
            os.path.join(ROOT, "src", "load_data.cpp"),
            {
                "safe option value parsing": r"require_value",
                "thread option": r"-threads",
                "window validation": r"overlap\s*>=\s*lengthE",
            },
        ),
        check_file(
            os.path.join(ROOT, "CMakeLists.txt"),
            {
                "Release O3": r"-O3",
                "OpenMP": r"OpenMP",
                "LTO": r"-flto",
                "BPmap target": r"add_executable\(BPmap\s+apps/BPmap\.cpp\)",
                "pkg-config link flags": r"RNALIB_LDFLAGS",
            },
        ),
        check_file(
            os.path.join(ROOT, "Makefile"),
            {
                "RNAlib2 pkg-config": r"RNALIB_PKG\s*\?=\s*RNAlib2",
                "manual RNAlib override": r"RNALIB_CFLAGS",
                "BPmap_test target": r"BPmap_test:",
                "new BPmap target": r"BPmap:",
            },
        ),
        check_file(
            os.path.join(ROOT, "apps", "BPmap.cpp"),
            {
                "contact map loader": r"load_data\(OPT\.filea",
                "alignment call": r"Ali\.mapalignment\(gap_b,\s*gap_a,\s*P_SCO\)",
                "empty sequence fallback": r"seq\.assign\(map\.mtx\.size\(\),\s*'N'\)",
            },
        ),
    ]

    print("\nSummary")
    if all(checks):
        print("All optimization checks passed.")
        return 0
    print("Some checks failed.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
