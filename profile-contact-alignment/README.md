# RESM-Search

This project has been organized into a standard source layout:

```text
apps/       Command-line entry points
src/        Core implementation files
include/    Public project headers
scripts/    Build, benchmark, and verification scripts
docs/       Original README and optimization notes
data/       Small generated/sample data files
legacy/     Old BPmap sources not wired into the main CMake build
artifacts/  Legacy binaries and generated tag files
bin/        CMake runtime output directory
```

Build from the project root:

```bash
cmake . -DCMAKE_BUILD_TYPE=Release
cmake --build . -j 8
```

Useful scripts:

```bash
scripts/build_and_test.sh
scripts/benchmark_performance.sh
python scripts/verify_optimizations.py
```

Align two precomputed contact maps:

```bash
./bin/BPmap -a query.map -b target.map
```

Build BPmap with Make and RNAlib2 from pkg-config:

```bash
make BPmap
make BPmap_test
```

If `pkg-config` cannot find RNAlib2, point it at ViennaRNA's pkgconfig
directory:

```bash
export PKG_CONFIG_PATH=/path/to/ViennaRNA/lib/pkgconfig:$PKG_CONFIG_PATH
make BPmap
```

Or override the flags directly:

```bash
make BPmap RNALIB_CFLAGS="-I/path/to/ViennaRNA/include" RNALIB_LIBS="-L/path/to/ViennaRNA/lib -lRNA"
```

Contact map files are read by `Map::load_data` and should contain `LEN`/`SIZE`
and `CON i j [score]` lines. Optional `PRF` lines provide sequence/profile
information; when they are absent, BPmap uses `N` placeholders for alignment
printing.
