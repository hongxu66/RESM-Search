# Contact-Map Alignment Optimization Summary

## Implemented Optimizations

### Option 1: OpenMP Parallel Parameter Sweep

Status: completed

File: `src/align.cpp`

Method: use `#pragma omp parallel for collapse(3)` to parallelize the three-level nested parameter sweep.

- Parameter combinations: 3 x 6 x 4 = 72 alignment runs
- Scheduling strategy: `schedule(dynamic)` for load balancing
- Default thread count: 8
- Thread safety: `#pragma omp critical` protects best-result updates

Code locations:

- `int Alignment::mapalignment(...stringstream &buf)`
- `double Alignment::mapalignment(...stringstream &buf)`
- `vec_int Alignment::mapalignment(...)`

Estimated speedup: 40-50%, depending on CPU core count and workload.

### Option 2/3: Memory Allocation Optimization

Status: completed

File: `src/align.cpp`

Method: remove variable-length arrays and use preallocated thread-local buffers.

```cpp
// Before: VLA on stack
double M[A[k] * B[k]];

// After: preallocated thread-local buffer
static thread_local std::vector<double> M_buffer(65536);
```

Benefits:

- Avoids stack overflow risks from large VLAs.
- Reuses allocated memory and reduces allocation/deallocation overhead.
- Gives each thread an independent buffer and avoids data races.

Estimated speedup: 8-15%.

### Option 5: Compiler Optimization Flags

Status: completed

File: `CMakeLists.txt`

Release flags:

```cmake
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}
    -O3
    -march=native
    -ffast-math
    -flto
    -fno-math-errno
    -finline-functions
    -foptimize-sibling-calls
")
```

Benefits:

- `-O3`: enables aggressive compiler optimization.
- `-march=native`: uses CPU-specific instructions such as SIMD and AVX.
- `-ffast-math`: allows faster floating-point optimization.
- `-flto`: enables link-time optimization across translation units.
- `-finline-functions`: increases function inlining.

Estimated speedup: 5-15%.

### Option 4: Early-Termination Strategy

Status: optional future work

Early termination can be useful in selected workloads, but the current OpenMP parallelization already provides the main speedup. The expected additional gain is limited unless the search space is much larger or pruning quality improves.

## Performance Summary

| Optimization | Estimated speedup | Status |
| --- | --- | --- |
| OpenMP parameter sweep | 40-50% | completed |
| Memory optimization | 8-15% | completed |
| Compiler optimization | 5-15% | completed |
| Total | 60-80% | completed |

The expected end-to-end acceleration is roughly 1.6x to 5x, depending on CPU core count and memory-access behavior.

## Build Instructions

Build a release version with full optimization:

```bash
cd /path/to/RESM-Search/profile-contact-alignment
cmake . -DCMAKE_BUILD_TYPE=Release
make -j8
```

Verify optimization flags:

```bash
cmake . -DCMAKE_BUILD_TYPE=Release
make VERBOSE=1 | grep -E "march=native|flto|-O3"
```

Check for OpenMP symbols:

```bash
objdump -t bin/RESM_Search | grep omp
```

## Key Changes

### Parallel Three-Level Sweep

```cpp
#pragma omp parallel for collapse(3) schedule(dynamic) num_threads(8)
for (int sx = 0; sx < sep_x_steps.size(); sx++) {
    for (int sy = 0; sy < sep_y_steps.size(); sy++) {
        for (int g_e = 0; g_e < gap_e_steps.size(); g_e++) {
            #pragma omp critical(alignment_update)
            {
                if (current_score > con_max + gap_max + prf_max) {
                    // Thread-safe update of the best result.
                }
            }
        }
    }
}
```

### Thread-Local Buffer

```cpp
static thread_local std::vector<double> M_buffer(65536);
```

Each thread owns its buffer and reuses it across calls.

### Compiler Flags

```cmake
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -ffast-math -flto ...")
```

## Modified Files

- `src/align.cpp`: parallelized map-alignment functions and optimized `ini_SCO`.
- `CMakeLists.txt`: added release optimization flags.
- `include/align.h`: confirmed inline-function optimization.

## Future Improvements

1. Profile-guided optimization with PGO.
2. Manual SIMD vectorization for matrix operations.
3. CUDA or OpenCL acceleration for compute-heavy kernels.
4. Improved parameter-search strategy to reduce unnecessary runs.
5. Better data-access patterns for higher cache hit rates.

## Performance Test

```bash
cmake . -DCMAKE_BUILD_TYPE=Release
make -j8
time ./bin/RESM_Search_align -a query.dbn -b target.fa

cmake . -DCMAKE_BUILD_TYPE=Debug
make -j8
time ./bin/RESM_Search_align -a query.dbn -b target.fa
```

Optimization completion date: 2026-05-11
