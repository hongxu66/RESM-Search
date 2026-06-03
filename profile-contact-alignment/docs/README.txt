RESM-Search
============

Version
-------
0.2.0 - 2026-03-23

Overview
--------
RESM-Search is an RNA sequence alignment and homologous sequence search toolkit.
The project provides three command-line programs:

1. `cmapbuild`: build a cmap model from RNA structure/sequence input.
2. `RESM_Search_align`: align RNA sequences against a query model.
3. `cmapsearch`: search homologous sequences from genome/database FASTA files.

Requirements
------------
- C++11 compiler
- CMake (>= 3.0)
- ViennaRNA / RNAlib2
- OpenMP (optional but recommended)

Build
-----
From project root:

1. `cmake .`
2. `cmake --build .`

Binaries are generated in `./bin`:
- `bin/cmapbuild`
- `bin/RESM_Search`
- `bin/RESM_Search_align`
- `bin/cmapsearch`
- `bin/RESM_Search_alignNCBI`
- `bin/BPmap`

Source Layout
-------------
- `apps/`: command-line program entry points
- `src/`: core implementation files
- `include/`: project headers
- `scripts/`: build, benchmark, and verification scripts
- `docs/`: project notes
- `legacy/`: older BPmap sources not wired into the main CMake build

Program Inputs and Outputs
--------------------------
`cmapbuild`
- Input:
  - Stockholm MSA file, or
  - DBN file + MSA file, or
  - DBN file only
- Output:
  - cmap model file

`RESM_Search_align`
- Input:
  - cmap model (or DBN/MSA-derived model)
  - target RNA sequence FASTA file
- Output:
  - alignment results (mapping + score)

`cmapsearch`
- Input:
  - query DBN file (optionally with MSA)
  - target genome/database FASTA file
- Output:
  - homologous sequence hits
  - pairwise alignment mapping
- score report

`BPmap`
- Input:
  - two precomputed contact map files containing `LEN`/`SIZE` and `CON` lines
  - optional `PRF` lines for sequence/profile-aware alignment output
- Output:
  - pairwise contact-map alignment and BestScore

Refactor Notes (2026-03-23)
---------------------------
- Project build is refactored with a shared core library target (`resm_core`) to reduce duplicated source lists.
- P-value calculation functions and data loading chain are removed.
- Alignment/search now report and filter with alignment score (`BestScore`) instead of P-value.
- Core alignment algorithm details (`ini_SCO`, `mod_SCO`, `chk`, dynamic programming flow) are unchanged.
