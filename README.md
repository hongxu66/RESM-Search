# RESM-Search

RESM-Search is an RNA homologous sequence search project. It combines RNA language-model inference, contact-map/profile generation, and contact-map alignment/search algorithms. The database storage layer is planned, but the database files are not included in this repository yet.

The typical workflow is:

1. Use `RNA-ESM2` to infer RNA profiles and contact maps from FASTA input.
2. Use `profile-contact-alignment` to build or align contact maps.
3. Rank candidate homologous sequences by alignment score or `BestScore`.

## Repository Layout

```text
RESM-Search/
+-- RNA-ESM2/                    # RNA-ESM2 model, inference scripts, training code, sample FASTA files
|   +-- rna_esm/                 # Local RNA-ESM/ESM model code
|   +-- evo/                     # MSA, tokenization, metrics, and helper modules
|   +-- fasta/                   # Example and test sequences
|   +-- ckpt/                    # Local checkpoints, ignored by Git
|   +-- data/                    # Local pretrained weights and training data, ignored by Git
+-- profile-contact-alignment/   # C++ contact-map alignment and search tools
    +-- apps/                    # Command-line entry points
    +-- src/                     # Core algorithm implementation
    +-- include/                 # Header files
    +-- scripts/                 # Build, verification, and benchmark scripts
    +-- docs/                    # Subproject notes
```

## Requirements

Python side:

- Python 3.8+
- PyTorch
- pytorch-lightning
- hydra-core / omegaconf
- biopython
- numpy
- scipy
- scikit-learn
- pandas
- tqdm
- numba
- transformers
- tape-proteins
- fair-esm, or another locally available `esm` package

C++ side:

- C++11 compiler
- CMake 3.0+
- pkg-config
- ViennaRNA / RNAlib2
- OpenMP, recommended for parallel alignment

If `pkg-config` cannot find RNAlib2, set ViennaRNA's pkgconfig path:

```bash
export PKG_CONFIG_PATH=/path/to/ViennaRNA/lib/pkgconfig:$PKG_CONFIG_PATH
```

## Model Weights and Data

Large model weights, checkpoints, generated tensors, and search databases are not suitable for direct GitHub storage. This repository's `.gitignore` excludes common large artifacts, including:

- `RNA-ESM2/ckpt/*.pt`
- `RNA-ESM2/ckpt/*.ckpt`
- `RNA-ESM2/data/*.pt`
- `*.npz`, `*.npy`, `*_allmaps.bpmap`
- C++ build products and Python cache files

Prepare the required weights locally before running inference. Common paths referenced by scripts include:

```text
RNA-ESM2/ckpt/RUNGPURNA-ESM2-trans-mappro-KDNY-epoch-30-step-15650.pt
RNA-ESM2/ckpt/RNA-ESM2-trans-2a100-mappro-KDNY-epoch=07-valid_F1=0.561-v1.ckpt
RNA-ESM2/data/esm2_t30_150M_UR50D.pt
RNA-ESM2/data/esm2_t33_650M_UR50D.pt
```

Search database files are not included. Put genome or database FASTA files in a local path and pass that path with `-b`, `-genome`, or `-database`.

## RNA-ESM2 Inference

Enter the model directory:

```bash
cd RNA-ESM2
```

Common contact-map/profile inference scripts:

- `esm2embatt_650M_con.py`: 650M configuration with contact-map output
- `esm2embatt_150M_con.py`: 150M configuration with contact-map output
- `esm_profile_test_batch.py`: batch profile/contact export
- `esm_profile_test_cpu.py`: CPU inference variant

Hydra parameters can be overridden from the command line:

```bash
CUDA_VISIBLE_DEVICES=0 python esm2embatt_650M_con.py \
  produce.seqfile=./fasta/testRNA42.fasta \
  produce.outname=testRNA42 \
  data.device=cuda:0 \
  data.model_path=./ckpt/RUNGPURNA-ESM2-trans-mappro-KDNY-epoch-30-step-15650.pt
```

Common outputs:

- `<outname>_allmaps.bpmap`: contact-map text for all sequences
- `<outname>_mapdict.json`: sequence ID to map-index metadata
- `<outname>_seqdict.json`: sequence ID to sequence mapping
- `<outname>_contacts.npz`: contact tensor
- `<outname>_features.npz`: profile or feature tensor

Some research scripts contain hard-coded historical cluster paths, such as `/mnt/remote_home/...`. Before running them on a new machine, check `results_path`, input paths, output paths, and checkpoint paths.

## Model Training

`RNA-ESM2/main.py` is the training entry point and uses Hydra configuration. Default training paths include:

- `data.ffindex_path=data/UniRef30_2020_02_a3m`
- `data.trrosetta_path=data/trrosetta`
- `data.trrosetta_train_split=train.txt`
- `data.trrosetta_valid_split=test.txt`

Example:

```bash
cd RNA-ESM2
python main.py \
  data.ffindex_path=/path/to/UniRef30_2020_02_a3m \
  data.trrosetta_path=/path/to/trrosetta \
  train.gpus=1 \
  train.max_epochs=100
```

Training datasets are usually large. Keep them in local storage, object storage, or institutional storage instead of committing them to GitHub.

## Build Contact-Map Alignment Tools

Enter the C++ subproject and build:

```bash
cd profile-contact-alignment
cmake . -DCMAKE_BUILD_TYPE=Release
cmake --build . -j 8
```

Binaries are generated under:

```text
profile-contact-alignment/bin/
```

Main commands:

- `cmapbuild`: build a contact-map model from DBN or structure input
- `cmapsearch`: search genome/database FASTA files with a query structure
- `RESM_Search`: search or align a query structure against target sequences
- `RESM_Search_align`: contact-map or sequence alignment entry point
- `RESM_Search_alignNCBI`: NCBI-oriented search entry point
- `BPmap`: directly align two precomputed contact maps

Useful helper scripts:

```bash
scripts/build_and_test.sh
scripts/benchmark_performance.sh
python scripts/verify_optimizations.py
```

## Input Formats

DBN files use three lines per record:

```text
>sequence_id
AUGCUAGCUA
(((....)))
```

FASTA files use standard FASTA format:

```text
>sequence_id
AUGCUAGCUA
```

Contact-map files support:

```text
LEN 100
CON 3 45 1.0
PRF 0 A N 0.85 0.05 0.05 0.05
```

Rules:

- `LEN` or `SIZE` defines map length.
- `CON i j score` defines a contact between positions `i` and `j`.
- `PRF` provides optional profile information.
- If `PRF` is missing, `BPmap` uses `N` as a fallback sequence/profile symbol.

## Usage Examples

Build a cmap model from DBN:

```bash
./bin/cmapbuild -dbn query.dbn -o query.cmp
```

Build a cmap model from DBN plus MSA:

```bash
./bin/cmapbuild -dbn query.dbn -a2m query.a2m -o query.cmp
./bin/cmapbuild -dbn query.dbn -sto query.sto -o query.cmp
```

Search a genome FASTA file:

```bash
./bin/cmapsearch -dbn query.dbn -genome genome.fa
```

Search a database FASTA file:

```bash
./bin/cmapsearch -dbn query.dbn -database database.fa
```

Run RESM search with a query DBN and target FASTA:

```bash
./bin/RESM_Search -a query.dbn -b targets.fa -threads 8 -sep_cut 1
```

Directly align two precomputed contact maps:

```bash
./bin/BPmap -a query.bpmap -b target.bpmap -threads 8 -use_prf
```

Common options:

- `-a`: query RNA secondary structure or contact-map file
- `-b`: target FASTA or contact-map file
- `-threads` / `-num_threads`: OpenMP thread count, default `8`
- `-gap_o`: gap opening penalty, default `-1`
- `-gap_e`: gap extension penalty, default `-0.01`
- `-sep_cut`: minimum sequence-separation filter threshold, default `1`
- `-iter`: iteration count, default `20`
- `-use_prf`: enable profile scoring
- `-prf_w`: profile score weight
- `-use_gap_ss`: enable secondary-structure-aware gap adjustment
- `-gap_ss_w`: secondary-structure gap penalty weight
- `-silent`: reduce log output

## Recommended Workflow

1. Prepare a query RNA sequence and structure, or prepare FASTA input.
2. Use `RNA-ESM2` to infer profiles and contact maps.
3. Use `cmapbuild` or generated `.bpmap` files to build the query map.
4. Prepare genome/database FASTA files for the target search space.
5. Run `cmapsearch`, `RESM_Search`, or `BPmap`.
6. Rank candidates by alignment score or `BestScore`.

## Notes

- `T` may be converted to `U` by some RNA input readers.
- Some scripts are research-stage utilities and may include historical cluster paths or PBS settings. Check input, output, and checkpoint paths before use.
- Model inference can require substantial GPU memory. Use the CPU script or reduce batch size/sequence length when GPU memory is limited.
- Store large weights and databases outside GitHub, for example in object storage, shared storage, Hugging Face, GitHub Releases, or institutional storage. Record download locations in configuration or documentation.

## Version-Control Rules

Recommended content to commit:

- Python and C++ source code
- Build scripts and lightweight verification scripts
- Small example FASTA/DBN files
- README and project documentation

Content that should not be committed:

- `.pt` / `.ckpt` model weights
- Large FASTA databases
- `.npz` / `.npy` inference outputs
- `__pycache__`
- CMake build directories, binaries, and benchmark artifacts
