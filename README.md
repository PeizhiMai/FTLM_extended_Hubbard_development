# FTLM Extended Hubbard `n(mu)`

This is a fresh C++ implementation of the FTLM thermodynamic core, built by learning from:

- `reference_ftlm_hub_cond/fltm_hub_cond/hubTri2Dcond_omp.f`
- `reference_ftlm_hub_cond/fltm_hub_cond/cond_spect_omp.f`

The conductivity-specific parts from the Fortran code were intentionally dropped. The target here is the grand-canonical density curve `n` versus chemical potential `mu`, plus the charge-fluctuation compressibility, at fixed inverse temperature `beta` for an extended Hubbard model on a periodic rectangular lattice.

## Model

The code currently includes:

- rectangular `lx x ly` lattice with periodic boundary conditions
- nearest-neighbor hopping `t_x`, `t_y`
- second-neighbor diagonal hopping `t'`
- homogeneous twists `phi_x`, `phi_y` in units of `2*pi`
- onsite interaction `U`
- nearest-neighbor density interaction `V`

## RAM strategy

To avoid the RAM-heavy pattern from older local code, the implementation:

- never builds the full Hamiltonian matrix
- stores only spin-sector bitstring lists, not full many-body objects
- applies `H` on the fly during Lanczos
- keeps only three Lanczos vectors in memory

## Build

```bash
cmake -S . -B build
cmake --build build
```

On Apple Silicon, a practical OpenMP path is Apple Clang plus a native `libomp` install. If `libomp` is installed under `/opt/homebrew/opt/libomp`, configure with:

```bash
cmake -S . -B build-clang-omp \
  -DCMAKE_C_COMPILER=/usr/bin/clang \
  -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
  -DFTLM_LIBOMP_ROOT=/opt/homebrew/opt/libomp
cmake --build build-clang-omp
```

If `FTLM_LIBOMP_ROOT` is not provided or no compatible OpenMP runtime is found, the FTLM executable still builds and runs in serial mode.

## Example

```bash
./build/ftlm_n_vs_mu \
  --lx 2 --ly 2 \
  --tx 1.0 --ty 1.0 --tp 0.2 \
  --phix 1.0 --phiy 1.0 \
  --u 4.0 --v 1.0 \
  --beta 1.5 \
  --samples 5 \
  --lanczos-steps 40 \
  --mu-min -2 --mu-max 6 --mu-count 33 \
  --output n_vs_mu.csv
```

To reuse one FTLM/Lanczos pass for several inverse temperatures, pass a
comma-separated beta list:

```bash
./build/ftlm_n_vs_mu \
  --lx 4 --ly 4 \
  --u 7.0 \
  --beta-list 2.857142857142857,12.5 \
  --samples 128 \
  --lanczos-steps 80 \
  --mu-min -3 --mu-max 4 --mu-count 281 \
  --output n_vs_mu_multi_beta.csv
```

With `--beta-list`, the output CSV contains one row per `(beta, mu)` and adds a
leading `beta` column. With a single `--beta`, the legacy output format is
unchanged.

The default Lanczos depth is `80` steps per random vector.
The default FTLM sample count is `5`.
By default, the FTLM executable diagonalizes small momentum blocks exactly when
their dimension is at most `256`; change this with `--exact-block-threshold N`,
or set it to `0` to force stochastic FTLM on every non-empty block.

If OpenMP is available in the compiler toolchain, FTLM can parallelize over random samples:

```bash
./build/ftlm_n_vs_mu --threads 4 ...
```

Using more threads should reduce runtime with only a modest RSS increase because each thread only needs its own Lanczos work vectors for the active block.

The output CSV contains:

- `mu`
- `n` as electrons per site
- `charge_correlation` as the connected total charge fluctuation per site:
  `(<N^2> - <N>^2) / N_sites`
- `compressibility = beta * charge_correlation`, equivalent to `d n / d mu`
  in the grand-canonical ensemble
- `partition_like` as the shifted grand-canonical FTLM sum used internally

For checkpointing long cluster runs, pass a checkpoint path:

```bash
./build/ftlm_n_vs_mu \
  --lx 4 --ly 4 \
  --u 7 \
  --beta-list 2.857142857142857,12.5 \
  --samples 128 \
  --lanczos-steps 80 \
  --checkpoint twist_000.ftlmcp \
  --output twist_000.csv
```

The checkpoint stores compact per-block spectra after each completed
`(N_up,N_down,k)` block. Rerunning the same command with the same checkpoint
skips completed blocks and resumes from the next missing block. The companion
`*.meta` file records the model, twist, sampling, seed, and exact-block
threshold; a mismatched command fails early instead of mixing incompatible
spectra.

## Exact Diagonalization

A separate exact-diagonalization executable is also available:

```bash
./build/ftlm_ed_n_vs_mu \
  --lx 4 --ly 2 \
  --beta 20 \
  --mu-min -5 --mu-max 15 --mu-count 81 \
  --output run_lx4_ly2_mu_m5_p15_ed.csv
```

This ED path reuses the same fixed-particle sectors and momentum-block construction as the FTLM code, but it diagonalizes each active `(N_up, N_down, k)` block exactly instead of using Lanczos sampling.

The twist parameters `--phix` and `--phiy` are homogeneous twists in units of `2*pi`. The default is `phix=phiy=1`, which is gauge-equivalent to no twist.

## 4x2 FTLM vs ED Comparison

For the `4x2` cluster at `beta=20`, `U=8`, `V=0`, `tx=ty=1`, the workspace now includes:

- FTLM result:
  - `run_lx4_ly2_mu_m5_p15_ftlm_compare.csv`
- ED result:
  - `run_lx4_ly2_mu_m5_p15_ed.csv`
- pointwise comparison:
  - `run_lx4_ly2_mu_m5_p15_ftlm_vs_ed.csv`
- overlay plot:
  - `run_lx4_ly2_ftlm_vs_ed.png`
- error plot (`FTLM - ED`):
  - `run_lx4_ly2_ftlm_minus_ed.png`

The pointwise comparison file contains:

- `mu`
- `n_ftlm`
- `n_ed`
- `diff_ftlm_minus_ed`

## Important current limitation

This version is intentionally clean and memory-aware, but it does not yet reproduce the translational symmetry reduction used in the Fortran parent-state code. That means the Hilbert-space vectors still scale with the full fixed-`(N_up,N_down)` sector size. A guard is included through `--max-sector-dim` so large sectors fail early instead of exhausting memory.

## Thread Handoff Notes

This section summarizes what was done in this workspace so a new Codex thread can pick up quickly.

### Original goal

Build a fresh C++ FTLM code for the extended Hubbard model on a periodic rectangular lattice, using the provided Fortran reference only as the algorithmic guide, with the immediate physics target:

- density `n` versus chemical potential `mu`
- at fixed inverse temperature `beta`

### Reference used

The implementation work was guided only by:

- `reference_ftlm_hub_cond/fltm_hub_cond/hubTri2Dcond_omp.f`
- `reference_ftlm_hub_cond/fltm_hub_cond/cond_spect_omp.f`

The main ideas taken from the reference were:

- FTLM over fixed particle-number sectors
- translational parent-state construction
- momentum-sector cancellation and normalization via orbit sums
- grand-canonical reweighting over `mu`

### What the C++ code does now

The current executable is:

- `build/ftlm_n_vs_mu`

and the source is:

- `src/main.cpp`

Current features:

- rectangular periodic lattice `lx x ly`
- nearest-neighbor hopping `tx`, `ty`
- second-neighbor diagonal hopping `tp`
- homogeneous twists `phix`, `phiy`
- onsite interaction `U`
- nearest-neighbor density interaction `V`
- FTLM sampling
- exact diagonalization fallback for small `(N_up,N_down,k)` blocks
- translation-symmetry parent-state reduction
- momentum-block (`k`) construction from translation orbits
- direct CSV output for `n(mu)` at fixed `beta`
- charge-fluctuation compressibility from the grand-canonical charge correlation

### Memory-optimization history

The implementation went through several rounds of memory work.

#### Stage 1

Initial C++ version:

- used fixed `(N_up,N_down)` sectors
- did not yet reduce into momentum blocks
- stored all FTLM spectra globally before the final `mu` loop

This was clean but not competitive in RAM use.

#### Stage 2

Added symmetry reduction:

- translational parent states
- momentum-sector cancellation and normalization
- explicit `k`-block Lanczos runs

This improved structure and runtime, but not RSS enough by itself.

#### Stage 3

First serious RAM pass:

- streamed each `k` block directly into the `n(mu)` accumulators
- stopped storing all spectra globally
- processed one block at a time
- replaced nested row storage with a compact CSR-like structure:
  - `row_ptr`
  - `col_idx`
  - `values`

#### Stage 4

Second RAM pass:

- compressed stored hops from full target states to compact orbit references:
  - target parent index
  - translation shift
  - fermionic sign
- kept the state-to-relation map only locally during sector construction

This gave another real RSS reduction.

#### Stage 5

Tried replacing the relation hash table with flat sorted state arrays and binary search.

Result:

- it did not help RSS on the tested `4x2` case
- it slightly increased memory

That change was reverted.

### Current best C++ memory result

The current code in `src/main.cpp` is the reverted-to-good state after the successful memory reductions.

Best measured C++ number on the `4x2` cluster:

- command:
  - `/usr/bin/time -l ./build/ftlm_n_vs_mu --lx 4 --ly 2 --beta 20 --mu-min -5 --mu-max 5 --mu-count 41 --output run_lx4_ly2_beta20_streamed2.csv`
- peak RSS:
  - `10,715,136` bytes
  - about `10.22 MiB`

Earlier C++ memory numbers on the same `4x2` case:

- original symmetry-reduced version:
  - `18,006,016` bytes
- after first streaming / CSR pass:
  - `12,271,616` bytes
- after compact-hop pass:
  - `10,715,136` bytes
- flat-array lookup experiment:
  - `10,813,440` bytes
  - reverted

### Reference Fortran memory number

A temporary gfortran-compatible copy of the reference `8`-site code was built and run in:

- `tmp_fortran_8/hubTri2Dcond_omp_8.f`

Measured reference result for the `8`-site run:

- peak RSS:
  - `7,467,008` bytes
  - about `7.12 MiB`

So the current C++ code is still above the reference by about:

- `3.25 MiB` on the tested `8`-site case

### Other measured C++ runs

Current C++ code was also run on:

- `3x2`, `beta=20`
- `4x2`, `beta=20`
- `3x3`, `beta=20`

Representative files:

- `run_lx3_ly2_beta20_kblocks.csv`
- `run_lx4_ly2_beta20_streamed2.csv`
- `run_lx3_ly3_beta20_kblocks.csv`

### Important caveats for a future thread

- The current code is still more memory-heavy than the Fortran reference.
- The largest remaining overhead is likely in sector construction:
  - `unordered_map<StateKey, Relation, ...>`
  - temporary orbit bookkeeping
  - repeated translation work
- The `3x3` run showed mild unphysical tail behavior at large `mu`, so numerical quality should still be treated as under active development.
- The Fortran reference computes much more than `n(mu)`:
  - optical conductivity workflow
  - second Lanczos
  - file output for postprocessing
  so direct runtime comparisons are not apples-to-apples.

### Best next steps

If a new thread continues optimization, the highest-value next steps are:

- remove or shrink the `unordered_map` orbit lookup without reintroducing the failed flat-array overhead
- reduce temporary allocation during sector construction
- reuse translation buffers and scratch vectors
- consider direct combinatorial ranking/unranking for state lookup
- validate `n(mu)` more carefully against exact diagonalization on very small clusters
- check whether the FTLM weighting can be made closer to the Fortran normalization conventions in the large-`mu` tails
