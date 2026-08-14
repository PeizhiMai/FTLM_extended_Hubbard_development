# FTLM Extended Hubbard `n(mu)`

This is a fresh C++ implementation of the FTLM thermodynamic core, built by learning from:

- `reference_ftlm_hub_cond/fltm_hub_cond/hubTri2Dcond_omp.f`
- `reference_ftlm_hub_cond/fltm_hub_cond/cond_spect_omp.f`

The first target was the grand-canonical density curve `n` versus chemical potential `mu`, plus the charge-fluctuation compressibility. The code now also includes exact-ED and two-Lanczos FTLM machinery for regular optical/DC conductivity.

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
- stores the active momentum block as sparse CSR and uses a compact 12-byte
  sector hopping record
- drops zero-amplitude hopping directions (notably all `t'=0` diagonals)
- keeps only three Lanczos vectors in memory
- loads checkpoint spectra in index-only mode during generation

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
  --samples 16 \
  --lanczos-steps 80 \
  --mu-min -3 --mu-max 4 --mu-count 281 \
  --output n_vs_mu_multi_beta.csv
```

With `--beta-list`, the output CSV contains one row per `(beta, mu)` and adds a
leading `beta` column. With a single `--beta`, the legacy output format is
unchanged.

The default Lanczos depth is `80` steps per random vector.
The default FTLM sample count is `5`. With an extensible v2 checkpoint,
`--samples R` means the target total number of permanent sample IDs per
stochastic momentum block.
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
- `log_partition` as the absolute `log(Z)` reconstructed from the shifted sum

For twist averaging, run independent jobs with different `--phix/--phiy`, then
combine their CSV files:

```bash
python3 scripts/average_twist_outputs.py \
  twist_*.csv \
  --mode observable-average \
  --expected-twists 16 \
  --output twist_average.csv
```

Rows are combined at fixed `(beta,mu)` and the output includes
`x = 1 - mean(n)`, twist standard deviations, and standard errors. Thus twists
that give different individual densities at the same `mu` do not require
interpolation before the thermodynamic average.

If symmetry-related twists are represented by one file, pass one positive
integer multiplicity per input. For example, square-lattice representatives
with diagonal weight one and off-diagonal weight two can reconstruct the full
grid with `--weights 1,2,... --expected-weight 16`. `--expected-twists` still
counts independently supplied files.

`observable-average` is the usual equal-weight average of measured observables
over twists. If you want to combine the twists by summing their partition
functions instead, use the absolute `log_partition` column:

```bash
python3 scripts/average_twist_outputs.py \
  twist_*.csv \
  --mode partition-sum \
  --sites 16 \
  --output twist_partition_sum.csv
```

For checkpointing long cluster runs, pass a checkpoint path:

```bash
./build/ftlm_n_vs_mu \
  --lx 4 --ly 4 \
  --u 7 \
  --beta-list 2.857142857142857,12.5 \
  --samples 16 \
  --lanczos-steps 80 \
  --checkpoint twist_000.ftlmcp \
  --max-runtime-minutes 210 \
  --progress-jsonl twist_000.progress.jsonl \
  --progress-state twist_000.state.json \
  --output twist_000.csv
```

New thermodynamic checkpoints use append-only format v2. Every completed
`(N_up,N_down,k,sample_id)` record is checksummed and fsynced independently, and
stores raw first-vector overlaps rather than weights already divided by `R`.
The RNG and seed derivation are versioned in `*.meta`. Rerunning at the same
target resumes only missing IDs; increasing `--samples 16` to `--samples 32`
generates only IDs `16-31`. Target `R` is deliberately not immutable metadata,
while lattice/model/twist/seed/Lanczos depth/exact threshold/RNG version are.

The process responds to `SIGUSR1`/`SIGTERM` by starting no more samples and
finishing durable active samples. Progress is printed after every checkpoint and
at the requested heartbeat interval. A partially written final record is
truncated to its last verified boundary on restart.

Reduce or inspect a checkpoint without reconstructing Hilbert-space bases:

```bash
./build/ftlm_reduce_checkpoint \
  --checkpoint twist_000.ftlmcp --samples 16 \
  --beta-list 2.857142857142857,12.5 \
  --mu-min -3 --mu-max 4 --mu-count 281 \
  --output twist_000_R16.csv

./build/ftlm_reduce_checkpoint \
  --checkpoint twist_000.ftlmcp --samples 32 --status
```

Status mode reports every missing sample ID and the next unit of work. Reduction
at `R=32` is refused until every stochastic block has IDs `0-31`; lower-R
reductions remain byte-reproducible after extension. Legacy v1 checkpoints are
still readable at their original `R`, but are read-only and cannot be extended.

The complete 4x4 CADES workflow, fixed twist manifest, resource gate, wave
submission, validation, and final plotting are documented in
`campaigns/lx4_ly4_u7/README.md`.

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

An exact optical/DC conductivity executable is available for small clusters and
benchmarking:

```bash
./build/ftlm_ed_conductivity \
  --lx 3 --ly 2 \
  --u 7 \
  --phix 0.2 --phiy 0.3 \
  --beta-list 2.857142857142857,12.5 \
  --mu-min -3 --mu-max 4 --mu-count 281 \
  --omega-min 0 --omega-max 12 --omega-count 1201 \
  --eta 0.05 \
  --conductivity-direction x \
  --output optical_conductivity.csv \
  --dc-output dc_conductivity.csv
```

The optical CSV contains:

- `beta`, `mu`, `omega`
- `sigma`, the Lorentzian-broadened regular optical conductivity
- `n`, `charge_correlation`, `compressibility`, and `log_partition`

The DC CSV reports `sigma_dc = sigma(omega=0)` using the same broadening
parameter `eta`. This implementation evaluates the Kubo formula exactly inside
each momentum block and is therefore intended for clusters/blocks that can be
fully diagonalized. Use `--max-conductivity-block-dim N` to guard memory.

The twist parameters `--phix` and `--phiy` are homogeneous twists in units of `2*pi`. The default is `phix=phiy=1`, which is gauge-equivalent to no twist.


## Two-Lanczos FTLM Conductivity

For larger clusters, use the two-Lanczos current-current executable:

```bash
./build/ftlm_conductivity \
  --lx 4 --ly 4 \
  --u 7 \
  --phix 0.2 --phiy 0.3 \
  --beta-list 2.857142857142857,12.5 \
  --mu-min -3 --mu-max 4 --mu-count 281 \
  --samples 128 \
  --lanczos-steps 80 \
  --omega-min 0 --omega-max 12 --omega-count 1201 \
  --eta 0.05 \
  --conductivity-direction x \
  --conductivity-checkpoint twist_000.condcp \
  --conductivity-memory-mode hybrid \
  --conductivity-basis-memory-mb 512 \
  --max-runtime-minutes 230 \
  --output optical_conductivity.csv \
  --dc-output dc_conductivity.csv
```

The executable uses exact ED conductivity for momentum blocks with
`basis_dim <= --exact-block-threshold` and two Lanczos runs per random vector for
larger blocks. In the first Lanczos run the starting vector is the random FTLM
state; in the second Lanczos run the starting vector is the normalized current
state `J|r>`. The compact per-sample records are appended to
`--conductivity-checkpoint`, so cluster jobs can be restarted without repeating
completed `(N_up,N_down,k,sample)` records.

The default memory mode is `hybrid`: if storing both Lanczos bases is estimated
to fit under `--conductivity-basis-memory-mb`, the code stores them for speed;
otherwise it recomputes basis vectors while forming the projected current matrix.
Use `--conductivity-memory-mode store-bases` or `recompute` to force either path.

If `--max-runtime-minutes` is set, the run exits cleanly after the next completed
sample record and can be resumed by rerunning the same command. Final optical and
DC CSVs are written only once all requested samples are present. The checkpoint
metadata validates generation-critical parameters such as lattice, model, twist,
direction, sample count, Lanczos depth, seed, and exact-block threshold; the raw
records remain reusable for different beta, mu, omega, or eta postprocessing
grids.

The optical CSV contains:

- `beta`, `mu`, `omega`
- `sigma`, the Lorentzian-broadened regular optical conductivity
- `n`, `charge_correlation`, `compressibility`, and `log_partition`

The DC CSV contains the same thermodynamic columns and `sigma_dc = sigma(0)`.
This first implementation reports only the regular broadened conductivity;
diamagnetic tau/Drude-weight diagnostics are intentionally left for a later
extension.

Twist-resolved optical or DC conductivity CSVs can be combined with the same
averaging script used for thermodynamics. For optical data, rows are grouped by
`(beta, mu, omega)`; for DC data they are grouped by `(beta, mu)`:

```bash
python3 scripts/average_twist_outputs.py \
  twist_*_optical.csv \
  --mode observable-average \
  --output optical_twist_average.csv
```

For partition-weighted conductivity averages, pass `--mode partition-sum --sites N`.

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

Lanczos vectors are reduced into `(N_up,N_down,k_x,k_y)` momentum blocks, but
sector construction still enumerates the full fixed-`(N_up,N_down)` state
space once to form translation orbits and compact hopping relations. For the
4x4 half-filled sector this means scanning `165,636,900` product states even
though the largest momentum block has dimension `10,353,252`. The
`--max-sector-dim` guard therefore remains important.

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

#### Stage 6

For the 4x4 campaign, the node-based state-relation hash was replaced by a
dense combination-rank lookup tailored to the small spin-bitmask address
space. Orbit relations are now stored in one flat Cartesian-rank array and
packed into at most eight bytes each. The 4x4 half-filled relation table is
therefore about `1.23 GiB`, with no per-state hash-node allocation. Checkpoints
and thermodynamic CSVs were verified byte-for-byte against the previous hash
implementation on 2x2, 3x3, and a selected 3x4 block.

### Small-cluster C++ memory history

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

At that historical stage the C++ code was above the reference by about:

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

- Small-cluster RSS comparisons predate the Stage 6 flat-relation change and
  should not be extrapolated to 4x4.
- The largest remaining costs are full-sector enumeration, repeated
  translations, the compact hopping table, and the momentum-block Lanczos
  vectors.
- The `3x3` run showed mild unphysical tail behavior at large `mu`, so numerical quality should still be treated as under active development.
- The Fortran reference computes much more than `n(mu)`:
  - optical conductivity workflow
  - second Lanczos
  - file output for postprocessing
  so direct runtime comparisons are not apples-to-apples.

### Best next steps

If a new thread continues optimization, the highest-value next steps are:

- avoid rescanning the full product sector where direct orbit generation is possible
- reduce the compact hopping-table footprint or generate rows on demand
- reuse translation buffers and scratch vectors
- benchmark direct combinatorial ranking against the current dense spin-rank table
- validate `n(mu)` more carefully against exact diagonalization on very small clusters
- check whether the FTLM weighting can be made closer to the Fortran normalization conventions in the large-`mu` tails
