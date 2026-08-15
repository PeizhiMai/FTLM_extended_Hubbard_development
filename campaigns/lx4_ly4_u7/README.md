# 4x4 U=7 compressibility campaign

This directory fixes the production parameters used for the CADES campaign.
It evaluates 66 symmetry-unique representatives of a 20x20 Gamma-centered
quadrature of the 4x4-cluster reduced Brillouin zone (rBZ), with physical
spacing `Delta kx = Delta ky = pi/40`:

- `lx=ly=4`, `U=7`, `tx=ty=1`, `tp=V=0`
- one recurrence to `m_max=300`, with reusable prefix spectra at
  `m=80,120,160,200,250,300`; exact momentum-block threshold `256`
- production target `R=36` (with retained `R=8,16,32` reductions)
- `beta=2.857142857142857,12.5` and 281 points on `mu=[-3,4]`
- physical twists in `0 <= kx <= ky <= pi/4`, listed in `twists.csv`
- periodic D4 multiplicities 1, 2, 4, or 8, giving 400 effective rBZ points
  from 66 jobs
- one exclusive `burst` node per twist, 36 OpenMP workers, **250 GiB**,
  four hours, using the `default` QOS with no Slurm dependencies or
  campaign-side concurrency throttle

The executable's `phix` and `phiy` are boundary-flux fractions, not physical
momenta.  For this 4x4 cluster,

```text
kx = 2*pi*phix/4,  ky = 2*pi*phiy/4.
```

The full periodic rBZ is `[-pi/4,pi/4) x [-pi/4,pi/4)`. Its one-dimensional
mesh is `kx/pi = -10/40,-9/40,...,9/40`; it contains zero and has 20 distinct
periodic points. Independent sign reflections and axis exchange reduce it to
the closed triangular representative wedge `0 <= kx <= ky <= pi/4`. The
`+pi/4` representatives stand for the single periodic boundary point
`-pi/4 == +pi/4 (mod pi/2)`, so that boundary is not double-counted. In CLI
units, `phix,phiy = 0,1/20,...,10/20`. `twist_quadrature.json` records the
convention and expected multiplicities.

Existing `twist_005`, `(phix,phiy)=(0.25,0.25)`, is exactly the production-grid
point `(kx,ky)=(5*pi/40,5*pi/40)`. Its seed and append-only checkpoint are
therefore retained. The zero-twist point keeps ID `000`; the other 64
representatives use IDs `100` through `163`, avoiding all prior checkpoint
directories. The canceled, superseded 8x8 pilot
`twist_804` is not in this manifest and must never be averaged into production.
The largest block remains `(Nup,Ndown,mx,my)=(8,8,0,0)`, with dimension
`10,353,252`.

The twist reduction does **not** remove many-body total-momentum blocks. For a
generic rBZ twist, every `(mx,my)` block must still be included; the D4 weights
only avoid separate jobs at symmetry-related twist vectors.

## Why OpenMP, MKL, and no MPI

The production binary uses OpenMP over random vectors so all workers share one
sector basis and one sparse momentum-block Hamiltonian. MPI ranks on the same
node would independently rebuild those dominant data structures and multiply
memory use, so MPI is intentionally not used. Intel MKL is linked in sequential
mode (`MKL_NUM_THREADS=1`) for LAPACK work; the large-block cost is the sparse
Lanczos matvec, not dense LAPACK. The CADES build gate can compare Intel and GCC,
but production defaults to Intel oneAPI 2024.1 plus MKL 2024.1.

## Deployment and gate

Deployment refuses a dirty git tree and installs a commit-addressed source and
binary snapshot:

```bash
campaigns/lx4_ly4_u7/cades/deploy_to_cades.sh intel
```

The original `m=80` pilot uses the legacy checkpoint and may remain active.
Submit the accuracy-first `m_max=300` largest-block probe as a separate job and
separate `mext` checkpoint:

```bash
CAMPAIGN_ROOT=/lustre/or-scratch24/scratch/9pm/ftlm_codex/lx4_ly4_u7_campaign \
  source/campaigns/lx4_ly4_u7/cades/submit_m300_resource_probe.sh
```

After the job finishes, apply the gate:

```bash
/usr/bin/python3.11 \
  source/campaigns/lx4_ly4_u7/cades/validate_m300_resource_probe.py
```

The validator requires all sample IDs `0-15`, dimension `10,353,252`, all six
prefix records for every sample, and at least 5 GiB RSS headroom inside the
250-GiB request. Sample time is reported rather than used as an arbitrary
accuracy cutoff. It writes `RESOURCE_GATE_M300_PASSED.json`, without which
multi-prefix production submission is refused.

## Production waves

The active twist-005 checkpoint is both the resource/convergence pilot and one
of the 66 production representatives. Complete and validate it before launching
the other 65:

```bash
/usr/bin/python3.11 source/campaigns/lx4_ly4_u7/cades/submit_wave.py \
  --samples 36 --threads 36 --twist-id 005 --checkpoint-series mext \
  --lanczos-max-steps 300 \
  --lanczos-save-steps 80,120,160,200,250,300
/usr/bin/python3.11 source/campaigns/lx4_ly4_u7/cades/campaign_status.py \
  --samples 36 --checkpoint-series mext --lanczos-steps 300
/usr/bin/python3.11 source/campaigns/lx4_ly4_u7/validate_thermo_csv.py \
  runs/twist_005/twist_005_thermo_m300_R036.csv \
  --output runs/twist_005/validation_m300_R036.json
```

Each four-hour job stops internally after 210 minutes, appends every completed
sample durably, and exits cleanly. Run the same pilot command again for another
wave until twist 005 is complete. Each sample performs 300 Hamiltonian
applications once; the six requested tridiagonal prefix spectra are then saved
and all six CSV reductions are generated by the same completed job. Only after
that pilot and the resource gate
pass, submit every remaining representative concurrently:

```bash
/usr/bin/python3.11 source/campaigns/lx4_ly4_u7/cades/submit_wave.py \
  --samples 36 --threads 36 --checkpoint-series mext --lanczos-max-steps 300 \
  --lanczos-save-steps 80,120,160,200,250,300
```

Production checkpoints remain fixed at target `R=36`; never run concurrent
writers against one checkpoint. The same checkpoint can be reduced exactly at
`R=8`, `R=16`, `R=32`, or `R=36`. A future aggregate `R=72` calculation must
use two deterministic, disjoint `R=36` checkpoint shards plus a validated
merge/reducer rather than setting one job to `--samples 72`.

## Reduction, averaging, and the x axis

Twists are combined at fixed `(beta,mu)`: average `n` and `kappa`, then set
`x=1-mean(n)`. Reflection maps `kx` to `-kx` (and independently for `ky`),
while square-lattice axis exchange maps `(kx,ky)` to `(ky,kx)`. We therefore
run only `0 <= kx <= ky <= pi/4`. Generic off-diagonal interior points have
weight 8; axes, diagonals, and points touching the periodic boundary have
smaller orbits and weights 1, 2, or 4. The weights in `twists.csv` sum to 400,
exactly reproducing the 20x20 Gamma-centered periodic quadrature. An unweighted
average over the 66 files would be incorrect.
Fixed-`mu` averaging also handles the fact that individual twists have
different densities without interpolating mismatched thermodynamic control
parameters.

Once all twist CSVs exist:

```bash
for m in 80 120 160 200 250 300; do
  /usr/bin/python3.11 source/campaigns/lx4_ly4_u7/finalize_campaign.py \
    --campaign-root "$CAMPAIGN_ROOT" --samples 36 \
    --checkpoint-series mext --lanczos-steps "$m" --skip-plot
done
```

The finalizer requires all 66 representatives (effective multiplicity 400), two
betas, 281 rows per beta, monotone density, nonnegative compressibility, and
complete `x=[0,0.35]` coverage. It writes the symmetry-weighted observable
average. CADES' system Python does not include
Matplotlib, so copy that average to the local workspace and make the
publication-style PNG/PDF with the requested environment:

```bash
source ~/.venvs/myenv/bin/activate
python scripts/plot_compressibility_vs_hole_doping.py \
  twist_average_rbz20_gamma_dk_pi40_m300_R036.csv \
  --output compressibility_vs_x_rbz20_gamma_dk_pi40_m300_R036.png
```

Because beta and mu are reducer inputs rather than checkpoint metadata,
retained checkpoints can produce new temperature or chemical-potential grids
without more Lanczos work.
